#include <intrinsicslim.hpp>
#include <distortion_energy.hpp>
#include <iostream>
#include <intrinsicgrad.hpp>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/writeOBJ.h>

double compute_symdir_energy(DataGeo &data_mesh, const Eigen::MatrixXd &UV){

  auto energy = symmetric_dirichlet;
  Eigen::SparseMatrix<double> Dx, Dy;
  Eigen::VectorXd areas;

  computeGrad_intrinsic(data_mesh, Dx, Dy, areas);

  Eigen::VectorXd Dxu = Dx * UV.col(0);
  Eigen::VectorXd Dxv = Dx * UV.col(1);
  Eigen::VectorXd Dyu = Dy * UV.col(0);
  Eigen::VectorXd Dyv = Dy * UV.col(1);

  double total_energy = 0;
  unsigned flipped_triangles = 0;
  double total_area = 0;
  for(int i=0; i<data_mesh.intTri->intrinsicMesh->nFaces(); ++i){
    Eigen::Matrix2d J;
		J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
    if(J.determinant()<0){
      flipped_triangles++;
    }
    double locen = energy(J)*areas(i);
    total_area+=areas(i);
    total_energy += locen;
  }
  return total_energy/total_area;
}

unsigned slim_tillconverges(DataGeo &data_mesh, igl::SLIMData& slimdata, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &UV_init, unsigned max_iterations, bool igrad){
  if(!slimdata.has_pre_calc){
    Eigen::MatrixXd fixed_UV_positions;
    Eigen::VectorXi fixed_UV_indices;
    igl::boundary_loop(F, fixed_UV_indices);
    igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
    igl::slim_precompute(V, F, UV_init, slimdata, igl::MappingEnergyType::SYMMETRIC_DIRICHLET,
        fixed_UV_indices, fixed_UV_positions, 0);
    //std::cout << "Initial energy: " << slimdata.energy <<
    //  "; intri en: " << compute_symdir_energy(data_mesh, UV_init)  << std::endl;
  }

  if(igrad){
    Eigen::VectorXd areas;
    Eigen::SparseMatrix<double> Dx, Dy;
    computeGrad_intrinsic(data_mesh, Dx, Dy, areas);

    slimdata.Dx = Dx;
    slimdata.Dy = Dy;
    slimdata.mesh_area = areas.sum();
    slimdata.M = areas*2; // in igl::slim this was the way
    slimdata.F = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
    //std::cout << "slim energy: " << slimdata.energy <<
    //  "; intri en: " << compute_symdir_energy(data_mesh, new_UV)  << std::endl;
    slimdata.energy = compute_symdir_energy(data_mesh, UV_init) * 2;
  }
   // std::cout << "Inital energy: " << slimdata.energy << std::endl;

  double tol = 1e-8;
  double past_energy = 0;
  double curr_energy = slimdata.energy;
  unsigned itr = 0;
  while(itr < max_iterations && std::abs(past_energy-curr_energy)>tol){
    igl::slim_solve(slimdata, 1);
    past_energy = curr_energy;
    curr_energy = slimdata.energy;
    ++itr;
   // std::cout << " Energy itr. " << itr << " : " << curr_energy << std::endl;
  }
  return itr;
}



unsigned intrinsicslim(DataGeo &data_mesh, Eigen::MatrixXd &UV_init, Eigen::MatrixXd &UV, unsigned slim_maxitr, unsigned intrinsic_maxitr, std::ostream &fout){
  igl::SLIMData slimdata;

  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;
  Eigen::MatrixXd V = data_mesh.V;
  Eigen::MatrixXi F = data_mesh.F;
  auto flip_func = edgeorder_flip;
  auto start = std::chrono::high_resolution_clock::now();
  if(!slimdata.has_pre_calc){
    Eigen::MatrixXd fixed_UV_positions;
    Eigen::VectorXi fixed_UV_indices;
    igl::boundary_loop(F, fixed_UV_indices);
    igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
    igl::slim_precompute(V, F, UV_init, slimdata, igl::MappingEnergyType::SYMMETRIC_DIRICHLET,
        fixed_UV_indices, fixed_UV_positions, 0);
    std::cout << "Initial energy: " << slimdata.energy <<
      "; intri en: " << compute_symdir_energy(data_mesh, UV_init)  << std::endl;
  }

  // first extrinsic
  unsigned total_iterations = slim_tillconverges(data_mesh, slimdata, V, F, UV_init, slim_maxitr, true);

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  double curr_energy =  slimdata.energy;
  double past_energy =  compute_symdir_energy(data_mesh, UV_init)*2;
  double tol = 1e-8;
  unsigned max_iterations = intrinsic_maxitr;
  unsigned itr = 0;

  // dividing by 2 because double_area usage in igl::slim
  fout << past_energy/2 << "," << total_iterations << "," << duration << "," << curr_energy/2 << ",";


  while(itr<max_iterations && std::abs(past_energy-curr_energy)>tol){
    UV = slimdata.V_o;

    //intrinsic flipping
    unsigned flips, del;
    unsigned total_flips = 0;
    unsigned total_del_flips = 0;

    start = std::chrono::high_resolution_clock::now();
    while(flips=flip_func(data_mesh, UV, del, EnergyType::SYMMETRIC_DIRICHLET)){
        total_flips+=flips;
        total_del_flips+=del;
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    fout << total_flips << "," << total_del_flips << "," << duration << ",";
    fout << compute_symdir_energy(data_mesh, UV) << ",";

    start = std::chrono::high_resolution_clock::now();
    
    total_iterations = slim_tillconverges(data_mesh, slimdata, V, F, UV, slim_maxitr, true);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    fout << total_iterations << "," << duration << "," << slimdata.energy/2 << ",";
    past_energy = curr_energy;
    curr_energy = slimdata.energy;
    ++itr;
    std::cout << "Intrinsic Itr " << itr << ": " << curr_energy << std::endl; 
  }
  UV = slimdata.V_o;
  return itr;
}


unsigned intrinsicslim(DataGeo &data_mesh, Eigen::MatrixXd &UV_init, Eigen::MatrixXd &UV, unsigned slim_maxitr, unsigned intrinsic_maxitr, std::ostream &fout, std::string path, std::string mesh_name){
  igl::SLIMData slimdata;

  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;
  Eigen::MatrixXd V = data_mesh.V;
  Eigen::MatrixXi F = data_mesh.F;
  auto flip_func = edgeorder_flip;
  if(!slimdata.has_pre_calc){
    Eigen::MatrixXd fixed_UV_positions;
    Eigen::VectorXi fixed_UV_indices;
    igl::boundary_loop(F, fixed_UV_indices);
    igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
    igl::slim_precompute(V, F, UV_init, slimdata, igl::MappingEnergyType::SYMMETRIC_DIRICHLET,
        fixed_UV_indices, fixed_UV_positions, 0);
    //std::cout << "Initial energy: " << slimdata.energy <<
    //  "; intri en: " << compute_symdir_energy(data_mesh, UV_init)  << std::endl;
  }

  double past_energy =  0;
  double curr_energy =  slimdata.energy;
  double tol = 1e-8;
  unsigned max_iterations = intrinsic_maxitr;
  unsigned itr = 0;
  
  while(itr<max_iterations && std::abs(past_energy-curr_energy)>tol){
    UV = slimdata.V_o;

    //intrinsic flipping
    unsigned flips, del;
    unsigned total_flips = 0;
    unsigned total_del_flips = 0;

    auto start = std::chrono::high_resolution_clock::now();
    while(flips=flip_func(data_mesh, UV, del, EnergyType::SYMMETRIC_DIRICHLET)){
        total_flips+=flips;
        total_del_flips+=del;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    fout << total_flips << "," << total_del_flips << "," << duration << ",";
    fout << compute_symdir_energy(data_mesh, UV) << ",";

    start = std::chrono::high_resolution_clock::now();
    
    unsigned total_iterations = slim_tillconverges(data_mesh, slimdata, V, F, UV, slim_maxitr, true);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    fout << total_iterations << "," << duration << "," << slimdata.energy/2 << ",";
    past_energy = curr_energy;
    curr_energy = slimdata.energy;
    ++itr;

    Eigen::MatrixXd CN;
    Eigen::MatrixXi FN;
    std::string to_store = path + "/" + mesh_name + std::to_string(itr);
    to_store += ".obj";
    igl::writeOBJ(to_store, data_mesh.V, data_mesh.F, CN, FN, slimdata.V_o, data_mesh.F);
  }
  UV = slimdata.V_o;
  return itr;
}

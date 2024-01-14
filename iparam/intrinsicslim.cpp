#include <intrinsicslim.hpp>
#include <distortion_energy.hpp>
#include <iostream>
#include <intrinsicgrad.hpp>
#include <igl/slim.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>

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


Eigen::MatrixXd intrinsicslim(DataGeo &data_mesh, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &UV_init,
    unsigned number_iterations, unsigned flip_granularity, const FlipType& ft){
  igl::SLIMData slimdata;

  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy; 
  
  auto flip_func = edgeorder_flip;
  if(ft==FlipType::GREEDY) flip_func = greedy_flip;
  if(ft==FlipType::HEURISTIC) flip_func = heuristic_flip;
  if(ft==FlipType::RANDOM) flip_func = random_flip;

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

  unsigned iterations = number_iterations / flip_granularity;
  unsigned rem = number_iterations % flip_granularity;
  
  while(iterations--){
    Eigen::MatrixXd new_UV = igl::slim_solve(slimdata, flip_granularity);
    /*
    //intrinsic flipping
    unsigned total_flips;
    while(total_flips=flip_func(data_mesh, new_UV, EnergyType::SYMMETRIC_DIRICHLET)) 
      std::cout << "flips: " << total_flips << std::endl; 
    computeGrad_intrinsic(data_mesh, Dx, Dy, areas);
    slimdata.Dx = Dx;
    slimdata.Dy = Dy;
    slimdata.mesh_area = areas.sum();
    slimdata.M = areas*2; // in igl::slim this was the way
    slimdata.F = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
    //std::cout << "slim energy: " << slimdata.energy << 
    //  "; intri en: " << compute_symdir_energy(data_mesh, new_UV)  << std::endl;
    slimdata.energy = compute_symdir_energy(data_mesh, new_UV) * 2;
    */
  }
  //std::cout << "slim energy: " << slimdata.energy << std::endl;
  return slimdata.V_o;
}

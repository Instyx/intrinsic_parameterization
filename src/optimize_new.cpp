#include "optimize_new.hpp"

#include <chrono>
#include <filesystem>

#include <igl/slim.h>
#include <igl/writeOBJ.h>

#include "energy.hpp"
#include "intrinsicslim.hpp"
#include "parameterization.hpp"
#include "store_intrinsic_mesh.hpp"
#include "harmonic.hpp"
#include "lscm.hpp"
#include "tutte.hpp"
#include "adjacency.hpp"
#include "distortion_energy.hpp"
#include "map_to_boundary.hpp"

#include <iostream>
#include <fstream>

inline void saveData(std::string fileName, Eigen::MatrixXd matrix)
{
	//https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
	const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

	std::ofstream file(fileName);
	if (file.is_open())
	{
		file << std::setprecision(128) << matrix.format(CSVFormat);
		file.close();
	}
}

Results optimize_single_new(Eigen::MatrixXd &V, Eigen::MatrixXi &F, EnergyType method, std::string dir, std::string mesh_name){
  return optimize_single(V, F, method, dir, mesh_name, false, false);
}

Results optimize_single_new(Eigen::MatrixXd &V, Eigen::MatrixXi &F, EnergyType method, std::string dir, std::string mesh_name, const bool init_with_intrinsic, const bool priority_queue_flips){
  size_t lastindex = mesh_name.find_last_of(".");
  std::string mesh_name_wo_extension =  mesh_name.substr(0, lastindex);

  // Results
  Results res;

  Eigen::MatrixXd CN;
  Eigen::MatrixXi FN;
  std::string to_store_dir = dir + "/" + mesh_name_wo_extension;
  std::filesystem::create_directory(to_store_dir);
  std::string readable_name;
  switch (method) {
    case EnergyType::DIRICHLET:
      readable_name = "dirichlet";
      break;
    case EnergyType::ASAP:
      readable_name = "asap";
      break;
    case EnergyType::ARAP:
      readable_name = "arap";
      break;
    case EnergyType::SYMMETRIC_DIRICHLET:
      readable_name = "symdirichlet";
      break;
    default:
      throw std::invalid_argument("Method unknown");
  }
  std::cout << "RUNNING MESH " << mesh_name << " WITH " << readable_name << "" <<std::endl;
  to_store_dir+= "/"+readable_name;

  std::filesystem::create_directory(to_store_dir);
  //std::string to_store_dir_all = to_store_dir + "/inbetween";
  //std::filesystem::create_directory(to_store_dir_all);

  // Datastructure construction (to be fair this is only necessary for iparam)
  auto start = std::chrono::high_resolution_clock::now();
  DataGeo data_mesh;
  data_mesh.V = V;
  data_mesh.F = F;
  data_mesh.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh.F));
  data_mesh.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh.inputMesh, data_mesh.V));
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));

  Eigen::MatrixXd UV_ext;
  int initial_flips = -1;
  if (init_with_intrinsic) {
    initial_flips = asIDTasPossible(data_mesh);
  }
  res.flips.push_back(initial_flips);
  res.flips_delaunay.push_back(initial_flips);

  Eigen::MatrixXd UV_iparam;
  // just for slim methods
  igl::SLIMData slimdata;

  double curr_energy;
  auto end = std::chrono::high_resolution_clock::now();
  auto duration_init = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  res.init_time = duration_init;

  double tol = 1e-8;
  unsigned max_iterations = 1000;
  unsigned itr = 0;

  double (*energyfunc)(const Eigen::Matrix2d &);
  Eigen::VectorXd tri_wise_E_before, tri_wise_E_after;

  // Extrinsic normal first run
  std::cout << "------------ EXTRINSIC --------------" <<std::endl;
  std::cout << "Initialization:" << std::endl;
  std::cout << "  constructing datastructure took " << duration_init << "ms" << std::endl;
  std::cout << "Interation: " << itr << std::endl;
  unsigned iterations = 1;


  Eigen::VectorXi VT, VTi;
  Eigen::Matrix<int, -1, 3> TT;
  Eigen::VectorXi VV, VVi;
  Eigen::VectorXi B;
  Eigen::SparseMatrix<double> L;

  Eigen::MatrixXi F_touse = init_with_intrinsic ? data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>() : F;

  start = std::chrono::high_resolution_clock::now();
  vt_adjacency(F_touse, V, VT, VTi);
  tt_adjacency(F_touse, VT, VTi, TT);
  vv_adjacency(F_touse, V, TT, VV, VVi);
  bdy_loop(F_touse, TT, VT, VTi, B);
  end = std::chrono::high_resolution_clock::now();
  auto boundary_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  if (method == EnergyType::DIRICHLET || method == EnergyType::ASAP){
    start = std::chrono::high_resolution_clock::now();
    if(init_with_intrinsic){
      data_mesh.intTri->requireCotanLaplacian();
      L = data_mesh.intTri->cotanLaplacian;
      data_mesh.intTri->unrequireCotanLaplacian();
    }
    else{
      igl::cotmatrix(V,F,L);
      L = -L;
    }
    end = std::chrono::high_resolution_clock::now();
  }
  auto L_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  unsigned fixed1, fixed2;
  if (method == EnergyType::ASAP){
    start = std::chrono::high_resolution_clock::now();
    fixed1 = B[0], fixed2 = B[0];
    for (unsigned i = 0; i < B.rows(); ++i) {
      for (unsigned j = 0; j < B.rows(); ++j) {
        if((V.row(B[i])-V.row(B[j])).norm()>(V.row(fixed1)-V.row(fixed2)).norm()){
          fixed1 = B[i];
          fixed2 = B[j];
        }
      }
    }
    end = std::chrono::high_resolution_clock::now();
  }
  auto fixedpoint_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  duration_init = 0;
  switch (method) {
    case EnergyType::DIRICHLET:{
      start = std::chrono::high_resolution_clock::now();
      circle_boundary_proportional(V, B, UV_ext);
      harmonic(L, B, UV_ext, 2);
      end = std::chrono::high_resolution_clock::now();
      curr_energy = compute_total_energy_localjacob(data_mesh, UV_ext, method);
      duration_init += boundary_time + L_time;
      energyfunc = &dirichlet;
      break;
    }
    case EnergyType::ASAP:{
      start = std::chrono::high_resolution_clock::now();
      lscm(F_touse, V, L, fixed1, fixed2, UV_ext);
      end = std::chrono::high_resolution_clock::now();
      curr_energy = compute_total_energy_localjacob(data_mesh, UV_ext, method);
      duration_init += boundary_time + L_time + fixedpoint_time;
      energyfunc = &asap;
      break;
    }
    case EnergyType::ARAP:{
      start = std::chrono::high_resolution_clock::now();
      Eigen::Matrix<double, -1, 2> UV_ext_init;
      circle_boundary_proportional(V, B, UV_ext_init);
      tutte(VV, VVi, B, UV_ext_init, 0);
      iterations = ARAP_tillconverges(data_mesh, UV_ext_init, UV_ext, max_iterations, true, init_with_intrinsic);
      end = std::chrono::high_resolution_clock::now();
      curr_energy = compute_total_energy_localjacob(data_mesh, UV_ext, method);
      energyfunc = &arap;
      duration_init += boundary_time; // needed for tutte
      break;
    }
    case EnergyType::SYMMETRIC_DIRICHLET:{
      start = std::chrono::high_resolution_clock::now();
      Eigen::Matrix<double, -1, 2> UV_ext_init;
      circle_boundary_proportional(V, B, UV_ext_init);
      tutte(VV, VVi, B, UV_ext_init, 0);
      iterations = slim_tillconverges(data_mesh, slimdata, V, F, UV_ext_init, 1000, init_with_intrinsic);
      end = std::chrono::high_resolution_clock::now();
      UV_ext = slimdata.V_o;
      curr_energy = slimdata.energy/2;
      energyfunc = &symmetric_dirichlet;
      duration_init += boundary_time; // needed for tutte
      break;
    }
    default:
      throw std::invalid_argument("Method unknown");
  }
  duration_init += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  // collecting results
  res.opt_durations.push_back(duration_init);
  res.opt_iterations.push_back(iterations);
  res.energies.push_back(curr_energy);

  std::cout << "  parameterization took " << duration_init << "ms in "<< iterations << " iterations" << std::endl;
  std::cout << "    energy: " << std::setprecision(32) << curr_energy << std::endl;
  std::string str = to_store_dir + "/" + mesh_name_wo_extension + "_ext" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_ext, F);

  // store energy per triangle
  tri_wise_energy(data_mesh, UV_ext, energyfunc, init_with_intrinsic, tri_wise_E_before);
  str = to_store_dir + "/" + mesh_name_wo_extension + "_ext" + ".energy";
  saveData(str, tri_wise_E_before);

  //IPARAM
  std::cout << "------------ IPARAM --------------" <<std::endl;
  double past_energy = 0;

  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;

  UV_iparam = UV_ext;
  while(itr<max_iterations && std::abs(past_energy-curr_energy)>tol){
    std::cout << "Interation: " << itr+1 << std::endl;
    past_energy = curr_energy;
    //intrinsic flipping
    unsigned total_flips = 0;
    unsigned total_del_flips = 0;

    // total_flips = queue_flip(data_mesh, UV_iparam, total_del_flips, method);
    unsigned flips, del;
    start = std::chrono::high_resolution_clock::now();
    if (priority_queue_flips) {
      total_flips = priority_queue_flip(data_mesh, UV_iparam, total_del_flips, method);
    } else {
      total_flips = queue_flip(data_mesh, UV_iparam, total_del_flips, method);
    }
    end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "  found " << total_flips << " improving flips, " << total_del_flips <<" were delaunay in "<< duration << "ms" << std::endl;
    res.flips.push_back(total_flips);
    res.flips_delaunay.push_back(total_del_flips);
    res.flip_durations.push_back(duration);
    if (total_flips == 0) {
      break;
    }
    curr_energy = compute_total_energy_localjacob(data_mesh, UV_iparam, method);
    res.energies.push_back(curr_energy);
    std::cout << "    energy: " << std::setprecision(32) << curr_energy << std::endl;
    // Parametrize using new connectivity + Jacobians
    iterations = 1;
    switch (method) {
      case EnergyType::DIRICHLET:{
        start = std::chrono::high_resolution_clock::now();
        data_mesh.intTri->requireCotanLaplacian();
        L = data_mesh.intTri->cotanLaplacian;
        data_mesh.intTri->unrequireCotanLaplacian();
        harmonic(L, B, UV_iparam, 2);
        end = std::chrono::high_resolution_clock::now();
        curr_energy = compute_total_energy_localjacob(data_mesh, UV_iparam, method);
        break;
      }
      case EnergyType::ASAP:{
        // compute with gauss seidel
        start = std::chrono::high_resolution_clock::now();
        data_mesh.intTri->requireCotanLaplacian();
        L = data_mesh.intTri->cotanLaplacian;
        data_mesh.intTri->unrequireCotanLaplacian();
        F_touse = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
        Eigen::MatrixXd UV_new;
        lscm(F_touse, V, L, UV_iparam, fixed1, fixed2, UV_new);
        end = std::chrono::high_resolution_clock::now();
        UV_iparam = UV_new;
        curr_energy = compute_total_energy_localjacob(data_mesh, UV_iparam, method);
        break;
      }
      case EnergyType::ARAP:{
        start = std::chrono::high_resolution_clock::now();
        Eigen::MatrixXd UV_new;
        iterations = ARAP_tillconverges(data_mesh, UV_iparam, UV_new, max_iterations, true, true);
        end = std::chrono::high_resolution_clock::now();
        UV_iparam = UV_new;
        curr_energy = compute_total_energy_localjacob(data_mesh, UV_iparam, method);
        break;
      }
      case EnergyType::SYMMETRIC_DIRICHLET:{
        start = std::chrono::high_resolution_clock::now();
        iterations = slim_tillconverges(data_mesh, slimdata, V, F, UV_iparam, max_iterations, true);
        end = std::chrono::high_resolution_clock::now();
        UV_iparam = slimdata.V_o;
        curr_energy = slimdata.energy/2;
        break;
      }
      default:
        throw std::invalid_argument("Method unknown");
    }
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    res.opt_durations.push_back(duration);
    res.opt_iterations.push_back(iterations);
    res.energies.push_back(curr_energy);

    std::cout << "  parameterization took " << duration << "ms in "<< iterations << " iterations" << std::endl;
    std::cout << "    energy: " << std::setprecision(32) << curr_energy << std::endl;

    ++itr;
    // std::string to_store = to_store_dir + "/inbetween/" + mesh_name +"_"+ std::to_string(itr);
    // to_store += ".obj";
    // igl::writeOBJ(to_store, data_mesh.V, data_mesh.F, CN, FN, UV_iparam, data_mesh.F);
    // try {
    //   store_intrinsic_edges(data_mesh, to_store_dir + "/inbetween/" + mesh_name_wo_extension +"_"+ std::to_string(itr));
    // } catch (...) {
    //   std::cout << "Error in mesh: " << mesh_name_wo_extension << std::endl;
    // }
    // try {
    //   store_intrinsic_mesh(data_mesh, UV_iparam, to_store_dir + "/inbetween/" + mesh_name_wo_extension +"_"+ std::to_string(itr));
    // } catch (...) {
    //   std::cout << "Error in mesh: " << mesh_name_wo_extension << std::endl;
    // }
  }

  str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam, F);

  // store energy per triangle
  tri_wise_energy(data_mesh, UV_iparam, energyfunc, true, tri_wise_E_after);
  str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam" + ".energy";
  saveData(str, tri_wise_E_after);

  // try {
  //   store_intrinsic_edges(data_mesh, to_store_dir + "/" + mesh_name_wo_extension);
  // } catch (...) {
  //   std::cout << "Error in mesh: " << mesh_name_wo_extension << std::endl;
  // }
  try {
    store_intrinsic_mesh(data_mesh, UV_iparam, to_store_dir + "/" + mesh_name_wo_extension, tri_wise_E_before, tri_wise_E_after, res);
  } catch (...) {
    std::cout << "Error in mesh: " << mesh_name_wo_extension << std::endl;
  }
  return res;
}

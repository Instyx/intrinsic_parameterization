#include "optimize.hpp"

#include <chrono>
#include <filesystem>

#include <igl/slim.h>
#include <igl/writeOBJ.h>

#include "intrinsicslim.hpp"
#include "parameterization.hpp"
#include "store_intrinsic_mesh.hpp"

Results optimize_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, EnergyType method, std::string dir, std::string mesh_name){
  return optimize_single(V, F, method, dir, mesh_name, false, false);
}

Results optimize_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, EnergyType method, std::string dir, std::string mesh_name, const bool init_with_intrinsic, const bool priority_queue_flips){
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
  std::string to_store_dir_all = to_store_dir + "/inbetween";
  std::filesystem::create_directory(to_store_dir_all);

  // Datastructure construction (to be fair this is only necessary for iparam)
  auto start = std::chrono::high_resolution_clock::now();
  DataGeo data_mesh;
  data_mesh.V = V;
  data_mesh.F = F;
  data_mesh.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh.F));
  data_mesh.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh.inputMesh, data_mesh.V));
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->requireFaceAreas();
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
  // Extrinsic normal first run
  std::cout << "------------ EXTRINSIC --------------" <<std::endl;
  std::cout << "Initialization:" << std::endl;
  std::cout << "  constructing datastructure took " << duration_init << "ms" << std::endl;
  std::cout << "Interation: " << itr << std::endl;
  unsigned iterations = 1;
  switch (method) {
    case EnergyType::DIRICHLET:{
      start = std::chrono::high_resolution_clock::now();
      UV_ext = harmonic(data_mesh, init_with_intrinsic);
      end = std::chrono::high_resolution_clock::now();
      curr_energy = compute_total_energy(data_mesh, UV_ext, method, init_with_intrinsic);
      break;
    }
    case EnergyType::ASAP:{
      start = std::chrono::high_resolution_clock::now();
      UV_ext = LSCM(data_mesh, true, init_with_intrinsic);
      end = std::chrono::high_resolution_clock::now();
      curr_energy = compute_total_energy(data_mesh, UV_ext, method, init_with_intrinsic);
      break;
    }
    case EnergyType::ARAP:{
      start = std::chrono::high_resolution_clock::now();
      Eigen::MatrixXd UV_ext_init = tutte(data_mesh, init_with_intrinsic);
      iterations = ARAP_tillconverges(data_mesh, UV_ext_init, UV_ext, max_iterations, true, init_with_intrinsic);
      end = std::chrono::high_resolution_clock::now();
      curr_energy = compute_total_energy(data_mesh, UV_ext, method, init_with_intrinsic);
      break;
    }
    case EnergyType::SYMMETRIC_DIRICHLET:{
      start = std::chrono::high_resolution_clock::now();
      Eigen::MatrixXd UV_ext_init = tutte(data_mesh, init_with_intrinsic);
      iterations = slim_tillconverges(data_mesh, slimdata, V, F, UV_ext_init, 1000, init_with_intrinsic);
      end = std::chrono::high_resolution_clock::now();
      UV_ext = slimdata.V_o;
      curr_energy = slimdata.energy/2;
      break;
    }
    default:
      throw std::invalid_argument("Method unknown");
  }
  duration_init = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  // collecting results
  res.opt_iterations.push_back(iterations);
  res.energies.push_back(curr_energy);

  std::cout << "  parameterization took " << duration_init << "ms in "<< iterations << " iterations" << std::endl;
  std::cout << "    energy: " << std::setprecision(32) << curr_energy << std::endl;
  std::string str = to_store_dir + "/" + mesh_name_wo_extension + "_ext" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_ext, F);

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

    auto start = std::chrono::high_resolution_clock::now();
    // total_flips = queue_flip(data_mesh, UV_iparam, total_del_flips, method);
    unsigned flips, del;
    start = std::chrono::high_resolution_clock::now();
    if (priority_queue_flips) {
      total_flips = priority_queue_flip(data_mesh, UV_iparam, total_del_flips, method);
    } else {
      total_flips = queue_flip(data_mesh, UV_iparam, total_del_flips, method);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "  found " << total_flips << " improving flips, " << total_del_flips <<" were delaunay in "<< duration << "ms" << std::endl;
    res.flips.push_back(total_flips);
    res.flips_delaunay.push_back(total_del_flips);
    res.flip_durations.push_back(duration);
    if (total_flips == 0) {
      break;
    }
    curr_energy = compute_total_energy(data_mesh, UV_iparam, method, true);
    res.energies.push_back(curr_energy);
    std::cout << "    energy: " << std::setprecision(32) << curr_energy << std::endl;
    // Parametrize using new connectivity + Jacobians
    iterations = 1;
    switch (method) {
      case EnergyType::DIRICHLET:{
        start = std::chrono::high_resolution_clock::now();
        UV_iparam = harmonic(data_mesh, true);
        end = std::chrono::high_resolution_clock::now();
        curr_energy = compute_total_energy(data_mesh, UV_iparam, method, true);
        break;
      }
      case EnergyType::ASAP:{
        start = std::chrono::high_resolution_clock::now();
        UV_iparam = LSCM(data_mesh, true, true);
        end = std::chrono::high_resolution_clock::now();
        curr_energy = compute_total_energy(data_mesh, UV_iparam, method, true);
        break;
      }
      case EnergyType::ARAP:{
        start = std::chrono::high_resolution_clock::now();
        Eigen::MatrixXd UV_new;
        iterations = ARAP_tillconverges(data_mesh, UV_iparam, UV_new, max_iterations, true, true);
        end = std::chrono::high_resolution_clock::now();
        UV_iparam = UV_new;
        curr_energy = compute_total_energy(data_mesh, UV_iparam, method, true);
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
    std::string to_store = to_store_dir + "/inbetween/" + mesh_name +"_"+ std::to_string(itr);
    to_store += ".obj";
    igl::writeOBJ(to_store, data_mesh.V, data_mesh.F, CN, FN, UV_iparam, data_mesh.F);
    try {
      store_intrinsic_edges(data_mesh, to_store_dir + "/inbetween/" + mesh_name_wo_extension +"_"+ std::to_string(itr));
    } catch (...) {
      std::cout << "Error in mesh: " << mesh_name_wo_extension << std::endl;
    }
    try {
      store_intrinsic_mesh(data_mesh, UV_iparam, to_store_dir + "/inbetween/" + mesh_name_wo_extension +"_"+ std::to_string(itr));
    } catch (...) {
      std::cout << "Error in mesh: " << mesh_name_wo_extension << std::endl;
    }
  }

  str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam, F);

  try {
    store_intrinsic_edges(data_mesh, to_store_dir + "/" + mesh_name_wo_extension);
  } catch (...) {
    std::cout << "Error in mesh: " << mesh_name_wo_extension << std::endl;
  }
  try {
    store_intrinsic_mesh(data_mesh, UV_iparam, to_store_dir + "/" + mesh_name_wo_extension);
  } catch (...) {
    std::cout << "Error in mesh: " << mesh_name_wo_extension << std::endl;
  }
  return res;
}

#include "datageo.hpp"
#include <intrinsicflip.hpp>
#include <parameterization.hpp>
#include <test.hpp>
#include <iostream>
#include <chrono>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>
#include <dirent.h>
#include <intrinsicslim.hpp>
#include <filesystem>
// to check if deluanay flips always decrease the energy
bool isDelaunayFlipBad(DataGeo &data_mesh, const Eigen::MatrixXd &UV){
  for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
    if(data_mesh.intTri->isDelaunay(e)) continue;
    if(e.isBoundary()) continue;
    double diff = flippeddiff(data_mesh, UV, e, EnergyType::ASAP);
    if(diff==0) std::cout << "Deluanay flip not possible" << std::endl;
    if(diff>0) return true;
  }
  return false;
}

DataGeo compareIDTvsGreedy(DataGeo &data_mesh){
  DataGeo data_mesh_idt;
  data_mesh_idt.V = data_mesh.V;
  data_mesh_idt.F = data_mesh.F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->flipToDelaunay();

  Eigen::MatrixXd UV, new_UV, UV_idt, new_UV_idt;

  computeParameterization(data_mesh_idt, data_mesh_idt.V, data_mesh_idt.F, UV_idt, new_UV_idt,
      false, true, '2');
  UV_idt = new_UV_idt;

  computeParameterization(data_mesh, data_mesh.V, data_mesh.F, UV, new_UV, false, false, '2');
  UV = new_UV;
  unsigned del;
  while(greedy_flip(data_mesh, UV, del, EnergyType::DIRICHLET));
  computeParameterization(data_mesh, data_mesh.V, data_mesh.F, UV, new_UV, false, true, '2');
  UV = new_UV;
  std::cout << "Intrinsic Energy: " << std::endl;
  std::cout << "  IDT energy: " << compute_total_energy(data_mesh_idt, UV_idt, EnergyType::DIRICHLET, true) <<
    " ;  Greedy energy: " << compute_total_energy(data_mesh, UV, EnergyType::DIRICHLET, true) << std::endl;

  std::cout << "Extrinsic Energy: " << std::endl;
  std::cout << "  IDT energy: " << compute_total_energy(data_mesh_idt, UV_idt, EnergyType::DIRICHLET, false) <<
    " ;  Greedy energy: " << compute_total_energy(data_mesh, UV, EnergyType::DIRICHLET, false) << std::endl;
  return data_mesh_idt;
}


unsigned flipEdgesifCoplanar(DataGeo &data_mesh, bool onlyDelaunay){
  Eigen::MatrixXd V = data_mesh.V;
  double tolarence = 1e-8;
  unsigned total_flips=0;
  for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
    if(e.isBoundary()) continue;
    if(onlyDelaunay && data_mesh.intTri->isDelaunay(e)) continue;

    double edge_len_b = data_mesh.intTri->edgeLengths[e];
    data_mesh.intTri->flipEdgeIfPossible(e);
    double edge_len = data_mesh.intTri->edgeLengths[e];

   // std::cout << edge_len << "   " << edge_len_b << std::endl;

    size_t v1 = data_mesh.intTri->vertexIndices[e.firstVertex()];
    size_t v2 = data_mesh.intTri->vertexIndices[e.secondVertex()];
    Eigen::RowVector3d vec = V.row(v1) - V.row(v2);
    double vec_norm = vec.norm();

    // flip back if not coplanar
    if(std::abs(vec_norm-edge_len)>tolarence) data_mesh.intTri->flipEdgeIfPossible(e);
    else total_flips++;
  }
  return total_flips;
}

/*
void compareIDTvsIPARAM(DataGeo &data_mesh, bool isFreeBoundary, const EnergyType &et, Eigen::MatrixXd &UV_idt, Eigen::MatrixXd &UV_iparam){
  DataGeo data_mesh_idt;
  data_mesh_idt.V = data_mesh.V;
  data_mesh_idt.F = data_mesh.F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->flipToDelaunay();

  Eigen::MatrixXd UV_idt_init;

  reset_constraints();
  computeParameterization(data_mesh_idt, data_mesh_idt.V, data_mesh_idt.F, UV_idt, UV_idt_init, false, true, '2');
  UV_idt.resize(data_mesh.V.rows(),2);
  UV_iparam.resize(data_mesh.V.rows(),2);
  if(et == EnergyType::ARAP){
    std::cout << "------ IDT --------" << std::endl;
    UV_idt = ARAP_tillconverges(data_mesh_idt, UV_idt_init, 500, isFreeBoundary, true);
    std::cout << "------ IPARAM --------" << std::endl;
    UV_iparam = intrinsic_ARAP(data_mesh, 150, 20, isFreeBoundary);
    std::cout << "last IDT energy int: " << compute_total_energy(data_mesh_idt, UV_idt, EnergyType::ARAP, true) << std::endl;
    std::cout << "last IDT energy ext: " << compute_total_energy(data_mesh_idt, UV_idt, EnergyType::ARAP, false) << std::endl;
    std::cout << "last IPARAM energy int: " << compute_total_energy(data_mesh, UV_iparam, EnergyType::ARAP, true) << std::endl;
    std::cout << "last IPARAM energy ext: " << compute_total_energy(data_mesh, UV_iparam, EnergyType::ARAP, false) << std::endl;

  }
  // symmetric Dirichlet
  else{
    UV_idt
  }
}
*/
void test_ARAP_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F,  std::string mesh_name, bool isFreeBoundary, std::ostream &fout){

  //IPARAM
  std::cout << "------------ IPARAM -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam" << ",";

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

  Eigen::MatrixXd UV_iparam;
  unsigned total_iterations = intrinsic_ARAP(data_mesh, UV_iparam, 1000, 50, isFreeBoundary, fout);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  fout << total_iterations << "," <<  duration << ","; // this is total duration
  for(int i=total_iterations*7;i<350;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  //extrinsic
  std::cout << "------------ EXTRINSIC -------------- " <<std::endl;
  fout << mesh_name << "," << "ext" << ",";
  start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXd UV_ext;
  Eigen::MatrixXd UV_ext_init = tutte(data_mesh, false);
  total_iterations = ARAP_tillconverges(data_mesh, UV_ext_init, UV_ext, 1000, isFreeBoundary, false);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << compute_total_energy(data_mesh, UV_ext_init, EnergyType::ARAP , false) << ",";
  fout << total_iterations <<  "," << duration << "," << compute_total_energy(data_mesh, UV_ext, EnergyType::ARAP , false) << ",";
  for(int i = 0; i<352;++i){
    fout << -1 << ",";
  }
  fout << "\n";


  // IDT
  std::cout << "------------ IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "idt" << ",";
  start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh_idt;
  data_mesh_idt.V = V;
  data_mesh_idt.F = F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->requireEdgeLengths();
  data_mesh_idt.intTri->requireVertexIndices();
  data_mesh_idt.intTri->requireFaceAreas();

  data_mesh_idt.intTri->flipToDelaunay();

  Eigen::MatrixXd UV_int;
  Eigen::MatrixXd UV_int_init = tutte(data_mesh_idt, true);

  total_iterations = ARAP_tillconverges(data_mesh_idt, UV_int_init, UV_int, 1000, isFreeBoundary, true);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << compute_total_energy(data_mesh, UV_int_init, EnergyType::ARAP , true) << ",";
  fout << total_iterations << "," << duration << "," << compute_total_energy(data_mesh_idt, UV_int, EnergyType::ARAP , true) << ",";
  for(int i = 0; i<352;++i){
    fout << -1 << ",";
  }
  fout << "\n";

}

void test_ARAP_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string dir, std::string mesh_name, bool isFreeBoundary, std::ostream &fout){
  size_t lastindex = mesh_name.find_last_of(".");
  std::string mesh_name_wo_extension =  mesh_name.substr(0, lastindex);

  Eigen::MatrixXd CN;
  Eigen::MatrixXi FN;
  std::string to_store_dir = dir + "/" + mesh_name_wo_extension;
  std::filesystem::create_directory(to_store_dir);
  to_store_dir += "/arap";
  std::filesystem::create_directory(to_store_dir);
  std::string to_store_dir_all = to_store_dir + "/inbetween";
  std::filesystem::create_directory(to_store_dir_all);

  //IPARAM
  std::cout << "------------ IPARAM -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam" << ",";

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
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  Eigen::MatrixXd UV_iparam;
  unsigned total_iterations = intrinsic_ARAP(data_mesh, UV_iparam, 1000, 50, isFreeBoundary, fout, to_store_dir_all, mesh_name_wo_extension);
  fout << total_iterations << "," <<  duration << ","; // this is init duration
  for(int i=total_iterations*7;i<350;++i){
    fout << -1 << ",";
  }
  fout << "\n";
  std::string str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam, F);

  //extrinsic
  std::cout << "------------ EXTRINSIC -------------- " <<std::endl;
  fout << mesh_name << "," << "ext" << ",";
  start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXd UV_ext;
  Eigen::MatrixXd UV_ext_init = tutte(data_mesh, false);
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  fout << duration << "," << compute_total_energy(data_mesh, UV_ext_init, EnergyType::ARAP , false) << ",";
  
  start = std::chrono::high_resolution_clock::now();
  total_iterations = ARAP_tillconverges(data_mesh, UV_ext_init, UV_ext, 1000, isFreeBoundary, false);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << total_iterations <<  "," << duration << "," << compute_total_energy(data_mesh, UV_ext, EnergyType::ARAP , false) << ",";
  for(int i = 0; i<352;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_ext" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_ext, F);

  // IDT
  std::cout << "------------ IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "idt" << ",";
  start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh_idt;
  data_mesh_idt.V = V;
  data_mesh_idt.F = F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->requireEdgeLengths();
  data_mesh_idt.intTri->requireVertexIndices();
  data_mesh_idt.intTri->requireFaceAreas();

  data_mesh_idt.intTri->flipToDelaunay();
  end = std::chrono::high_resolution_clock::now();
  auto duration_idt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  
  start = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd UV_int;
  Eigen::MatrixXd UV_int_init = tutte(data_mesh_idt, true);
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << duration_idt << "," << duration << "," << compute_total_energy(data_mesh, UV_int_init, EnergyType::ARAP , true) << ",";

  start = std::chrono::high_resolution_clock::now();
  total_iterations = ARAP_tillconverges(data_mesh_idt, UV_int_init, UV_int, 1000, isFreeBoundary, true);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout  << total_iterations << "," << duration << "," << compute_total_energy(data_mesh_idt, UV_int, EnergyType::ARAP , true) << ",";
  for(int i = 0; i<351;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_idt" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_int, F);


  //IDT_IPARAM
  std::cout << "------------ IPARAM_IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam_idt" << ",";
  to_store_dir_all = to_store_dir + "/inbetween_idt";
  std::filesystem::create_directory(to_store_dir_all);

  Eigen::MatrixXd UV_iparam_idt;
  total_iterations = intrinsic_ARAP(data_mesh_idt, UV_iparam_idt, 1000, 50, isFreeBoundary, fout, to_store_dir_all, mesh_name_wo_extension);
  fout << total_iterations << "," << duration_idt << ","; // this is init duration
  for(int i=total_iterations*7;i<350;++i){
    fout << -1 << ",";
  }
  fout << "\n";
  str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam_idt" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam_idt, F);

}


void test_Dirichlet_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F,  std::string mesh_name,  std::ostream &fout){

  std::cout << "------------ IPARAM -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam" << ",";

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

  Eigen::MatrixXd UV_iparam;

  unsigned total_iterations = intrinsic_harmonic(data_mesh, UV_iparam, 50, fout);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  fout << total_iterations << "," << duration << ","; // this is total duration
  for(int i=total_iterations*7;i<300;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  //extrinsic
  std::cout << "------------ EXTRINSIC -------------- " <<std::endl;
  fout << mesh_name << "," << "ext" << ",";
  start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXd UV_ext = harmonic(data_mesh, false);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  fout << duration << "," << compute_total_energy(data_mesh, UV_ext, EnergyType::DIRICHLET , false) << ",";
  for(int i = 0; i<302; ++i){
    fout << -1 << ",";
  }
  fout << "\n";

  // IDT
  std::cout << "------------ IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "idt" << ",";
  start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh_idt;
  data_mesh_idt.V = V;
  data_mesh_idt.F = F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->requireEdgeLengths();
  data_mesh_idt.intTri->requireVertexIndices();
  data_mesh_idt.intTri->requireFaceAreas();

  data_mesh_idt.intTri->flipToDelaunay();

  Eigen::MatrixXd UV_int = harmonic(data_mesh_idt, true);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << duration << "," << compute_total_energy(data_mesh_idt, UV_int, EnergyType::DIRICHLET, true) << ",";
  for(int i = 0; i<302;++i){
    fout << -1 << ",";
  }
  fout << "\n";

}

void test_Dirichlet_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string dir, std::string mesh_name,  std::ostream &fout){

  size_t lastindex = mesh_name.find_last_of(".");
  std::string mesh_name_wo_extension =  mesh_name.substr(0, lastindex);

  Eigen::MatrixXd CN;
  Eigen::MatrixXi FN;
  std::string to_store_dir = dir + "/" + mesh_name_wo_extension;
  std::filesystem::create_directory(to_store_dir);
  to_store_dir += "/dirichlet";
  std::filesystem::create_directory(to_store_dir);
  std::string to_store_dir_all = to_store_dir + "/inbetween";
  std::filesystem::create_directory(to_store_dir_all);

  std::cout << "------------ IPARAM -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam" << ",";

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
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  Eigen::MatrixXd UV_iparam;

  unsigned total_iterations = intrinsic_harmonic(data_mesh, UV_iparam, 50, fout, to_store_dir_all, mesh_name_wo_extension);
  fout << total_iterations << "," << duration << ","; // this is init duration
  for(int i=total_iterations*7;i<300;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  std::string str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam, F);

  //extrinsic
  std::cout << "------------ EXTRINSIC -------------- " <<std::endl;
  fout << mesh_name << "," << "ext" << ",";
  start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXd UV_ext = harmonic(data_mesh, false);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  fout << duration << "," << compute_total_energy(data_mesh, UV_ext, EnergyType::DIRICHLET , false) << ",";
  for(int i = 0; i<302; ++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_ext" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_ext, F);

  // IDT
  std::cout << "------------ IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "idt" << ",";
  start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh_idt;
  data_mesh_idt.V = V;
  data_mesh_idt.F = F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->requireEdgeLengths();
  data_mesh_idt.intTri->requireVertexIndices();
  data_mesh_idt.intTri->requireFaceAreas();

  data_mesh_idt.intTri->flipToDelaunay();
  end = std::chrono::high_resolution_clock::now();
  auto duration_idt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  start = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd UV_int = harmonic(data_mesh_idt, true);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << duration_idt << "," << duration << "," << compute_total_energy(data_mesh_idt, UV_int, EnergyType::DIRICHLET, true) << ",";
  for(int i = 0; i<301;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_idt" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_int, F);

  //IDT_IPARAM
  std::cout << "------------ IPARAM_IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam_idt" << ",";
  to_store_dir_all = to_store_dir + "/inbetween_idt";
  std::filesystem::create_directory(to_store_dir_all);

  Eigen::MatrixXd UV_iparam_idt;
  
  total_iterations = intrinsic_harmonic(data_mesh_idt, UV_iparam_idt, 50, fout, to_store_dir_all, mesh_name_wo_extension);
  fout << total_iterations << "," <<  duration_idt << ","; 
  for(int i=total_iterations*7;i<300;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam_idt" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam_idt, F);

}



void test_ASAP_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string mesh_name, bool isFreeBoundary, std::ostream &fout){


  std::cout << "------------ IPARAM -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam" << ",";

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

  Eigen::MatrixXd UV_iparam;

  unsigned total_iterations = intrinsic_LSCM(data_mesh, UV_iparam, 50, true, fout);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  fout << total_iterations << "," << duration << ","; // this is total duration
  for(int i=total_iterations*7;i<300;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  //extrinsic
  std::cout << "------------ EXTRINSIC -------------- " <<std::endl;
  fout << mesh_name << "," << "ext" << ",";
  start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXd UV_ext = LSCM(data_mesh, isFreeBoundary, false);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  fout << duration << "," << compute_total_energy(data_mesh, UV_ext, EnergyType::ASAP, false) << ",";
  for(int i = 0; i<302; ++i){
    fout << -1 << ",";
  }
  fout << "\n";

  // IDT
  std::cout << "------------ IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "idt" << ",";
  start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh_idt;
  data_mesh_idt.V = V;
  data_mesh_idt.F = F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->requireEdgeLengths();
  data_mesh_idt.intTri->requireVertexIndices();
  data_mesh_idt.intTri->requireFaceAreas();

  data_mesh_idt.intTri->flipToDelaunay();

  Eigen::MatrixXd UV_int = LSCM(data_mesh_idt, isFreeBoundary, true);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << duration << "," << compute_total_energy(data_mesh_idt, UV_int, EnergyType::ASAP, true) << ",";
  for(int i = 0; i<302;++i){
    fout << -1 << ",";
  }
  fout << "\n";

}

void test_ASAP_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string dir, std::string mesh_name, bool isFreeBoundary, std::ostream &fout){

  size_t lastindex = mesh_name.find_last_of(".");
  std::string mesh_name_wo_extension =  mesh_name.substr(0, lastindex);

  Eigen::MatrixXd CN;
  Eigen::MatrixXi FN;
  std::string to_store_dir = dir + "/" + mesh_name_wo_extension;
  std::filesystem::create_directory(to_store_dir);
  to_store_dir += "/asap";
  std::filesystem::create_directory(to_store_dir);
  std::string to_store_dir_all = to_store_dir + "/inbetween";
  std::filesystem::create_directory(to_store_dir_all);


  std::cout << "------------ IPARAM -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam" << ",";

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
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  Eigen::MatrixXd UV_iparam;

  unsigned total_iterations = intrinsic_LSCM(data_mesh, UV_iparam, 50, true, fout, to_store_dir_all, mesh_name_wo_extension);
  fout << total_iterations << "," << duration << ","; // this is init duration
  for(int i=total_iterations*7;i<300;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  std::string str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam, F);

  //extrinsic
  std::cout << "------------ EXTRINSIC -------------- " <<std::endl;
  fout << mesh_name << "," << "ext" << ",";
  start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXd UV_ext = LSCM(data_mesh, isFreeBoundary, false);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  fout << duration << "," << compute_total_energy(data_mesh, UV_ext, EnergyType::ASAP, false) << ",";
  for(int i = 0; i<302; ++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_ext" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_ext, F);


  // IDT
  std::cout << "------------ IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "idt" << ",";
  start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh_idt;
  data_mesh_idt.V = V;
  data_mesh_idt.F = F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->requireEdgeLengths();
  data_mesh_idt.intTri->requireVertexIndices();
  data_mesh_idt.intTri->requireFaceAreas();

  data_mesh_idt.intTri->flipToDelaunay();
  end = std::chrono::high_resolution_clock::now();
  auto duration_idt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  start = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd UV_int = LSCM(data_mesh_idt, isFreeBoundary, true);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << duration_idt << "," <<duration << "," << compute_total_energy(data_mesh_idt, UV_int, EnergyType::ASAP, true) << ",";
  for(int i = 0; i<301;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_int, F);

  //IDT_IPARAM
  std::cout << "------------ IPARAM_IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam_idt" << ",";
  to_store_dir_all = to_store_dir + "/inbetween_idt";
  std::filesystem::create_directory(to_store_dir_all);

  Eigen::MatrixXd UV_iparam_idt;

  total_iterations = intrinsic_LSCM(data_mesh_idt, UV_iparam_idt, 50, true, fout, to_store_dir_all, mesh_name_wo_extension);
  fout << total_iterations << "," << duration_idt << ","; // this is init duration
  for(int i=total_iterations*7;i<300;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam_idt" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam_idt, F);

}


void test_SymDirichlet_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string mesh_name, std::ostream &fout){

  //IPARAM
  std::cout << "------------ IPARAM -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam" << ",";

  auto start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh;
  igl::SLIMData slimdata;
  data_mesh.V = V;
  data_mesh.F = F;
  data_mesh.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh.F));
  data_mesh.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh.inputMesh, data_mesh.V));
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->requireFaceAreas();

  Eigen::MatrixXd UV_iparam;
  Eigen::MatrixXd UV_iparam_init = tutte(data_mesh, false);

  unsigned total_iterations = intrinsicslim(data_mesh, UV_iparam_init, UV_iparam, 1000, 50, fout);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  fout << total_iterations << "," <<  duration << ","; // this is total duration
  for(int i=total_iterations*7;i<350;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  //extrinsic
  std::cout << "------------ EXTRINSIC -------------- " <<std::endl;
  fout << mesh_name << "," << "ext" << ",";
  start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXd UV_ext;
  Eigen::MatrixXd UV_ext_init = tutte(data_mesh, false);
  igl::SLIMData slimdata_ext;
  total_iterations = slim_tillconverges(data_mesh, slimdata_ext, V, F, UV_ext_init, 1000, false);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << compute_total_energy(data_mesh, UV_ext_init, EnergyType::SYMMETRIC_DIRICHLET , false) << ",";
  fout << total_iterations <<  "," << duration << "," << slimdata_ext.energy/2 << ",";
  for(int i = 0; i<352;++i){
    fout << -1 << ",";
  }
  fout << "\n";


  // IDT
  std::cout << "------------ IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "idt" << ",";
  start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh_idt;
  igl::SLIMData slimdata_idt;
  data_mesh_idt.V = V;
  data_mesh_idt.F = F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->requireEdgeLengths();
  data_mesh_idt.intTri->requireVertexIndices();
  data_mesh_idt.intTri->requireFaceAreas();

  data_mesh_idt.intTri->flipToDelaunay();

  Eigen::MatrixXd UV_int;
  Eigen::MatrixXd UV_int_init = tutte(data_mesh_idt, true);

  total_iterations = slim_tillconverges(data_mesh_idt, slimdata_idt, V, F, UV_int_init, 1000, true);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << compute_total_energy(data_mesh, UV_int_init, EnergyType::SYMMETRIC_DIRICHLET , true) << ",";
  fout << total_iterations << "," << duration << "," << slimdata_idt.energy/2 << ",";
  for(int i = 0; i<352;++i){
    fout << -1 << ",";
  }
  fout << "\n";

}

void test_SymDirichlet_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string dir, std::string mesh_name, std::ostream &fout){

  size_t lastindex = mesh_name.find_last_of(".");
  std::string mesh_name_wo_extension =  mesh_name.substr(0, lastindex);

  Eigen::MatrixXd CN;
  Eigen::MatrixXi FN;
  std::string to_store_dir = dir + "/" + mesh_name_wo_extension;
  std::filesystem::create_directory(to_store_dir);
  to_store_dir += "/symdirichlet";
  std::filesystem::create_directory(to_store_dir);
  std::string to_store_dir_all = to_store_dir + "/inbetween";
  std::filesystem::create_directory(to_store_dir_all);

  //IPARAM
  std::cout << "------------ IPARAM -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam" << ",";

  auto start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh;
  igl::SLIMData slimdata;
  data_mesh.V = V;
  data_mesh.F = F;
  data_mesh.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh.F));
  data_mesh.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh.inputMesh, data_mesh.V));
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->requireFaceAreas();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration_init = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


  start = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd UV_iparam;
  Eigen::MatrixXd UV_iparam_init = tutte(data_mesh, false);
  end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << duration << ",";
  unsigned total_iterations = intrinsicslim(data_mesh, UV_iparam_init, UV_iparam, 1000, 50, fout, to_store_dir_all, mesh_name_wo_extension);
  fout << total_iterations << "," << duration_init << ","; // this is init duration
  for(int i=total_iterations*7;i<350;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  std::string str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam, F);

  //extrinsic
  std::cout << "------------ EXTRINSIC -------------- " <<std::endl;
  fout << mesh_name << "," << "ext" << ",";
  start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXd UV_ext;
  Eigen::MatrixXd UV_ext_init = tutte(data_mesh, false);
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << duration << "," << compute_total_energy(data_mesh, UV_ext_init, EnergyType::SYMMETRIC_DIRICHLET , false) << ",";
  start = std::chrono::high_resolution_clock::now();
  igl::SLIMData slimdata_ext;
  total_iterations = slim_tillconverges(data_mesh, slimdata_ext, V, F, UV_ext_init, 1000, false);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << total_iterations <<  "," << duration << "," << slimdata_ext.energy/2 << ",";
  for(int i = 0; i<352;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_ext" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_ext, F);


  // IDT
  std::cout << "------------ IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "idt" << ",";
  start = std::chrono::high_resolution_clock::now();

  DataGeo data_mesh_idt;
  igl::SLIMData slimdata_idt;
  data_mesh_idt.V = V;
  data_mesh_idt.F = F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh_idt.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh_idt.V));
  data_mesh_idt.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh_idt.inputMesh, *data_mesh_idt.inputGeometry));
  data_mesh_idt.intTri->requireEdgeLengths();
  data_mesh_idt.intTri->requireVertexIndices();
  data_mesh_idt.intTri->requireFaceAreas();

  data_mesh_idt.intTri->flipToDelaunay();
  end = std::chrono::high_resolution_clock::now();
  auto duration_idt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  start = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd UV_int;
  Eigen::MatrixXd UV_int_init = tutte(data_mesh_idt, true);
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << duration_idt << "," << duration << "," << compute_total_energy(data_mesh, UV_int_init, EnergyType::SYMMETRIC_DIRICHLET , true) << ",";

  start = std::chrono::high_resolution_clock::now();
  total_iterations = slim_tillconverges(data_mesh_idt, slimdata_idt, V, F, UV_int_init, 1000, true);

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  fout << total_iterations << "," << duration << "," << slimdata_idt.energy/2 << ",";
  for(int i = 0; i<351;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_int" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_int, F);

  //IDT_IPARAM
  std::cout << "------------ IPARAM_IDT -------------- " <<std::endl;
  fout << mesh_name << "," << "iparam_idt" << ",";
  to_store_dir_all = to_store_dir + "/inbetween_idt";
  std::filesystem::create_directory(to_store_dir_all);

  Eigen::MatrixXd UV_iparam_idt;
  total_iterations = intrinsicslim(data_mesh_idt, UV_int_init, UV_iparam_idt, 1000, 50, fout, to_store_dir_all, mesh_name_wo_extension);
  fout << total_iterations << "," <<  duration_idt << ","; // this is init duration
  for(int i=total_iterations*7;i<350;++i){
    fout << -1 << ",";
  }
  fout << "\n";

  str = to_store_dir + "/" + mesh_name_wo_extension + "_iparam_idt" + ".obj";
  igl::writeOBJ(str, V, F, CN, FN, UV_iparam_idt, F);


}


void test_all(){
  const char* folderPath = "../data";
  DIR* directory = opendir(folderPath);
  if (directory == NULL) {
    std::cerr << "Failed to open directory." << std::endl;
    return;
  }
  std::fstream fout_dirichlet, fout_asap, fout_arap, fout_symdirichlet;
  // opens an existing csv file or creates a new file.
  fout_arap.open("results_arap.csv", std::ios::out | std::ios::app);
  fout_asap.open("results_asap.csv", std::ios::out | std::ios::app);
  fout_dirichlet.open("results_dirichlet.csv", std::ios::out | std::ios::app);
  fout_symdirichlet.open("results_symdirichlet.csv", std::ios::out | std::ios::app);

  // Read directory entries
  struct dirent* entry;
  while ((entry = readdir(directory)) != NULL) {
    // Skip "." and ".." entries
    if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
        continue;
    }
    std::string filePath = std::string(folderPath) + "/" + std::string(entry->d_name);
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    DataGeo data_mesh, data_mesh_idt;
    igl::read_triangle_mesh(filePath,V,F);
    test_ARAP_single(V, F, std::string(entry->d_name), true, fout_arap);
    test_Dirichlet_single(V, F, std::string(entry->d_name), fout_dirichlet);
    test_ASAP_single(V, F, std::string(entry->d_name), true, fout_asap);
    test_SymDirichlet_single(V, F, std::string(entry->d_name), fout_symdirichlet);
  }
  fout_arap.close();
  fout_asap.close();
  fout_dirichlet.close();
  fout_symdirichlet.close();
  // Close the directory
  closedir(directory);
}

void test_all_withtextures(){
  const char* folderPath = "../test_data";
  DIR* directory = opendir(folderPath);
  if (directory == NULL) {
    std::cerr << "Failed to open directory." << std::endl;
    return;
  }
  std::fstream fout_dirichlet, fout_asap, fout_arap, fout_symdirichlet;
  // opens an existing csv file or creates a new file.
  fout_arap.open("results_arap.csv", std::ios::out | std::ios::app);
  fout_asap.open("results_asap.csv", std::ios::out | std::ios::app);
  fout_dirichlet.open("results_dirichlet.csv", std::ios::out | std::ios::app);
  fout_symdirichlet.open("results_symdirichlet.csv", std::ios::out | std::ios::app);

  // Read directory entries
  struct dirent* entry;
  while ((entry = readdir(directory)) != NULL) {
    // Skip "." and ".." entries
    if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
        continue;
    }
    std::string filePath = std::string(folderPath) + "/" + std::string(entry->d_name);
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    DataGeo data_mesh, data_mesh_idt;
    igl::read_triangle_mesh(filePath,V,F);
    test_ARAP_single(V, F, folderPath, std::string(entry->d_name), true, fout_arap);
    test_Dirichlet_single(V, F, folderPath, std::string(entry->d_name),  fout_dirichlet);
    test_ASAP_single(V, F, folderPath, std::string(entry->d_name), true, fout_asap);
    test_SymDirichlet_single(V, F, folderPath, std::string(entry->d_name),  fout_symdirichlet);
  }
  fout_arap.close();
  fout_asap.close();
  fout_dirichlet.close();
  fout_symdirichlet.close();

  // Close the directory
  closedir(directory);
}

#include "datageo.hpp"
#include <intrinsicflip.hpp>
#include <parameterization.hpp>
#include <test.hpp>
#include <iostream>

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
  while(greedy_flip(data_mesh, UV, EnergyType::DIRICHLET));
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
  /*
  else{
    UV_idt
  }
  */
}


void test_ARAP(DataGeo &data_mesh, DataGeo &data_mesh_idt, bool isFreeBoundary){
  
  Eigen::MatrixXd UV_idt_init = tutte(data_mesh_idt, true);
  Eigen::MatrixXd UV_idt = ARAP_tillconverges(data_mesh_idt, UV_idt_init, 1000, isFreeBoundary, true);
  
}

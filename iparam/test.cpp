#include <intrinsicflip.hpp>
#include <parameterization.hpp>
#include <test.hpp>

// FOR DIRICHLET

bool isDelaunayFlipBad(DataGeo &data_mesh, const Eigen::MatrixXd &UV){
  for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
    if(data_mesh.intTri->isDelaunay(e)) continue;
    if(flippeddiff(data_mesh, UV, e, EnergyType::DIRICHLET)>0) return true;
  }
  return false;
}

DataGeo compareIDTvsGreedy(DataGeo &data_mesh){
  DataGeo data_mesh_idt;
  data_mesh_idt.V = data_mesh.V;
  data_mesh_idt.F = data_mesh.F;
  data_mesh_idt.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh.F));
  data_mesh_idt.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh_idt.inputMesh, data_mesh.V));
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




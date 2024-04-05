#include "parameterization.hpp"

#include <string>
#include <filesystem>
#include "read_mesh.hpp"

bool isFlipped(DataGeo &data_mesh, Eigen::MatrixXd UV){
  for(gcs::Face f : data_mesh.intTri->intrinsicMesh->faces()){
    Eigen::Matrix2d J;
    faceJacobian(data_mesh, UV, f, J);
    if(J.determinant()<0) {
      return true;
    }
  }
  return false;
}

bool toFilter(Eigen::MatrixXd &V, Eigen::MatrixXi &F){
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

  std::cout << "------ HARMONIC --------" <<  std::endl;
  UV_ext = harmonic(data_mesh, false);
  if(isFlipped(data_mesh, UV_ext)){
    return true;
  }

  std::cout << "------ CONFORMAL --------" << std::endl;
  UV_ext = LSCM(data_mesh, true, false);
  if(isFlipped(data_mesh, UV_ext)){
    return true;
  }

  std::cout << "------ TUTTE --------" <<  std::endl;
  UV_ext = tutte_ext(data_mesh);
  if(isFlipped(data_mesh, UV_ext)){
    return true;
  }
  
  std::cout << "------ PASSED THE FILTER --------" <<  std::endl;
  return false;

}

bool toFilterSingle(Eigen::MatrixXd &V, Eigen::MatrixXi F, std::string method){
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

  if(method == "harmonic"){
    UV_ext = harmonic(data_mesh, false);
  }
  else if(method == "conformal"){
    UV_ext = LSCM(data_mesh, true, false);
  }
  else if(method == "tutte"){
    UV_ext = tutte_ext(data_mesh);
  }
  else{
    std::cout << "wrong method name" << std::endl;
    return true;
  }

  if(isFlipped(data_mesh, UV_ext)){
    return true;
  }
  
  return false;
}

int main(int argc, char *argv[]) {
  if(argc < 3){
    std::cout << "Specify mesh and output folder" << std::endl;
    std::cout << "./filter <mesh> <output> [method]" << std::endl;
    std::cout << "<method> in [harmonic, conformal, tutte]" << std::endl;
    std::cout << "if no method is specified, runs on all methods" << std::endl;
    std::cout << "<mesh>   as filepath" << std::endl;
    std::cout << "<output> as directory" << std::endl;
  }

  std::filesystem::path mesh_path = std::string(argv[1]);
  std::filesystem::path folder_path = std::string(argv[2]);
  auto target = folder_path / mesh_path.filename();
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  read_mesh(mesh_path, V, F);
  if(argc == 4){
    std::string method = argv[3];
    if(!toFilterSingle(V, F, method))
      std::filesystem::copy_file(mesh_path, target, std::filesystem::copy_options::overwrite_existing);
    return 0;
  }
  if(!toFilter(V,F)){
    std::filesystem::copy_file(mesh_path, target, std::filesystem::copy_options::overwrite_existing);
  }

  return 0;
}

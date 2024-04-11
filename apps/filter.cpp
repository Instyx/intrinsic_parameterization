#include "parameterization.hpp"
#include "lscm.hpp"
#include "harmonic.hpp"
#include "adjacency.hpp"
#include "map_to_boundary.hpp"

#include <string>
#include <filesystem>
#include "read_mesh.hpp"
#include <igl/copyleft/cgal/orient2D.h>
#include <igl/cotmatrix.h>


void flipped_elements(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::VectorXi& I
){
  I.setZero(F.rows());
  for(int i=0;i<F.rows();i++){
    double a[2] = {V(F(i,0),0),V(F(i,0),1)};
    double b[2] = {V(F(i,1),0),V(F(i,1),1)};
    double c[2] = {V(F(i,2),0),V(F(i,2),1)};
    if(igl::copyleft::cgal::orient2D(a,b,c) <= 0) {
      I(i) = 1; // if cw or collinear, it's flipped
    }else
      I(i) = 0; // if ccw, it's not flipped
  }
}

bool isFlipped(DataGeo &data_mesh, Eigen::MatrixXd UV){
  auto F = data_mesh.F;
  for(int i=0;i<F.rows();i++){
    double a[2] = {UV(F(i,0),0),UV(F(i,0),1)};
    double b[2] = {UV(F(i,1),0),UV(F(i,1),1)};
    double c[2] = {UV(F(i,2),0),UV(F(i,2),1)};
    if(igl::copyleft::cgal::orient2D(a,b,c) <= 0) {
      return true;
    }
  }
  return false;
}

int toFilter(Eigen::MatrixXd &V, Eigen::MatrixXi &F){
  DataGeo data_mesh;
  data_mesh.V = V;
  data_mesh.F = F;
  data_mesh.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh.F));
  data_mesh.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh.inputMesh, data_mesh.V));
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->requireFaceAreas();
  Eigen::Matrix<double, -1, 2> UV_ext;
  Eigen::VectorXi B;
  bdy_loop(F, V, B);

  std::cout << "------ HARMONIC --------" <<  std::endl;
  // UV_ext = harmonic(data_mesh, false);
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  circle_boundary_proportional(V, B, UV_ext);
  harmonic(-L,B,UV_ext,2);
  if(isFlipped(data_mesh, UV_ext)){
    return 1;
  }

  std::cout << "------ CONFORMAL --------" << std::endl;
  lscm(F, V, -L, UV_ext);
  // UV_ext = LSCM(data_mesh, true, false);
  if(isFlipped(data_mesh, UV_ext)){
    return 2;
  }

  std::cout << "------ TUTTE --------" <<  std::endl;
  UV_ext = tutte_ext(data_mesh);
  if(isFlipped(data_mesh, UV_ext)){
    return 3;
  }

  std::cout << "------ PASSED THE FILTER --------" <<  std::endl;
  return 0;

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
  if(argc < 2){
    std::cout << "Specify mesh and output folder" << std::endl;
    std::cout << "./filter <mesh> [method]" << std::endl;
    std::cout << "<method> in [harmonic, conformal, tutte]" << std::endl;
    std::cout << "if no method is specified, runs on all methods" << std::endl;
    std::cout << "<mesh>   as filepath" << std::endl;
    // std::cout << "<output> as directory" << std::endl;
    return 0;
  }

  std::filesystem::path mesh_path = std::string(argv[1]);
  // std::filesystem::path folder_path = std::string(argv[2]);
  // auto target = folder_path / mesh_path.filename();
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  read_mesh(mesh_path, V, F);
  if(argc == 3){
    std::string method = argv[2];
    if(!toFilterSingle(V, F, method))
      std::cout << "Passed" << std::endl;
    else
      std::cout << "Failed: " << method << std::endl;
    return 0;
  }
  int f = toFilter(V,F);
  std::string method[4] = {"", "harmonic", "conformal", "tutte"};
  if(f == 0){
    std::cout << "Passed" << std::endl;
  } else {
    std::cout << "Failed: " << method[f] << std::endl;
  }

  return 0;
}

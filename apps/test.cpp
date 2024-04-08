#include <string>
#include <filesystem>
#include <iostream>
#include <Eigen/Sparse>

#include <igl/cotmatrix.h>
#include <igl/map_vertices_to_circle.h>

#include "read_mesh.hpp"
#include "adjacency.hpp"
#include "time.hpp"

#include "harmonic.hpp"
#include "tutte.hpp"
#include "lscm.hpp"
#include "solver.hpp"
#include "map_to_boundary.hpp"
#include "parameterization.hpp"
#include <igl/opengl/glfw/Viewer.h>

void uniform_laplace(Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& L){
  typedef Eigen::Triplet<double> T;
  std::vector<T> ijv;
  ijv.reserve(F.size()*2);
  for(int i = 0;i<F.rows();i++)
  {
    // Loop over this **simplex**
    for(int j = 0;j<F.cols();j++)
    for(int k = j+1;k<F.cols();k++)
    {
      // Get indices of edge: s --> d
      int s = F(i,j);
      int d = F(i,k);
      ijv.push_back(T(s,s,-1)); // insert any diagonal elements (it just reserves the space)
      ijv.push_back(T(d,d,-1));
      ijv.push_back(T(s,d,1)); // Every hinge edge is considered twice, need to clip to one in next step
      ijv.push_back(T(d,s,1));
    }
  }
  const int n = F.maxCoeff()+1;
  L.resize(n,n);
  L.setFromTriplets(ijv.begin(),ijv.end());
  // force all non-zeros to be negative one
  for(int k=0; k<L.outerSize(); ++k) {
    for(typename Eigen::SparseMatrix<double>::InnerIterator it(L,k); it; ++it) {
      assert(it.value() != 0);
      L.coeffRef(it.row(),it.col()) = -1;
    }
  }
  // Set diagonal
  int sum;
  for(int k=0; k<L.outerSize(); ++k) {
    sum = 0;
    for(typename Eigen::SparseMatrix<double>::InnerIterator it(L,k); it; ++it) {
      ++sum;
    }
    L.coeffRef(k,k) = (sum-1);
  }
}

int main(int argc, char *argv[]) {
  std::filesystem::path mesh_path = std::string(argv[1]);
  // std::filesystem::path folder_path = std::string(argv[2]);
  // auto target = folder_path / mesh_path.filename();
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  read_mesh(mesh_path, V, F);
  Eigen::VectorXi VT, VTi;
  Eigen::Matrix<int, -1, 3> TT;
  Eigen::VectorXi VV, VVi;
  Eigen::VectorXi B;
  Eigen::SparseMatrix<double> L;
  Eigen::Matrix<double, -1, 2> UV;

  DataGeo data_mesh;
  Eigen::Matrix<double, -1, 2> X;
  data_mesh.V = V;
  data_mesh.F = F;
  data_mesh.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh.F));
  data_mesh.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh.inputMesh, data_mesh.V));
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->requireFaceAreas();
  // TIME_BLOCK("old",
  //   // UV = harmonic(data_mesh, false);
  //   UV = LSCM(data_mesh, true, false);
  // )

  TIME_BLOCK("new",
    // vt_adjacency(F, V, VT, VTi);
    // tt_adjacency(F, VT, VTi, TT);
    // vv_adjacency(F, V, TT, VV, VVi);
    // bdy_loop(F, TT, VT, VTi, B);
    igl::cotmatrix(V,F,L);

    // circle_boundary_proportional(V, B, UV);
    // harmonic(-L, B, UV, 2);
    // solve_with_known_cholmod(-L, B, UV);
    lscm(F, V, -L, UV);
  )

  Eigen::SparseMatrix<double> A;

  // circle_boundary_proportional(V, B, UV);
  // uniform_laplace(F, A);
  // TIME_BLOCK("tutte_h",
  //   harmonic(A, B, UV, 2);
  // )
  // circle_boundary_proportional(V, B, UV);
  // TIME_BLOCK("tutte_d",
  //   tutte(VV, VVi, B, UV, 2);
  // )

  // igl::opengl::glfw::Viewer viewer;
  // viewer.data().set_mesh(UV, F);
  // viewer.launch();

  // igl::opengl::glfw::Viewer viewer2;
  // viewer2.data().set_mesh(UV, F);
  // viewer2.launch();
  return 0;
}

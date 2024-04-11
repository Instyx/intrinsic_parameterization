#include <string>
#include <filesystem>
#include <iostream>
#include <Eigen/Sparse>

#include <igl/cotmatrix.h>
#include <igl/map_vertices_to_circle.h>

#include "Eigen/src/Core/Matrix.h"
#include "distortion_energy.hpp"
#include "energy.hpp"
#include "read_mesh.hpp"
#include "adjacency.hpp"
#include "time.hpp"

#include "harmonic.hpp"
#include "tutte.hpp"
#include "lscm.hpp"
#include "solver.hpp"
#include "map_to_boundary.hpp"
#include "parameterization.hpp"
#include "intrinsicflip.hpp"
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
  Eigen::Matrix<double, -1, 2> UV, UV_new;
  Eigen::MatrixXd UV_old, UV_old2;

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
  Eigen::VectorXd E;
  TIME_BLOCK("ext",
    // vt_adjacency(F, V, VT, VTi);
    // tt_adjacency(F, VT, VTi, TT);
    // vv_adjacency(F, V, TT, VV, VVi);
    // bdy_loop(F, TT, VT, VTi, B);
    igl::cotmatrix(V,F,L);

    // circle_boundary_proportional(V, B, UV);
    // harmonic(-L, B, UV, 2);
    // solve_with_known_cholmod(-L, B, UV);
    lscm(F, V, -L, UV);
    //UV_old = LSCM(data_mesh, true, false);    
    //computeParameterization(data_mesh, V, F, UV_old, UV_old2, true, false, '3');

    //tri_wise_energy(data_mesh, UV, asap, false, E);
    //std::cout << "energy with tri_wise: " << E.sum() << std::endl; 
    //std::cout << "energy with old lscm impl: " << compute_total_energy_localjacob(data_mesh, UV_old, EnergyType::ASAP) << std::endl;
    //std::cout << "energy with old lscm impl w ext grad: " << compute_total_energy(data_mesh, UV_old, EnergyType::ASAP, false) << std::endl;
    //std::cout << "energy with old2 lscm impl w ext grad: " << compute_total_energy(data_mesh, UV_old2, EnergyType::ASAP, false) << std::endl;
    //std::cout << "energy with old lscm impl w int grad: " << compute_total_energy(data_mesh, UV_old, EnergyType::ASAP, true) << std::endl;
    //std::cout << "energy with old2 lscm impl w int grad: " << compute_total_energy(data_mesh, UV_old2, EnergyType::ASAP, true) << std::endl;
    std::cout << "energy: " << compute_total_energy_localjacob(data_mesh, UV, EnergyType::ASAP) << std::endl;
    //std::cout << "energy w ext grad: " << compute_total_energy(data_mesh, UV, EnergyType::ASAP, false) << std::endl;
    //std::cout << "energy w intt grad: " << compute_total_energy(data_mesh, UV, EnergyType::ASAP, true) << std::endl;
  )
  unsigned total_flips;
  TIME_BLOCK("int flip",
    unsigned total_del_flips = 0;
    total_flips = edgeorder_flip(data_mesh, UV, total_del_flips, EnergyType::ASAP);
    //std::cout << "energy w grad: " << compute_total_energy(data_mesh, UV, EnergyType::ASAP, true) << std::endl;
  )
  std::cout << "total fips: " << total_flips << std::endl;
  std::cout << "energy: " << compute_total_energy_localjacob(data_mesh, UV, EnergyType::ASAP) << std::endl;
  Eigen::Matrix<double, -1, 2> UV_gs_maxitr, UV_gs_convergence;
  Eigen::MatrixXi F_new;
  TIME_BLOCK("laplacian",
    F_new = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();

    data_mesh.intTri->requireCotanLaplacian();
    L =  data_mesh.intTri->cotanLaplacian;
  )

  TIME_BLOCK("iparam", 
    lscm(F_new, V, L, UV_new);
   )

  std::cout << "energy: " << compute_total_energy(data_mesh, UV_new, EnergyType::ASAP, true) << std::endl;
  double energy_convergence, energy_maxitr; 
  TIME_BLOCK("gauss conv",
    //lscm(F_new, V, L, UV, UV_gs_convergence);
    //std::cout << "energy w grad: " << compute_total_energy(data_mesh, UV_gs_convergence, EnergyType::ASAP, true) << std::endl;
  )
  energy_convergence = compute_total_energy_localjacob(data_mesh, UV_gs_convergence, EnergyType::ASAP);
  
  TIME_BLOCK("gauss 5 itr",
    lscm(F_new, V, L, UV, 5, UV_gs_maxitr);
    //std::cout << "energy w grad: " << compute_total_energy(data_mesh, UV_gs_maxitr, EnergyType::ASAP, true) << std::endl;
  )
  energy_maxitr = compute_total_energy_localjacob(data_mesh, UV_gs_maxitr, EnergyType::ASAP);

  std::cout << "energy after convergence: " << energy_convergence << std::endl;
  std::cout << "energy after 5 itr: " << energy_maxitr << std::endl;
  std::cout << "convergence energy decrease: " << energy_convergence/energy_maxitr << std::endl;

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

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(UV, F);
  viewer.launch();

  igl::opengl::glfw::Viewer viewer2;
  viewer2.data().set_mesh(UV_old, F);
  viewer2.launch();

  igl::opengl::glfw::Viewer viewer3;
  viewer3.data().set_mesh(UV_old2, F);
  viewer3.launch();

  return 0;
}

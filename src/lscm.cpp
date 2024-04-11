#include "lscm.hpp"

#include <igl/vector_area_matrix.h>
#include <igl/repdiag.h>
#include <igl/boundary_loop.h>

#include "solver.hpp"
#include "time.hpp"
#include "adjacency.hpp"


void vector_area_matrix(const Eigen::MatrixXd& V, const Eigen::VectorXi& B, Eigen::SparseMatrix<double>& A){
  //Prepare a vector of triplets to set the matrix
  typedef Eigen::Triplet<double> T;
  std::vector<T> ijv;
  const int n = V.rows();
  ijv.reserve(2*n);

  for(int k = 0; k < B.rows(); k++)
  {
    int i = B(k);
    int j = B((k+1)%B.rows());
        ijv.push_back(T(i+n, j, 0.25));
        ijv.push_back(T(j, i+n, 0.25));
        ijv.push_back(T(i, j+n, -0.25));
        ijv.push_back(T(j+n, i, -0.25));
  }

  //Set A from triplets (Eigen will sum triplets with same coordinates)
  A.resize(n * 2, n * 2);
  A.setFromTriplets(ijv.begin(), ijv.end());
}

void lscm(
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& V,
  const Eigen::SparseMatrix<double>& L,
  Eigen::Matrix<double, -1,2>& UV
){
  const int nV = V.rows();

  Eigen::VectorXi boundary;
  bdy_loop(F, V, boundary);

  Eigen::SparseMatrix<double> A;
  vector_area_matrix(V, boundary, A);

  // Assemble the cotan laplacian matrix
  Eigen::SparseMatrix<double> L_flat;
  igl::repdiag(L,2,L_flat);
  Eigen::SparseMatrix<double> Q;
  Q = L_flat - 2.*A;

  Eigen::VectorXi B;
  B.resize(4);

  Eigen::MatrixXd X;

  X.resize(nV*2,1);
  // igl::boundary_loop(F, boundary);
  unsigned fixed1 = boundary[0], fixed2 = boundary[0];
  for (unsigned i = 0; i < boundary.rows(); ++i) {
    for (unsigned j = 0; j < boundary.rows(); ++j) {
      if((V.row(boundary[i])-V.row(boundary[j])).norm()>(V.row(fixed1)-V.row(fixed2)).norm()){
        fixed1 = boundary[i];
        fixed2 = boundary[j];
      }
    }
  }
  B(0) = fixed1;
  B(1) = fixed2;
  B(2) = fixed1+nV;
  B(3) = fixed2+nV;
  X(B(0)) = 1;
  X(B(0)+nV) = 0;
  X(B(1)) = 0;
  X(B(1)+nV) = 1;
  solve_with_known_try_cholmod(Q, B, X);
  UV.resize(nV,2);
  for (unsigned i=0;i<UV.cols();++i)
  {
    UV.col(i) = X.block(nV*i,0,nV,1);
  }
}

void lscm(
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& V,
  const Eigen::SparseMatrix<double>& L,
  const unsigned fixed1,
  const unsigned fixed2,
  Eigen::MatrixXd& UV
){
  const int nV = V.rows();

  Eigen::VectorXi boundary;
  bdy_loop(F, V, boundary);

  Eigen::SparseMatrix<double> A;
  vector_area_matrix(V, boundary, A);

  // Assemble the cotan laplacian matrix
  Eigen::SparseMatrix<double> L_flat;
  igl::repdiag(L,2,L_flat);
  Eigen::SparseMatrix<double> Q;
  Q = L_flat - 2.*A;

  Eigen::VectorXi B;
  B.resize(4);

  Eigen::MatrixXd X;

  X.resize(nV*2,1);
  // igl::boundary_loop(F, boundary);
  B(0) = fixed1;
  B(1) = fixed2;
  B(2) = fixed1+nV;
  B(3) = fixed2+nV;
  X(B(0)) = 1;
  X(B(0)+nV) = 0;
  X(B(1)) = 0;
  X(B(1)+nV) = 1;
  solve_with_known_try_cholmod(Q, B, X);
  UV.resize(nV,2);
  for (unsigned i=0;i<UV.cols();++i)
  {
    UV.col(i) = X.block(nV*i,0,nV,1);
  }
}


void lscm(
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& V,
  const Eigen::SparseMatrix<double>& L,
  const Eigen::MatrixXd& UV_init,
  const unsigned fixed1,
  const unsigned fixed2,
  Eigen::MatrixXd& UV
){
  const int nV = V.rows();

  Eigen::VectorXi boundary;
  bdy_loop(F, V, boundary);

  Eigen::SparseMatrix<double> A;
  vector_area_matrix(V, boundary, A);

  // Assemble the cotan laplacian matrix
  Eigen::SparseMatrix<double> L_flat;
  igl::repdiag(L,2,L_flat);
  Eigen::SparseMatrix<double> Q;
  Q = L_flat - 2.*A;

  Eigen::VectorXi B;
  B.resize(4);

  Eigen::MatrixXd X, X_init;

  X.resize(nV*2,1);

  X << UV_init.col(0), UV_init.col(1);
  
  B(0) = fixed1;
  B(1) = fixed2;
  B(2) = fixed1+nV;
  B(3) = fixed2+nV;
  X(B(0)) = 1;
  X(B(0)+nV) = 0;
  X(B(1)) = 0;
  X(B(1)+nV) = 1;
  solve_with_known_gs(Q, B, X);
  UV.resize(nV,2);
  for (unsigned i=0;i<UV.cols();++i)
  {
    UV.col(i) = X.block(nV*i,0,nV,1);
  }
}

void lscm(
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& V,
  const Eigen::SparseMatrix<double>& L,
  const Eigen::Matrix<double,-1,2>& UV_init,
  const int maxIter,
  Eigen::Matrix<double, -1,2>& UV
){
  const int nV = V.rows();

  Eigen::VectorXi boundary;
  bdy_loop(F, V, boundary);

  Eigen::SparseMatrix<double> A;
  vector_area_matrix(V, boundary, A);

  // Assemble the cotan laplacian matrix
  Eigen::SparseMatrix<double> L_flat;
  igl::repdiag(L,2,L_flat);
  Eigen::SparseMatrix<double> Q;
  Q = L_flat - 2.*A;

  Eigen::VectorXi B;
  B.resize(4);

  Eigen::MatrixXd X, X_init;

  X.resize(nV*2,1);

  X << UV_init.col(0), UV_init.col(1);
  // igl::boundary_loop(F, boundary);
  unsigned fixed1 = boundary[0], fixed2 = boundary[0];
  for (unsigned i = 0; i < boundary.rows(); ++i) {
    for (unsigned j = 0; j < boundary.rows(); ++j) {
      if((V.row(boundary[i])-V.row(boundary[j])).norm()>(V.row(fixed1)-V.row(fixed2)).norm()){
        fixed1 = boundary[i];
        fixed2 = boundary[j];
      }
    }
  }
  B(0) = fixed1;
  B(1) = fixed2;
  B(2) = fixed1+nV;
  B(3) = fixed2+nV;
  X(B(0)) = 1;
  X(B(0)+nV) = 0;
  X(B(1)) = 0;
  X(B(1)+nV) = 1;
  solve_with_known_gs(Q, B, maxIter, X);
  UV.resize(nV,2);
  for (unsigned i=0;i<UV.cols();++i)
  {
    UV.col(i) = X.block(nV*i,0,nV,1);
  }
}

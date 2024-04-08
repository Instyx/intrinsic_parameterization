#include "harmonic.hpp"

#include <Eigen/CholmodSupport>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>


#include <type_traits>
#include <iostream>
#include <math.h>
#include "map_to_boundary.hpp"

// Solves  L*x = 0 s.t. X_known stays X_known (constraints need to be set in result X)
template <typename Derived, typename Solver>
void solve_with_known(
      const Eigen::SparseMatrix<typename Derived::Scalar>& L,
      const Eigen::Ref<const Eigen::VectorXi> known,
      Solver& solver,
      Eigen::MatrixBase<Derived>& X)
{
    typedef typename Derived::Scalar Scalar;
    const int nv = L.rows();
    Eigen::VectorXi idx = Eigen::VectorXi::Zero(nv);
    int bc = known.size();
    if (bc == nv){
      return;
    }
    for (int i = 0; i < bc; i++){
      idx[known[i]] = -i-1;
    }
    int row = 0;
    for (int i = 0; i < nv; i++) if (idx[i]==0) idx[i] = row++;

    std::vector<Eigen::Triplet<Scalar>> coef;
    coef.reserve( (nv-known.rows())*6 );
    Eigen::Matrix<Scalar,-1,-1> b = Eigen::Matrix<Scalar,-1,-1>::Zero(nv-known.rows(),X.cols());
    for (int i = 0; i < L.outerSize(); ++i) {
      for (typename Eigen::SparseMatrix<Scalar,Eigen::ColMajor>::InnerIterator it(L,i); it; ++it) {
        if (idx[it.row()] >= 0){
          if (idx[it.col()] >= 0){
            coef.push_back(Eigen::Triplet<Scalar>(idx[it.row()], idx[it.col()], it.value()));
          } else {
            b.row(idx[it.row()]) -= X.row(it.col())*it.value();
          }
        }
      }
    }

    Eigen::SparseMatrix<Scalar> A(nv-known.rows(),nv-known.rows());
    A.setFromTriplets(coef.begin(), coef.end());
    Eigen::Matrix<Scalar,-1,-1> x;

    solver.compute(A);
    if(solver.info() != Eigen::Success){
      return;
    }

    // divide into two functions?
    x = solver.solve(b);

    row = 0;
    for (int i = 0; i < nv; i++)
        if (idx[i] >= 0)
            X.row(i) = x.row(row++);
}

void solve_with_known_cholmod(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::MatrixXd &X){
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver;
  solve_with_known(L, known, solver, X);
}

void solve_with_known_cholmod(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::Matrix<double,-1,2> &X){
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver;
  solve_with_known(L, known, solver, X);
}

void solve_with_known_lu(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::MatrixXd &X){
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solve_with_known(L, known, solver, X);
}

void solve_with_known_lu(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::Matrix<double,-1,2> &X){
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solve_with_known(L, known, solver, X);
}

void solve_with_known_try_cholmod(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::MatrixXd &X){
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver;
  solve_with_known(L, known, solver, X);
  if(solver.info() != Eigen::Success){
    Eigen::SparseLU<Eigen::SparseMatrix<double>> backup;
    solve_with_known(L, known, backup, X);
  }
}

void solve_with_known_try_cholmod(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::Matrix<double,-1,2> &X){
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver;
  solve_with_known(L, known, solver, X);
  if(solver.info() != Eigen::Success){
    Eigen::SparseLU<Eigen::SparseMatrix<double>> backup;
    solve_with_known(L, known, backup, X);
  }
}

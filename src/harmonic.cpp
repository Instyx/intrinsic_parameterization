#include "harmonic.hpp"

#include <Eigen/CholmodSupport>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>


#include <type_traits>
#include <iostream>
#include <math.h>
#include "map_to_boundary.hpp"


template <typename Derived>
void harmonic(const Eigen::SparseMatrix<typename Derived::Scalar>& L,
              const Eigen::Ref<const Eigen::VectorXi> B,
              Eigen::MatrixBase<Derived>& X,
              int boundary_type)
{
    typedef typename Derived::Scalar Scalar;

    const int nv = L.rows();

    X.derived().resize(nv,2);
    switch (boundary_type) {
      case 2:
        break;
      case 1:
        square_boundary(B, X);
        break;
      case 0:
      default:
        circle_boundary(B, X);
    }
    Eigen::VectorXi idx = Eigen::VectorXi::Zero(nv);
    int bc = B.size();
    if (bc == nv){
      return;
    }
    for (int i = 0; i < bc; i++){
      idx[B[i]] = -i-1;
    }
    int row = 0;
    for (int i = 0; i < nv; i++) if (idx[i]==0) idx[i] = row++;

    std::vector<Eigen::Triplet<Scalar>> coef;
    coef.reserve( (nv-B.rows())*6 );
    Eigen::Matrix<Scalar,Eigen::Dynamic,2> b = Eigen::Matrix<Scalar,Eigen::Dynamic,2>::Zero(nv-B.rows(),2);
    for (int i = 0; i < L.outerSize(); ++i) {
      for (typename Eigen::SparseMatrix<Scalar,Eigen::ColMajor>::InnerIterator it(L,i); it; ++it) {
        if (idx[it.row()] >= 0){
          if (idx[it.col()] >= 0){
            coef.push_back(Eigen::Triplet<Scalar>(idx[it.row()], idx[it.col()], it.value()));
          } else {
            b(idx[it.row()],0) -= X(it.col(),0)*it.value();
            b(idx[it.row()],1) -= X(it.col(),1)*it.value();
          }
        }
      }
    }

    Eigen::SparseMatrix<Scalar> A(nv-B.rows(),nv-B.rows());
    A.setFromTriplets(coef.begin(), coef.end());
    Eigen::Matrix<Scalar,Eigen::Dynamic,2> x;
    if constexpr(std::is_same_v<Scalar, double>)
    {
      Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<Scalar> > solver;
      solver.compute(A);
      x = solver.solve(b);
    } else {
      Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar> > solver;
      solver.compute(A);
      x = solver.solve(b);
    }

    row = 0;
    for (int i = 0; i < nv; i++)
        if (idx[i] >= 0)
            X.row(i) = x.row(row++);
}


template void harmonic<Eigen::Matrix<double,-1,2>>(
            const Eigen::SparseMatrix<double>& L,
           const Eigen::Ref<const Eigen::VectorXi> B,
           Eigen::MatrixBase<Eigen::Matrix<double,-1,2>> & X,
           int boundary_type=0);

template void harmonic<Eigen::Matrix<float,-1,2>>(
          const Eigen::SparseMatrix<float>& L,
          const Eigen::Ref<const Eigen::VectorXi> B,
          Eigen::MatrixBase<Eigen::Matrix<float,-1,2>> & X,
          int boundary_type=0);

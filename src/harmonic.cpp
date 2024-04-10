#include "harmonic.hpp"

#include <Eigen/CholmodSupport>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>


#include <type_traits>
#include <iostream>
#include <math.h>
#include "map_to_boundary.hpp"
#include "solver.hpp"


template <typename Derived>
void harmonic(const Eigen::SparseMatrix<typename Derived::Scalar>& L,
              const Eigen::Ref<const Eigen::VectorXi> B,
              Eigen::MatrixBase<Derived>& X,
              int boundary_type)
{
    // typedef typename Derived::Scalar Scalar;

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
    //  std::cout << X << std::endl;
    solve_with_known_cholmod(L,B,X.derived());
}

template <typename Derived>
void harmonic(const Eigen::SparseMatrix<typename Derived::Scalar>& L,
              const Eigen::Ref<const Eigen::VectorXi> B,
              const Eigen::MatrixBase<Derived>& X_init,
              Eigen::MatrixBase<Derived>& X,
              int boundary_type)
{
    // typedef typename Derived::Scalar Scalar;

    const int nv = L.rows();

    X.derived().resize(nv,2);
    X = X_init;
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
    // std::cout << X << std::endl;
    solve_with_known_gs(L,B,X.derived());
}



template void harmonic<Eigen::Matrix<double,-1,2>>(
            const Eigen::SparseMatrix<double>& L,
           const Eigen::Ref<const Eigen::VectorXi> B,
           Eigen::MatrixBase<Eigen::Matrix<double,-1,2>> & X,
           int boundary_type=0);

template void harmonic<Eigen::MatrixXd>(
            const Eigen::SparseMatrix<double>& L,
           const Eigen::Ref<const Eigen::VectorXi> B,
           Eigen::MatrixBase<Eigen::MatrixXd> & X,
           int boundary_type=0);


template void harmonic<Eigen::Matrix<double,-1,2>>(
            const Eigen::SparseMatrix<double>& L,
           const Eigen::Ref<const Eigen::VectorXi> B,
           const Eigen::MatrixBase<Eigen::Matrix<double,-1,2>> & X_init,
           Eigen::MatrixBase<Eigen::Matrix<double,-1,2>> & X,
           int boundary_type=0);

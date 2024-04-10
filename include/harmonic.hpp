#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

template <typename Derived>
void harmonic(
           const Eigen::SparseMatrix<typename Derived::Scalar>& L,
           const Eigen::Ref<const Eigen::VectorXi> B,
           Eigen::MatrixBase<Derived> & X,
           int boundary_type);

template <typename Derived>
void harmonic(const Eigen::SparseMatrix<typename Derived::Scalar>& L,
              const Eigen::Ref<const Eigen::VectorXi> B,
              const Eigen::MatrixBase<Derived>& X_init,
              Eigen::MatrixBase<Derived>& X,
              int boundary_type);

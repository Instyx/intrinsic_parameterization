#pragma once

#include <Eigen/Core>


template <typename DerivedIn,typename DerivedOut>
void circle_boundary_proportional(
  const Eigen::MatrixBase<DerivedIn>& V,
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<DerivedOut>& X);


template <typename Derived>
void circle_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Derived>& X);


template <typename Derived>
void square_boundary(
  const Eigen::Ref<const Eigen::VectorXi> B,
  Eigen::MatrixBase<Derived>& X);

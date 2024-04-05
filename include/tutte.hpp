#pragma once

#include <Eigen/Core>

template <typename Derived>
void tutte(const Eigen::Ref<const Eigen::VectorXi> VV,
           const Eigen::Ref<const Eigen::VectorXi> VVi,
           const Eigen::Ref<const Eigen::VectorXi> B,
           //Eigen::Ref<Eigen::MatrixXd> X)
           Eigen::MatrixBase<Derived> & X,
           int boundary_type);

#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

void lscm(
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& V,
  const Eigen::SparseMatrix<double>& L,
  Eigen::Matrix<double, -1,2>& UV
);

#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

void lscm(
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& V,
  const Eigen::SparseMatrix<double>& L,
  Eigen::Matrix<double, -1,2>& UV
);

void lscm(
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& V,
  const Eigen::SparseMatrix<double>& L,
  const Eigen::Matrix<double,-1,2>& UV_init,
  Eigen::Matrix<double, -1,2>& UV
);

void lscm(
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& V,
  const Eigen::SparseMatrix<double>& L,
  const Eigen::Matrix<double,-1,2>& UV_init,
  const int maxIter,
  Eigen::Matrix<double, -1,2>& UV
);

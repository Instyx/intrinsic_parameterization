#pragma once

#include <Eigen/Sparse>

Eigen::SparseMatrix<double> eiLaplace(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
Eigen::SparseMatrix<double> asym_but_pos_eiLaplace(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
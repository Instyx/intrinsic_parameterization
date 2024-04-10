#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>


void gauss_seidel(const Eigen::SparseMatrix<double> &A, const Eigen::Matrix<double, -1, -1> &B, Eigen::Matrix<double, -1, -1> &X);

void gauss_seidel(const Eigen::SparseMatrix<double> &A, const Eigen::Matrix<double, -1, -1> &B, const int maxIter, Eigen::Matrix<double, -1, -1> &X);

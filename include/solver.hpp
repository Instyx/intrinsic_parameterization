#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

void solve_with_known_cholmod(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::Matrix<double,-1,2> &X);
void solve_with_known_cholmod(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::MatrixXd &X);


void solve_with_known_lu(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::MatrixXd &X);
void solve_with_known_lu(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::Matrix<double,-1,2> &X);

void solve_with_known_try_cholmod(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::MatrixXd &X);
void solve_with_known_try_cholmod(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::Matrix<double,-1,2> &X);

void solve_with_known_gs(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::MatrixXd &X);
void solve_with_known_gs(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, Eigen::Matrix<double, -1, 2> &X);

void solve_with_known_gs(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, const int maxIter, Eigen::MatrixXd &X);
void solve_with_known_gs(const Eigen::SparseMatrix<double> &L, const Eigen::Ref<const Eigen::VectorXi> known, const int maxIter, Eigen::Matrix<double, -1, 2> &X);

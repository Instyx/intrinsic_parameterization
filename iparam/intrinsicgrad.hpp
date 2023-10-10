#pragma once

#include <Eigen/Core>
#include <datageo.hpp>

void computeSurfaceGradientMatrix(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::SparseMatrix<double> & D1, Eigen::SparseMatrix<double> & D2);

void computeGrad_intrinsic(DataGeo &data_mesh, Eigen::SparseMatrix<double> & Dx, 
    Eigen::SparseMatrix<double> & Dy, Eigen::VectorXd &areas);


#pragma once

#include <Eigen/Core>
#include <datageo.hpp>
void computeGrad_intrinsic(DataGeo &data_mesh, Eigen::SparseMatrix<double> & Dx, 
    Eigen::SparseMatrix<double> & Dy, Eigen::VectorXd &areas);


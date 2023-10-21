#pragma once
#include <Eigen/Core>
#include <datageo.hpp>
#include <igl/slim.h>

void slim_parameterization(DataGeo &data_mesh, igl::SLIMData &slimdata, Eigen::MatrixXd &UV, bool igrad, bool isFreeBoundary);
void boundary(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, bool isFreeBoundary, Eigen::VectorXi &fixed_UV_indices, Eigen::MatrixXd &fixed_UV_positions); 

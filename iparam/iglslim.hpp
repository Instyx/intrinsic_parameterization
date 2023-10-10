#pragma once
#include <Eigen/Core>
#include <datageo.hpp>
#include <intrinsicflip.hpp>
#include <igl/slim.h>
void slim_parameterization(DataGeo &data_mesh, igl::SLIMData &slimdata, Eigen::MatrixXd &UV, bool igrad, bool isFreeBoundary, const FlipType &ft);

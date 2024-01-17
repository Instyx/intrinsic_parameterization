#pragma once

#include <Eigen/Core>
#include <datageo.hpp>
#include <intrinsicflip.hpp>
#include <igl/slim.h>



void slim_tillconverges(DataGeo &data_mesh, igl::SLIMData& slimdata, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &UV_init, unsigned max_iterations, bool igrad);

Eigen::MatrixXd intrinsicslim(DataGeo &data_mesh, Eigen::MatrixXd &UV_init);



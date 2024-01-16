#pragma once

#include <Eigen/Core>
#include <datageo.hpp>
#include <intrinsicflip.hpp>

Eigen::MatrixXd intrinsicslim(DataGeo &data_mesh, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &UV_init,
    unsigned number_iterations, unsigned flip_granularity, const FlipType& ft);

Eigen::MatrixXd intrinsicslim_tillconverges(DataGeo &data_mesh, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &UV_init,
    unsigned max_iterations, const FlipType& ft);

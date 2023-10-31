#pragma once

#include <Eigen/Core>
#include <datageo.hpp>


void minmax_distortions(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &UV, std::vector<double> &res);

void minmax_distortions_intri(DataGeo &data_mesh, const Eigen::MatrixXd &UV, std::vector<double> &res);

void compute_metrics(DataGeo &data_mesh, const Eigen::MatrixXd &UV_o, std::vector<double> &res);

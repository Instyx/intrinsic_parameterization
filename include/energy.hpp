#pragma once

#include <Eigen/Core>
#include "datageo.hpp"

void tri_wise_energy_int(DataGeo& mesh_data, const Eigen::MatrixXd& UV, double (*energy)(const Eigen::Matrix2d &), Eigen::VectorXd& E);

void tri_wise_energy_ext(DataGeo& mesh_data, const Eigen::MatrixXd& UV, double (*energy)(const Eigen::Matrix2d &), Eigen::VectorXd& E);

void tri_wise_energy(DataGeo& mesh_data, const Eigen::MatrixXd& UV, double (*energy)(const Eigen::Matrix2d &), const bool intrinsic, Eigen::VectorXd& E);

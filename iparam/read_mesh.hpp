#pragma once

#include <Eigen/Core>

void read_mesh(const std::string path, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

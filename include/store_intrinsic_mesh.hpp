#pragma once

#include <datageo.hpp>
#include <optimize.hpp>

void store_intrinsic_edges(const DataGeo &data_mesh, const std::string filepath);

void store_intrinsic_mesh(const DataGeo &data_mesh, const Eigen::MatrixXd &UV, const std::string filename);
void store_intrinsic_mesh(const DataGeo &data_mesh, const Eigen::MatrixXd &UV, const std::string filename, Results &res);
void store_intrinsic_mesh(const DataGeo &data_mesh, const Eigen::MatrixXd &UV, const std::string filename, const Eigen::VectorXd& beforeTriEnergy, const Eigen::VectorXd& afterTriEnergy);
void store_intrinsic_mesh(const DataGeo &data_mesh, const Eigen::MatrixXd &UV, const std::string filename, const Eigen::VectorXd& beforeTriEnergy, const Eigen::VectorXd& afterTriEnergy, Results &res);

#pragma once

#include <datageo.hpp>
#include <distortion_energy.hpp>
#include <optimize.hpp>

Results optimize_single_new(Eigen::MatrixXd &V, Eigen::MatrixXi &F, EnergyType method, std::string dir, std::string mesh_name);

Results optimize_single_new(Eigen::MatrixXd &V, Eigen::MatrixXi &F, EnergyType method, std::string dir, std::string mesh_name, const bool init_with_intrinsic, const bool priority_queue_flips);


Results optimize_single_new(Eigen::MatrixXd &V, Eigen::MatrixXi &F, EnergyType method, std::string dir, std::string mesh_name, const bool init_with_intrinsic, const bool priority_queue_flips, int random_runs);

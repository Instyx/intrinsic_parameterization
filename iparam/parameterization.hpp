#pragma once
#include "datageo.hpp"
#include <Eigen/Core>
#include "intrinsicflip.hpp"
#include "distortion_energy.hpp"
void computeParameterization(DataGeo &data_mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &UV, Eigen::MatrixXd &new_UV,
    bool isFreeBoundary, bool igrad, const FlipType &ft, const EnergyType &et, int type);

double compute_total_energy(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et);


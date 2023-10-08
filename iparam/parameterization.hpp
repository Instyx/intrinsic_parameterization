#pragma once
#include "datageo.hpp"
#include <Eigen/Core>
#include "intrinsicflip.hpp"

enum class EnergyType{
  DIRICHLET,
  ASAP,
  ARAP,
  SYMMETRIC_DIRICHLET
};

void computeParameterization(DataGeo &data_mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &UV, Eigen::MatrixXd &new_UV,
    bool isFreeBoundary, bool igrad, const FlipType &ft, const EnergyType &et, int type);




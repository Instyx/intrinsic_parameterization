#pragma once
#include <Eigen/Core>

enum class EnergyType{
  DIRICHLET,
  ASAP,
  ARAP,
  SYMMETRIC_DIRICHLET
};

// distortion energies
double dirichlet(const Eigen::Matrix2d &J);
double asap(const Eigen::Matrix2d &J);
double arap(const Eigen::Matrix2d &J);
double symmetric_dirichlet(const Eigen::Matrix2d &J);
double symmetric_dirichlet_alt(const Eigen::Matrix2d &J);

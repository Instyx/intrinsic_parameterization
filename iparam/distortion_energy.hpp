#pragma once
#include <Eigen/Core>

enum class EnergyType{
  DIRICHLET,
  ASAP,
  ARAP,
  SYMMETRIC_DIRICHLET
};


double dirichlet(Eigen::Matrix2d J);
double asap(Eigen::Matrix2d J);
double arap(Eigen::Matrix2d J);
double symmetric_dirichlet(Eigen::Matrix2d J);

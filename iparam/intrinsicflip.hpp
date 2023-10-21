#pragma once

#include <Eigen/Core>
#include "datageo.hpp"
#include <distortion_energy.hpp>
enum class FlipType { 
  EDGEORDER, 
  GREEDY,
  RANDOM,
  HEURISTIC 
};

double flippeddiff(DataGeo &data_mesh, const Eigen::MatrixXd &UV, gcs::Edge e, const EnergyType &et);

unsigned greedy_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et);

unsigned heuristic_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et);

unsigned random_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et);

unsigned edgeorder_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et);

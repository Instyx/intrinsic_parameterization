#pragma once

#include <Eigen/Core>
#include "datageo.hpp"

enum class FlipType { 
  EDGEORDER, 
  GREEDY,
  RANDOM,
  HEURISTIC 
};

double flippeddiff(DataGeo &data_mesh, const Eigen::MatrixXd &UV, gcs::Edge e, double (*energy)(Eigen::Matrix2d));

unsigned greedy_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, double (*energy)(Eigen::Matrix2d));

unsigned heuristic_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, double (*energy)(Eigen::Matrix2d));

unsigned random_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, double (*energy)(Eigen::Matrix2d));

unsigned edgeorder_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, double (*energy)(Eigen::Matrix2d));

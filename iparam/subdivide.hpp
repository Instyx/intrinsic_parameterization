#pragma once
#include "distortion_energy.hpp"
#include <Eigen/Core>
#include "datageo.hpp"

double calc_energy(DataGeo &data_mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXd &UV, gcs::Halfedge he, const EnergyType &et);

unsigned flipThroughEdges(DataGeo &data_mesh, const Eigen::MatrixXd V, const Eigen::MatrixXd UV, const EnergyType &et);

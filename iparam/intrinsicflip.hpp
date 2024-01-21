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



// computes the energy difference when the edge e flipped
// Input : data_mesh : data structure for intrinsic triangulation
//         UV : #vertices x 2 : UV vertex positions
//         e : edge to be flipped
//         et : energy type that is considered for the flip
// Output: returns the energy diff
double flippeddiff(DataGeo &data_mesh, const Eigen::MatrixXd &UV, gcs::Edge e, const EnergyType &et);


// one iteration of intrinsic flipping (goes through all edges once)
// Input : data_mesh : data structure for intrinsic triangulation
//         UV : #vertices x 2 : UV vertex positions
//         e : edge to be flipped
//         et : energy type that is considered for the flip
// Output: updates the data_mesh intrinsic triangulation
//         and returns the number of flips
unsigned greedy_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et);

unsigned heuristic_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et);

unsigned random_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et);

unsigned edgeorder_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et);

unsigned delaunay_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et);

unsigned queue_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et);

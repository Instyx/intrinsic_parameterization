#pragma once
#include "distortion_energy.hpp"
#include <Eigen/Core>
#include "datageo.hpp"
// calculates the energy of the faces in the local diamond boundary of an halfedge by subdividing the mesh 
// by calculating the common subdivision of the extrinsic and intrinsic triangulation
// Input : data_mesh : data structure for intrinsic triangulation
//         V : #vertices x 3 : mesh input 
//         UV : #vertices x 2 : UV vertex positions
//         ge : input halfedge 
//         et : energy type that is considered for the flip
// Output: returns the energy in of the faces by weighting them with face areas
double calc_energy(DataGeo &data_mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXd &UV, gcs::Halfedge he, const EnergyType &et);

// one iteration of intrinsic flipping (goes through all edges once)
// Input : data_mesh : data structure for intrinsic triangulation
//         V : #vertices x 3 : mesh input 
//         UV : #vertices x 2 : UV vertex positions
//         et : energy type that is considered for the flip
// Output: updates the data_mesh intrinsic triangulation 
//         and returns the number of flips
unsigned flipThroughEdges(DataGeo &data_mesh, const Eigen::MatrixXd V, const Eigen::MatrixXd UV, const EnergyType &et);

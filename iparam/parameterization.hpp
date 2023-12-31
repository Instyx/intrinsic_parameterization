#pragma once
#include "datageo.hpp"
#include <Eigen/Core>
#include "distortion_energy.hpp"

// Input : data_mesh : data structure for intrinsic triangulation
//         V : #vertices x 3 : mesh input 
//         UV : #vertices x 2 : initial UV vertex positions (needed for ARAP)
//         isFreeBoundary : true if free boundary, false if fixed boundary
//         igrad : true for using with intrinsic triangulation 
//                 false for using with extrinsic input triangulation
//         type : '1' : uniform laplacian
//                '2' : cotangent laplacian
//                '3' : lscm
//                '4' : ARAP
// Output: new_UV : #vertices x 2 : new UV vertex positions 
void computeParameterization(DataGeo &data_mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &UV, Eigen::MatrixXd &new_UV,
    bool isFreeBoundary, bool igrad, int type);

// computes the total distortion energy of the mesh normalized by the total area of the mesh 
// Input : data_mesh : data structure for intrinsic triangulation
//         V : #vertices x 3 : mesh input 
//         UV : #vertices x 2 : initial UV vertex positions (needed for ARAP)
//         et : energy type that is considered for the flip
//         igrad : true for using with intrinsic triangulation 
//                 false for using with extrinsic input triangulation
// Output: returns the total energy normalized by the area of the mesh 

double compute_total_energy(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et, bool igrad);

// boundary constraints are stored as static variables, for the purpose of testing this resets the static variables
void reset_constraints();

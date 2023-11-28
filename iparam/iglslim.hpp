#pragma once
#include <Eigen/Core>
#include <datageo.hpp>
#include <igl/slim.h>


// 1 iteration of SLIM method using intrinsic gradients
// Input: data structures for DataGeo and SLIMdata, 
//        UV : #vertices x 2 initial UV parameterization 
//        igrad: whether to use intrinsic properties
//        isFreeBoundary: not used (always boundary-free)
// Output: UV : new UV vertex positions
void slim_parameterization(DataGeo &data_mesh, igl::SLIMData &slimdata, Eigen::MatrixXd &UV, bool igrad, bool isFreeBoundary);

// mapping boundary vertices to UV domain
// Input: V : #vertices x 3 vertex positions 
//        F : #faces x 3 face - vertex index matrix 
//        isFreeBoundary : wheter to use free or fixed boundary 
// Output: fixed_UV_indices: boundary conditions x 1   indices of the fixed vertices 
//         fixed_UV_positions : boudnary conditions x 2   position of the fixed vertices in UV domain
void boundary(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, bool isFreeBoundary, Eigen::VectorXi &fixed_UV_indices, Eigen::MatrixXd &fixed_UV_positions); 

#pragma once

#include <Eigen/Core>

// computes the SVD of the matrix J=USV^T  
// Input : J : 2 x 2 
// Output: U : 2 x 2
//         S : 2 x 2 : diagonal includes singular values
//         V : 2 x 2 
void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V);


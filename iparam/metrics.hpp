#pragma once

#include <Eigen/Core>
#include <datageo.hpp>

// computes max/min conformal, isometric, authalic distortions using extrinsic triangulation
// Input : V : #vertices x 3 : mesh input
//         F : #faces x 3 : face - vertex index matrix
//         UV : #vertices x 2 : UV vertex positions
// Output : res[0], res[1] : min, max conformal distortion 
//          res[2], res[3] : min, max isometric distortion 
//          res[4], res[5] : min, max authalic distortion 
void minmax_distortions(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &UV, std::vector<double> &res);

// computes max/min conformal, isometric, authalic distortions using intrinsic triangulation
// Input : data_mesh : data structure for intrinsic triangulation
//         UV : #vertices x 2 : UV vertex positions
// Output : res[0], res[1] : min, max conformal distortion (intrinsic) 
//          res[2], res[3] : min, max isometric distortion (intrinsic) 
//          res[4], res[5] : min, max authalic distortion (intrinsic) 
void minmax_distortions_intri(DataGeo &data_mesh, const Eigen::MatrixXd &UV, std::vector<double> &res);



// Input : data_mesh : data structure for intrinsic triangulation
//         UV : #vertices x 2 : UV vertex positions
// Output : res[0] - res[4] extrinsic metrics 
//          res[5] - res[9] intrinsic metrics (has bugs and were not included in the thesis report) 
//          each of them consists {flipped percentage, max area distortion, average area error, max angle distortion, average angle error}
void compute_metrics(DataGeo &data_mesh, const Eigen::MatrixXd &UV_o, std::vector<double> &res);

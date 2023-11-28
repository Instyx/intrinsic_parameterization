#pragma once

#include <Eigen/Core>
#include <datageo.hpp>


// computes the extrinsic gradient operators w.r.t. local frame
// Input : V : #vertices x 3 : mesh vertex positions
//         F : #faces x 3 : face - vertex index matrix  
// Output : D1 : #faces x #vertices : gradient operator w.r.t. the local axis that runs along the 23 edge
//          D2 : #faces x #vertices : gradient operator w.r.t. the local axis in the tangent plane that is perpendicular to the 23 edge
void computeSurfaceGradientMatrix(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::SparseMatrix<double> & D1, Eigen::SparseMatrix<double> & D2);

// computes the intrinsic gradient operators w.r.t. local frame
// Input : data_mesh : data structure for intrinsic triangulation
// Output : Dx : #faces x #vertices : gradient operator w.r.t. the local axis that runs along the 23 edge
//          Dy : #faces x #vertices : gradient operator w.r.t. the local axis in the tangent plane that is perpendicular to the 23 edge
//          areas : #faces : face areas 
void computeGrad_intrinsic(DataGeo &data_mesh, Eigen::SparseMatrix<double> & Dx, 
    Eigen::SparseMatrix<double> & Dy, Eigen::VectorXd &areas);


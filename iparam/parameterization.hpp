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

//double compute_energy_ext(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &UV, const EnergyType &et);

double compute_total_energy_fast(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const Eigen::SparseMatrix<double> &Dx,
    const Eigen::SparseMatrix<double> &Dy, const Eigen::VectorXd &areas, const EnergyType &et);

double compute_total_energy_localjacob(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et);

// boundary constraints are stored as static variables, for the purpose of testing this resets the static variables
void reset_constraints();


unsigned ARAP_tillconverges(DataGeo &data_mesh, Eigen::MatrixXd &UV_init, Eigen::MatrixXd &UV, unsigned max_iterations, bool isFreeBoundary, bool igrad);

unsigned intrinsic_ARAP(DataGeo &data_mesh, Eigen::MatrixXd &UV, unsigned ARAP_maxitr, unsigned intrinsic_maxitr, bool isFreeBoundary, std::fstream &fout);

Eigen::MatrixXd LSCM(DataGeo &data_mesh, bool isFreeBoundary, bool igrad);

unsigned intrinsic_LSCM(DataGeo &data_mesh, Eigen::MatrixXd &UV, unsigned max_iterations, bool isFreeBoundary, std::fstream &fout);

Eigen::MatrixXd harmonic(DataGeo &data_mesh, bool igrad);

unsigned intrinsic_harmonic(DataGeo &data_mesh, Eigen::MatrixXd &UV, unsigned max_iterations, std::fstream &fout);


Eigen::MatrixXd tutte(DataGeo &data_mesh, bool igrad);

// these also saves inbetween meshes with textures 

unsigned intrinsic_ARAP(DataGeo &data_mesh, Eigen::MatrixXd &UV, unsigned ARAP_maxitr, unsigned intrinsic_maxitr, bool isFreeBoundary, std::fstream &fout, std::string path, std::string mesh_name);

unsigned intrinsic_harmonic(DataGeo &data_mesh, Eigen::MatrixXd &UV, unsigned max_iterations, std::fstream &fout, std::string path, std::string mesh_name);

unsigned intrinsic_LSCM(DataGeo &data_mesh, Eigen::MatrixXd &UV, unsigned max_iterations, bool isFreeBoundary, std::fstream &fout, std::string path, std::string mesh_name);



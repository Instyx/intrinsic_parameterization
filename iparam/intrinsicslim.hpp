#pragma once

#include <Eigen/Core>
#include <datageo.hpp>
#include <intrinsicflip.hpp>
#include <igl/slim.h>



unsigned slim_tillconverges(DataGeo &data_mesh, igl::SLIMData& slimdata, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, 
    const Eigen::MatrixXd &UV_init, unsigned max_iterations, bool igrad);



unsigned intrinsicslim(DataGeo &data_mesh, Eigen::MatrixXd &UV_init, Eigen::MatrixXd &UV, unsigned slim_maxitr,
    unsigned intrinsic_maxitr, std::fstream &fout);

// this also saves the inbetween mesh with textures
unsigned intrinsicslim(DataGeo &data_mesh, Eigen::MatrixXd &UV_init, Eigen::MatrixXd &UV, unsigned slim_maxitr, unsigned intrinsic_maxitr, std::fstream &fout, std::string path, std::string mesh_name);

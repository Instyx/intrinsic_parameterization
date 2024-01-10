#pragma once

#include <datageo.hpp>
#include <distortion_energy.hpp>

bool isDelaunayFlipBad(DataGeo &data_mesh, const Eigen::MatrixXd &UV);
DataGeo compareIDTvsGreedy(DataGeo &data_mesh);
unsigned flipEdgesifCoplanar(DataGeo &data_mesh, bool onlyDelaunay);

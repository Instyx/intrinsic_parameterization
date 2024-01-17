#pragma once

#include <datageo.hpp>
#include <distortion_energy.hpp>

bool isDelaunayFlipBad(DataGeo &data_mesh, const Eigen::MatrixXd &UV);
DataGeo compareIDTvsGreedy(DataGeo &data_mesh);
unsigned flipEdgesifCoplanar(DataGeo &data_mesh, bool onlyDelaunay);

void compareIDTvsIPARAM(DataGeo &data_mesh, bool isFreeBoundary, const EnergyType &et, Eigen::MatrixXd &UV_idt, Eigen::MatrixXd &UV_iparam);

bool load_mesh_test(DataGeo &datageo, Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::string filename);

void test_ARAP_single(DataGeo &data_mesh,  std::string mesh_name, bool isFreeBoundary, std::ostream &fout);
void test_ARAP();

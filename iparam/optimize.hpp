#pragma once

#include <datageo.hpp>
#include <distortion_energy.hpp>


struct Results{
  double init_time;
  double cs_time;
  double cm_time;
  int ext_vertex_count;
  int cs_vertex_count;
  std::vector<double> energies;
  std::vector<int> flips;
  std::vector<int> flips_delaunay;
  std::vector<double> flip_durations;
  std::vector<double> opt_durations;
  std::vector<int> opt_iterations;
};

Results optimize_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, EnergyType method, std::string dir, std::string mesh_name);

Results optimize_single(Eigen::MatrixXd &V, Eigen::MatrixXi &F, EnergyType method, std::string dir, std::string mesh_name, const bool init_with_intrinsic, const bool priority_queue_flips);

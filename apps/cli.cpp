#include <string>
#include <ostream>
#include <iostream>
#include "read_mesh.hpp"
#include "test.hpp"
#include <filesystem>
#include <igl/face_areas.h>
#include "optimize_new.hpp"

int main(int argc, char *argv[]) {
  // evaluation
if(argc < 4) {
    std::cout << "You have to specify the parameterization" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "./cli <method> <mesh> <output> [start] [queue]" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "<method> in [dirichlet, arap, asap, symdirichlet]" << std::endl;
    std::cout << "<mesh>   as filepath" << std::endl;
    std::cout << "<output> as directory" << std::endl;
    std::cout << "[start] in [normal (dc), delaunay]" << std::endl;
    std::cout << "[queue] in [queue (dc), prio]" << std::endl;
    return 0;
  }
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXd A;
  DataGeo data_mesh;
  std::string method = std::string(argv[1]);
  std::string mesh_path = std::string(argv[2]);
  std::string outdir = std::string(argv[3]);
  bool start = argc > 4 && (std::string(argv[4]) == "delaunay");
  // std::cout << start << "\n";
  bool queue = argc > 5 && (std::string(argv[5]) == "prio");
  // std::cout << queue << "\n";
  read_mesh(mesh_path,V,F);
  igl::doublearea(V,F,A);
  V /=sqrt(A.sum());

  std::fstream log;
  std::string filename = mesh_path.substr(mesh_path.find_last_of("/\\") + 1);
  size_t lastindex = filename.find_last_of(".");
  std::string mesh_name_wo_extension =  filename.substr(0, lastindex);
  std::string to_store_dir = outdir + "/" + mesh_name_wo_extension;
  Results res;
  if (method == "dirichlet"){
    // test_Dirichlet_single(V, F, outdir, filename, std::cout);
    res = optimize_single_new(V, F, EnergyType::DIRICHLET, outdir, filename, start, queue);
  } else if (method == "arap"){
    // test_ARAP_single(V, F, outdir, filename, true, log);
    res = optimize_single_new(V, F, EnergyType::ARAP, outdir, filename, start, queue);
  } else if (method == "asap"){
    // test_ASAP_single(V, F, outdir, filename, true, log);
    res = optimize_single_new(V, F, EnergyType::ASAP, outdir, filename, start, queue);
  } else if (method == "symdirichlet"){
    // test_SymDirichlet_single(V, F, outdir, filename, log);
    res = optimize_single_new(V, F, EnergyType::SYMMETRIC_DIRICHLET, outdir, filename, start, queue);
  } else {
    std::cout << "Method '" << method << "' is unknown." << std::endl;
    return 1;
  }
  std::filesystem::create_directory(to_store_dir);

  log.open(to_store_dir+"/results_"+method+".log", std::ios::out);
  double sum = res.init_time;
  // print res
  log << filename << "," << method << "," << (start ? "delaunay_init" : "normal_init") << "," << (queue ? "priority_queue" : "normal_queue") << "\n";

  log << "init_time," << res.init_time << "\n";
  log << "energies";
  for(double i : res.energies)
    log << "," << std::setprecision(32) << i;
  log << "\n";

  log << "flips";
  for(int i : res.flips)
    log << "," << i;
  log << "\n";

  log << "flips_delaunay";
  for(int i : res.flips_delaunay)
    log << "," << i;
  log << "\n";

  sum = 0;
  log << "flip_durations";
  for(double i : res.flip_durations){
    sum += i;
    log << "," << i;
  }
  log << "\n";

  log << "opt_durations";

  for(double i : res.opt_durations){
    sum += i;
    log << "," << i;
  }
  log << "\n";

  log << "opt_iterations";
  for(int i : res.opt_iterations)
    log << "," << i;
  log << "\n";
  log << "total_duration_iparam," << sum << "\n";
  log << "total_duration_normal," << res.init_time+res.opt_durations.front() << "\n";
  log << "subdivision_time," << res.cs_time << "\n";
  log << "mesh_construction_time," << res.cm_time << "\n";
  log << "ext_vertices," << res.ext_vertex_count << "\n";
  log << "cs_vertices," << res.cs_vertex_count << "\n";
	return 0;
}

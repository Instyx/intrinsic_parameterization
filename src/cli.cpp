#include <string>
#include <ostream>
#include <iostream>
#include "read_mesh.hpp"
#include "test.hpp"
#include <filesystem>

int main(int argc, char *argv[]) {
  // evaluation
  if(argc != 4) {
    std::cout << "You have to specify the parameterization" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "./cli <method> <mesh> <output>" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "<method> in [dirichlet, arap, asap, symdirichlet]" << std::endl;
    std::cout << "<mesh>   as filepath" << std::endl;
    std::cout << "<output> as directory" << std::endl;
    return 0;
  }
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  DataGeo data_mesh;
  read_mesh(std::string(argv[2]),V,F);
  std::fstream log;
  std::string filename = std::string(argv[2]).substr(std::string(argv[2]).find_last_of("/\\") + 1);
  size_t lastindex = filename.find_last_of(".");
  std::string mesh_name_wo_extension =  filename.substr(0, lastindex);
  std::string to_store_dir = std::string(argv[3]) + "/" + mesh_name_wo_extension;

  if (std::string(argv[1]) == "dirichlet"){
    std::filesystem::create_directory(to_store_dir);
    log.open(to_store_dir+"/results_"+std::string(argv[1])+".log", std::ios::out);
    test_Dirichlet_single(V, F, std::string(argv[3]), filename, log);
  } else if (std::string(argv[1]) == "arap"){
    std::filesystem::create_directory(to_store_dir);
    log.open(to_store_dir+"/results_"+std::string(argv[1])+".log", std::ios::out);
    test_ARAP_single(V, F, std::string(argv[3]), filename, true, log);
  } else if (std::string(argv[1]) == "asap"){
    std::filesystem::create_directory(to_store_dir);
    log.open(to_store_dir+"/results_"+std::string(argv[1])+".log", std::ios::out);
    test_ASAP_single(V, F, std::string(argv[3]), filename, true, log);
  } else if (std::string(argv[1]) == "symdirichlet"){
    std::filesystem::create_directory(to_store_dir);
    log.open(to_store_dir+"/results_"+std::string(argv[1])+".log", std::ios::out);
    test_SymDirichlet_single(V, F, std::string(argv[3]), filename, log);
  } else {
    std::cout << "Method '" << argv[1] << "' is unknown." << std::endl;
  }

	return 0;
}

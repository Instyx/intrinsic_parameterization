#include <string>
#include <ostream>
#include <iostream>
#include "read_mesh.hpp"
#include "test.hpp"

int main(int argc, char *argv[]) {
  // evaluation
  if(argc != 4) {
    std::cout << "You have to specify the parameterization" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "./cli <method> <mesh> <output>" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "<method> in [Dirichlet, ARAP, Conformal, SymDirichlet]" << std::endl;
    std::cout << "<mesh>   as filepath" << std::endl;
    std::cout << "<output> as directory" << std::endl;
    return 0;
  }
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  DataGeo data_mesh;
  read_mesh(std::string(argv[2]),V,F);
  std::fstream log;

  size_t lastindex = mesh_name.find_last_of(".");
  std::string mesh_name_wo_extension =  mesh_name.substr(0, lastindex);
  std::string to_store_dir = dir + "/" + mesh_name_wo_extension;
  std::filesystem::create_directory(to_store_dir);

  // opens an existing csv file or creates a new file.
  log.open(to_store_dir+"/results_"+std::string(argv[1])+".log", std::ios::out);

  if (std::string(argv[1]) == "Dirichlet"){

  } else if (std::string(argv[1]) == "ARAP"){
    test_ARAP_single(V, F, std::string(argv[3]), std::string(argv[2]).substr(std::string(argv[2]).find_last_of("/\\") + 1), true, log);
  } else if (std::string(argv[1]) == "Conformal"){

  } else if (std::string(argv[1]) == "SymDirichlet"){

  } else {
    std::cout << "Method '" << argv[1] << "' is unknown." << std::endl;
  }

	return 0;
}

#include <string>
#include <ostream>
#include <iostream>

#include "test.hpp"

int main(int argc, char *argv[]) {
  // evaluation
  if(argc != 3) {
    std::cout << "You have to specify the parameterization" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "./cli <method> <mesh> <output>" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "<method> in [Dirichlet, ARAP, Conformal, SymDirichlet]" << std::endl;
    std::cout << "<mesh>   as filepath" << std::endl;
    std::cout << "<output> as filepath" << std::endl;
    return 0;
  }
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  DataGeo data_mesh;
  load_mesh_test(data_mesh, V, F, argv[2]);
  if (std::string(argv[1]) == "Dirichlet"){

  } else if (std::string(argv[1]) == "ARAP"){
    test_ARAP_single(data_mesh, std::string(argv[2]), true, std::cout);
  } else if (std::string(argv[1]) == "Conformal"){

  } else if (std::string(argv[1]) == "SymDirichlet"){

  } else {
    std::cout << "Method '" << argv[1] << "' is unknown." << std::endl;
  }

	return 0;
}

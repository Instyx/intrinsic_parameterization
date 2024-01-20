#include "read_mesh.hpp"
#include <igl/read_triangle_mesh.h>

void read_mesh(const std::string path, Eigen::MatrixXd& V, Eigen::MatrixXi& F){
  igl::read_triangle_mesh(path,V,F);
}

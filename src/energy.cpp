#pragma cling add_include_path("")
#pragma cling add_library_path("")
#include "energy.hpp"
#include <iostream>

inline Eigen::Matrix2d get_J(const Eigen::MatrixXd& UV, const double lens[3], const int indices[3]){
  const double len2sq = lens[2]*lens[2];
  Eigen::Matrix2d J1, J2inv;
  const double x = (lens[0]*lens[0] - lens[1]*lens[1] + len2sq)/(2*lens[0]);
  const double y = sqrt(len2sq - x*x);
  // J2 << len0, x, 0 , y;
  J2inv << 1/lens[0], -x/(lens[0]*y), 0 , 1/y;
  J1 << UV(indices[1],0) - UV(indices[0],0), UV(indices[2],0) - UV(indices[0],0),
        UV(indices[1],1) - UV(indices[0],1), UV(indices[2],1) - UV(indices[0],1);
  const Eigen::Matrix2d J = J1*J2inv;

  return J;
};

// void tri_wise_dblareas(const Eigen::MatrixXi& F, const Eigen::MatrixXd& V, Eigen::VectorXd& A){
//   A.resize(F.rows());
//   for (size_t i = 0; i < F.rows(); i++) {
//     const int indices[3] = { F(i,0), F(i,1), F(i,2)};
//     const Eigen::Vector3d v0 = V.row(indices[0]);
//     const Eigen::Vector3d v1 = V.row(indices[1])-v0.transpose();
//     const Eigen::Vector3d v2 = V.row(indices[2])-v0.transpose();
//     A(i) = v1.cross(v2).norm()/2;
//     // const double a = (v1-v0).norm();
//     // const double b = (v2-v1).norm();
//     // const double c = (v0-v2).norm();
//     // const double s = 0.5*(a+b+c);
//     // A(i) = sqrt(s*(s-a)*(s-b)*(s-c));
//   }
// }

void tri_wise_energy_int(DataGeo& mesh_data, const Eigen::MatrixXd& UV, double (*energy)(const Eigen::Matrix2d &), Eigen::VectorXd& E){
    double areassum = 0;
    E.resize( mesh_data.intTri->intrinsicMesh->nFaces());

    mesh_data.intTri->requireEdgeLengths();
    size_t i = 0;
    for(auto f : mesh_data.intTri->intrinsicMesh->faces()) {
      double lens[3];
      int indices[3];
      size_t k = 0;
      for(auto he : f.adjacentHalfedges()) {
        lens[k] = mesh_data.intTri->edgeLengths[he.edge()];
        indices[k] = he.tailVertex().getIndex();
        ++k;
      }
      const double area = mesh_data.intTri->faceArea(f);
      areassum += area;
      E(i++) = energy(get_J(UV, lens, indices))*area;
    }
    E /= areassum;
}

void tri_wise_energy_ext(DataGeo& mesh_data, const Eigen::MatrixXd& UV, double (*energy)(const Eigen::Matrix2d &), Eigen::VectorXd& E){
    double areassum = 0;
    const Eigen::MatrixXi& F = mesh_data.F;
    const Eigen::MatrixXd& V = mesh_data.V;
    E.resize(F.rows());

    for (size_t i = 0; i < F.rows(); i++) {
      const int indices[3] = { F(i,0), F(i,1), F(i,2)};
      const Eigen::Vector3d v0 = V.row(indices[0]);
      const Eigen::Vector3d v1 = V.row(indices[1]);
      const Eigen::Vector3d v2 = V.row(indices[2]);
      const double lens[3] = {(v1-v0).norm(), (v2-v1).norm(), (v0-v2).norm()};
      const double area = (v1-v0).cross((v2-v0)).norm()/2;

      areassum += area;
      E(i) = energy(get_J(UV, lens, indices))*area;
    }
    E /= areassum;
}

void tri_wise_energy(DataGeo& mesh_data, const Eigen::MatrixXd& UV, double (*energy)(const Eigen::Matrix2d &), const bool intrinsic, Eigen::VectorXd& E){
  if (intrinsic)
    return tri_wise_energy_int(mesh_data, UV, energy, E);
  return tri_wise_energy_ext(mesh_data, UV, energy, E);
}

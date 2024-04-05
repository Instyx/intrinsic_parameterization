#include <nanobind/nanobind.h>
#include <datageo.hpp>
#include <nanobind/eigen/dense.h>

namespace nb = nanobind;

DataGeo create_datastructure(const nb::DRef<Eigen::MatrixXd> &V, const nb::DRef<Eigen::MatrixXi> &F){
  DataGeo data_mesh;
  data_mesh.V = V;
  data_mesh.F = F;
  data_mesh.inputMesh.reset(new gcs::ManifoldSurfaceMesh(data_mesh.F));
  data_mesh.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh.inputMesh, data_mesh.V));
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->requireFaceAreas();
  return data_mesh;
}

void init_create_datastructure(nb::module_ &m) {
  m.def("create_datastructure", &create_datastructure);
}

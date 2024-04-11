#pragma once

#include <Eigen/Core>
#include <vector>
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

namespace gc = geometrycentral;
namespace gcs = gc::surface;



// data structure for holding intrinsic triangulation and input mesh 

struct DataGeo{
  // input
  Eigen::Matrix<double, -1, 3> V;
  Eigen::Matrix<int, -1, 3> F;
  
  // intrinsic triangulation
  std::unique_ptr<gcs::SignpostIntrinsicTriangulation> intTri;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> inputMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> inputGeometry;

};

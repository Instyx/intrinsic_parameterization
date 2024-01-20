#include "store_intrinsic_mesh.hpp"
#include "happly.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/rich_surface_mesh_data.h"
#include <iomanip>

namespace gc = geometrycentral;

Eigen::VectorXd interpolate(const gcs::SurfacePoint& pt, const Eigen::MatrixXd& V){
  if (pt.type == gcs::SurfacePointType::Vertex){
    return V.row(pt.vertex.getIndex());
  } else {
    return V.row(pt.edge.firstVertex().getIndex())*(1-pt.tEdge) + V.row(pt.edge.secondVertex().getIndex())*pt.tEdge;
  }
}

void store_intrinsic_edges(const DataGeo &data_mesh, std::string filepath) {
  std::ofstream out(filepath+".int.edges");
  std::vector<Eigen::VectorXd> points;
  std::vector<std::array<int, 2>> edges, og_edges;
  edges.clear();
  points.clear();
  for(gcs::Edge e : data_mesh.intTri->intrinsicMesh->edges()) {
    std::vector<gcs::SurfacePoint> pointVec = data_mesh.intTri->traceIntrinsicHalfedgeAlongInput(e.halfedge());

    for (int k = 0; k < pointVec.size(); k++) {
      points.push_back(interpolate(pointVec[k], data_mesh.V));
      // auto p = interpolate(pointVec[k], data_mesh.V);
      if (k == 0) continue;
      int t = points.size();
      std::array<int, 2> tmp{t-2, t-1};
      if (pointVec.size() == 2) {
        og_edges.push_back(tmp);
      }else{
        edges.push_back(tmp);
      }

    }
  }
  out << points.size() << " " << og_edges.size() << " " << edges.size() << std::endl;
  for(auto p : points)
    out << std::fixed << std::setprecision(32) << p[0] << "," << p[1] << "," << p[2] << std::endl;
  for(auto e : og_edges)
    out << e[0] << "," << e[1] << std::endl;
  out << std::endl;
  for(auto e : edges)
    out << e[0] << "," << e[1] << std::endl;
}


void store_intrinsic_mesh(DataGeo &data_mesh, std::string filename){
  gcs::CommonSubdivision& cs = data_mesh.intTri->getCommonSubdivision();
  cs.constructMesh();

  gcs::FaceData<double> faceIDs = niceColors(*data_mesh.intTri->intrinsicMesh, 7);
  auto data = cs.copyFromB(faceIDs);
  // auto data = cs.copyFromA(faceIDs);
  gcs::VertexData<gc::Vector3> csPositions = cs.interpolateAcrossA(data_mesh.inputGeometry->vertexPositions);
  // separate x/y/z coordinates
  gcs::VertexData<double> x(*cs.mesh);
  gcs::VertexData<double> y(*cs.mesh);
  gcs::VertexData<double> z(*cs.mesh);

  for (auto v : cs.mesh->vertices()) {
    gc::Vector3 p = csPositions[v];
    x[v] = p.x;
    y[v] = p.y;
    z[v] = p.z;
  }
  gcs::RichSurfaceMeshData richData(*cs.mesh);
  richData.outputFormat = happly::DataFormat::ASCII;
  richData.addMeshConnectivity();
  richData.addVertexProperty("x", x);
  richData.addVertexProperty("y", y);
  richData.addVertexProperty("z", z);
  richData.addFaceProperty("OGFaceIds", data);
  richData.write(filename+".int.ply");
}

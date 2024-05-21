#include "store_intrinsic_mesh.hpp"
#include "happly.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/rich_surface_mesh_data.h"
#include <iomanip>
#include <chrono>

namespace gc = geometrycentral;

Eigen::VectorXd interpolate(const gcs::SurfacePoint& pt, const Eigen::MatrixXd& V){
  if (pt.type == gcs::SurfacePointType::Vertex){
    return V.row(pt.vertex.getIndex());
  } else {
    return V.row(pt.edge.firstVertex().getIndex())*(1-pt.tEdge) + V.row(pt.edge.secondVertex().getIndex())*pt.tEdge;
  }
}

void store_intrinsic_edges(const DataGeo &data_mesh, const std::string filepath) {
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

void store_intrinsic_mesh(const DataGeo &data_mesh, const Eigen::MatrixXd &UV, const std::string filename){
  Results res;
  return store_intrinsic_mesh(data_mesh, UV, filename, res);
}

void store_intrinsic_mesh(const DataGeo &data_mesh, const Eigen::MatrixXd &UV, const std::string filename, Results &res){
  const int n = data_mesh.intTri->intrinsicMesh->nFaces();
  const Eigen::VectorXd afterTriEnergy(n), beforeTriEnergy(n);
  return store_intrinsic_mesh(data_mesh, UV, filename, beforeTriEnergy, afterTriEnergy, res);
}

void store_intrinsic_mesh(const DataGeo &data_mesh, const Eigen::MatrixXd &UV, const std::string filename, const Eigen::VectorXd& beforeTriEnergy, const Eigen::VectorXd& afterTriEnergy){
  Results res;
  return store_intrinsic_mesh(data_mesh, UV, filename, beforeTriEnergy, afterTriEnergy, res);
}

void store_intrinsic_mesh(const DataGeo &data_mesh, const Eigen::MatrixXd &UV, const std::string filename, const Eigen::VectorXd& beforeTriEnergy, const Eigen::VectorXd& afterTriEnergy, Results &res){
  data_mesh.intTri->requireFaceAreas();
  data_mesh.intTri->unrequireFaceAreas();

  auto start = std::chrono::high_resolution_clock::now();
  gcs::CommonSubdivision& cs = data_mesh.intTri->getCommonSubdivision();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  res.cs_time = duration;
  start = std::chrono::high_resolution_clock::now();
  cs.constructMesh();
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  res.cm_time = duration;
  res.ext_vertex_count = data_mesh.intTri->intrinsicMesh->nVertices();
  res.cs_vertex_count = cs.mesh->nVertices();
  // return;

  gcs::FaceData<double> colorIDs = niceColors(*data_mesh.intTri->intrinsicMesh, 7);
  gcs::FaceData<double> beforeTriEnergyData = gcs::FaceData<double>(*data_mesh.inputMesh, beforeTriEnergy);
  gcs::FaceData<double> afterTriEnergyData = gcs::FaceData<double>(*data_mesh.intTri->intrinsicMesh, afterTriEnergy);
  gcs::FaceData<double> extrId(*cs.mesh), intrId(*cs.mesh);
  gcs::FaceData<double> intrArea = data_mesh.intTri->faceAreas;
  auto color = cs.copyFromB(colorIDs);
  auto extr = cs.copyFromA(beforeTriEnergyData);
  auto intr = cs.copyFromB(afterTriEnergyData);
  auto area = cs.copyFromB(intrArea);
  for (auto f : cs.mesh->faces()) {
    extrId[f] = cs.sourceFaceA[f].getIndex();
    intrId[f] = cs.sourceFaceB[f].getIndex();
  }
  // auto data = cs.copyFromA(colorIDs);
  gcs::VertexData<gc::Vector3> csPositions = cs.interpolateAcrossA(data_mesh.inputGeometry->vertexPositions);
  gcs::VertexData<double> u(*data_mesh.intTri->intrinsicMesh);
  gcs::VertexData<double> v(*data_mesh.intTri->intrinsicMesh);
  for (auto ver : data_mesh.intTri->intrinsicMesh->vertices()) {
    u[ver] = UV(ver.getIndex(),0);
    v[ver] = UV(ver.getIndex(),1);
  }
  gcs::VertexData<double> uPositions = cs.interpolateAcrossB(u);
  gcs::VertexData<double> vPositions = cs.interpolateAcrossB(v);
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

  std::vector<std::vector<size_t>> faceIndices = cs.mesh->getFaceVertexList();
  richData.plyData.addFaceIndices(faceIndices);

  richData.addVertexProperty("x", x);
  richData.addVertexProperty("y", y);
  richData.addVertexProperty("z", z);
  richData.addVertexProperty("u", uPositions);
  richData.addVertexProperty("v", vPositions);
  richData.addFaceProperty("ExtrID", extrId);
  richData.addFaceProperty("IntrColor", color);
  richData.addFaceProperty("IntrID", intrId);
  richData.addFaceProperty("IntrArea", area);
  richData.addFaceProperty("ExtrEnergy", extr);
  richData.addFaceProperty("IntrEnergy", intr);
  richData.write(filename+".int.ply");
}

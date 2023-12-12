#include "intrinsicflip.hpp"
#include "distortion_energy.hpp"
#include <random>
#include <algorithm>
#include <math.h>

double cross2d(Eigen::Vector2d &v1, Eigen::Vector2d &v2){
  return v1(0)*v2(1)-v1(1)*v2(0);
}

// the points are in order with the polygon rotation
bool isConcave(std::vector<Eigen::Vector2d>& points){
  Eigen::Vector2d e1 = points[1]-points[0];
  Eigen::Vector2d e2 = points[2]-points[1];
  Eigen::Vector2d e3 = points[3]-points[2];
  Eigen::Vector2d e4 = points[0]-points[3];
  
  double cross1 = cross2d(e1,e2);
  double cross2 = cross2d(e2,e3);
  double cross3 = cross2d(e3,e4);
  double cross4 = cross2d(e4,e1);

  if (cross1 * cross2 < 0 || cross2 * cross3 < 0 || cross3 * cross4 < 0) {
    return true;
  }
    
  return false;
}

bool diamondJacobians(DataGeo &data_mesh, const Eigen::MatrixXd &UV, gcs::Edge &e, Eigen::Matrix2d &J1, Eigen::Matrix2d &J2){
  std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
  double fi_len1 = data_mesh.intTri->edgeLengths[e];
  double fi_len2 = data_mesh.intTri->edgeLengths[halfedges[0].edge()];
  double fi_len3 = data_mesh.intTri->edgeLengths[halfedges[1].edge()];

  double se_len1 = data_mesh.intTri->edgeLengths[e];
  double se_len2 = data_mesh.intTri->edgeLengths[halfedges[2].edge()];
  double se_len3 = data_mesh.intTri->edgeLengths[halfedges[3].edge()];

  Eigen::Matrix2d E1, E2, E1_tilde, E2_tilde;
  double temp = (fi_len2*fi_len2 - fi_len1*fi_len1 - fi_len3*fi_len3)/(-2*fi_len1); 
  E1_tilde << fi_len1, temp, 0 , sqrt(fi_len3*fi_len3 - temp*temp); 
  temp = (se_len2*se_len2 - se_len1*se_len1 - se_len3*se_len3)/(-2*se_len1); 
  E2_tilde << se_len1, temp, 0 , sqrt(se_len3*se_len3 - temp*temp);
 
  data_mesh.inputGeometry->requireVertexPositions();

  //        v3 /\
  //          /  \
  //         /    \
  //     v1 / _ e _\ v2
  //        \     /
  //         \   /  
  //          \ /
  //          v4
  size_t v1 = data_mesh.intTri->vertexIndices[halfedges[1].tipVertex()];
  size_t v2 = data_mesh.intTri->vertexIndices[halfedges[0].tailVertex()];
  size_t v3 = data_mesh.intTri->vertexIndices[halfedges[0].tipVertex()];
  size_t v4 = data_mesh.intTri->vertexIndices[halfedges[2].tipVertex()];
 
  std::vector<Eigen::Vector2d> points(4);
  points[0] = UV.row(v1).transpose();
  points[1] = UV.row(v4).transpose();
  points[2] = UV.row(v2).transpose();
  points[3] = UV.row(v3).transpose();
  //cout << "intri:   "<< data_mesh.inputGeometry->vertexPositions[v1] << ";  start mesh" << V.row(v1) << endl;

  // if concave the intrinsic flip is not possible
  if(isConcave(points)) return false;

  E1 << UV(v2,0) - UV(v1,0), UV(v3,0) - UV(v1,0),
        UV(v2,1) - UV(v1,1), UV(v3,1) - UV(v1,1);
  E2 << UV(v1,0) - UV(v2,0), UV(v4,0) - UV(v2,0),
        UV(v1,1) - UV(v2,1), UV(v4,1) - UV(v2,1);
  J1 = E1 * E1_tilde.inverse();
  J2 = E2 * E2_tilde.inverse();
  
  return true;
}

double flippeddiff(DataGeo &data_mesh, const Eigen::MatrixXd &UV, gcs::Edge e, const EnergyType &et){
  if(e.isBoundary()) return 0;

  auto energy = dirichlet;
  if(et==EnergyType::DIRICHLET) energy = dirichlet;
  if(et==EnergyType::ASAP) energy = asap;
  if(et==EnergyType::ARAP) energy = arap;
  if(et==EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;
  
  gcs::Face f1 = e.halfedge().face(); 
  gcs::Face f2 = e.halfedge().twin().face();
  Eigen::Matrix2d J1, J2, J1_prime, J2_prime;

  // if flip is not possible return 0
  if(!diamondJacobians(data_mesh, UV, e, J1, J2)){
    return 0;
  }
  double before = energy(J1) * data_mesh.intTri->faceArea(f1) +
                  energy(J2) * data_mesh.intTri->faceArea(f2);
  data_mesh.intTri->flipEdgeIfPossible(e);
  gcs::Edge flipped = e;

  diamondJacobians(data_mesh, UV, flipped, J1_prime, J2_prime);

  double after = energy(J1_prime) * data_mesh.intTri->faceArea(flipped.halfedge().face()) + energy(J2_prime) * data_mesh.intTri->faceArea(flipped.halfedge().twin().face());

  double tolerance = 1e-6; // set tolerance to 1e-6 
  data_mesh.intTri->flipEdgeIfPossible(flipped);
  if (fabs(before - after) / std::max(fabs(before), fabs(after)) > tolerance) {
    return after - before;
  }
  else {
    return 0; 
  }
}

unsigned greedy_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et){
  
  auto energy = dirichlet;
  if(et==EnergyType::DIRICHLET) energy = dirichlet;
  if(et==EnergyType::ASAP) energy = asap;
  if(et==EnergyType::ARAP) energy = arap;
  if(et==EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;

  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceAreas();

  unsigned totalflips=0;
  gcs::EdgeData<double> diffs(* (data_mesh.intTri->intrinsicMesh));
  std::vector<size_t> indices(diffs.size());
  std::vector<int> visited(diffs.size());
  int i=0;
  for(gcs::Edge e : data_mesh.intTri->intrinsicMesh->edges()){
    double energydiff = flippeddiff(data_mesh, UV, e, et);
    diffs[e] = energydiff;
    indices[i]= e.getIndex();
    ++i;
  }

  unsigned delaunay = 0;
  double tolerance = 1e-6; // set tolerance to 1e-6 
  for (size_t i = 0; i < indices.size(); i++) {
    size_t idx = indices[i];
    if(diffs[idx]>=0   || visited[idx]) continue;
    gcs::Edge e = data_mesh.intTri->intrinsicMesh->edge(idx);
    if(e.isBoundary()) continue;
    data_mesh.intTri->flipEdgeIfPossible(e);
    if(data_mesh.intTri->isDelaunay(e)) delaunay++;
    std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
    visited[idx]=1;
    visited[halfedges[0].edge().getIndex()]=1;
    visited[halfedges[1].edge().getIndex()]=1;
    visited[halfedges[2].edge().getIndex()]=1;
    visited[halfedges[3].edge().getIndex()]=1;
    ++totalflips;
  }

  data_mesh.intTri->refreshQuantities();
  std::cout << " delaunay flips: " << delaunay << std::endl;
  return totalflips;
}

unsigned heuristic_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et){
  
  auto energy = dirichlet;
  if(et==EnergyType::DIRICHLET) energy = dirichlet;
  if(et==EnergyType::ASAP) energy = asap;
  if(et==EnergyType::ARAP) energy = arap;
  if(et==EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;


  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceAreas();

  unsigned totalflips=0;
  gcs::EdgeData<double> diffs(* (data_mesh.intTri->intrinsicMesh));
  std::vector<size_t> indices(diffs.size());
  std::vector<int> visited(diffs.size());
  int i=0;
  for(gcs::Edge e : data_mesh.intTri->intrinsicMesh->edges()){
    double energydiff = flippeddiff(data_mesh, UV, e, et);
    diffs[e] = energydiff;
    indices[i]= e.getIndex();
    ++i;
  }
  gcs::EdgeData<double> heuristic(* (data_mesh.intTri->intrinsicMesh));
  
  for(gcs::Edge e : data_mesh.intTri->intrinsicMesh->edges()){
    if(e.isBoundary()){
      heuristic[e] = 0;
      continue;
    }
    std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
    heuristic[e] = diffs[e] - (diffs[halfedges[0].edge()]+diffs[halfedges[1].edge()]+diffs[halfedges[2].edge()]+diffs[halfedges[3].edge()]);
  }  
  unsigned delaunay = 0;
  for (size_t i = 0; i < indices.size(); i++) {
    size_t idx = indices[i];
    if(diffs[idx]>=0) continue;
    gcs::Edge e = data_mesh.intTri->intrinsicMesh->edge(idx);
    if(visited[idx] || e.isBoundary()) continue;
    std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
    data_mesh.intTri->flipEdgeIfPossible(e);
    visited[idx]=1;
    if(data_mesh.intTri->isDelaunay(e)) delaunay++;
    visited[halfedges[0].edge().getIndex()]=1;
    visited[halfedges[1].edge().getIndex()]=1;
    visited[halfedges[2].edge().getIndex()]=1;
    visited[halfedges[3].edge().getIndex()]=1;
    ++totalflips;
  }

  data_mesh.intTri->refreshQuantities();
  std::cout << " delaunay flips: " << delaunay << std::endl;
  return totalflips;
}

unsigned random_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et){
  
  auto energy = dirichlet;
  if(et==EnergyType::DIRICHLET) energy = dirichlet;
  if(et==EnergyType::ASAP) energy = asap;
  if(et==EnergyType::ARAP) energy = arap;
  if(et==EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;

  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceAreas();
  unsigned totalflips = 0;
  unsigned total_concaves =0;

  std::vector<size_t> indices(data_mesh.intTri->intrinsicMesh->nEdges());
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), std::mt19937 {std::random_device{}()});
  unsigned delaunay = 0;
  for(size_t idx : indices){
    gcs::Edge e = data_mesh.intTri->intrinsicMesh->edge(idx);
    if(e.isBoundary()) continue;
    gcs::Face f1 = e.halfedge().face(); 
    gcs::Face f2 = e.halfedge().twin().face();
    Eigen::Matrix2d J1, J2, J1_prime, J2_prime;
    if(!diamondJacobians(data_mesh, UV, e, J1, J2)) continue;

    double before = energy(J1) * data_mesh.intTri->faceArea(f1) +
                    energy(J2) * data_mesh.intTri->faceArea(f2);
    data_mesh.intTri->flipEdgeIfPossible(e);
    gcs::Edge flipped = e;
    diamondJacobians(data_mesh, UV, flipped, J1_prime, J2_prime);
    
    double after = energy(J1_prime) * data_mesh.intTri->faceArea(flipped.halfedge().face()) + energy(J2_prime) * data_mesh.intTri->faceArea(flipped.halfedge().twin().face());
    double tolerance = 1e-6; // set tolerance to 1e-6 

    if (fabs(before - after) / std::max(fabs(before), fabs(after)) > tolerance) {
      if (before > after) {
        totalflips++;
      if(data_mesh.intTri->isDelaunay(flipped)) delaunay++;
      }
      else {
        data_mesh.intTri->flipEdgeIfPossible(flipped);
      }
    }
    else {
      data_mesh.intTri->flipEdgeIfPossible(flipped);
    }
  }
  data_mesh.intTri->refreshQuantities();
  std::cout << " delaunay flips: " << delaunay << std::endl;
  return totalflips;
}

unsigned edgeorder_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et){
  
  auto energy = dirichlet;
  if(et==EnergyType::DIRICHLET) energy = dirichlet;
  if(et==EnergyType::ASAP) energy = asap;
  if(et==EnergyType::ARAP) energy = arap;
  if(et==EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;

  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceAreas();
  unsigned totalflips = 0;
  unsigned delaunay = 0;
  for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
    if(e.isBoundary()) continue;
    gcs::Face f1 = e.halfedge().face(); 
    gcs::Face f2 = e.halfedge().twin().face();
    Eigen::Matrix2d J1, J2, J1_prime, J2_prime;
    if(!diamondJacobians(data_mesh, UV, e, J1, J2)) continue;
    double before = energy(J1) * data_mesh.intTri->faceArea(f1) +
                    energy(J2) * data_mesh.intTri->faceArea(f2);
    data_mesh.intTri->flipEdgeIfPossible(e);
    gcs::Edge flipped = e;
    diamondJacobians(data_mesh, UV, flipped, J1_prime, J2_prime);
    
    double after = energy(J1_prime) * data_mesh.intTri->faceArea(flipped.halfedge().face()) + energy(J2_prime) * data_mesh.intTri->faceArea(flipped.halfedge().twin().face());
    double tolerance = 1e-6; // set tolerance to 1e-6 

    if (fabs(before - after) / std::max(fabs(before), fabs(after)) > tolerance) {
      if (before > after) {
        totalflips++;
        if(data_mesh.intTri->isDelaunay(e)) delaunay++;
      }
      else {
        data_mesh.intTri->flipEdgeIfPossible(flipped);
      }
    }
    else {
      data_mesh.intTri->flipEdgeIfPossible(flipped);
    }

  }
  data_mesh.intTri->refreshQuantities();
  std::cout << " delaunay flips: " << delaunay << std::endl;
  return totalflips;
}


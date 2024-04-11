#include "intrinsicflip.hpp"
#include "distortion_energy.hpp"
#include "parameterization.hpp"
#include <cmath>
#include <ostream>
#include <random>
#include <algorithm>
#include <math.h>
#include <queue>
#include <tuple>

extern "C" {
void exactinit(void);
double orient2d(const double *pa, const double *pb, const double *pc);
}
// helper to call Shewchuk predicate with eigen types transparently
double orient2d(const Eigen::Vector2d &pa,
                const Eigen::Vector2d &pb,
                const Eigen::Vector2d &pc)
{
    return orient2d(pa.data(),pb.data(),pc.data());
}

double cross2d(Eigen::Vector2d &v1, Eigen::Vector2d &v2){
  return v1(0)*v2(1)-v1(1)*v2(0);
}


bool isConcave_robust(std::vector<Eigen::Vector2d>& points){
  Eigen::Vector2d p1 = points[0];
  Eigen::Vector2d p2 = points[1];
  Eigen::Vector2d p3 = points[2];
  Eigen::Vector2d p4 = points[3];
  double v1 = orient2d(p1, p2, p3);
  double v2 = orient2d(p2, p3, p4);
  double v3 = orient2d(p3, p4, p1);
  double v4 = orient2d(p4, p1, p2);

  if (v1 * v2 < 0 || v2 * v3 < 0 || v3 * v4 < 0) {
    return true;
  }
  return false;
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


bool isFlipPossible(DataGeo &data_mesh, const Eigen::MatrixXd &UV, gcs::Edge &e){
  std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
  
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


  // if flip causes self loop, don't flip
  if(v3==v4) return false;

  // if one of the end vertices has degree smaller than 3, don't flip
  // if the vertex is in the boundary, check if smaller than 2
  if(halfedges[1].tipVertex().isBoundary()){
    if(halfedges[1].tipVertex().degree()<=2) return false;
  }
  else if(halfedges[1].tipVertex().degree()<=3) return false;

  if(halfedges[0].tailVertex().isBoundary()){
    if(halfedges[0].tailVertex().degree()<=2) return false;
  }
  else if(halfedges[0].tailVertex().degree()<=3) return false;

  // if concave the intrinsic flip is not possible
  if(isConcave_robust(points)) return false;
  // if(isConcave(points)) return false;

  return true;

}

// DO NOT USE
bool diamondJacobians(DataGeo &data_mesh, const Eigen::MatrixXd &UV, gcs::Edge &e, Eigen::Matrix2d &J1, Eigen::Matrix2d &J2){
  std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
  
  size_t v1 = data_mesh.intTri->vertexIndices[halfedges[1].tipVertex()];
  size_t v2 = data_mesh.intTri->vertexIndices[halfedges[0].tailVertex()];
  size_t v3 = data_mesh.intTri->vertexIndices[halfedges[0].tipVertex()];
  size_t v4 = data_mesh.intTri->vertexIndices[halfedges[2].tipVertex()];

  std::vector<Eigen::Vector2d> points(4);
  points[0] = UV.row(v1).transpose();
  points[1] = UV.row(v4).transpose();
  points[2] = UV.row(v2).transpose();
  points[3] = UV.row(v3).transpose();


  // if flip causes self loop, don't flip
  if(v3==v4) return false;
  // if one of the end vertices has degree smaller than 3, don't flip
  if(halfedges[1].tipVertex().degree()<=3 || halfedges[0].tailVertex().degree()<=3) return false;
  // if concave the intrinsic flip is not possible
  if(isConcave_robust(points)) return false;
  // if(isConcave(points)) return false;

  faceJacobian(data_mesh, UV, e.halfedge().face(), J1);
  faceJacobian(data_mesh, UV, e.halfedge().twin().face(), J2);
  return true;
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

 // data_mesh.inputGeometry->requireVertexPositions();

  //        v3 /\
  //          /  \
  //         /    \
  //     v1 / _ e _\ v2
  //        \     /
  //         \   /
  //          \ /
  //          v4
  
    //cout << "intri:   "<< data_mesh.inputGeometry->vertexPositions[v1] << ";  start mesh" << V.row(v1) << endl;


  
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
  if(!data_mesh.intTri->flipEdgeIfPossible(e)) return 0;
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

unsigned greedy_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et){

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

  std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
    return diffs[i] < diffs[j];
    });
  delaunay_flips = 0;
  double tolerance = 1e-6; // set tolerance to 1e-6
  for (size_t i = 0; i < indices.size(); i++) {
    size_t idx = indices[i];
    if(diffs[idx]>=0   || visited[idx]) continue;
    gcs::Edge e = data_mesh.intTri->intrinsicMesh->edge(idx);
    if(e.isBoundary()) continue;
    data_mesh.intTri->flipEdgeIfPossible(e);
    if(data_mesh.intTri->isDelaunay(e)) delaunay_flips++;
    std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
    visited[idx]=1;
    visited[halfedges[0].edge().getIndex()]=1;
    visited[halfedges[1].edge().getIndex()]=1;
    visited[halfedges[2].edge().getIndex()]=1;
    visited[halfedges[3].edge().getIndex()]=1;
    ++totalflips;
  }

  data_mesh.intTri->refreshQuantities();
  return totalflips;
}

unsigned heuristic_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et){

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
  delaunay_flips = 0;
  for (size_t i = 0; i < indices.size(); i++) {
    size_t idx = indices[i];
    if(diffs[idx]>=0) continue;
    gcs::Edge e = data_mesh.intTri->intrinsicMesh->edge(idx);
    if(visited[idx] || e.isBoundary()) continue;
    std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
    data_mesh.intTri->flipEdgeIfPossible(e);
    visited[idx]=1;
    if(data_mesh.intTri->isDelaunay(e)) delaunay_flips++;
    visited[halfedges[0].edge().getIndex()]=1;
    visited[halfedges[1].edge().getIndex()]=1;
    visited[halfedges[2].edge().getIndex()]=1;
    visited[halfedges[3].edge().getIndex()]=1;
    ++totalflips;
  }

  data_mesh.intTri->refreshQuantities();
  return totalflips;
}

unsigned random_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et){

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
  delaunay_flips = 0;
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
      if(data_mesh.intTri->isDelaunay(flipped)) delaunay_flips;
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
  return totalflips;
}


namespace std{
    namespace
    {
        // Code from boost
        // Reciprocal of the golden ratio helps spread entropy
        //     and handles duplicates.
        // See Mike Seymour in magic-numbers-in-boosthash-combine:
        //     http://stackoverflow.com/questions/4948780
        template <class T>
        inline void hash_combine(std::size_t& seed, T const& v)
        {
            seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }

        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
            hash_combine(seed, std::get<Index>(tuple));
          }
        };

        template <class Tuple>
        struct HashValueImpl<Tuple,0>
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            hash_combine(seed, std::get<0>(tuple));
          }
        };
    }

    template <typename ... TT>
    struct hash<std::tuple<TT...>>
    {
        size_t
        operator()(std::tuple<TT...> const& tt) const
        {
            size_t seed = 0;
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
            return seed;
        }

    };
}

unsigned queue_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et){
  auto energy = dirichlet;
  if(et==EnergyType::DIRICHLET) energy = dirichlet;
  else if(et==EnergyType::ASAP) energy = asap;
  else if(et==EnergyType::ARAP) energy = arap;
  else if(et==EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->unrequireEdgeLengths();
  unsigned totalflips = 0;
  delaunay_flips = 0;
  std::unordered_map<std::tuple<int, int, int, int>, bool> checked_diamonds;
  std::queue<gcs::Edge> q;
  for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
    q.push(e);
  }
  while (!q.empty()){
    gcs::Edge e = q.front();q.pop();
    if(e.isBoundary()) continue;
    gcs::Face f1 = e.halfedge().face();
    gcs::Face f2 = e.halfedge().twin().face();
    std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
    std::array<size_t, 4> quad = {halfedges[0].vertex().getIndex(), halfedges[1].vertex().getIndex(), halfedges[2].vertex().getIndex(), halfedges[3].vertex().getIndex()};
    std::sort(quad.begin(), quad.end());
    if (checked_diamonds.find(std::make_tuple(quad[0],quad[1],quad[2],quad[3])) != checked_diamonds.end()) continue;
    else checked_diamonds[std::make_tuple(quad[0],quad[1],quad[2],quad[3])] = true;
    Eigen::Matrix2d J1, J2, J1_prime, J2_prime;
    if(!isFlipPossible(data_mesh, UV, e)) continue;
    faceJacobian(data_mesh, UV, f1, J1);
    faceJacobian(data_mesh, UV, f2, J2);
    double before = energy(J1) * data_mesh.intTri->faceArea(f1) +
                    energy(J2) * data_mesh.intTri->faceArea(f2);
    // store edge quantaties before flipping
    double length = data_mesh.intTri->edgeLengths[e];
    double forwardangle = data_mesh.intTri->signpostAngle[e.halfedge()];
    double reverseangle = data_mesh.intTri->signpostAngle[e.halfedge().twin()];
    bool isOrig = data_mesh.intTri->edgeIsOriginal[e];

    //std::cout << "areas before: " << data_mesh.intTri->faceArea(f1) << "  " << data_mesh.intTri->faceArea(f2) << std::endl;
    //std::cout << "edge before: " << data_mesh.intTri->edgeLengths[f1.halfedge().edge()] << "  " << data_mesh.intTri->edgeLengths[f1.halfedge().next().edge()] << "  " <<data_mesh.intTri->edgeLengths[f1.halfedge().next().next().edge()] << std::endl;
    if (!data_mesh.intTri->flipEdgeIfPossible(e)) continue;
    data_mesh.intTri->refreshQuantities();
    gcs::Edge flipped = e;
    f1 = flipped.halfedge().face();
    f2 = flipped.halfedge().twin().face();
    faceJacobian(data_mesh, UV, f1, J1_prime);
    faceJacobian(data_mesh, UV, f2, J2_prime);

    double after = energy(J1_prime) * data_mesh.intTri->faceArea(f1) +
                   energy(J2_prime) * data_mesh.intTri->faceArea(f2);

    double tolerance = 1e-4;
    double abs_tolerance = 1e-10;
    //std::cout << "global: " << before-after << std::endl; 
    //std::cout << "local: " << before - after << std::endl; 
    // std::cout << "+0.3 local: " << (before_local+0.3) - (after_local+0.3) << std::endl << std::endl;

    double diff = fabs(before-after);

    if (before > after && (diff / std::max(fabs(before), fabs(after)) > tolerance && diff>abs_tolerance )) {
      //std::cout << "diff " << diff << std::endl;
      //std::cout << "global diff: " << before_global - after_global << std::endl << std::endl;
      totalflips++;
      if(data_mesh.intTri->isDelaunay(e)) delaunay_flips++;
      for (size_t i = 0; i < 4; i++) {
        q.push(halfedges[i].edge());
      }
      //std::cout << "FLIPPED" << std::endl;
    } else {
      //data_mesh.intTri->flipEdgeIfPossible(flipped);
      data_mesh.intTri->flipEdgeManual(flipped, length, forwardangle, reverseangle, isOrig, true);
      data_mesh.intTri->refreshQuantities();
     // std::cout << "areas before: " << data_mesh.intTri->faceArea(flipped.halfedge().face()) << "  " << data_mesh.intTri->faceArea(flipped.halfedge().twin().face()) << std::endl;
      //std::cout << "edge before: " << data_mesh.intTri->edgeLengths[f1.halfedge().edge()] << "  " << data_mesh.intTri->edgeLengths[f1.halfedge().next().edge()] << "  " <<data_mesh.intTri->edgeLengths[f1.halfedge().next().next().edge()] << std::endl;
      //double en = compute_total_energy_localjacob(data_mesh, UV, et);
      //std::cout << "NOT FLIPPED: " << en << std::endl;
    }
  }
  return totalflips;
}

unsigned priority_queue_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et){
  auto energy = dirichlet;
  if(et==EnergyType::DIRICHLET) energy = dirichlet;
  else if(et==EnergyType::ASAP) energy = asap;
  else if(et==EnergyType::ARAP) energy = arap;
  else if(et==EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->unrequireEdgeLengths();

  unsigned totalflips = 0;
  delaunay_flips = 0;
  std::unordered_map<std::tuple<int, int, int, int>, bool> checked_diamonds;
  std::priority_queue<std::pair<double, gcs::Edge> > q;
  std::vector<double> actual_energies(data_mesh.intTri->intrinsicMesh->nEdges(),  0);
  for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
    double energy_diff =  flippeddiff(data_mesh, UV, e, et);
    if(energy_diff<0)
      q.push(std::make_pair(-1*energy_diff, e));
    actual_energies[e.getIndex()] = energy_diff;
  }
  while (!q.empty()){
    gcs::Edge e = q.top().second;
    double energy_diff = -1*q.top().first;
    if(energy_diff>=0) break; // there is no decreasing flip
    q.pop();
    if(e.isBoundary()) continue;
    if (energy_diff!=actual_energies[e.getIndex()]) continue;
    std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
    std::array<size_t, 4> quad = {halfedges[0].vertex().getIndex(), halfedges[1].vertex().getIndex(), halfedges[2].vertex().getIndex(), halfedges[3].vertex().getIndex()};
    std::sort(quad.begin(), quad.end());
    if (checked_diamonds.find(std::make_tuple(quad[0],quad[1],quad[2],quad[3])) != checked_diamonds.end()) continue;
    else checked_diamonds[std::make_tuple(quad[0],quad[1],quad[2],quad[3])] = true;

    data_mesh.intTri->flipEdgeIfPossible(e);
    totalflips++;
    if(data_mesh.intTri->isDelaunay(e)) delaunay_flips++;
    for (size_t i = 0; i < 4; i++) {
      gcs::Edge next_e = halfedges[i].edge();
      double diff = flippeddiff(data_mesh, UV, next_e, et);
      if(diff<0)
        q.push(std::make_pair(-1*diff, next_e));
      actual_energies[next_e.getIndex()] = diff;
    }
  }
  return totalflips;
}


unsigned edgeorder_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, unsigned &delaunay_flips, const EnergyType &et){

  auto energy = dirichlet;
  if(et==EnergyType::DIRICHLET) energy = dirichlet;
  if(et==EnergyType::ASAP) energy = asap;
  if(et==EnergyType::ARAP) energy = arap;
  if(et==EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;

  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceAreas();
  unsigned totalflips = 0;
  delaunay_flips = 0;
  for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
    if(e.isBoundary()) continue;
    gcs::Face f1 = e.halfedge().face();
    gcs::Face f2 = e.halfedge().twin().face();
    Eigen::Matrix2d J1, J2, J1_prime, J2_prime;
    if(!diamondJacobians(data_mesh, UV, e, J1, J2)) continue;
    //double before = energy(J1) * data_mesh.intTri->faceArea(f1) +
    //                energy(J2) * data_mesh.intTri->faceArea(f2);

    double before = compute_total_energy_localjacob(data_mesh, UV, et);
    data_mesh.intTri->flipEdgeIfPossible(e);
    gcs::Edge flipped = e;
    diamondJacobians(data_mesh, UV, flipped, J1_prime, J2_prime);

    double after = compute_total_energy_localjacob(data_mesh, UV, et);

    //double after = energy(J1_prime) * data_mesh.intTri->faceArea(flipped.halfedge().face()) + energy(J2_prime) * data_mesh.intTri->faceArea(flipped.halfedge().twin().face());
    double tolerance = 1e-6; // set tolerance to 1e-6

    if (fabs(before - after) / std::max(fabs(before), fabs(after)) > tolerance) {
      if (before > after) {
        totalflips++;
        if(data_mesh.intTri->isDelaunay(e)) delaunay_flips++;
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
  return totalflips;
}

unsigned delaunay_flip(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et){
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
    if(data_mesh.intTri->isDelaunay(e)) continue;
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
        std::cout << "Delaunay flip increased the energy" << std::endl;
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



unsigned asIDTasPossible(DataGeo &data_mesh){
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->unrequireVertexIndices();
  unsigned totalflips = 0;
  unsigned flips = 1;
  while(flips){
    unsigned local_flips = 0;
    for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
      if(e.isBoundary()) continue;
      std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();

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
      // if flip causes self loop, don't flip
      if(v3==v4) {
        //std::cout << "self loop" << std::endl;
        continue;
      }
      // if one of the end vertices has degree smaller than 3, don't flip
      if(halfedges[1].tipVertex().degree()<=3 || halfedges[0].tailVertex().degree()<=3) {
        //std::cout << "degree 3" << std::endl;
        continue;
      }
      if(data_mesh.intTri->flipEdgeIfNotDelaunay(e)){
          local_flips++;
          data_mesh.intTri->refreshQuantities();

      }

    }
    //std::cout << "flipso: " << local_flips << std::endl;
    flips = local_flips;
    totalflips += flips;
  }

  return totalflips;
}

#include "subdivide.hpp"
#include <math.h>

Eigen::SparseVector<double> b(const gcs::SurfacePoint& pt, unsigned nvertices){
  Eigen::SparseVector<double> result;
  result.resize(nvertices);
  if (pt.type == gcs::SurfacePointType::Vertex){
    result.insert(pt.vertex.getIndex()) = 1;
  } 
  else { 
    result.insert(pt.edge.firstVertex().getIndex()) = (1-pt.tEdge);
    result.insert(pt.edge.secondVertex().getIndex()) = pt.tEdge;
  }
  return result;
}

double subdivide(const Eigen::MatrixXd &V, const Eigen::MatrixXd &UV, std::vector<gcs::SurfacePoint>& vec1, std::vector<gcs::SurfacePoint>& vec2, 
    double len1, double len2, double len3,bool complete, Eigen::MatrixXd &J, Eigen::VectorXd &areas){

  unsigned nvertices = V.rows();
  Eigen::Matrix2d E, E_tilde;
  double temp = (len3*len3 - len1*len1 - len2*len2)/(-2*len1); 
  E_tilde << len1, temp, 0 , sqrt(len2*len2 - temp*temp);
  
  //subdivided triangle size
  int face_size = 1 + (vec2.size()-2)*2;
  if(complete){
    face_size += vec1.size() - vec2.size(); 
  }

  J.resize(face_size*2,2); 
  areas.resize(face_size);

  int i=0;
  int j=0;
  double seg_len1 = 0;
  double seg_len2 = 0;
  Eigen::Vector2d edgevec1;
  Eigen::Vector2d edgevec2;
  Eigen::Matrix2d EE_tilde, EE;
  seg_len1 += (V.transpose()*(b(vec1[i+1], nvertices)-b(vec1[i], nvertices))).norm();
  seg_len2 += (V.transpose()*(b(vec2[j+1], nvertices)-b(vec2[j], nvertices))).norm();
  Eigen::Vector2d edgevec_before1 = seg_len1/len1 * E_tilde.col(0);
  Eigen::Vector2d edgevec_before2 = seg_len2/len2 * E_tilde.col(1);

  // add first triangle 
  EE_tilde << edgevec_before1, edgevec_before2;
  EE << UV.transpose()*(b(vec1[i+1], nvertices)-b(vec1[i], nvertices)), UV.transpose()*(b(vec2[j+1], nvertices)-b(vec1[i], nvertices));
  J.block(0,0,2,2) << EE * EE_tilde.inverse();
  areas(0) = abs(EE_tilde.determinant())/2;
  int idx = 1;
  ++i;
  ++j;
  while(j<vec2.size()-1){
    seg_len1 += (V.transpose()*(b(vec1[i+1], nvertices)-b(vec1[i], nvertices))).norm();
    seg_len2 += (V.transpose()*(b(vec2[j+1], nvertices)-b(vec2[j], nvertices))).norm();

    edgevec1 = seg_len1/len1 * E_tilde.col(0);
    edgevec2 = seg_len2/len2 * E_tilde.col(1);

    EE_tilde << (edgevec1-edgevec_before1), (edgevec_before2-edgevec_before1);
    EE << UV.transpose()*(b(vec1[i+1], nvertices)-b(vec1[i], nvertices)), UV.transpose()*(b(vec2[j], nvertices)-b(vec1[i], nvertices));
    J.block(idx*2,0,2,2) << EE * EE_tilde.inverse();
    areas(idx) = abs(EE_tilde.determinant())/2;

    EE_tilde << (edgevec_before2-edgevec2), (edgevec1-edgevec2);
    EE << UV.transpose()*(b(vec2[j], nvertices)-b(vec2[j+1], nvertices)), UV.transpose()*(b(vec1[i+1], nvertices)-b(vec2[j+1], nvertices));
    J.block(idx*2+2,0,2,2) << EE * EE_tilde.inverse();
    areas(idx+1) = abs(EE_tilde.determinant())/2;

    ++i;
    ++j;
    idx+=2;
    edgevec_before1 = edgevec1;
    edgevec_before2 = edgevec2;
  }

  if(complete){
    while(i<vec1.size()-1){
      seg_len1 += (V.transpose()*(b(vec1[i+1], nvertices)-b(vec1[i], nvertices))).norm();
      edgevec1 = seg_len1/len1 * E_tilde.col(0);
      EE_tilde << (edgevec1-edgevec_before1), (edgevec_before2-edgevec_before1);
      EE << UV.transpose()*(b(vec1[i+1], nvertices)-b(vec1[i], nvertices)), UV.transpose()*(b(vec2[j], nvertices)-b(vec1[i], nvertices));
      J.block(idx*2,0,2,2) << EE * EE_tilde.inverse();
      areas(idx) = abs(EE_tilde.determinant())/2;

      ++i;
      ++idx;
      edgevec_before1=edgevec1;
    }
  }
  return len1 - seg_len1;
}

bool compareVectors(const std::pair<std::vector<gcs::SurfacePoint>, double >& a, const std::pair<std::vector<gcs::SurfacePoint>, double >& b) {
    return a.first.size() < b.first.size();
}
  
double calc_energy(DataGeo &data_mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXd &UV, gcs::Halfedge he, const EnergyType &et){
  
  double (*energy)(Eigen::Matrix2d);

  if(et == EnergyType::DIRICHLET) energy = dirichlet;
  if(et == EnergyType::ASAP) energy = asap;
  if(et == EnergyType::ARAP) energy = arap;
  if(et == EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;


  gcs::Halfedge he1 = he;
  gcs::Halfedge he2 = he1.next();
  gcs::Halfedge he3 = he2.next();
  std::vector<std::pair<std::vector<gcs::SurfacePoint>, double > > vec(3);
  vec[0].first = data_mesh.intTri->traceIntrinsicHalfedgeAlongInput(he1);
  vec[1].first = data_mesh.intTri->traceIntrinsicHalfedgeAlongInput(he2);
  vec[2].first = data_mesh.intTri->traceIntrinsicHalfedgeAlongInput(he3);
  vec[0].second = data_mesh.intTri->edgeLengths[he1.edge()];
  vec[1].second = data_mesh.intTri->edgeLengths[he2.edge()];
  vec[2].second = data_mesh.intTri->edgeLengths[he3.edge()];
  
  // sort the intrinsic edges w.r.t. their number of middlepoints
  std::sort(vec.begin(), vec.end(), compareVectors);

  int i = 0;
  int j = 0 ;

  if(vec[2].first[0] == vec[1].first.back()){
    std::reverse(vec[1].first.begin(),vec[1].first.end());
  }
  else{
    std::reverse(vec[2].first.begin(),vec[2].first.end());
  }

  Eigen::MatrixXd J;
  Eigen::VectorXd areas;
  double new_seglen = subdivide(V, UV, vec[2].first, vec[1].first, vec[2].second, vec[1].second, vec[0].second, false, J, areas);
  double total_energy = 0;
  for(int i = 0; i<areas.size(); ++i){
    Eigen::Matrix2d temp = J.block(i*2,0,2,2);
    total_energy += areas(i) * energy(temp);
  }
  if(vec[2].first.back() == vec[0].first.back()){
    std::reverse(vec[0].first.begin(), vec[0].first.end());
  }
  double len3 = (V.transpose()*(b(vec[1].first.back(), V.rows())-b(vec[2].first[vec[1].first.size()-1], V.rows()))).norm();
  std::reverse(vec[2].first.begin(), vec[2].first.end());
  if(vec[2].first.size()!=vec[1].first.size()){
    std::vector<gcs::SurfacePoint> new_seg(vec[2].first.begin(),vec[2].first.begin()+vec[2].first.size()-vec[1].first.size()+1);
  
    if(new_seg.size()>vec[0].first.size()){ 
      subdivide(V, UV, new_seg, vec[0].first, new_seglen, vec[0].second, len3, true, J, areas);
    }
    else{
      subdivide(V, UV, vec[0].first, new_seg, vec[0].second, new_seglen, len3, true, J, areas);
    }
    for(int i = 0; i<areas.size(); ++i){
      Eigen::Matrix2d temp = J.block(i*2,0,2,2);
      total_energy += areas(i) * energy(temp);
    }
  }
  return total_energy;
}

unsigned flipThroughEdges(DataGeo &data_mesh, const Eigen::MatrixXd V, const Eigen::MatrixXd UV, const EnergyType &et){
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceAreas();
  unsigned totalflips = 0;
  for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
    if(e.isBoundary()) continue;
    gcs::Halfedge h1 = e.halfedge();
    gcs::Halfedge h2 = h1.twin();
    double before = calc_energy(data_mesh, V, UV, h1, et) + calc_energy(data_mesh, V, UV, h2, et);
    data_mesh.intTri->flipEdgeIfPossible(e);
    h1= e.halfedge();
    h2= h1.twin();
    double after = calc_energy(data_mesh, V, UV, h1, et)+ calc_energy(data_mesh, V, UV, h2, et); 
    double tolerance = 1e-6; // set tolerance to 1e-6 

    if (fabs(before - after) / std::max(fabs(before), fabs(after)) > tolerance) {
      if (before > after) {
        totalflips++;
      }
      else {
        data_mesh.intTri->flipEdgeIfPossible(e);
      }
    }
    else {
      data_mesh.intTri->flipEdgeIfPossible(e);
    }
  }
    
  return totalflips;
}


#include "laplace.hpp"

#include "iarap.hpp"
#include "segments.hpp"

typedef Eigen::SparseVector<double>::InnerIterator SVIter;

Eigen::SparseMatrix<double> eiLaplace(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
  iARAPData data;
  initMesh(V, F, data);
  
  std::vector<Segment> all_segs;
  calc_segments(data.V, data.intTri, all_segs);
  int nSegments = all_segs.size();
  
  typedef Eigen::Triplet<double> T;
  std::vector<T> BarysList;
  // conservative estimation
  BarysList.reserve(nSegments*4);
  
  Eigen::VectorXd W;
  W.resize(nSegments);
  
  for (size_t i = 0; i < nSegments; i++) {
    const int u = all_segs[i].startV;
    const int v = all_segs[i].endV;
    W[i] = all_segs[i].weight;
    for (SVIter iter(all_segs[i].bary); iter; ++iter){
      BarysList.push_back(T(i,iter.index(), iter.value()));
    }
  }
  
  Eigen::SparseMatrix<double> Barys;
  Barys.resize(nSegments, V.rows());
  Barys.setZero();
  Barys.setFromTriplets(BarysList.begin(), BarysList.end());
  return Barys.transpose()*W.asDiagonal()*Barys;
}


Eigen::SparseMatrix<double> asym_but_pos_eiLaplace(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
  iARAPData data;
  initMesh(V, F, data);
  
  std::vector<Segment> all_segs;
  calc_segments(data.V, data.intTri, all_segs);
  int nSegments = all_segs.size();
  
  typedef Eigen::Triplet<double> T;
  std::vector<T> BarysList;
  // conservative estimation
  const int nEdges = data.intTri->intrinsicMesh->nEdges();
  BarysList.reserve(nEdges*3);
  
  Eigen::VectorXd W;
  W.resize(nEdges*2);
  int prevV = -1;
  int j = 0;
  for (size_t i = 0; i < nSegments; i++) {
    const int u = all_segs[i].startV;
    const int v = all_segs[i].endV;
    if (prevV != u){
      // last segment of prev vertex star
      if (i > 0) {
        W[j] = all_segs[i-1].weight;
        for (SVIter iter(all_segs[i-1].bary); iter; ++iter){
          BarysList.push_back(T(j,iter.index(), iter.value()));
        }
        j++;
      }
      // first segment of current vertex star
      W[j] = all_segs[i].weight;
      for (SVIter iter(all_segs[i].bary); iter; ++iter){
        BarysList.push_back(T(j,iter.index(), iter.value()));
      }
      j++;
    }
    prevV = u;
  }
  W[j] = all_segs[nSegments-1].weight;
  for (SVIter iter(all_segs[nSegments-1].bary); iter; ++iter){
    BarysList.push_back(T(j,iter.index(), iter.value()));
  }
  
  Eigen::SparseMatrix<double> Barys;
  Barys.resize(nSegments, V.rows());
  Barys.setZero();
  Barys.setFromTriplets(BarysList.begin(), BarysList.end());
  return Barys.transpose()*W.asDiagonal()*Barys;
}

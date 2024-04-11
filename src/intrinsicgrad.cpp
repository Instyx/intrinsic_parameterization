#include "intrinsicgrad.hpp"
#include <igl/grad_intrinsic.h>
#include <igl/local_basis.h>
#include <igl/grad.h>

void computeSurfaceGradientMatrix(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::SparseMatrix<double> & D1, Eigen::SparseMatrix<double> & D2) {
  Eigen::MatrixXd F1, F2, F3;
  Eigen::SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}


void computeGrad_intrinsic(DataGeo &data_mesh, Eigen::SparseMatrix<double> & Dx,
    Eigen::SparseMatrix<double> & Dy, Eigen::VectorXd &areas){
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceIndices();

  Eigen::SparseMatrix<double>G;
  Eigen::SparseMatrix<double>OG;
  unsigned nfaces = data_mesh.intTri->intrinsicMesh->nFaces();
  unsigned nvertices = data_mesh.intTri->intrinsicMesh->nVertices();
  areas.resize(nfaces);
  Eigen::MatrixXd lengths(nfaces,3);

  // ordering the lenthgs correct order for igl::grad_intrinsic
  Eigen::MatrixXi F_new = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
  for(gcs::Face face : data_mesh.intTri->intrinsicMesh->faces()){
    size_t idx = data_mesh.intTri->faceIndices[face];
    uint e_iter = 0;
    for(gcs::Edge e : face.adjacentEdges()){
      std::array<gcs::Vertex, 2> verts = e.adjacentVertices();
      lengths(idx, (e_iter+2)%3) = data_mesh.intTri->edgeLengths[e];
      // can be removed for efficiency if we trust geometry central.
      if (!(
        (F_new(idx, e_iter) == verts[0].getIndex() && F_new(idx, (e_iter+1)%3) == verts[1].getIndex()) ||
        (F_new(idx, e_iter) == verts[1].getIndex() && F_new(idx, (e_iter+1)%3) == verts[0].getIndex()))
      ) {
        // This should never happen:
        std::cout << "This should never happen: Order fucked in triangle " << idx << "\n";
      }
      e_iter++;
    }
    areas(idx) = data_mesh.intTri->faceArea(face);
  }
  data_mesh.intTri->unrequireEdgeLengths();
  data_mesh.intTri->unrequireFaceIndices();
  // igl::grad(data_mesh.intTri->inputGeom->vertexPositions, data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>(), OG);
  igl::grad_intrinsic(lengths, F_new, G);
  Dx.resize(nfaces, nvertices);
  Dy.resize(nfaces, nvertices);
  Dx=G.block(0, 0, nfaces, nvertices);
  Dy=G.block(nfaces, 0, nfaces, nvertices);
}

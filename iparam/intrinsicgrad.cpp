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
  data_mesh.intTri->requireFaceIndices();
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceAreas();
 
  Eigen::SparseMatrix<double>G;
  unsigned nfaces = data_mesh.intTri->intrinsicMesh->nFaces();
  unsigned nvertices = data_mesh.intTri->intrinsicMesh->nVertices();
  areas.resize(nfaces);
  Eigen::MatrixXd lengths(nfaces,3);

  Eigen::MatrixXi F_new = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
  for(gcs::Face face : data_mesh.intTri->intrinsicMesh->faces()){
    size_t idx = data_mesh.intTri->faceIndices[face];
    size_t v0 = F_new(idx,0); 
    size_t v1 = F_new(idx,1); 
    size_t v2 = F_new(idx,2); 
    for(gcs::Edge e : face.adjacentEdges()){
      std::array<gcs::Vertex, 2> verts = e.adjacentVertices();
      if(data_mesh.intTri->vertexIndices[verts[0]]==v0){
        if(data_mesh.intTri->vertexIndices[verts[1]]==v1){ // v0 - v1
          lengths(idx, 2) = data_mesh.intTri->edgeLengths[e];
        }
        else{ // v0 - v2
          lengths(idx, 1) = data_mesh.intTri->edgeLengths[e];
        }
      }
      else if(data_mesh.intTri->vertexIndices[verts[0]]==v1){
        if(data_mesh.intTri->vertexIndices[verts[1]]==v0){ // v1 - v0
          lengths(idx, 2) = data_mesh.intTri->edgeLengths[e];
        }
        else{ // v1 - v2
          lengths(idx, 0) = data_mesh.intTri->edgeLengths[e];
        }
      }
      else{
        if(data_mesh.intTri->vertexIndices[verts[1]]==v0){ // v2 - v0
          lengths(idx, 1) = data_mesh.intTri->edgeLengths[e];
        }
        else{ // v2 - v1
          lengths(idx, 0) = data_mesh.intTri->edgeLengths[e];
        }
      }
    }
    areas(idx) = data_mesh.intTri->faceArea(face);
  }
  igl::grad_intrinsic(lengths, F_new, G);
  Dx.resize(nfaces, nvertices);
  Dy.resize(nfaces, nvertices);
  Dx=G.block(0, 0, nfaces, nvertices);
  Dy=G.block(nfaces, 0, nfaces, nvertices);
}


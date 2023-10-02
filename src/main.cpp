#include <asm-generic/errno-base.h>
#include <cinttypes>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>


/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/adjacency_list.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/grad_intrinsic.h>
#include <igl/edge_lengths.h>
#include <igl/intrinsic_delaunay_triangulation.h>

#include <laplace.hpp>
#include <datageo.hpp>
#include "geometrycentral/surface/edge_length_geometry.h"
#include <fstream>

using namespace std;
using namespace Eigen;

ofstream outfile("output.txt");

using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

DataGeo data_mesh;

int option_en;

bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

bool reset=false;
bool prevFreeBoundary=false;
VectorXi fixed_UV_indices;
MatrixXd fixed_UV_positions;


// void (*gradp)(SparseMatrix<double, 0, int>&, SparseMatrix<double, 0, int>&) ;
bool igrad;

MatrixXd colors;
MatrixXd C1, P1;
MatrixXd C2, P2;
vector<Eigen::Vector3d> colored_points;
unsigned iterations=0;
void Redraw()
{
	viewer.data().clear();

	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);

		viewer.data().set_face_based(false);

    if(UV.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
    }
	}
	else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);
	}
	if(colors.size()!=0&&!reset) viewer.data().set_colors(colors);
}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}

double sq(double val) {return val*val;}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}


static void computeGrad_intrinsic(SparseMatrix<double> & Dx, SparseMatrix<double> & Dy, VectorXd &areas){
    data_mesh.intTri->requireFaceIndices();
    
    SparseMatrix<double>G;
    areas.resize(F.rows());
    MatrixXd lengths(F.rows(),3);


    // TODO: check if the vertex indices are same as in V 
    MatrixXi F_new = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
   /* data_mesh.inputGeometry->requireVertexIndices();
    for(gcs::Vertex v : data_mesh.intTri->intrinsicMesh->vertices()){
      int i = v.getIndex();
      cout << data_mesh.inputGeometry->vertexPositions[v] << endl;
      cout << "V:   " << V.row(i) << endl; 
    }*/
    for(gcs::Face face : data_mesh.intTri->intrinsicMesh->faces()){
      size_t idx = data_mesh.intTri->faceIndices[face];
      size_t v0 = F_new(idx,0); 
      size_t v1 = F_new(idx,1); 
      size_t v2 = F_new(idx,2); 
      //cout << v0 << "  " << v1 << "  " << v2 << endl;
      //cout << " - edges - " << endl; 
      for(gcs::Edge e : face.adjacentEdges()){
        array<gcs::Vertex, 2> verts = e.adjacentVertices();
        //cout << data_mesh.intTri->vertexIndices[verts[0]] << "    " << data_mesh.intTri->vertexIndices[verts[1]] << endl;
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
    Dx.resize(F.rows(),V.rows());
    Dy.resize(F.rows(),V.rows());
    Dx=G.block(0,0,F.rows(), V.rows());
    Dy=G.block(F.rows(),0,F.rows(), V.rows());
}


static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}

Eigen::SparseMatrix<double> compute_L_uniform(){
    Eigen::SparseMatrix<double> Laplacian(V.rows(),V.rows());
    std::vector<Eigen::Triplet<double> > tripletList;
		std::vector<std::vector<int> > VV;
    igl::adjacency_list(F, VV);
    for (unsigned i = 0; i < VV.size(); ++i) {
        for (unsigned j = 0; j < VV[i].size(); ++j) {
            tripletList.push_back(Eigen::Triplet<double>(i, i, -1));
            tripletList.push_back(Eigen::Triplet<double>(i, VV[i][j], 1));
        }
    }
    Laplacian.setFromTriplets(tripletList.begin(), tripletList.end());
    return Laplacian;
}

SparseMatrix<double> compute_L_intrinsic(){
 // TODO check why this and cotmatrix are same ?
  SparseMatrix<double> L_int, L_cot;
  L_int = eiLaplace(V,F);
  igl::cotmatrix(V,F,L_cot);
  cout << "intrinsic - cot laplacian diff: " << (L_int + L_cot).norm() << endl; // we add because igl::cotmatrix is negated
  return L_int;

}

void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double> &C, VectorXd &d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.
	
	// size of C : #constraints x 2#V
	// C contains one 1 per row and it corresponds to one u or v coordinate
	
	// u coordinates
	vector<Triplet<double> > list;
	d = VectorXd::Zero(2*indices.size());
	for (int i = 0; i < indices.size(); ++i) {
		list.push_back(Triplet<double>(i,indices(i),1)); 
		d(i) = positions(i,0);
	}
	// v coordinates
	for (int i = 0; i < indices.size(); ++i) {
		list.push_back(Triplet<double>(i+indices.size(),indices(i)+V.rows(),1));
		d(i+indices.size()) = positions(i,1);
	}
	C.resize(2*indices.size(),2*V.rows());
	C.setFromTriplets(list.begin(), list.end());

}

void testt(){
  //data_mesh.intTri->flipToDelaunay(); 
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireEdgeCotanWeights();
  data_mesh.inputGeometry->requireVertexIndices();
  vector<gcs::Edge> v;
  double sum_nondeluanay = 0;
  double sum_deluanay = 0;
  for(gcs::Edge curr: data_mesh.intTri->intrinsicMesh->edges()) {
    if(data_mesh.intTri->isDelaunay(curr)){
      continue;
    }
    cout << "found non delaunay edge" << endl;
    v.push_back(curr);
  }
  for (gcs::Edge e : v) {
    std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
    double fi_len1 = data_mesh.intTri->edgeLengths[e];
    double fi_len2 = data_mesh.intTri->edgeLengths[halfedges[0].edge()];
    double fi_len3 = data_mesh.intTri->edgeLengths[halfedges[1].edge()];

    double se_len1 = data_mesh.intTri->edgeLengths[e];
    double se_len2 = data_mesh.intTri->edgeLengths[halfedges[2].edge()];
    double se_len3 = data_mesh.intTri->edgeLengths[halfedges[3].edge()];


    Matrix2d J1, J2;
    double temp = (fi_len2*fi_len2 - fi_len1*fi_len1 - fi_len3*fi_len3)/(-2*fi_len1); 
    J1 << fi_len1, temp, 0 , sqrt(fi_len3*fi_len3 - temp*temp); 
    temp = (se_len2*se_len2 - se_len1*se_len1 - se_len3*se_len3)/(-2*se_len1); 
    J2 << se_len1, temp, 0 , sqrt(se_len3*se_len3 - temp*temp);

    Matrix2d U1, U2, V1, V2, S1, S2;
    SSVD2x2(J1, U1, S1, V1); 
    SSVD2x2(J2, U2, S2, V2);

    sum_nondeluanay+=  S1(0,0)-S1(1,1) + S2(0,0)-S2(1,1); 
    cout << S1(0,0) << "   " << S1(1,1) << endl;
    cout << S2(0,0) << "   " << S2(1,1) << endl;
    cout << "-------------------" << endl;

 //   if(isDelaunay(e)) continue;
    data_mesh.intTri->flipEdgeIfPossible(e);

    halfedges = e.diamondBoundary();
    fi_len1 = data_mesh.intTri->edgeLengths[e];
    fi_len2 = data_mesh.intTri->edgeLengths[halfedges[0].edge()];
    fi_len3 = data_mesh.intTri->edgeLengths[halfedges[1].edge()];

    se_len1 = data_mesh.intTri->edgeLengths[e];
    se_len2 = data_mesh.intTri->edgeLengths[halfedges[2].edge()];
    se_len3 = data_mesh.intTri->edgeLengths[halfedges[3].edge()];


    Matrix2d J1_prime, J2_prime;
    temp = (fi_len2*fi_len2 - fi_len1*fi_len1 - fi_len3*fi_len3)/(-2*fi_len1); 
    J1_prime << fi_len1, temp, 0 , sqrt(fi_len3*fi_len3 - temp*temp); 
      temp = (se_len2*se_len2 - se_len1*se_len1 - se_len3*se_len3)/(-2*se_len1); 
    J2_prime << se_len1, temp, 0 , sqrt(se_len3*se_len3 - temp*temp);

    Matrix2d S1_prime, S2_prime;
    SSVD2x2(J1_prime, U1, S1_prime, V1); 
    SSVD2x2(J2_prime, U2, S2_prime, V2);
    sum_deluanay+=  S1_prime(0,0) - S1_prime(1,1) + S2_prime(0,0) - S2_prime(1,1); 
    cout  << S1_prime(0,0) << "   " <<  S1_prime(1,1) << endl;
    cout  << S2_prime(0,0) << "   " << S2_prime(1,1) << endl;
   // assert(sq(S1_prime(0,0)-1)+sq(S1_prime(1,1)-1) + sq(S2_prime(0,0)-1)+sq(S2_prime(1,1)-1) <= sq(S1(0,0)-1)+sq(S1(1,1)-1) + sq(S2(0,0)-1)+sq(S2(1,1)-1));
    cout << "****** new edge ********" << endl;
    cout << "****** new edge ********" << endl;

  }
  cout << sqrt(sum_nondeluanay) << endl;
  cout << sqrt(sum_deluanay) << endl;
  
 

}

double (*energy)(Matrix2d &J);

double drichlet(Matrix2d &J){
  return J.norm()*J.norm();
}

double symmetricDrichlet(Matrix2d &J){
  if(J.determinant()>=0){
    return J.norm() * J.norm() + J.inverse().norm()*J.inverse().norm();
  }
  else{
    cout << "determinant negative" << endl;
    return 1e7; // should be inf
  }

}

double arap(Matrix2d &J){
  Matrix2d U, S, VV;
  SSVD2x2(J, U, S, VV);
  Matrix2d R = U*VV.transpose();
  return (J-R).norm() * (J-R).norm();
   
}

double asap(Matrix2d &J){
  return (J(0,0)-J(1,1)) *  (J(0,0)-J(1,1)) +  (J(1,0)+J(0,1)) * (J(1,0)+J(0,1)); 

}

Eigen::SparseVector<double> b(const gcs::SurfacePoint& pt){
  Eigen::SparseVector<double> result;
  result.resize(V.rows());
  if (pt.type == gcs::SurfacePointType::Vertex){
    result.insert(pt.vertex.getIndex()) = 1;
  } 
  else if (pt.type == gcs::SurfacePointType::Edge) {
    result.insert(pt.edge.firstVertex().getIndex()) = (1-pt.tEdge);
    result.insert(pt.edge.secondVertex().getIndex()) = pt.tEdge;
  }
  else{
    int i = 0;
    for(gcs::Vertex v : pt.face.adjacentVertices()){
      result.insert(v.getIndex())= pt.faceCoords[i];
      ++i;
    }
    colored_points.push_back(V.transpose() * result);
    cout << "FACEEE in b" << endl;
    cout << pt.faceCoords << endl;
  }
  return result;
}

double cross2d(Vector2d &v1, Vector2d &v2){
  return v1(0)*v2(1)-v1(1)*v2(0);
}

// the points are in order with the polygon rotation
bool isConcave(vector<Vector2d>& points){
  Vector2d e1 = points[1]-points[0];
  Vector2d e2 = points[2]-points[1];
  Vector2d e3 = points[3]-points[2];
  Vector2d e4 = points[0]-points[3];
  
  double cross1 = cross2d(e1,e2);
  double cross2 = cross2d(e2,e3);
  double cross3 = cross2d(e3,e4);
  double cross4 = cross2d(e4,e1);

  if (cross1 * cross2 > 0 || cross2 * cross3 > 0 || cross3 * cross4 > 0) {
    return false;
  }
    
    return true;
}



// TODO check whether the order of edges change the Jacobian
bool diamondJacobians(gcs::Edge &e, Matrix2d &J1, Matrix2d &J2){
  std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
  double fi_len1 = data_mesh.intTri->edgeLengths[e];
  double fi_len2 = data_mesh.intTri->edgeLengths[halfedges[0].edge()];
  double fi_len3 = data_mesh.intTri->edgeLengths[halfedges[1].edge()];

  double se_len1 = data_mesh.intTri->edgeLengths[e];
  double se_len2 = data_mesh.intTri->edgeLengths[halfedges[2].edge()];
  double se_len3 = data_mesh.intTri->edgeLengths[halfedges[3].edge()];

  Matrix2d E1, E2, E1_tilde, E2_tilde;
  double temp = (fi_len2*fi_len2 - fi_len1*fi_len1 - fi_len3*fi_len3)/(-2*fi_len1); 
  E1_tilde << fi_len1, temp, 0 , sqrt(fi_len3*fi_len3 - temp*temp); 
  temp = (se_len2*se_len2 - se_len1*se_len1 - se_len3*se_len3)/(-2*se_len1); 
  E2_tilde << se_len1, temp, 0 , sqrt(se_len3*se_len3 - temp*temp);
 
  data_mesh.inputGeometry->requireVertexPositions();

  size_t v1 = data_mesh.intTri->vertexIndices[halfedges[1].tipVertex()];
  size_t v2 = data_mesh.intTri->vertexIndices[halfedges[0].tailVertex()];
  size_t v3 = data_mesh.intTri->vertexIndices[halfedges[0].tipVertex()];
  size_t v4 = data_mesh.intTri->vertexIndices[halfedges[2].tipVertex()];
 
  vector<Vector2d> points(4);
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

void faceJacobian(gcs::Face &f, Matrix2d &J){
  int i = 0;
  std::array<gcs::Halfedge, 3> halfedges;

  for(gcs::Halfedge he : f.adjacentHalfedges()){
    if(i==3) break; // assumed triangle faces
    halfedges[i]=he;
    ++i;
  } 
  double fi_len1 = data_mesh.intTri->edgeLengths[halfedges[0].edge()];
  double fi_len2 = data_mesh.intTri->edgeLengths[halfedges[1].edge()];
  double fi_len3 = data_mesh.intTri->edgeLengths[halfedges[2].edge()];

  Matrix2d E, E_tilde;
  double temp = (fi_len2*fi_len2 - fi_len1*fi_len1 - fi_len3*fi_len3)/(-2*fi_len1); 
  E_tilde << fi_len1, temp, 0 , sqrt(fi_len3*fi_len3 - temp*temp); 
  

  size_t v1 = data_mesh.intTri->vertexIndices[halfedges[1].tipVertex()];
  size_t v2 = data_mesh.intTri->vertexIndices[halfedges[0].tailVertex()];
  size_t v3 = data_mesh.intTri->vertexIndices[halfedges[0].tipVertex()];
  E << UV(v3,0) - UV(v2,0), UV(v1,0) - UV(v2,0),
        UV(v3,1) - UV(v2,1), UV(v1,1) - UV(v2,1);
  J = E * E_tilde.inverse();
  
}

void triangleJacobian(std::vector<gcs::SurfacePoint> &vec, Matrix2d &J){
  double len1 = (V.transpose()*(b(vec[1])-b(vec[0]))).norm();
  double len2 = (V.transpose()*(b(vec[2])-b(vec[0]))).norm();
  double len3 = (V.transpose()*(b(vec[2])-b(vec[1]))).norm();
  Matrix2d E, E_tilde;
  double temp = (len2*len2 - len1*len1 - len3*len3)/(-2*len1); 
  E_tilde << len1, temp, 0 , sqrt(len3*len3 - temp*temp);
  Vector2d seg1 = UV.transpose() * (b(vec[1])-b(vec[0]));
  Vector2d seg2 = UV.transpose() * (b(vec[2])-b(vec[0]));
  E << seg1, seg2;
  J = E * E_tilde.inverse();

  
}

Vector2d c(vector<gcs::SurfacePoint> &vec, double t){
  return UV.transpose()*b(vec[0])*(1-t) + UV.transpose()*b(vec.back())*t;
}


double triangle(vector<gcs::SurfacePoint>& vec1, vector<gcs::SurfacePoint>& vec2, double len1, double len2, double len3, bool complete, MatrixXd &J, VectorXd &areas){
 // cout << "lengts: " << len1 << "  " << len2 << "  " << len3 << endl;
  Matrix2d E, E_tilde;
  double temp = (len3*len3 - len1*len1 - len2*len2)/(-2*len1); 
  E_tilde << len1, temp, 0 , sqrt(len2*len2 - temp*temp);
  if(len2*len2 < temp* temp) cout <<  " ABOOOOOO " << endl;
  int face_size = 1 + (vec2.size()-2)*2;
  if(complete){
    face_size += vec1.size() - vec2.size(); 
  }
  J.resize(face_size*2,2); 
  areas.resize(face_size);
  P1.resize(face_size,2);
  P2.resize(face_size,2);

  int i=0;
  int j=0;
  double seg_len1 = 0;
  double seg_len2 = 0;
  Vector2d edgevec1;
  Vector2d edgevec2;
  Matrix2d EE_tilde, EE;
  seg_len1 += (V.transpose()*(b(vec1[i+1])-b(vec1[i]))).norm();
  seg_len2 += (V.transpose()*(b(vec2[j+1])-b(vec2[j]))).norm();
  Vector2d edgevec_before1 = seg_len1/len1 * E_tilde.col(0);
  Vector2d edgevec_before2 = seg_len2/len2 * E_tilde.col(1);
  // add first triangle 
  EE_tilde << edgevec_before1, edgevec_before2;
  EE << UV.transpose()*(b(vec1[i+1])-b(vec1[i])), UV.transpose()*(b(vec2[j+1])-b(vec1[i]));
  J.block(0,0,2,2) << EE * EE_tilde.inverse();
  areas(0) = abs(EE_tilde.determinant())/2;
  P1.row(0) = b(vec1[i+1]).transpose() * UV; 
  P2.row(0) = b(vec2[j+1]).transpose() * UV; 
  int idx = 1;
  ++i;
  ++j;
  //cout << " 3D flat: " << endl; 
  //cout << EE_tilde << endl;
  //cout << " UV " << endl;
  //cout << EE << endl;
  while(j<vec2.size()-1){
    seg_len1 += (V.transpose()*(b(vec1[i+1])-b(vec1[i]))).norm();
    seg_len2 += (V.transpose()*(b(vec2[j+1])-b(vec2[j]))).norm();
    edgevec1 = seg_len1/len1 * E_tilde.col(0);
    edgevec2 = seg_len2/len2 * E_tilde.col(1);
    EE_tilde << (edgevec1-edgevec_before1), (edgevec_before2-edgevec_before1);
    EE << UV.transpose()*(b(vec1[i+1])-b(vec1[i])), UV.transpose()*(b(vec2[j])-b(vec1[i]));
    J.block(idx*2,0,2,2) << EE * EE_tilde.inverse();
    areas(idx) = abs(EE_tilde.determinant())/2;
    EE_tilde << (edgevec_before2-edgevec2), (edgevec1-edgevec2);
    EE << UV.transpose()*(b(vec2[j])-b(vec2[j+1])), UV.transpose()*(b(vec1[i+1])-b(vec2[j+1]));
    J.block(idx*2+2,0,2,2) << EE * EE_tilde.inverse();
    areas(idx+1) = abs(EE_tilde.determinant())/2;

    P1.row(idx) = b(vec1[i]).transpose()*UV;
    P2.row(idx) = b(vec2[j+1]).transpose()*UV;
    P1.row(idx+1) = b(vec1[i+1]).transpose()*UV;
    P2.row(idx+1) = b(vec2[j+1]).transpose()*UV;

    ++i;
    ++j;
    idx+=2;
    edgevec_before1 = edgevec1;
    edgevec_before2 = edgevec2;
  }
  if(complete){
    while(i<vec1.size()-1){
      seg_len1 += (V.transpose()*(b(vec1[i+1])-b(vec1[i]))).norm();
      edgevec1 = seg_len1/len1 * E_tilde.col(0);
      EE_tilde << (edgevec1-edgevec_before1), (edgevec_before2-edgevec_before1);
      EE << UV.transpose()*(b(vec1[i+1])-b(vec1[i])), UV.transpose()*(b(vec2[j])-b(vec1[i]));
      J.block(idx*2,0,2,2) << EE * EE_tilde.inverse();
      areas(idx) = abs(EE_tilde.determinant())/2;
      P1.row(idx) = UV.transpose()*b(vec1[i+1]);
      P2.row(idx) = UV.transpose()*b(vec2[j]);

      ++i;
      ++idx;
      edgevec_before1=edgevec1;
    }
  }
  return len1 - seg_len1;
}

double triangle_c(vector<gcs::SurfacePoint>& vec1, vector<gcs::SurfacePoint>& vec2, double len1, double len2, double len3, bool complete, MatrixXd &J, VectorXd &areas){
  vector<double> tvalues1 = data_mesh.intTri->recoverTraceTValues(vec1);
  vector<double> tvalues2 = data_mesh.intTri->recoverTraceTValues(vec2);
  Matrix2d E, E_tilde;
  double temp = (len3*len3 - len1*len1 - len2*len2)/(-2*len1); 
  E_tilde << len1, temp, 0 , sqrt(len2*len2 - temp*temp);
  if(len2*len2 < temp* temp) cout <<  " ABOOOOOO " << endl;
  int face_size = 1 + (vec2.size()-2)*2;
  if(complete){
    face_size += vec1.size() - vec2.size(); 
  }
  J.resize(face_size*2,2); 
  areas.resize(face_size);
  P1.resize(face_size,2);
  P2.resize(face_size,2);

  int i=0;
  int j=0;

  Vector2d edgevec1;
  Vector2d edgevec2;
  Matrix2d EE_tilde, EE;
  Vector2d edgevec_before1 = tvalues1[i+1] * E_tilde.col(0);
  Vector2d edgevec_before2 = tvalues2[j+1] * E_tilde.col(1);
  // add first triangle 
  EE_tilde << edgevec_before1, edgevec_before2;
  EE << c(vec1, tvalues1[i+1]) - c(vec1,tvalues1[i]), c(vec2,tvalues2[j+1]) - c(vec1,tvalues1[i]);
  J.block(0,0,2,2) << EE * EE_tilde.inverse();
  areas(0) = abs(EE_tilde.determinant())/2;
  P1.row(0) = b(vec1[i+1]).transpose() * UV; 
  P2.row(0) = b(vec2[j+1]).transpose() * UV; 
  int idx = 1;
  ++i;
  ++j;
  while(j<vec2.size()-1){
    edgevec1 = tvalues1[i+1] * E_tilde.col(0);
    edgevec2 = tvalues2[j+1] * E_tilde.col(1);
    EE_tilde << (edgevec1-edgevec_before1), (edgevec_before2-edgevec_before1);
    EE << c(vec1, tvalues1[i+1]) - c(vec1,tvalues1[i]), c(vec2,tvalues2[j]) - c(vec1,tvalues1[i]);

    J.block(idx*2,0,2,2) << EE * EE_tilde.inverse();
    areas(idx) = abs(EE_tilde.determinant())/2;
    EE_tilde << (edgevec_before2-edgevec2), (edgevec1-edgevec2);
    EE << c(vec2, tvalues2[j]) - c(vec2,tvalues2[j+1]), c(vec1,tvalues1[i+1]) - c(vec2,tvalues2[j+1]);

    J.block(idx*2+2,0,2,2) << EE * EE_tilde.inverse();
    areas(idx+1) = abs(EE_tilde.determinant())/2;

    P1.row(idx) = b(vec1[i]).transpose()*UV;
    P2.row(idx) = b(vec2[j+1]).transpose()*UV;
    P1.row(idx+1) = b(vec1[i+1]).transpose()*UV;
    P2.row(idx+1) = b(vec2[j+1]).transpose()*UV;

    ++i;
    ++j;
    idx+=2;
    edgevec_before1 = edgevec1;
    edgevec_before2 = edgevec2;
  }
  if(complete){
    while(i<vec1.size()-1){
      edgevec1 = tvalues1[i+1] * E_tilde.col(0);
      EE_tilde << (edgevec1-edgevec_before1), (edgevec_before2-edgevec_before1);
      EE << c(vec1, tvalues1[i+1]) - c(vec1,tvalues1[i]), c(vec2,tvalues2[j]) - c(vec1,tvalues1[i]);
      J.block(idx*2,0,2,2) << EE * EE_tilde.inverse();
      areas(idx) = abs(EE_tilde.determinant())/2;
      P1.row(idx) = UV.transpose()*b(vec1[i+1]);
      P2.row(idx) = UV.transpose()*b(vec2[j]);

      ++i;
      ++idx;
      edgevec_before1=edgevec1;
    }
  }
  return len1 * (1-tvalues1[i]);
}



bool compareVectors(const pair<vector<gcs::SurfacePoint>, double >& a, const pair<vector<gcs::SurfacePoint>, double >& b) {
    return a.first.size() < b.first.size();
}
  
// TODO:: think about the orientation when reversed and so, make sure they have the same orientation maybe??
double calc_energy(gcs::Halfedge he){
  gcs::Halfedge he1 = he;
  gcs::Halfedge he2 = he1.next();
  gcs::Halfedge he3 = he2.next();
  vector<pair<vector<gcs::SurfacePoint>, double > > vec(3);
  vec[0].first = data_mesh.intTri->traceIntrinsicHalfedgeAlongInput(he1);
  vec[1].first = data_mesh.intTri->traceIntrinsicHalfedgeAlongInput(he2);
  vec[2].first = data_mesh.intTri->traceIntrinsicHalfedgeAlongInput(he3);
  vec[0].second = data_mesh.intTri->edgeLengths[he1.edge()];
  vec[1].second = data_mesh.intTri->edgeLengths[he2.edge()];
  vec[2].second = data_mesh.intTri->edgeLengths[he3.edge()];
  sort(vec.begin(), vec.end(), compareVectors);
  //cout << " intrinsic sizes " << endl;
  //cout << vec[0].first.size() << "  " << vec[1].first.size() << "  " << vec[2].first.size() << endl;
  int i = 0;
  int j =0 ;
  /*vector<double> ts = data_mesh.intTri->recoverTraceTValues(vec[2].first);
  for(double d : ts){
    cout << d << "  ";
  }
  cout << endl;*/
  if(vec[2].first[0] == vec[1].first.back()){
    reverse(vec[1].first.begin(),vec[1].first.end());
  }
  else{
    reverse(vec[2].first.begin(),vec[2].first.end());
  }
  MatrixXd J;
  VectorXd areas;
  MatrixXd temp1, temp2, temp3, temp4;
  // cout << "start vertices: " <<vec[1].first[0] << "  " << vec[2].first[0] << endl; 
  double new_seglen = triangle(vec[2].first, vec[1].first, vec[2].second, vec[1].second, vec[0].second, false, J, areas);
  double total_energy = 0;
  for(int i = 0; i<areas.size(); ++i){
    Matrix2d temp = J.block(i*2,0,2,2);
    total_energy += areas(i) * energy(temp);
  }
  temp1=P1;
  temp2=P2;
  //cout << " first done " << endl;
  // cout << total_energy << endl;
  if(vec[2].first.back() == vec[0].first.back()){
    reverse(vec[0].first.begin(), vec[0].first.end());
  }
  double len3 = (V.transpose()*(b(vec[1].first.back())-b(vec[2].first[vec[1].first.size()-1]))).norm();
  reverse(vec[2].first.begin(), vec[2].first.end());
  if(vec[2].first.size()!=vec[1].first.size()){
   // cout << " second in " << endl;
    vector<gcs::SurfacePoint> new_seg(vec[2].first.begin(),vec[2].first.begin()+vec[2].first.size()-vec[1].first.size()+1);
    //cout << new_seg.size() << "   " << vec[0].first.size() << endl;
    // cout << "start vertices 2: " <<vec[0].first[0] << "  " << new_seg[0] << endl; 
    // cout << new_seglen << "  " << vec[0].second << "  " << len3 << endl;
  
    if(new_seg.size()>vec[0].first.size()){ 
      triangle(new_seg, vec[0].first, new_seglen, vec[0].second, len3, true, J, areas);
    }
    else{
      triangle(vec[0].first, new_seg, vec[0].second, new_seglen, len3, true, J, areas);
    }
    temp3 = P1;
    temp4 = P2;
    P1.resize(temp1.rows()+temp3.rows(),2);
    P2.resize(temp2.rows()+temp4.rows(),2);
    P1 << temp1, temp3;
    P2 << temp2, temp4;
   // cout << " second done " << endl;
    for(int i = 0; i<areas.size(); ++i){
      Matrix2d temp = J.block(i*2,0,2,2);
      total_energy += areas(i) * energy(temp);
    }
  //  cout << total_energy << endl;
  }
  return total_energy;
}


// intrinsic edges are straight in the UV mapping
double calc_energy_intri(gcs::Halfedge &he){
  gcs::Halfedge he1 = he;
  gcs::Halfedge he2 = he1.next();
  gcs::Halfedge he3 = he2.next();
  vector<pair<vector<gcs::SurfacePoint>, double > > vec(3);
  vec[0].first = data_mesh.intTri->traceInputHalfedgeAlongIntrinsic(he1);
  vec[1].first = data_mesh.intTri->traceInputHalfedgeAlongIntrinsic(he2);
  vec[2].first = data_mesh.intTri->traceInputHalfedgeAlongIntrinsic(he3);
  vec[0].second = data_mesh.inputGeometry->edgeLength(he1.edge());
  vec[1].second = data_mesh.inputGeometry->edgeLength(he2.edge());
  vec[2].second = data_mesh.inputGeometry->edgeLength(he3.edge());
  sort(vec.begin(), vec.end(), compareVectors);
  //cout << " intrinsic sizes " << endl;
  //cout << vec[0].first.size() << "  " << vec[1].first.size() << "  " << vec[2].first.size() << endl;
  int i = 0;
  int j = 0 ;
  if(vec[2].first[0] == vec[1].first.back()){
    reverse(vec[1].first.begin(),vec[1].first.end());
  }
  else{
    reverse(vec[2].first.begin(),vec[2].first.end());
  }
  MatrixXd J;
  VectorXd areas;
  MatrixXd temp1, temp2, temp3, temp4;
  // cout << "start vertices: " <<vec[1].first[0] << "  " << vec[2].first[0] << endl; 
  double new_seglen = triangle_c(vec[2].first, vec[1].first, vec[2].second, vec[1].second, vec[0].second, false, J, areas);
  double total_energy = 0;
  for(int i = 0; i<areas.size(); ++i){
    Matrix2d temp = J.block(i*2,0,2,2);
    total_energy += areas(i) * energy(temp);
  }
  temp1=P1;
  temp2=P2;
  //cout << " first done " << endl;
  // cout << total_energy << endl;
  if(vec[2].first.back() == vec[0].first.back()){
    reverse(vec[0].first.begin(), vec[0].first.end());
  }
  double len3 = (V.transpose()*(b(vec[1].first.back())-b(vec[2].first[vec[1].first.size()-1]))).norm();
  reverse(vec[2].first.begin(), vec[2].first.end());
  if(vec[2].first.size()!=vec[1].first.size()){
   // cout << " second in " << endl;
    vector<gcs::SurfacePoint> new_seg(vec[2].first.begin(),vec[2].first.begin()+vec[2].first.size()-vec[1].first.size()+1);
    //cout << new_seg.size() << "   " << vec[0].first.size() << endl;
    // cout << "start vertices 2: " <<vec[0].first[0] << "  " << new_seg[0] << endl; 
    // cout << new_seglen << "  " << vec[0].second << "  " << len3 << endl;
  
    if(new_seg.size()>vec[0].first.size()){
      triangle_c(new_seg, vec[0].first, new_seglen, vec[0].second, len3, true, J, areas);
    }
    else{
      triangle_c(vec[0].first, new_seg, vec[0].second, new_seglen, len3, true, J, areas);
    }
    temp3 = P1;
    temp4 = P2;
    P1.resize(temp1.rows()+temp3.rows(),2);
    P2.resize(temp2.rows()+temp4.rows(),2);
    P1 << temp1, temp3;
    P2 << temp2, temp4;
   // cout << " second done " << endl;
    for(int i = 0; i<areas.size(); ++i){
      Matrix2d temp = J.block(i*2,0,2,2);
      total_energy += areas(i) * energy(temp);
    }
  //  cout << total_energy << endl;
  }
  return total_energy;
}


/*
void diamondJacobians_uv(gcs::Edge &e, Matrix2d &J1, Matrix2d &J2){
  double len1, len2, len3;
 MatrixXi F_new = data_mesh.inputMesh->getFaceVertexMatrix<int>();
  for(gcs::Face face : e.adjacentFaces()){
    size_t idx = data_mesh.intTri->faceIndices[face];
    size_t v0 = F_new(idx,0); 
    size_t v1 = F_new(idx,1); 
    size_t v2 = F_new(idx,2); 
    for(gcs::Edge e : face.adjacentEdges()){
      array<gcs::Vertex, 2> verts = e.adjacentVertices();
      if(data_mesh.intTri->vertexIndices[verts[0]]==v0){
        if(data_mesh.intTri->vertexIndices[verts[1]]==v1){ // v0 - v1
          len3 = data_mesh.intTri->edgeLengths[e];
        }
        else{ // v0 - v2
          len2 = data_mesh.intTri->edgeLengths[e];
        }
      }
      else if(data_mesh.intTri->vertexIndices[verts[0]]==v1){
        if(data_mesh.intTri->vertexIndices[verts[1]]==v0){ // v1 - v0
          len3 = data_mesh.intTri->edgeLengths[e];
        }
        else{ // v1 - v2
          len1 = data_mesh.intTri->edgeLengths[e];
        }
      }
      else{
        if(data_mesh.intTri->vertexIndices[verts[1]]==v0){ // v2 - v0
          len2 = data_mesh.intTri->edgeLengths[e];
        }
        else{ // v2 - v1
          len1 = data_mesh.intTri->edgeLengths[e];
        }
      }
    }
  Matrix2d E, E_tilde;
  double temp = (len2*len2 - len1*len1 - len3*len3)/(-2*len1); 
  E_tilde << len1, temp, 0 , sqrt(len3*len3 - temp*temp); 
  

  size_t v1 = data_mesh.intTri->vertexIndices[halfedges[1].tipVertex()];
  size_t v2 = data_mesh.intTri->vertexIndices[halfedges[0].tailVertex()];
  size_t v3 = data_mesh.intTri->vertexIndices[halfedges[0].tipVertex()];
  
  E1 << UV(v2,0) - UV(v1,0), UV(v3,0) - UV(v1,0),
        UV(v2,1) - UV(v2,1), UV(v3,1) - UV(v1,1);
  J1 = E1 * E1_tilde.inverse();
  
    areas(idx) = data_mesh.intTri->faceAreas[face];
  }
  
}
*/

double flippeddiff(gcs::Edge e){
  if(e.isBoundary()) return 0;
  gcs::Face f1 = e.halfedge().face(); 
  gcs::Face f2 = e.halfedge().twin().face();
  Matrix2d J1, J2, J1_prime, J2_prime;

  diamondJacobians(e, J1, J2);
  double before = energy(J1) * data_mesh.intTri->faceArea(f1) +
                  energy(J2) * data_mesh.intTri->faceArea(f2);
  //cout << " before flip " << data_mesh.intTri->edgeLengths[e] << endl;
  data_mesh.intTri->flipEdgeIfPossible(e);
  //cout << " after flip " << endl;
  //cout << " after flip " << data_mesh.intTri->edgeLengths[e] << endl;
  gcs::Edge flipped = e;
  //cout << " after finding the flipped edge" << endl;

  // if flip is not possible flip it back and return 0
  if(!diamondJacobians(flipped, J1_prime, J2_prime)){
    data_mesh.intTri->flipEdgeIfPossible(flipped);
    return 0;
  }

  double after = energy(J1_prime) * data_mesh.intTri->faceArea(flipped.halfedge().face()) + energy(J2_prime) * data_mesh.intTri->faceArea(flipped.halfedge().twin().face());
  data_mesh.intTri->flipEdgeIfPossible(flipped);
  return after-before;
}


unsigned greedyflip(){
  unsigned totalflips=0;
  gcs::EdgeData<double> diffs(* (data_mesh.intTri->intrinsicMesh));
  vector<size_t> indices(diffs.size());
  vector<int> visited(diffs.size());
  int i=0;
  for(gcs::Edge e : data_mesh.intTri->intrinsicMesh->edges()){
    double energydiff = flippeddiff(e);
    diffs[e] = energydiff;
    indices[i]= e.getIndex();
    ++i;
  }
  sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
    return diffs[i] < diffs[j];
    });
  for (size_t i = 0; i < diffs.size(); i++) {
    size_t idx = indices[i];
    if(diffs[idx]>=0) break;
    if(visited[idx]) continue;
    gcs::Edge e = data_mesh.intTri->intrinsicMesh->edge(idx);
    data_mesh.intTri->flipEdgeIfPossible(e);
    std::array<gcs::Halfedge, 4> halfedges = e.diamondBoundary();
    visited[idx]=1;
    visited[halfedges[0].edge().getIndex()]=1;
    visited[halfedges[1].edge().getIndex()]=1;
    visited[halfedges[2].edge().getIndex()]=1;
    visited[halfedges[3].edge().getIndex()]=1;
    ++totalflips;
  }
  cout << " totalflips: " << totalflips << endl;
  return totalflips;
}


unsigned flipThroughEdges(){
  assert(UV.size()!=0 && "Computer parameterization first!!");
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceAreas();
  //while(true){
  // TODO: compute the before and the after energy per face and not in this function per edge
    unsigned totalflips = 0;
    /*double energy_before=0, energy_after=0;
    for (gcs::Face face : data_mesh.intTri->intrinsicMesh->faces() ){
      Matrix2d J;
      faceJacobian(face,J);
      energy_before+= energy(J) * data_mesh.intTri->faceArea(face);
    }*/
    for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
      if(e.isBoundary()) continue;
      gcs::Face f1 = e.halfedge().face(); 
      gcs::Face f2 = e.halfedge().twin().face();
      Matrix2d J1, J2, J1_prime, J2_prime;
      diamondJacobians(e, J1, J2);
      double before = energy(J1) * data_mesh.intTri->faceArea(f1) +
                      energy(J2) * data_mesh.intTri->faceArea(f2);
      //cout << " before flip " << data_mesh.intTri->edgeLengths[e] << endl;
      data_mesh.intTri->flipEdgeIfPossible(e);
      //cout << " after flip " << endl;
      //cout << " after flip " << data_mesh.intTri->edgeLengths[e] << endl;
      gcs::Edge flipped = e;
      //cout << " after finding the flipped edge" << endl;
      diamondJacobians(flipped, J1_prime, J2_prime);
      // if one of the determinants is negative, then the axis of the triangle got different direction
      if(J1_prime.determinant()<0 || J2_prime.determinant()<0){
          data_mesh.intTri->flipEdgeIfPossible(flipped);
          continue;
      }
      double after = energy(J1_prime) * data_mesh.intTri->faceArea(flipped.halfedge().face()) + energy(J2_prime) * data_mesh.intTri->faceArea(flipped.halfedge().twin().face());
      //cout << " energies: " << before << " -> " << after << endl;
      double tolerance = 1e-6; // set tolerance to 1e-6 

      if (fabs(before - after) / max(fabs(before), fabs(after)) > tolerance) {
        if (before > after) {
          totalflips++;
          // cout << "flipped diff:  " << before - after << endl;
        }
        else {
          data_mesh.intTri->flipEdgeIfPossible(flipped);
        }
      }
      else {
        data_mesh.intTri->flipEdgeIfPossible(flipped);
      }

    //if(totalflips==0) break;
    }
    data_mesh.intTri->refreshQuantities();
    /*
    for (gcs::Face face : data_mesh.intTri->intrinsicMesh->faces() ){
      Matrix2d J;
      faceJacobian(face,J);
      energy_after+= energy(J) * data_mesh.intTri->faceArea(face);
    }
  cout << " ENERGY: " << energy_before << "  -->  " << energy_after << endl;
  */
  cout << " totalflips: " << totalflips << endl;
 /* for(int i : visited){
    cout << i << endl;
  } */
  return totalflips;
}
unsigned flipThroughEdges_new(){
  assert(UV.size()!=0 && "Computer parameterization first!!");
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireFaceAreas();
  //while(true){
    unsigned totalflips = 0;
    double energy_before=0, energy_after=0;
    for (gcs::Face face : data_mesh.intTri->intrinsicMesh->faces() ){
      energy_before+= calc_energy(face.halfedge());
    }
    for(gcs::Edge e: data_mesh.intTri->intrinsicMesh->edges()) {
      if(e.isBoundary()) continue;
      gcs::Halfedge h1 = e.halfedge();
      gcs::Halfedge h2 = h1.twin();
      double before = calc_energy(h1) + calc_energy(h2);
      cout << " before flip " << data_mesh.intTri->edgeLengths[e] << endl;
      data_mesh.intTri->flipEdgeIfPossible(e);
      //cout << " after flip " << endl;
      cout << " after flip " << data_mesh.intTri->edgeLengths[e] << endl;
      h1= e.halfedge();
      h2= h1.twin();
      //cout << " after finding the flipped edge" << endl;
      double after = calc_energy(h1)+ calc_energy(h2); 
      cout << " energies: " << before << " -> " << after << endl;
      double tolerance = 1e-6; // set tolerance to 1e-6 

      if (fabs(before - after) / max(fabs(before), fabs(after)) > tolerance) {
        if (before > after) {
          totalflips++;
          cout << "flipped diff:  " << before - after << endl;
        }
        else {
          data_mesh.intTri->flipEdgeIfPossible(e);
        }
      }
      else {
        data_mesh.intTri->flipEdgeIfPossible(e);
        

      }

    //if(totalflips==0) break;
    }
    
    for (gcs::Face face : data_mesh.intTri->intrinsicMesh->faces() ){
      energy_after+= calc_energy(face.halfedge()); 
    }
  cout << " totalflips: " << totalflips << endl;
  cout << " ENERGY: " << energy_before << "  -->  " << energy_after << endl;
 /* for(int i : visited){
    cout << i << endl;
  } */
  return totalflips;
}

// TODO:: for the face jacobian, try to use a function from igl and check if the determinants are
// still negative
double compute_total_energy(){
  double total_energy = 0;
  int flipped_triangles = 0;
  for(gcs::Face f : data_mesh.intTri->intrinsicMesh->faces()){
    Matrix2d J;
    faceJacobian(f, J);
    if(J.determinant()<0){
      //cout << " determinant negative " << endl;
      flipped_triangles++;
    }
    total_energy += energy(J)*data_mesh.intTri->faceAreas[f];
  }
  cout << " number flipped triangles: " << flipped_triangles << endl;
  cout << " total triangles: " << data_mesh.intTri->intrinsicMesh->nFaces() << endl;
 return total_energy;
}

void computeParameterizationIntrinsic(int type){
  cout << "in parameterization" << endl;
	SparseMatrix<double> A;
	VectorXd b;
	Eigen::SparseMatrix<double> C;
	VectorXd d;
	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices

	// prevFreeBoundary is used so that contraints do not need to be computed each time since they are the same
	// and brute force is expensive
	if(prevFreeBoundary!=freeBoundary||fixed_UV_indices.size()==0){
		if (!freeBoundary)
		{
			// The boundary vertices should be fixed to positions on the unit disc. Find these position and
			// save them in the #V x 2 matrix fixed_UV_position.
			igl::boundary_loop(F, fixed_UV_indices);
			igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
		}
		else{
			// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.
			
			// brute force to find most distant 2 vertices
			unsigned fixed1 = 0, fixed2 = 0;
			for (unsigned i = 0; i < V.rows(); ++i) {
				for (unsigned j = 0; j < V.rows(); ++j) {
					if((V.row(i)-V.row(j)).norm()>(V.row(fixed1)-V.row(fixed2)).norm()){
						fixed1 = i;
						fixed2 = j;
					}
				}
			}
			fixed_UV_indices.resize(2);
			fixed_UV_positions.resize(2,2);
			fixed_UV_indices << fixed1, fixed2;
			fixed_UV_positions << 1,0,0,1;
			cout << fixed_UV_indices.size() << endl;
			cout << fixed_UV_positions << endl;

	  }
	}
	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d); 


	A.resize(2*V.rows(),2*V.rows());
  b.resize(2*V.rows());
  if(type=='1'){
    data_mesh.intTri->requireCotanLaplacian();
    SparseMatrix<double> L = data_mesh.intTri->cotanLaplacian;
    vector<Triplet<double> > tlist;
		for (int i = 0; i < L.outerSize(); ++i) {
			for (SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(L,i); it; ++it) {
				tlist.push_back(Triplet<double>(it.row(), it.col(), it.value())); 
				tlist.push_back(Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
		A.setFromTriplets(tlist.begin(), tlist.end());		
		b = VectorXd::Zero(2*V.rows());
  }

  if(type=='2'){
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!
		
		SparseMatrix<double> Dx, Dy;
		VectorXd areas(F.rows());
		computeGrad_intrinsic(Dx,Dy,areas);
		// A = (DxADx + DyADy    DxADy - DyADx)
		//     (-DxADy + DyADx   DxADx + DyADy)
		SparseMatrix<double> B1 = Dx.transpose()*areas.asDiagonal()*Dx + Dy.transpose()*areas.asDiagonal()*Dy; 
		SparseMatrix<double> B2 = Dx.transpose()*areas.asDiagonal()*Dy - Dy.transpose()*areas.asDiagonal()*Dx;
    cout << "as" << endl;
		vector<Triplet<double> > tlist;
		for (int i = 0; i < B1.outerSize(); ++i) {
			for (SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(B1,i); it; ++it) {
				tlist.push_back(Triplet<double>(it.row(), it.col(), it.value()));
				tlist.push_back(Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
  

		for (int i = 0; i < B2.outerSize(); ++i) {
			for (SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(B2,i); it; ++it) {
				tlist.push_back(Triplet<double>(it.row()+V.rows(), it.col(), it.value()));
				tlist.push_back(Triplet<double>(it.row(), it.col()+V.rows(), -it.value()));
			}
		}		
		A.setFromTriplets(tlist.begin(), tlist.end());
		b = VectorXd::Zero(2*V.rows());
  }


  if(type=='3'){
    if(UV.size()==0) computeParameterizationIntrinsic('1');
    data_mesh.intTri->requireFaceIndices();
    
    SparseMatrix<double> Dx, Dy, G;
    VectorXd areas(F.rows());
    MatrixXd lengths(F.rows(),3);


    //while(flipThroughEdges());

    unsigned before= 0;
    unsigned curr = flipThroughEdges();
    while(before!=curr){
      before = curr;
      curr = flipThroughEdges();
    }

    computeGrad_intrinsic(Dx, Dy, areas);
    VectorXd Dxu = Dx * UV.col(0);		
		VectorXd Dxv = Dx * UV.col(1);		
		VectorXd Dyu = Dy * UV.col(0);		
		VectorXd Dyv = Dy * UV.col(1);
		MatrixXd RR(F.rows(),4); // each row is the flattened closest rotation matrix 
    for (int i = 0; i < F.rows(); ++i) {
			Matrix2d J, U, S, VV;
			J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
			SSVD2x2(J, U, S, VV);
			Matrix2d R = U*VV.transpose();
			RR(i,0) = R(0,0);
			RR(i,1) = R(0,1);
			RR(i,2) = R(1,0);
			RR(i,3) = R(1,1);
		}		
		b << Dx.transpose() * (areas.asDiagonal() * RR.col(0)) + Dy.transpose() * (areas.asDiagonal() * RR.col(1)),
		Dx.transpose() * (areas.asDiagonal() * RR.col(2)) + Dy.transpose() * (areas.asDiagonal() * RR.col(3));

		SparseMatrix<double> B1 = Dx.transpose()*areas.asDiagonal()*Dx + Dy.transpose()*areas.asDiagonal()*Dy; 

		vector<Triplet<double> > tlist;
		for (int i = 0; i < B1.outerSize(); ++i) {
			for (SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(B1,i); it; ++it) {
				tlist.push_back(Triplet<double>(it.row(), it.col(), it.value()));
				tlist.push_back(Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
		A.setFromTriplets(tlist.begin(), tlist.end());
  }
	// build (A C^t; C 0)
	SparseMatrix<double> Ct, temp1, temp2, res;
	Ct = C.transpose();
	igl::cat(2, A, Ct, temp1);
	C.conservativeResize(C.rows(), C.rows()+C.cols());
	igl::cat(1, temp1, C, res);

	VectorXd rhs(2*V.rows()+C.rows());
	rhs << b,d;
	SparseLU<SparseMatrix<double> > solver;
	solver.analyzePattern(res);
	solver.factorize(res);
	cout << "solver info: " << solver.info() << endl;
	Eigen::VectorXd x = solver.solve(rhs);

	UV.resize(V.rows(),2);
	UV.col(0) = x.segment(0,V.rows());
 	UV.col(1) = x.segment(V.rows(),V.rows());
  cout << "end" << endl;
 	prevFreeBoundary = freeBoundary;

}


void computeParameterization(int type)
{
  cout << "in parameterization" << endl;
  if(igrad && UV.size()!=0) while(greedyflip());
	SparseMatrix<double> A;
	VectorXd b;
	Eigen::SparseMatrix<double> C;
	VectorXd d;
	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices

	// prevFreeBoundary is used so that contraints do not need to be computed each time since they are the same
	// and brute force is expensive
	if(prevFreeBoundary!=freeBoundary||fixed_UV_indices.size()==0){
		if (!freeBoundary)
		{
			// The boundary vertices should be fixed to positions on the unit disc. Find these position and
			// save them in the #V x 2 matrix fixed_UV_position.
			igl::boundary_loop(F, fixed_UV_indices);
			igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
		}
		else{
			// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.
			
			// brute force to find most distant 2 vertices
			unsigned fixed1 = 0, fixed2 = 0;
			for (unsigned i = 0; i < V.rows(); ++i) {
				for (unsigned j = 0; j < V.rows(); ++j) {
					if((V.row(i)-V.row(j)).norm()>(V.row(fixed1)-V.row(fixed2)).norm()){
						fixed1 = i;
						fixed2 = j;
					}
				}
			}
			fixed_UV_indices.resize(2);
			fixed_UV_positions.resize(2,2);
			fixed_UV_indices << fixed1, fixed2;
			fixed_UV_positions << 1,0,0,1;
			cout << fixed_UV_indices.size() << endl;
			cout << fixed_UV_positions << endl;

			//non-deterministic

			// srand(time(0));
			// int iterations = 5;
			// unsigned fixed1 = 0, fixed2 = 0;
			// while(iterations--){
			// 	unsigned index = rand() % V.rows();
			// 	// cout << index << " | " << V.rows() << endl;  
			// 	RowVector3d selected = V.row(index);
			// 	unsigned curr = index;
			// 	for(unsigned i=0; i<V.rows(); ++i){
			// 		if((V.row(i)-selected).norm()>(V.row(curr)-selected).norm()) curr = i;
			// 	}
			// 	// cout << "after for: " << curr << ", fixed1: " << fixed1 << ", fixed2: " << fixed2  << endl;
			// 	if((selected-V.row(curr)).norm()>(V.row(fixed1)-V.row(fixed2)).norm()){
			// 		fixed1 = index;
			// 		fixed2 = curr;
			// 	}
			// 	// cout << "after if"  << endl;
			// }
			// fixed_UV_indices.resize(2);
			// fixed_UV_positions.resize(2,2);
			// fixed_UV_indices << fixed1, fixed2;
			// fixed_UV_positions << 1,0,0,1;
		}
	}
	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d); 


	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	A.resize(2*V.rows(),2*V.rows());
	b.resize(2*V.rows());
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh
		// A = (L 0)
		//     (0 L)
		SparseMatrix<double> L = compute_L_uniform();
		vector<Triplet<double> > tlist;
		cout << "L: " << L.nonZeros() << endl;
		for (int i = 0; i < L.outerSize(); ++i) {
			for (SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(L,i); it; ++it) {
				tlist.push_back(Triplet<double>(it.row(), it.col(), it.value()));
				tlist.push_back(Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
		A.setFromTriplets(tlist.begin(), tlist.end());
		b = VectorXd::Zero(2*V.rows());
	}

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~

		//A = (L 0)

		//    (0 L)  
    cout << " in 2 " << endl;
		SparseMatrix<double> L;
     if(igrad){
      data_mesh.intTri->requireCotanLaplacian();
      L=data_mesh.intTri->cotanLaplacian;
    }
	  else	igl::cotmatrix(V,F,L);
  //  cout << "2: " << L << endl;
		vector<Triplet<double> > tlist;
		for (int i = 0; i < L.outerSize(); ++i) {
			for (SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(L,i); it; ++it) {
				tlist.push_back(Triplet<double>(it.row(), it.col(), -it.value())); // - comes from the warning above 
				tlist.push_back(Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), -it.value()));
			}
		}
		A.setFromTriplets(tlist.begin(), tlist.end());		
		b = VectorXd::Zero(2*V.rows());

	}
	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!
		
		SparseMatrix<double> Dx, Dy;
		VectorXd areas;
    if(igrad){
      computeGrad_intrinsic(Dx, Dy, areas);
    }
    else{
		  computeSurfaceGradientMatrix(Dx,Dy);
	  	igl::doublearea(V,F,areas);
      areas/=2;
    }
		
		// A = (DxADx + DyADy    DxADy - DyADx)
		//     (-DxADy + DyADx   DxADx + DyADy)
		SparseMatrix<double> B1 = Dx.transpose()*areas.asDiagonal()*Dx + Dy.transpose()*areas.asDiagonal()*Dy; 
		SparseMatrix<double> B2 = Dx.transpose()*areas.asDiagonal()*Dy - Dy.transpose()*areas.asDiagonal()*Dx; 
		vector<Triplet<double> > tlist;
		for (int i = 0; i < B1.outerSize(); ++i) {
			for (SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(B1,i); it; ++it) {
				tlist.push_back(Triplet<double>(it.row(), it.col(), it.value()));
				tlist.push_back(Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
  

		for (int i = 0; i < B2.outerSize(); ++i) {
			for (SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(B2,i); it; ++it) {
				tlist.push_back(Triplet<double>(it.row()+V.rows(), it.col(), it.value()));
				tlist.push_back(Triplet<double>(it.row(), it.col()+V.rows(), -it.value()));
			}
		}		
		A.setFromTriplets(tlist.begin(), tlist.end());
		b = VectorXd::Zero(2*V.rows());
	}

	if (type == '4') {
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices
		
		if(UV.size()==0) {
      cout <<  " start UV " << endl;
			computeParameterization('3');
      cout << " total energy: " << compute_total_energy() << endl;
		}
		VectorXd areas;
		SparseMatrix<double> Dx, Dy;
    if(igrad){
      computeGrad_intrinsic(Dx, Dy, areas);
    }
    else{
		  computeSurfaceGradientMatrix(Dx,Dy);
      igl::doublearea(V,F,areas);
		  areas/=2;
    }
		// to compute the SVD from the previous iteration
		VectorXd Dxu = Dx * UV.col(0);		
		VectorXd Dxv = Dx * UV.col(1);		
		VectorXd Dyu = Dy * UV.col(0);		
		VectorXd Dyv = Dy * UV.col(1);
		MatrixXd RR(F.rows(),4); // each row is the flattened closest rotation matrix
    int flipped_triangles = 0;
		for (int i = 0; i < F.rows(); ++i) {
			Matrix2d J, U, S, VV;
			J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
      if(J.determinant()<0) ++flipped_triangles;
			SSVD2x2(J, U, S, VV);
			Matrix2d R = U*VV.transpose();
			RR(i,0) = R(0,0);
			RR(i,1) = R(0,1);
			RR(i,2) = R(1,0);
			RR(i,3) = R(1,1);
		}		
    cout << " flippo: " << flipped_triangles << endl;
		b << Dx.transpose() * (areas.asDiagonal() * RR.col(0)) + Dy.transpose() * (areas.asDiagonal() * RR.col(1)),
		Dx.transpose() * (areas.asDiagonal() * RR.col(2)) + Dy.transpose() * (areas.asDiagonal() * RR.col(3));

		SparseMatrix<double> B1 = Dx.transpose()*areas.asDiagonal()*Dx + Dy.transpose()*areas.asDiagonal()*Dy; 

		vector<Triplet<double> > tlist;
		for (int i = 0; i < B1.outerSize(); ++i) {
			for (SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(B1,i); it; ++it) {
				tlist.push_back(Triplet<double>(it.row(), it.col(), it.value()));
				tlist.push_back(Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
		A.setFromTriplets(tlist.begin(), tlist.end());

	}

  // Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail

	// build (A C^t; C 0)
	SparseMatrix<double> Ct, temp1, temp2, res;
	Ct = C.transpose();
	igl::cat(2, A, Ct, temp1);
	C.conservativeResize(C.rows(), C.rows()+C.cols());
	igl::cat(1, temp1, C, res);

	VectorXd rhs(2*V.rows()+C.rows());
	rhs << b,d;
	SparseLU<SparseMatrix<double> > solver;
	solver.analyzePattern(res);
	solver.factorize(res);
	cout << "solver info: " << solver.info() << endl;
	Eigen::VectorXd x = solver.solve(rhs);

	UV.resize(V.rows(),2);
	UV.col(0) = x.segment(0,V.rows());
 	UV.col(1) = x.segment(V.rows(),V.rows());
  double en = compute_total_energy();
  cout << " total energy: " << en << endl;
  cout << "end" << endl;
 	prevFreeBoundary = freeBoundary;

}

void intrinsicUV(const std::unique_ptr<gcs::IntrinsicTriangulation>& intTri){
  vector<Eigen::Vector2d> points;
  vector<array<int,2> > edges;
  for(gcs::Edge e : intTri->intrinsicMesh->edges()){
    points.push_back(UV.row(e.firstVertex().getIndex()).transpose());
    points.push_back(UV.row(e.secondVertex().getIndex()).transpose());
    int t = points.size();
    std::array<int, 2> tmp{t-2, t-1};
    edges.push_back(tmp);
    
  }
  cout << "before resizing p1 p2" << endl;
  P1.resize(edges.size(),2);
  P2.resize(edges.size(),2);
  int i = 0;
  for(auto e : edges){
    P1.row(i) = points[e[0]].transpose();
    P2.row(i) = points[e[1]].transpose();
    i++;
  }
 
}

void intrinsicEdges(const std::unique_ptr<gcs::IntrinsicTriangulation>& intTri, const Eigen::MatrixXd& V, Eigen::MatrixXd& P1, Eigen::MatrixXd& P2, MatrixXd& colors) {
  vector<Eigen::Vector3d> points;
  vector<array<int,2> > edges;
  /*
  for(gcs::Edge e : data_mesh.inputMesh->edges()){
    points.push_back(V.row(e.firstVertex().getIndex()).transpose());
    points.push_back(V.row(e.secondVertex().getIndex()).transpose());
    int t = points.size();
    std::array<int, 2> tmp{t-2, t-1};
    edges.push_back(tmp);
    
  }
  int ne = edges.size();
  */
  for(gcs::Edge e : intTri->intrinsicMesh->edges()) {
    std::vector<gcs::SurfacePoint> pointVec = intTri->traceIntrinsicHalfedgeAlongInput(e.halfedge());
    //if(pointVec.size()==2) continue;
    for (int k = 0; k < pointVec.size(); k++) {
      points.push_back(V.transpose()*b(pointVec[k]));
      if (k == 0) continue;
      int t = points.size();
      std::array<int, 2> tmp{t-2, t-1};
      edges.push_back(tmp);
    }
  }
  P1.resize(edges.size(),3);
  P2.resize(edges.size(),3);
  /*
  for(int i = 0; i< points.size();++i){
    P.row(i) = points[i].transpose();
  }
  E.resize(edges.size(),2);
  */
  int i = 0;
  for(auto e : edges){
    P1.row(i) = points[e[0]].transpose();
    P2.row(i) = points[e[1]].transpose();
    i++;
  }
  colors.resize(edges.size(),3);
  /*
  for(int j = 0; j < ne; ++j){
    colors.row(j) = Eigen::RowVector3d(0.5,0.5,0);
  }*/
  for(int j = 0; j < edges.size(); ++j){
    colors.row(j) = Eigen::RowVector3d(1,0,0);
  }
}

void intrinsicEdgesUV(const std::unique_ptr<gcs::IntrinsicTriangulation>& intTri, const Eigen::MatrixXd& UV, Eigen::MatrixXd& P1, Eigen::MatrixXd& P2, MatrixXd& colors) {
  vector<Eigen::Vector2d> points;
  vector<array<int,2> > edges;
  /*
  for(gcs::Edge e : data_mesh.inputMesh->edges()){
    points.push_back(V.row(e.firstVertex().getIndex()).transpose());
    points.push_back(V.row(e.secondVertex().getIndex()).transpose());
    int t = points.size();
    std::array<int, 2> tmp{t-2, t-1};
    edges.push_back(tmp);
    
  }
  int ne = edges.size();
  */
  for(gcs::Edge e : intTri->intrinsicMesh->edges()) {
    std::vector<gcs::SurfacePoint> pointVec = intTri->traceIntrinsicHalfedgeAlongInput(e.halfedge());
    //if(pointVec.size()==2) continue;
    for (int k = 0; k < pointVec.size(); k++) {
      points.push_back(UV.transpose()*b(pointVec[k]));
      if (k == 0) continue;
      int t = points.size();
      std::array<int, 2> tmp{t-2, t-1};
      edges.push_back(tmp);
    }
  }
  P1.resize(edges.size(),2);
  P2.resize(edges.size(),2);
  /*
  for(int i = 0; i< points.size();++i){
    P.row(i) = points[i].transpose();
  }
  E.resize(edges.size(),2);
  */
  int i = 0;
  for(auto e : edges){
    P1.row(i) = points[e[0]].transpose();
    P2.row(i) = points[e[1]].transpose();
    i++;
  }
  colors.resize(edges.size(),3);
  /*
  for(int j = 0; j < ne; ++j){
    colors.row(j) = Eigen::RowVector3d(0.5,0.5,0);
  }*/
  for(int j = 0; j < edges.size(); ++j){
    colors.row(j) = Eigen::RowVector3d(1,0,0);
  }
}
bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
  if(option_en==0) energy=drichlet;
  if(option_en==1) energy=symmetricDrichlet;
  if(option_en==2) energy=asap;
  if(option_en==3) energy=arap;

	switch (key) {
  case '1':
	case '2':
	case '3':
	case '4':
	case '5':
		reset=true;
		if(key=='4'){
			unsigned its = iterations;
			while(its--) 
        computeParameterization(key);
		}
    computeParameterization(key);
			// Add your code for detecting and displaying flipped triangles in the
			// UV domain here
		break;
	
	case '6': // conformal
		{
		reset=false;
		cout << "Conformal distortion" << endl;
		SparseMatrix<double> Dx, Dy;
    VectorXd areas;
    if(igrad) computeGrad_intrinsic(Dx,Dy, areas);
    else	computeSurfaceGradientMatrix(Dx,Dy);	
		VectorXd Dxu = Dx * UV.col(0);
		VectorXd Dxv = Dx * UV.col(1);
		VectorXd Dyu = Dy * UV.col(0);
		VectorXd Dyv = Dy * UV.col(1);
		VectorXd distortion(F.rows());
		for (unsigned i = 0; i < F.rows(); ++i) {
			Matrix2d J;
			J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
			distortion(i) = (J + J.transpose() - J.trace()*Matrix2d::Identity()).norm();
			// cout << distortion(i) << endl;
		}
		double minn = distortion.minCoeff();
		double maxx = distortion.maxCoeff();
		cout << "Min distort.: "<<minn << ", Max distort.:" << maxx << endl;
		// linear interpolation for colors
		colors.resize(F.rows(),3);
		for (unsigned i = 0; i < F.rows(); ++i) {
			colors.row(i) = RowVector3d(0,1,1)*(maxx-distortion(i))/(maxx-minn) + RowVector3d(1,0,0);
			// cout << i << ": " << colors.row(i) <<endl;
		}
		viewer.data().set_colors(colors);

		}
		break;
	case '7': // isometric
		{
		reset=false;
		cout << "Isometric distortion" << endl;			
		SparseMatrix<double> Dx, Dy;
    VectorXd areas;
    if(igrad) computeGrad_intrinsic(Dx,Dy, areas);
    else	computeSurfaceGradientMatrix(Dx,Dy);
		VectorXd Dxu = Dx * UV.col(0);
		VectorXd Dxv = Dx * UV.col(1);
		VectorXd Dyu = Dy * UV.col(0);
		VectorXd Dyv = Dy * UV.col(1);
		VectorXd distortion(F.rows());
		for (unsigned i = 0; i < F.rows(); ++i) {
			Matrix2d J, U, S, VV;
			J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
			SSVD2x2(J,U,S,VV);
			distortion(i) = (J - U*VV.transpose()).norm();
		}
		double minn = distortion.minCoeff();
		double maxx = distortion.maxCoeff();
		cout << "Min distort.: "<<minn << ", Max distort.:" << maxx << endl;
		// linear interpolation for colors		
		colors.resize(F.rows(),3);
		for (unsigned i = 0; i < F.rows(); ++i) {
			colors.row(i) = RowVector3d(0,1,1)*(maxx-distortion(i))/(maxx-minn) + RowVector3d(1,0,0);
		}
		viewer.data().set_colors(colors);

		}
		break;
	case '8': // authalic
		{
		reset=false;
		cout << "Authalic distortion" << endl;			
		SparseMatrix<double> Dx, Dy;
    VectorXd areas;
    if(igrad) computeGrad_intrinsic(Dx,Dy, areas);
    else	computeSurfaceGradientMatrix(Dx,Dy);
		VectorXd Dxu = Dx * UV.col(0);
		VectorXd Dxv = Dx * UV.col(1);
		VectorXd Dyu = Dy * UV.col(0);
		VectorXd Dyv = Dy * UV.col(1);
		VectorXd distortion(F.rows());
		for (unsigned i = 0; i < F.rows(); ++i) {
			Matrix2d J;
			J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
			double temp = J(0,0)*J(1,1)-J(1,0)*J(0,1) - 1;
			distortion(i) = temp*temp;
		}
		double minn = distortion.minCoeff();
		double maxx = distortion.maxCoeff();
		cout << "Min distort.: "<<minn << ", Max distort.:" << maxx << endl;
		// linear interpolation for colors		
		colors.resize(F.rows(),3);
		for (unsigned i = 0; i < F.rows(); ++i) {
			colors.row(i) = RowVector3d(0,1,1)*(maxx-distortion(i))/(maxx-minn) + RowVector3d(1,0,0);
		}
		viewer.data().set_colors(colors);
		}
		break;
	case '9': // to reset from the distortion view
		reset=true;
		break;
  case 'f':
    {
      cout << " f " << endl; 
      if(UV.size()==0)
        computeParameterization('2');
      cout << " starting flipping " << endl;
      double total_energy_b = 0;
      for(gcs::Face f : data_mesh.intTri->intrinsicMesh->faces()){
        Matrix2d J;
        faceJacobian(f, J);
        if(J.determinant()<0) continue;
        total_energy_b += energy(J)*data_mesh.intTri->faceAreas[f];
      }
      flipThroughEdges();
      //data_mesh.intTri->flipToDelaunay();
      data_mesh.intTri->refreshQuantities();
      double total_energy = 0;
      for(gcs::Face f : data_mesh.intTri->intrinsicMesh->faces()){
        Matrix2d J;
        faceJacobian(f, J);
        if(J.determinant()<0) continue;
        total_energy += energy(J)*data_mesh.intTri->faceAreas[f];
      }
      cout << "before: " << total_energy_b << endl;
      cout << "after: " << total_energy << endl;
      
    }
    break;
  case 'g':
    {
      cout << " g " << endl; 
      if(UV.size()==0)
        computeParameterization('2');
      cout << " starting flipping " << endl;
      double total_energy_b = 0;
      for(gcs::Face f : data_mesh.intTri->intrinsicMesh->faces()){
        Matrix2d J;
        faceJacobian(f, J);
        if(J.determinant()<0) continue;
        total_energy_b += energy(J)*data_mesh.intTri->faceAreas[f];
      }
      greedyflip();
      data_mesh.intTri->refreshQuantities();
      double total_energy = 0;
      for(gcs::Face f : data_mesh.intTri->intrinsicMesh->faces()){
        Matrix2d J;
        faceJacobian(f, J);
        if(J.determinant()<0) continue;
        total_energy += energy(J)*data_mesh.intTri->faceAreas[f];
      }
      cout << "before: " << total_energy_b << endl;
      cout << "after: " << total_energy << endl;
      
    }
    break;

  case 'v':
    {
      Eigen::MatrixXd C;
      Eigen::MatrixXi E;
      Eigen::MatrixXd colors;
      intrinsicEdges(data_mesh.intTri, V, P1, P2, colors);
      viewer.data().clear();
      viewer.data().set_mesh(V,F);
      C.resize(colored_points.size(),3);
      for (int i=0; i < colored_points.size(); ++i) {
        C.row(i) = colored_points[i].transpose();
      }
      viewer.data().point_size = 12;    
      //viewer.data().add_points(C, Eigen::RowVector3d(0, 0, 1));
      viewer.data().add_edges(P1, P2, colors);
		  viewer.data().set_face_based(false);
      viewer.data().show_faces=true;

      break;
    }
  case 'b':
    {
      if(showingUV){
        viewer.data().clear();
        viewer.data().set_mesh(UV,F);
        cout << " in b " << endl;
        intrinsicUV(data_mesh.intTri);
        viewer.data().add_edges(P1, P2, RowVector3d(1,0,0));
        viewer.data().set_face_based(false);
        viewer.data().show_faces=true;
      }
      break;
    }


  case 't' :
    {
    if(UV.size()==0)
      computeParameterization('2');
    
    double energy_total_before = 0;
    for(gcs::Face f: data_mesh.intTri->intrinsicMesh->faces()){
      gcs::Halfedge he = f.halfedge();
      energy_total_before+= calc_energy(he);
    }
    //flipThroughEdges();
    data_mesh.intTri->flipToDelaunay();

    double energy_total = 0;
    for(gcs::Face f: data_mesh.intTri->intrinsicMesh->faces()){
      gcs::Halfedge he = f.halfedge();
      energy_total += calc_energy(he);
    }
    cout << " before: " << energy_total_before << "  ->   after: " << energy_total << endl;
    break;
  }
  case 'c':
    {
    // compute the total energy w.r.t to the actual UV 

      double total_energy = 0;
      for(gcs::Face f : data_mesh.intTri->intrinsicMesh->faces()){
        Matrix2d J;
        faceJacobian(f, J);
        if(J.determinant()<0){
          //cout << " determinant negative " << endl;
          total_energy += calc_energy(f.halfedge());
        }
        else{
          total_energy += energy(J)*data_mesh.intTri->faceAreas[f];
        }
      }
      cout << " total_energy = " << total_energy << endl;
      break;
    }
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core();
      viewer.core() = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core();
        viewer.core() = temp2D;
        viewer.data().show_texture = false;
		    viewer.data().set_mesh(UV, F);
        showingUV = true;
        viewer.data().set_face_based(false);
        viewer.data().show_faces=true;
        Eigen::MatrixXi E;
        intrinsicEdgesUV(data_mesh.intTri, UV, P1, P2, colors);
        viewer.data().add_edges(P1, P2, RowVector3d(1,0,0));
        for(gcs::Face f: data_mesh.intTri->intrinsicMesh->faces()){
          gcs::Halfedge he = f.halfedge();
          calc_energy(he);
          viewer.data().add_edges(P1,P2,RowVector3d(0,0,1));
        }
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	if(key!='v' && key!=' ' && key!= 'b' ) Redraw();	
	return true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  data_mesh.V=V;
  data_mesh.F=F;
  data_mesh.inputMesh.reset(new gcs::ManifoldSurfaceMesh(F));
  data_mesh.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh.inputMesh, V));
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->requireFaceAreas();
  cout << " Edges: " << data_mesh.intTri->intrinsicMesh->nEdges() << endl;
 
//  testt();
  
  Redraw();
  viewer.core().align_camera_center(V);
  showingUV = false;

  return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core();
	temp2D = viewer.core();
	temp2D.orthographic = true;

	return false;
}

int main(int argc,char *argv[]) {
  if(argc != 2) {
    int reps = 10;
    outfile << "test mode" << endl;
    outfile << " ------------------------- " << endl;
    outfile << " ENERGY: drichlet " << endl;
    energy = drichlet;
    load_mesh("../data/cathead.obj");
    outfile << "****** cathead" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/bunny.off");
    outfile << "****** bunny" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/Octo_cut2.obj");
    outfile << "****** octo" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/stripe1.obj");
    outfile << "****** stripe" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);
    
    outfile << " ENERGY: symmetricDrichlet" << endl;
    energy = symmetricDrichlet;
    load_mesh("../data/cathead.obj");
    outfile << "****** cathead" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/bunny.off");
    outfile << "****** bunny" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/Octo_cut2.obj");
    outfile << "****** octo" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/stripe1.obj");
    outfile << "****** stripe" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    outfile << " ENERGY: arap" << endl;
    energy = arap;
    load_mesh("../data/cathead.obj");
    outfile << "****** cathead" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/bunny.off");
    outfile << "****** bunny" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/Octo_cut2.obj");
    outfile << "****** octo" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/stripe1.obj");
    outfile << "****** stripe" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    outfile << " ENERGY: asap" << endl;
    energy = asap;
    load_mesh("../data/cathead.obj");
    outfile << "****** cathead" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/bunny.off");
    outfile << "****** bunny" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    load_mesh("../data/Octo_cut2.obj");
    outfile << "****** octo" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

   load_mesh("../data/stripe1.obj");
    outfile << "****** stripe" << endl;
    computeParameterization('3');
    for(int i=0;i<reps;i++) flipThroughEdges();
    UV.resize(0,0);

    outfile.close();
    cout << "DONE!!!" << endl;
    return 0;
  }
 
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }

  energy = drichlet;
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  viewer.plugins.push_back(&plugin);
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  plugin.widgets.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parmaterization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);
      ImGui::InputScalar("ARAP iterations", ImGuiDataType_U32, &iterations, 0, 0);
      ImGui::Checkbox("intrinsic grad", &igrad);
      ImGui::RadioButton("drichlet", &option_en, 0); 
      ImGui::RadioButton("symmetric drichlet", &option_en, 1); 
      ImGui::RadioButton("asap", &option_en, 2); 
      ImGui::RadioButton("arap", &option_en, 3);
			// TODO: Add more parameters to tweak here...
		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_init = callback_init;

  viewer.launch();
}

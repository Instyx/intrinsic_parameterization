#include "datageo.hpp"
#include <algorithm>
#include <metrics.hpp>
#include <intrinsicgrad.hpp>
#include <svd.hpp>
#include <igl/internal_angles.h>
#include <igl/internal_angles_intrinsic.h>
#include <igl/doublearea.h>
#include <iostream>
void minmax_distortions(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &UV, std::vector<double> &res){
    Eigen::SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(V,F,Dx,Dy);
    Eigen::VectorXd Dxu = Dx * UV.col(0);
    Eigen::VectorXd Dxv = Dx * UV.col(1);
    Eigen::VectorXd Dyu = Dy * UV.col(0);
    Eigen::VectorXd Dyv = Dy * UV.col(1);
    Eigen::VectorXd distortion_conformal(F.rows());
    Eigen::VectorXd distortion_isometric(F.rows());
    Eigen::VectorXd distortion_authalic(F.rows());
		for (unsigned i = 0; i < F.rows(); ++i) {
      Eigen::Matrix2d J, U, S, VV;
			J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
			distortion_conformal(i) = (J + J.transpose() - J.trace()*Eigen::Matrix2d::Identity()).norm();
      SSVD2x2(J,U,S,VV);
			distortion_isometric(i) = (J - U*VV.transpose()).norm();
      double temp = J(0,0)*J(1,1)-J(1,0)*J(0,1) - 1;
			distortion_authalic(i) = temp*temp;
		}
		double min_conformal = distortion_conformal.minCoeff();
		double max_conformal = distortion_conformal.maxCoeff();
    double min_isometric = distortion_isometric.minCoeff();
		double max_isometric = distortion_isometric.maxCoeff();
    double min_authalic = distortion_authalic.minCoeff();
		double max_authalic = distortion_authalic.maxCoeff();

		// linear interpolation for colors
		/*
    colors.resize(F.rows()*3,3);
		for (unsigned i = 0; i < F.rows(); ++i) {
			colors.row(i) = Eigen::RowVector3d(0,1,1)*(max_conformal - distortion_conformal(i))/(max_conformal - min_conformal) + Eigen::RowVector3d(1,0,0);
			colors.row(i+F.rows()) = Eigen::RowVector3d(0,1,1)*(max_isometric - distortion_isometric(i))/(max_isometric- min_isometric) + Eigen::RowVector3d(1,0,0);
			colors.row(i+2*F.rows()) = Eigen::RowVector3d(0,1,1)*(max_authalic - distortion_authalic(i))/(max_authalic - min_authalic) + Eigen::RowVector3d(1,0,0); 
    }
*/
    res.resize(6);
    res[0]=min_conformal;
    res[1]=max_conformal;
    res[2]=min_isometric;
    res[3]=max_isometric;
    res[4]=min_authalic;
    res[5]=max_authalic;
}

void minmax_distortions_intri(DataGeo &data_mesh, const Eigen::MatrixXd &UV, std::vector<double> &res){
    Eigen::VectorXd areas;
    Eigen::SparseMatrix<double> Dx, Dy;
    Eigen::MatrixXi F = data_mesh.F; 
		computeGrad_intrinsic(data_mesh, Dx, Dy, areas);
    Eigen::VectorXd Dxu = Dx * UV.col(0);
    Eigen::VectorXd Dxv = Dx * UV.col(1);
    Eigen::VectorXd Dyu = Dy * UV.col(0);
    Eigen::VectorXd Dyv = Dy * UV.col(1);
    Eigen::VectorXd distortion_conformal(F.rows());
    Eigen::VectorXd distortion_isometric(F.rows());
    Eigen::VectorXd distortion_authalic(F.rows());
		for (unsigned i = 0; i < F.rows(); ++i) {
      Eigen::Matrix2d J, U, S, VV;
			J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
			distortion_conformal(i) = (J + J.transpose() - J.trace()*Eigen::Matrix2d::Identity()).norm();
      SSVD2x2(J,U,S,VV);
			distortion_isometric(i) = (J - U*VV.transpose()).norm();
      double temp = J(0,0)*J(1,1)-J(1,0)*J(0,1) - 1;
			distortion_authalic(i) = temp*temp;
		}
		double min_conformal = distortion_conformal.minCoeff();
		double max_conformal = distortion_conformal.maxCoeff();
    double min_isometric = distortion_isometric.minCoeff();
		double max_isometric = distortion_isometric.maxCoeff();
    double min_authalic = distortion_authalic.minCoeff();
		double max_authalic = distortion_authalic.maxCoeff();

    res.resize(6);
    res[0]=min_conformal;
    res[1]=max_conformal;
    res[2]=min_isometric;
    res[3]=max_isometric;
    res[4]=min_authalic;
    res[5]=max_authalic;
}




void angle_distortion(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &UV, Eigen::VectorXd &angle_errors, Eigen::VectorXd &angle_dists){
  Eigen::MatrixXd K_V, K_UV;
  Eigen::SparseMatrix<double> Dx, Dy;
  //Eigen::VectorXd areas;
	
  computeSurfaceGradientMatrix(V,F,Dx,Dy);
  //igl::doublearea(V,F,areas);
  //areas/=2;

  igl::internal_angles(V,F,K_V);
  igl::internal_angles(UV,F,K_UV);
  
  Eigen::MatrixXd angle_diffs = (K_UV - K_V).cwiseAbs();
  angle_errors.resize(F.rows());
  angle_errors << angle_diffs.col(0) + angle_diffs.col(1) + angle_diffs.col(2);
  
  //double average_angle_error = areas.cwiseProduct(angle_errors).sum();

  Eigen::VectorXd Dxu = Dx * UV.col(0);
  Eigen::VectorXd Dxv = Dx * UV.col(1);
  Eigen::VectorXd Dyu = Dy * UV.col(0);
  Eigen::VectorXd Dyv = Dy * UV.col(1);
  Eigen::VectorXd angle_dist1(F.rows());
  Eigen::VectorXd angle_dist2(F.rows());
  for (unsigned i = 0; i < F.rows(); ++i) {
    Eigen::Matrix2d J, U, S, VV;
    //J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
    //std::cout << J << std::endl;
    SSVD2x2(J,U,S,VV);
    if(S(1,1)!=0)
      angle_dist1(i) = S(0,0) / S(1,1); 
    else 
      angle_dist1(i) = 1e7;
    if(S(0,0)!=0)
      angle_dist2(i) = S(1,1) / S(0,0);
    else
      angle_dist2(i) = 1e7;
  }
  angle_dists.resize(F.rows());
  angle_dists << angle_dist1.cwiseAbs() + angle_dist2.cwiseAbs();
  //std::cout << angle_dists << std::endl;
}

// is bugged
void angle_distortion_intri(DataGeo &data_mesh, const Eigen::MatrixXd &UV, Eigen::VectorXd &angle_errors, Eigen::VectorXd &angle_dists){

  Eigen::MatrixXi F = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
  Eigen::MatrixXd lengths(F.rows(),3);
  for(gcs::Face face : data_mesh.intTri->intrinsicMesh->faces()){
    size_t idx = data_mesh.intTri->faceIndices[face];
    size_t v0 = F(idx,0); 
    size_t v1 = F(idx,1); 
    size_t v2 = F(idx,2); 
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
  }
  Eigen::MatrixXd sq_lengths;
  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;
	computeGrad_intrinsic(data_mesh,Dx,Dy, areas);
 
  sq_lengths = lengths.cwiseProduct(lengths);
  Eigen::MatrixXd K_V, K_UV;
  igl::internal_angles_intrinsic(sq_lengths,K_V);
  igl::internal_angles(UV, F, K_UV);
  Eigen::MatrixXd angle_diffs = (K_UV - K_V).cwiseAbs();
  angle_errors.resize(F.rows());
  angle_errors << angle_diffs.col(0) + angle_diffs.col(1) + angle_diffs.col(2);

  Eigen::VectorXd Dxu = Dx * UV.col(0);
  Eigen::VectorXd Dxv = Dx * UV.col(1);
  Eigen::VectorXd Dyu = Dy * UV.col(0);
  Eigen::VectorXd Dyv = Dy * UV.col(1);
  Eigen::VectorXd angle_dist1(F.rows());
  Eigen::VectorXd angle_dist2(F.rows());
  for (unsigned i = 0; i < F.rows(); ++i) {
    Eigen::Matrix2d J, U, S, VV;
    J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
    SSVD2x2(J,U,S,VV);
    if(S(1,1)!=0)
      angle_dist1(i) = S(0,0) / S(1,1); 
    else 
      angle_dist1(i) = 1e7;
    if(S(0,0)!=0)
      angle_dist2(i) = S(1,1) / S(0,0);
    else
      angle_dist2(i) = 1e7;
   // std::cout << S(0) << "  " << S(3) << std::endl; 
  }
  angle_dists.resize(F.rows());
  angle_dists << angle_dist1.cwiseAbs() + angle_dist2.cwiseAbs();
 // std::cout << angle_dists << std:endl;
}


// returns percentage flipped
double flipped_triangles(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &UV){
  Eigen::SparseMatrix<double> Dx, Dy;
	computeSurfaceGradientMatrix(V,F,Dx,Dy);
  Eigen::VectorXd Dxu = Dx * UV.col(0);
  Eigen::VectorXd Dxv = Dx * UV.col(1);
  Eigen::VectorXd Dyu = Dy * UV.col(0);
  Eigen::VectorXd Dyv = Dy * UV.col(1);
  unsigned flipped = 0;
  
  for (unsigned i = 0; i < F.rows(); ++i) {
    Eigen::Matrix2d J, U, S, VV;
    J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
    if(J.determinant() < 0){
      ++flipped;
    }
  }
  return (double)flipped/(double)F.rows();
}

double flipped_triangles_intri(DataGeo &data_mesh, const Eigen::MatrixXd &UV){
  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;
	computeGrad_intrinsic(data_mesh,Dx,Dy, areas);

  Eigen::VectorXd Dxu = Dx * UV.col(0);
  Eigen::VectorXd Dxv = Dx * UV.col(1);
  Eigen::VectorXd Dyu = Dy * UV.col(0);
  Eigen::VectorXd Dyv = Dy * UV.col(1);
  unsigned flipped = 0;
  
  for (unsigned i = 0; i < data_mesh.F.rows(); ++i) {
    Eigen::Matrix2d J, U, S, VV;
    J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
    if(J.determinant() < 0){
      ++flipped;
    }
  }
  return (double)flipped/(double)data_mesh.F.rows();
}



void compute_metrics(DataGeo &data_mesh, const Eigen::MatrixXd &UV_o, std::vector<double> &res){
 
  // extrinsic metrics
  Eigen::MatrixXd V = data_mesh.V;
  Eigen::MatrixXi F = data_mesh.F;
  Eigen::MatrixXd UV = UV_o;

  Eigen::VectorXd extrinsic_areas, extrinsic_areasUV;
  igl::doublearea(V,F,extrinsic_areas);
  extrinsic_areas/=2;
  igl::doublearea(UV,F,extrinsic_areasUV);
  extrinsic_areasUV = extrinsic_areasUV.cwiseAbs();
  extrinsic_areasUV/=2;
  // std::cout << extrinsic_areasUV << std::endl;
  // normalize
  V = V / std::sqrt(extrinsic_areas.sum());
  UV = UV / std::sqrt(extrinsic_areasUV.sum());
  extrinsic_areas/=extrinsic_areas.sum(); 
  extrinsic_areasUV/=extrinsic_areasUV.sum(); 
  
  // compute metrics
  double extrinsic_flipped = flipped_triangles(V,F,UV);
  
  Eigen::VectorXd extrinsic_angle_dists, extrinsic_angle_errors;
  angle_distortion(V, F, UV, extrinsic_angle_errors, extrinsic_angle_dists);
  double extrinsic_average_angle_error = extrinsic_areas.cwiseProduct(extrinsic_angle_errors).sum();
  double extrinsic_max_angle_dist = extrinsic_angle_dists.maxCoeff()-2;

  Eigen::VectorXd extrinsic_area_dists = extrinsic_areasUV.cwiseQuotient(extrinsic_areas) + extrinsic_areas.cwiseQuotient(extrinsic_areasUV);;
  double extrinsic_max_area_dist = extrinsic_area_dists.maxCoeff() - 2;
  Eigen::VectorXd extrinsic_area_errors = (extrinsic_areasUV - extrinsic_areas).cwiseAbs();
  double extrinsic_average_area_error = extrinsic_area_errors.sum();

  // intrinsic metrics
  F = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>(); 
  UV = UV_o;
  Eigen::VectorXd intrinsic_areas(F.rows());
  int i = 0;
  for(gcs::Face f : data_mesh.intTri->intrinsicMesh->faces()){
    intrinsic_areas(i) =  data_mesh.intTri->faceArea(f);
    ++i;
  }
  Eigen::VectorXd intrinsic_areasUV;
  igl::doublearea(UV,F,intrinsic_areasUV);
  intrinsic_areasUV = intrinsic_areasUV.cwiseAbs();
  intrinsic_areasUV/=2;
  UV = UV / std::sqrt(intrinsic_areasUV.sum());

  intrinsic_areas /= intrinsic_areas.sum();
  intrinsic_areasUV /= intrinsic_areasUV.sum();

  double intrinsic_flipped = flipped_triangles_intri(data_mesh, UV);
  Eigen::VectorXd intrinsic_angle_dists, intrinsic_angle_errors;
  angle_distortion_intri(data_mesh, UV, intrinsic_angle_errors, intrinsic_angle_dists);
  double intrinsic_average_angle_error = intrinsic_areas.cwiseProduct(intrinsic_angle_errors).sum();
  double intrinsic_max_angle_dist = intrinsic_angle_dists.maxCoeff() - 2;

  Eigen::VectorXd intrinsic_area_dists = intrinsic_areasUV.cwiseQuotient(intrinsic_areas) + intrinsic_areas.cwiseQuotient(intrinsic_areasUV);;
  double intrinsic_max_area_dist = intrinsic_area_dists.maxCoeff() - 2;
  Eigen::VectorXd intrinsic_area_errors = (extrinsic_areasUV - extrinsic_areas).cwiseAbs();
  double intrinsic_average_area_error = intrinsic_area_errors.sum();
 
  res = {extrinsic_flipped, extrinsic_max_area_dist, extrinsic_average_area_error, extrinsic_max_angle_dist, extrinsic_average_angle_error,
          intrinsic_flipped, intrinsic_max_area_dist, intrinsic_average_area_error, intrinsic_max_angle_dist, intrinsic_average_angle_error}; 
}

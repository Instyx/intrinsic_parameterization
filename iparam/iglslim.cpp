#include <iglslim.hpp>
#include <intrinsicgrad.hpp>
#include <igl/doublearea.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <distortion_energy.hpp>
void boundary(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, bool isFreeBoundary, Eigen::VectorXi &fixed_UV_indices, Eigen::MatrixXd &fixed_UV_positions) {
  if (!isFreeBoundary) {
    // map boundary to the unit circle
    igl::boundary_loop(F, fixed_UV_indices);
    igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
  }
  else {
    // fix the first vertex 
    fixed_UV_indices.resize(1);
    fixed_UV_positions.resize(1,2);
    fixed_UV_indices << 0;
    fixed_UV_positions << 0,0;
  }

}
double compute_total_energy(DataGeo &data_mesh, const Eigen::MatrixXd &UV){
  
  auto energy = symmetric_dirichlet;
  Eigen::SparseMatrix<double> Dx, Dy;
  Eigen::VectorXd areas;
  if(false){
      computeGrad_intrinsic(data_mesh, Dx, Dy, areas);
    }
  else{
    computeSurfaceGradientMatrix(data_mesh.V, data_mesh.F, Dx,Dy);
    igl::doublearea(data_mesh.V,data_mesh.F,areas);
    areas/=2;
  }
  Eigen::VectorXd Dxu = Dx * UV.col(0);		
  Eigen::VectorXd Dxv = Dx * UV.col(1);		
  Eigen::VectorXd Dyu = Dy * UV.col(0);		
  Eigen::VectorXd Dyv = Dy * UV.col(1);

  double total_energy = 0;
  unsigned flipped_triangles = 0;
  double total_area = 0;
  for(int i=0; i<data_mesh.intTri->intrinsicMesh->nFaces(); ++i){
    Eigen::Matrix2d J;
		J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
    if(J.determinant()<0){
      flipped_triangles++;
    }
    double locen = energy(J)*areas(i);
    total_area+=areas(i);
    total_energy += locen; 
  }
  return total_energy/total_area;
}

void slim_parameterization(DataGeo &data_mesh, igl::SLIMData &slimdata, Eigen::MatrixXd &UV, bool igrad, bool isFreeBoundary){

  Eigen::MatrixXd V = data_mesh.V;
  Eigen::MatrixXi F = data_mesh.F;
 
  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;

  if(igrad){
    computeGrad_intrinsic(data_mesh, Dx, Dy, areas);
  }
  else{
    computeSurfaceGradientMatrix(V,F,Dx,Dy);
    igl::doublearea(V,F,areas);
    areas/=2;
  }
  slimdata.Dx = Dx;
  slimdata.Dy = Dy;
  slimdata.M = areas;
  slimdata.mesh_area = slimdata.M.sum();

  if(igrad) slimdata.F = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>(); 

  if(!slimdata.has_pre_calc){
    slimdata.V = V;
    if(!igrad) slimdata.F = F;

    slimdata.V_o = UV;

    slimdata.v_num = V.rows();
    slimdata.f_num = F.rows(); 
    slimdata.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;


    Eigen::VectorXi fixed_UV_indices;
    Eigen::MatrixXd fixed_UV_positions;

    boundary(V, F, isFreeBoundary, fixed_UV_indices, fixed_UV_positions);

    slimdata.b = fixed_UV_indices;
    slimdata.bc = fixed_UV_positions;
    slimdata.soft_const_p = 0;

    slimdata.proximal_p = 0.0001;
    slimdata.mesh_improvement_3d = false; // whether to use a jacobian derived from a real mesh or an abstract regular mesh (used for mesh improvement)

    slimdata.dim = 2;
    slimdata.v_n = slimdata.v_num;
    slimdata.f_n = slimdata.f_num;

    slimdata.W.resize(slimdata.f_n, slimdata.dim * slimdata.dim);
    slimdata.Ri.resize(slimdata.f_n, slimdata.dim * slimdata.dim);
    slimdata.Ji.resize(slimdata.f_n, slimdata.dim * slimdata.dim);
    slimdata.rhs.resize(slimdata.dim * slimdata.v_num);
    // flattened weight matrix
    slimdata.WGL_M.resize(slimdata.dim * slimdata.dim * slimdata.f_n);
    for (int i = 0; i < slimdata.dim * slimdata.dim; i++)
      for (int j = 0; j < slimdata.f_n; j++)
        slimdata.WGL_M(i * slimdata.f_n + j) = slimdata.M(j);
    slimdata.energy = compute_total_energy(data_mesh, UV);
    slimdata.first_solve = true;
    slimdata.has_pre_calc = true;
  }
  
  Eigen::MatrixXd newUV = igl::slim_solve(slimdata,1);
  UV = newUV;
}

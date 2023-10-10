#include <iglslim.hpp>
#include <intrinsicgrad.hpp>
#include <igl/doublearea.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <distortion_energy.hpp>

void boundary(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, bool isFreeBoundary, Eigen::VectorXi &fixed_UV_indices, Eigen::MatrixXd &fixed_UV_positions) {
  if (!isFreeBoundary) {
    // The boundary vertices should be fixed to positions on the unit disc. Find these position and
    // save them in the #V x 2 matrix fixed_UV_position.
    igl::boundary_loop(F, fixed_UV_indices);
    igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
  }
  else {
    
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
  }

}

void slim_parameterization(DataGeo &data_mesh, igl::SLIMData &slimdata, Eigen::MatrixXd &UV, bool igrad, bool isFreeBoundary, const FlipType &ft){

  auto flip_func = edgeorder_flip;
  if(ft == FlipType::EDGEORDER) flip_func = edgeorder_flip;
  if(ft == FlipType::GREEDY) flip_func = greedy_flip;
  if(ft == FlipType::RANDOM) flip_func = random_flip;
  if(ft == FlipType::HEURISTIC) flip_func = heuristic_flip;
  
  Eigen::MatrixXd V = data_mesh.V;
  Eigen::MatrixXi F = data_mesh.F;
 
  if(igrad) {
    unsigned flips = flip_func(data_mesh, UV, symmetric_dirichlet);
    while(flips){
      flips = flip_func(data_mesh, UV, symmetric_dirichlet);
    }
  }


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


  if(!slimdata.has_pre_calc){
    slimdata.V = V;
    if(igrad)
      slimdata.F = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
    else
      slimdata.F = F;

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

    slimdata.first_solve = true;
    slimdata.has_pre_calc = true;
  }
  
  
  
  Eigen::MatrixXd newUV = igl::slim_solve(slimdata,1);
  UV = newUV;
}

#include "datageo.hpp"
#include <parameterization.hpp>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/cat.h>
#include <igl/adjacency_list.h>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <math.h>
#include <Eigen/SparseLU>
#include <intrinsicgrad.hpp>
#include <intrinsicflip.hpp>
#include <distortion_energy.hpp>
#include <svd.hpp>
#include <stdio.h>
#include <chrono>

auto start = std::chrono::high_resolution_clock::now();

auto end = std::chrono::high_resolution_clock::now();

// global variables for the boundary conditions 
Eigen::SparseMatrix<double> C;
Eigen::VectorXd d;
  
void reset_constraints(){
  C.resize(0,0);
  d.resize(0);
}


Eigen::SparseMatrix<double> compute_L_uniform(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F){
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


void ConvertConstraintsToMatrixForm(Eigen::VectorXi indices, Eigen::MatrixXd positions, unsigned nvertices, Eigen::SparseMatrix<double> &C, Eigen::VectorXd &d) {
	
	// u coordinates
  std::vector<Eigen::Triplet<double> > list;
	d = Eigen::VectorXd::Zero(2*indices.size());
	for (int i = 0; i < indices.size(); ++i) {
		list.push_back(Eigen::Triplet<double>(i,indices(i),1)); 
		d(i) = positions(i,0);
	}
	// v coordinates
	for (int i = 0; i < indices.size(); ++i) {
		list.push_back(Eigen::Triplet<double>(i+indices.size(),indices(i)+nvertices,1));
		d(i+indices.size()) = positions(i,1);
	}
	C.resize(2*indices.size(), 2*nvertices);
	C.setFromTriplets(list.begin(), list.end());

}

void computeConstraints(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, bool isFreeBoundary, int type, Eigen::SparseMatrix<double> &C, Eigen::VectorXd &d) {
  Eigen::VectorXi fixed_UV_indices;
  Eigen::MatrixXd fixed_UV_positions;
  if (!isFreeBoundary) {
    // the boundary is fixed to unit circle 
    igl::boundary_loop(F, fixed_UV_indices);
    igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
  }
  else {
    // brute force to find most distant 2 vertices on the boundary
    if(type=='3'){
      igl::boundary_loop(F, fixed_UV_indices);

      unsigned fixed1 = fixed_UV_indices[0], fixed2 = fixed_UV_indices[0];
      for (unsigned i = 0; i < fixed_UV_indices.rows(); ++i) {
        for (unsigned j = 0; j < fixed_UV_indices.rows(); ++j) {
          if((V.row(fixed_UV_indices[i])-V.row(fixed_UV_indices[j])).norm()>(V.row(fixed1)-V.row(fixed2)).norm()){
            fixed1 = fixed_UV_indices[i];
            fixed2 = fixed_UV_indices[j];
            //fixed1 = i;
            //fixed2 = j;
          }
        }
      }
      fixed_UV_indices.resize(2);
      fixed_UV_positions.resize(2,2);
      fixed_UV_indices << fixed1, fixed2;
      fixed_UV_positions << 1,0,0,1;
    }
    else if(type=='4'){
      // fix the first vertex
      fixed_UV_indices.resize(1);
      fixed_UV_positions.resize(1,2);
      fixed_UV_indices << 0;
      fixed_UV_positions << 0,0;
    }
    else {
      std::cout << "Free boudnary not possible" << std::endl;
    }
  }
	
	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions,V.rows(), C, d); 
}

void computeParameterization(DataGeo &data_mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &UV, Eigen::MatrixXd &new_UV,
    bool isFreeBoundary, bool igrad, int type) {

  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
	if(d.size()==0)  computeConstraints(V, F, isFreeBoundary, type, C, d);

	A.resize(2*V.rows(),2*V.rows());
	b.resize(2*V.rows());
  // Tutte - uniform Laplacian
	if (type == '1') {
		// A = (L 0)
		//     (0 L)
    Eigen::SparseMatrix<double> L;
    if(igrad){
      Eigen::MatrixXi F_new = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
      L = compute_L_uniform(V, F_new);
    }
    else L = compute_L_uniform(V, F);
    std::vector<Eigen::Triplet<double> > tlist;
		for (int i = 0; i < L.outerSize(); ++i) {
			for (Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(L,i); it; ++it) {
				tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
				tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
		A.setFromTriplets(tlist.begin(), tlist.end());
		b = Eigen::VectorXd::Zero(2*V.rows());
	}

  // Harmonic - cotangent Laplacian
	if (type == '2') {
		//A = (L 0)
		//    (0 L)  
    Eigen::SparseMatrix<double> L;
    if(igrad){
      data_mesh.intTri->requireCotanLaplacian();
      L=data_mesh.intTri->cotanLaplacian;
    }
	  else igl::cotmatrix(V,F,L);

    std::vector<Eigen::Triplet<double> > tlist;
		for (int i = 0; i < L.outerSize(); ++i) {
			for (Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(L,i); it; ++it) {
				tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), -it.value())); // - comes from the warning above 
				tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), -it.value()));
			}
		}
		A.setFromTriplets(tlist.begin(), tlist.end());		
		b = Eigen::VectorXd::Zero(2*V.rows());

	}
  // LSCM
	if (type == '3') {
    Eigen::SparseMatrix<double> Dx, Dy;
    Eigen::VectorXd areas;
    if(igrad){
      computeGrad_intrinsic(data_mesh, Dx, Dy, areas);
    }
    else{
		  computeSurfaceGradientMatrix(V, F, Dx, Dy);
	  	igl::doublearea(V,F,areas);
      areas/=2;
    }
		
		// A = (DxADx + DyADy  -DxADy + DyADx)
		//     (DxADy - DyADx   DxADx + DyADy)
    Eigen::SparseMatrix<double> B1 = Dx.transpose()*areas.asDiagonal()*Dx + Dy.transpose()*areas.asDiagonal()*Dy; 
    Eigen::SparseMatrix<double> B2 = -Dx.transpose()*areas.asDiagonal()*Dy + Dy.transpose()*areas.asDiagonal()*Dx; 
    std::vector<Eigen::Triplet<double> > tlist;
		for (int i = 0; i < B1.outerSize(); ++i) {
			for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(B1,i); it; ++it) {
				tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
				tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
  
		for (int i = 0; i < B2.outerSize(); ++i) {
			for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(B2,i); it; ++it) {
				tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col(), it.value()));
				tlist.push_back(Eigen::Triplet<double>(it.row(), it.col()+V.rows(), -it.value()));
			}
		}		
		A.setFromTriplets(tlist.begin(), tlist.end());
		b = Eigen::VectorXd::Zero(2*V.rows());
	}
  // ARAP
	if (type == '4') {
    Eigen::VectorXd areas;
    Eigen::SparseMatrix<double> Dx, Dy;
    if(igrad){
      computeGrad_intrinsic(data_mesh, Dx, Dy, areas);
    }
    else{
		  computeSurfaceGradientMatrix(V, F, Dx,Dy);
      igl::doublearea(V,F,areas);
		  areas/=2;
    }
		// to compute the SVD from the previous iteration
    Eigen::VectorXd Dxu = Dx * UV.col(0);		
    Eigen::VectorXd Dxv = Dx * UV.col(1);		
    Eigen::VectorXd Dyu = Dy * UV.col(0);		
    Eigen::VectorXd Dyv = Dy * UV.col(1);
    Eigen::MatrixXd RR(F.rows(),4); // each row is the flattened closest rotation matrix
    int flipped_triangles = 0;
	  // local step
    for (int i = 0; i < F.rows(); ++i) {
      Eigen::Matrix2d J, U, S, VV;
			J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
      if(J.determinant()<0) ++flipped_triangles;
			SSVD2x2(J, U, S, VV);
      Eigen::Matrix2d R = U*VV.transpose();
			RR(i,0) = R(0,0);
			RR(i,1) = R(0,1);
			RR(i,2) = R(1,0);
			RR(i,3) = R(1,1);
		}		
		b << Dx.transpose() * (areas.asDiagonal() * RR.col(0)) + Dy.transpose() * (areas.asDiagonal() * RR.col(1)),
		Dx.transpose() * (areas.asDiagonal() * RR.col(2)) + Dy.transpose() * (areas.asDiagonal() * RR.col(3));

    Eigen::SparseMatrix<double> B1 = Dx.transpose()*areas.asDiagonal()*Dx + Dy.transpose()*areas.asDiagonal()*Dy; 

    std::vector<Eigen::Triplet<double> > tlist;
		for (int i = 0; i < B1.outerSize(); ++i) {
			for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(B1,i); it; ++it) {
				tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
				tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
			}
		}
		A.setFromTriplets(tlist.begin(), tlist.end());

	}

	// build (A C^t; C 0)
  Eigen::SparseMatrix<double> Ct, temp1, temp2, res;
  Eigen::SparseMatrix<double> zeros(C.rows(), C.rows()); 
	Ct = C.transpose();
	igl::cat(2, A, Ct, temp1);
	igl::cat(2, C, zeros, temp2);
	igl::cat(1, temp1, temp2, res);

  Eigen::VectorXd rhs(2*V.rows()+C.rows());
	rhs << b,d;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(res);
	solver.factorize(res);
	Eigen::VectorXd x = solver.solve(rhs);
  std::cout << " solver " << solver.info() << std::endl;
	new_UV.resize(V.rows(),2);
	new_UV.col(0) = x.segment(0,V.rows());
 	new_UV.col(1) = x.segment(V.rows(),V.rows());

}


unsigned ARAP_tillconverges(DataGeo &data_mesh, Eigen::MatrixXd &UV_init, Eigen::MatrixXd &UV, unsigned max_iterations, bool isFreeBoundary, bool igrad){
  Eigen::MatrixXd V = data_mesh.V;
  Eigen::MatrixXi F = data_mesh.F;

  computeConstraints(V, F, isFreeBoundary, '4', C, d);

  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;

	A.resize(2*V.rows(),2*V.rows());
	b.resize(2*V.rows());
  
  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;
  if(igrad){
    computeGrad_intrinsic(data_mesh, Dx, Dy, areas);
  }
  else{
    computeSurfaceGradientMatrix(V, F, Dx,Dy);
    igl::doublearea(V,F,areas);
    areas/=2;
  }

  Eigen::SparseMatrix<double> B1 = Dx.transpose()*areas.asDiagonal()*Dx + Dy.transpose()*areas.asDiagonal()*Dy; 

  std::vector<Eigen::Triplet<double> > tlist;
  for (int i = 0; i < B1.outerSize(); ++i) {
    for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(B1,i); it; ++it) {
      tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
      tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
    }
  }
  A.setFromTriplets(tlist.begin(), tlist.end());

  Eigen::SparseMatrix<double> Ct, temp1, temp2, res;
  Eigen::SparseMatrix<double> zeros(C.rows(), C.rows()); 
	Ct = C.transpose();
	igl::cat(2, A, Ct, temp1);
	igl::cat(2, C, zeros, temp2);
	igl::cat(1, temp1, temp2, res);

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(res);
	solver.factorize(res);

  UV = UV_init;

  double tol = 1e-8;
  double past_energy = 0;
  double curr_energy = compute_total_energy(data_mesh, UV, EnergyType::ARAP, igrad);
  std::cout << " Initial energy: " << curr_energy << std::endl;
  unsigned itr = 0;
  while(itr<max_iterations && std::abs(past_energy-curr_energy)>tol){
    // to compute the SVD from the previous iteration
    Eigen::VectorXd Dxu = Dx * UV.col(0);		
    Eigen::VectorXd Dxv = Dx * UV.col(1);		
    Eigen::VectorXd Dyu = Dy * UV.col(0);		
    Eigen::VectorXd Dyv = Dy * UV.col(1);
    Eigen::MatrixXd RR(F.rows(),4); // each row is the flattened closest rotation matrix
    // local step
    for (int i = 0; i < F.rows(); ++i) {
      Eigen::Matrix2d J, U, S, VV;
      J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
      SSVD2x2(J, U, S, VV);
      Eigen::Matrix2d R = U*VV.transpose();
      RR(i,0) = R(0,0);
      RR(i,1) = R(0,1);
      RR(i,2) = R(1,0);
      RR(i,3) = R(1,1);
    }		
    b << Dx.transpose() * (areas.asDiagonal() * RR.col(0)) + Dy.transpose() * (areas.asDiagonal() * RR.col(1)),
    Dx.transpose() * (areas.asDiagonal() * RR.col(2)) + Dy.transpose() * (areas.asDiagonal() * RR.col(3));
   
    //global step
    Eigen::VectorXd rhs(2*V.rows()+C.rows());
  	rhs << b,d;
    Eigen::VectorXd x = solver.solve(rhs);
    std::cout << " solver " << solver.info() << std::endl;
	  UV.col(0) = x.segment(0,V.rows());
 	  UV.col(1) = x.segment(V.rows(),V.rows());

    past_energy = curr_energy;
    curr_energy = compute_total_energy_fast(data_mesh, UV, Dx, Dy, areas, EnergyType::ARAP);
    ++itr;
    std::cout << " Energy itr. " << itr << " : " << curr_energy << std::endl;
  }
  return itr;
}


unsigned intrinsic_ARAP(DataGeo &data_mesh, Eigen::MatrixXd &UV, unsigned ARAP_maxitr, unsigned intrinsic_maxitr, bool isFreeBoundary, std::fstream &fout){
  Eigen::MatrixXd UV_init;
  
  reset_constraints();
  UV_init = tutte(data_mesh, false);
  
  double past_energy = compute_total_energy(data_mesh, UV_init, EnergyType::ARAP, false);
  fout << past_energy;

  double tol = 1e-8;
  
  unsigned max_iterations = intrinsic_maxitr;
  unsigned itr = 0;
  unsigned total_iterations;
  // first with extrinsic geometry
  start = std::chrono::high_resolution_clock::now();
  total_iterations = ARAP_tillconverges(data_mesh, UV_init, UV, ARAP_maxitr, isFreeBoundary, false);
  end = std::chrono::high_resolution_clock::now();
  
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  
  fout << total_iterations << duration;  

  double curr_energy = compute_total_energy(data_mesh, UV, EnergyType::ARAP, false);
  fout << curr_energy;


  while(itr < max_iterations && std::abs(past_energy-curr_energy)>tol){
    std::cout << "Intrinsic itr. " << itr << ":" << std::endl; 
    // intrinsic flipping
    unsigned flips, del_flips;
    unsigned total_flips = 0;
    unsigned total_del_flips = 0;

    start = std::chrono::high_resolution_clock::now();
    while(flips = edgeorder_flip(data_mesh, UV, del_flips, EnergyType::ARAP)){
      total_flips += flips;
      total_del_flips += del_flips;
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    fout << total_flips << total_del_flips << duration;
    fout << compute_total_energy(data_mesh, UV, EnergyType::ARAP, true) << std::endl; 

    // compute parameterization with new intrinsic geometry
    Eigen::MatrixXd new_UV;

    start = std::chrono::high_resolution_clock::now();
    total_iterations = ARAP_tillconverges(data_mesh, UV, new_UV, ARAP_maxitr ,isFreeBoundary, true);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    fout << total_iterations << duration;

    UV = new_UV;
    past_energy = curr_energy;
    curr_energy = compute_total_energy(data_mesh, UV, EnergyType::ARAP, true);
    fout << curr_energy;
    ++itr;
  }
  return itr;
}

Eigen::MatrixXd LSCM(DataGeo &data_mesh, bool isFreeBoundary, bool igrad){
  Eigen::MatrixXd V = data_mesh.V;
  Eigen::MatrixXi F = data_mesh.F;

  computeConstraints(V, F, isFreeBoundary, '3', C, d);

  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;

	A.resize(2*V.rows(),2*V.rows());
	b.resize(2*V.rows());
  
  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;
  if(igrad){
    computeGrad_intrinsic(data_mesh, Dx, Dy, areas);
  }
  else{
    computeSurfaceGradientMatrix(V, F, Dx,Dy);
    igl::doublearea(V,F,areas);
    areas/=2;
  }

  Eigen::SparseMatrix<double> B1 = Dx.transpose()*areas.asDiagonal()*Dx + Dy.transpose()*areas.asDiagonal()*Dy; 
  Eigen::SparseMatrix<double> B2 = -Dx.transpose()*areas.asDiagonal()*Dy + Dy.transpose()*areas.asDiagonal()*Dx; 
  std::vector<Eigen::Triplet<double> > tlist;
  for (int i = 0; i < B1.outerSize(); ++i) {
    for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(B1,i); it; ++it) {
      tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
      tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
    }
  }

  for (int i = 0; i < B2.outerSize(); ++i) {
    for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(B2,i); it; ++it) {
      tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col(), it.value()));
      tlist.push_back(Eigen::Triplet<double>(it.row(), it.col()+V.rows(), -it.value()));
    }
  }		
  A.setFromTriplets(tlist.begin(), tlist.end());
  b = Eigen::VectorXd::Zero(2*V.rows());

  Eigen::SparseMatrix<double> Ct, temp1, temp2, res;
  Eigen::SparseMatrix<double> zeros(C.rows(), C.rows()); 
	Ct = C.transpose();
	igl::cat(2, A, Ct, temp1);
	igl::cat(2, C, zeros, temp2);
	igl::cat(1, temp1, temp2, res);

  Eigen::VectorXd rhs(2*V.rows()+C.rows());
	rhs << b,d;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(res);
	solver.factorize(res);
	Eigen::VectorXd x = solver.solve(rhs);
  std::cout << " solver " << solver.info() << std::endl;

  Eigen::MatrixXd UV;
	UV.resize(V.rows(),2);
	UV.col(0) = x.segment(0,V.rows());
 	UV.col(1) = x.segment(V.rows(),V.rows());
  return UV;
}

Eigen::MatrixXd intrinsic_LSCM(DataGeo &data_mesh, unsigned max_iterations, bool isFreeBoundary){
  Eigen::MatrixXd UV_init = LSCM(data_mesh, isFreeBoundary, false);
  double past_energy = 0;
  double curr_energy = compute_total_energy(data_mesh, UV_init, EnergyType::ASAP, false);
  double tol = 1e-8;
  unsigned itr = 0;
  Eigen::MatrixXd UV = UV_init;
  while(itr<max_iterations && std::abs(past_energy-curr_energy)>tol){
    unsigned flips;
    while(flips = greedy_flip(data_mesh, UV, EnergyType::ASAP));
    UV = LSCM(data_mesh, isFreeBoundary, true);
    past_energy = curr_energy;
    curr_energy = compute_total_energy(data_mesh, UV, EnergyType::ASAP, true);
    ++itr;
  }
  return UV;
}

Eigen::MatrixXd harmonic(DataGeo &data_mesh, bool igrad){
  Eigen::MatrixXd V = data_mesh.V;
  Eigen::MatrixXi F = data_mesh.F;
  
  computeConstraints(V, F, false, '2', C, d);
 
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;

	A.resize(2*V.rows(),2*V.rows());
	b.resize(2*V.rows());
  
  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;

  Eigen::SparseMatrix<double> L;
  if(igrad){
    data_mesh.intTri->requireCotanLaplacian();
    L=data_mesh.intTri->cotanLaplacian;
  }
  else igl::cotmatrix(V,F,L);

  std::vector<Eigen::Triplet<double> > tlist;
  for (int i = 0; i < L.outerSize(); ++i) {
    for (Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(L,i); it; ++it) {
      tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), -it.value())); // - comes from the warning above 
      tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), -it.value()));
    }
  }
  A.setFromTriplets(tlist.begin(), tlist.end());		
  b = Eigen::VectorXd::Zero(2*V.rows());

  Eigen::SparseMatrix<double> Ct, temp1, temp2, res;
  Eigen::SparseMatrix<double> zeros(C.rows(), C.rows()); 
	Ct = C.transpose();
	igl::cat(2, A, Ct, temp1);
	igl::cat(2, C, zeros, temp2);
	igl::cat(1, temp1, temp2, res);

  Eigen::VectorXd rhs(2*V.rows()+C.rows());
	rhs << b,d;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(res);
	solver.factorize(res);
	Eigen::VectorXd x = solver.solve(rhs);
  std::cout << " solver " << solver.info() << std::endl;

  Eigen::MatrixXd UV;
	UV.resize(V.rows(),2);
	UV.col(0) = x.segment(0,V.rows());
 	UV.col(1) = x.segment(V.rows(),V.rows());
  return UV;
}

Eigen::MatrixXd intrinsic_harmonic(DataGeo &data_mesh, unsigned max_iterations){
  Eigen::MatrixXd UV_init = harmonic(data_mesh, false);
  double past_energy = 0;
  double curr_energy = compute_total_energy(data_mesh, UV_init, EnergyType::DIRICHLET, false);
  double tol = 1e-8;
  unsigned itr = 0;
  Eigen::MatrixXd UV = UV_init;
  while(itr<max_iterations && std::abs(past_energy-curr_energy)>tol){
    unsigned flips;
    while(flips = greedy_flip(data_mesh, UV, EnergyType::DIRICHLET));
    UV = harmonic(data_mesh, true); 
    past_energy = curr_energy;
    curr_energy = compute_total_energy(data_mesh, UV, EnergyType::DIRICHLET, true);
    ++itr;
  }

  return UV;
}


Eigen::MatrixXd tutte(DataGeo &data_mesh, bool igrad){
  Eigen::MatrixXd V = data_mesh.V;
  Eigen::MatrixXi F = data_mesh.F;
  
  computeConstraints(V, F, false, '1', C, d);
 
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;

	A.resize(2*V.rows(),2*V.rows());
	b.resize(2*V.rows());

  Eigen::SparseMatrix<double> L;
  
  if(igrad){
    Eigen::MatrixXi F_new = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
    L = compute_L_uniform(V, F_new);
  }
  else L = compute_L_uniform(V, F);

  std::vector<Eigen::Triplet<double> > tlist;
  for (int i = 0; i < L.outerSize(); ++i) {
    for (Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator it(L,i); it; ++it) {
      tlist.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
      tlist.push_back(Eigen::Triplet<double>(it.row()+V.rows(), it.col()+V.rows(), it.value()));
    }
  }
  A.setFromTriplets(tlist.begin(), tlist.end());
  b = Eigen::VectorXd::Zero(2*V.rows());

  Eigen::SparseMatrix<double> Ct, temp1, temp2, res;
  Eigen::SparseMatrix<double> zeros(C.rows(), C.rows()); 
	Ct = C.transpose();
	igl::cat(2, A, Ct, temp1);
	igl::cat(2, C, zeros, temp2);
	igl::cat(1, temp1, temp2, res);

  Eigen::VectorXd rhs(2*V.rows()+C.rows());
	rhs << b,d;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(res);
	solver.factorize(res);
	Eigen::VectorXd x = solver.solve(rhs);
  std::cout << " solver " << solver.info() << std::endl;

  Eigen::MatrixXd UV;
	UV.resize(V.rows(),2);
	UV.col(0) = x.segment(0,V.rows());
 	UV.col(1) = x.segment(V.rows(),V.rows());

  return UV;
}


void faceJacobian(DataGeo &data_mesh, const Eigen::MatrixXd &UV, gcs::Face f, Eigen::Matrix2d &J){
  int i = 0;
  std::array<gcs::Halfedge, 3> halfedges;
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();

  for(gcs::Halfedge he : f.adjacentHalfedges()){
    if(i==3) break; // assumed triangle faces
    halfedges[i]=he;
    ++i;
  } 
 
  // flattening the intrinsic triangle using edge lengths
  double fi_len1 = data_mesh.intTri->edgeLengths[halfedges[0].edge()];
  double fi_len2 = data_mesh.intTri->edgeLengths[halfedges[1].edge()];
  double fi_len3 = data_mesh.intTri->edgeLengths[halfedges[2].edge()];

  Eigen::Matrix2d E, E_tilde;
  double temp = (fi_len2*fi_len2 - fi_len1*fi_len1 - fi_len3*fi_len3)/(-2*fi_len1); 
  E_tilde << fi_len1, temp, 0 , sqrt(fi_len3*fi_len3 - temp*temp); 
 
  // computing the Jacobian
  size_t v1 = data_mesh.intTri->vertexIndices[halfedges[1].tipVertex()];
  size_t v2 = data_mesh.intTri->vertexIndices[halfedges[0].tailVertex()];
  size_t v3 = data_mesh.intTri->vertexIndices[halfedges[0].tipVertex()];
  E << UV(v3,0) - UV(v2,0), UV(v1,0) - UV(v2,0),
        UV(v3,1) - UV(v2,1), UV(v1,1) - UV(v2,1);
  Eigen::Matrix2d E_tilde_inverse;
  E_tilde_inverse << E_tilde(1,1), -E_tilde(0,1), -E_tilde(1,0), E_tilde(0,0);
  E_tilde_inverse = E_tilde_inverse/E_tilde.determinant();
  J = E * E_tilde.inverse();
  
}

double compute_total_energy_fast(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const Eigen::SparseMatrix<double> &Dx,
    const Eigen::SparseMatrix<double> &Dy, const Eigen::VectorXd &areas, const EnergyType &et){
  
  double (*energy)(Eigen::Matrix2d);

  if(et == EnergyType::DIRICHLET) energy = dirichlet;
  if(et == EnergyType::ASAP) energy = asap;
  if(et == EnergyType::ARAP) energy = arap;
  if(et == EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet_alt;
    
  Eigen::VectorXd Dxu = Dx * UV.col(0);		
  Eigen::VectorXd Dxv = Dx * UV.col(1);		
  Eigen::VectorXd Dyu = Dy * UV.col(0);		
  Eigen::VectorXd Dyv = Dy * UV.col(1);
  double total_energy = 0;
  for(int i=0;i<data_mesh.intTri->intrinsicMesh->nFaces();++i){
    Eigen::Matrix2d J;
		J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
    double temp = energy(J)*areas(i);
    if(std::isnan(temp))
      std::cout << " Nan found in Jacobian: " << J << std::endl;
    total_energy += energy(J)*areas(i);
  }
  return total_energy/areas.sum();
}


double compute_total_energy(DataGeo &data_mesh, const Eigen::MatrixXd &UV, const EnergyType &et, bool igrad){
  
  double (*energy)(Eigen::Matrix2d);

  if(et == EnergyType::DIRICHLET) energy = dirichlet;
  if(et == EnergyType::ASAP) energy = asap;
  if(et == EnergyType::ARAP) energy = arap;
  if(et == EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet_alt;
  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;
  if(igrad){
    computeGrad_intrinsic(data_mesh, Dx, Dy, areas);
  }
  else{
    computeSurfaceGradientMatrix(data_mesh.V, data_mesh.F, Dx, Dy);
    igl::doublearea(data_mesh.V, data_mesh.F, areas);
    areas/=2;
  }

  Eigen::VectorXd Dxu = Dx * UV.col(0);		
  Eigen::VectorXd Dxv = Dx * UV.col(1);		
  Eigen::VectorXd Dyu = Dy * UV.col(0);		
  Eigen::VectorXd Dyv = Dy * UV.col(1);
/*
  bool hasNaN = false;
  for (int k = 0; k < Dx.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, k); it; ++it) {
          if (std::isnan(it.value())) {
              hasNaN = true;
              break;
          }
      }
      if (hasNaN) {
          break;
      }
  }
  if(hasNaN) std::cout << " Dx problem " << std::endl;
  */
  double total_energy = 0;
  for(int i=0;i<data_mesh.intTri->intrinsicMesh->nFaces();++i){
    Eigen::Matrix2d J;
		J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
    double temp = energy(J)*areas(i);
    if(std::isnan(temp))
      std::cout << " Nan found in Jacobian: " << J << std::endl;
    total_energy += energy(J)*areas(i);
  }
  return total_energy/areas.sum();
}

double compute_energy_ext(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &UV, const EnergyType &et){

  double (*energy)(Eigen::Matrix2d);

  if(et == EnergyType::DIRICHLET) energy = dirichlet;
  if(et == EnergyType::ASAP) energy = asap;
  if(et == EnergyType::ARAP) energy = arap;
  if(et == EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet_alt;

  Eigen::VectorXd areas;
  Eigen::SparseMatrix<double> Dx, Dy;

  computeSurfaceGradientMatrix(V, F, Dx, Dy);
  igl::doublearea(V, F, areas);
	areas/=2;
  
  Eigen::VectorXd Dxu = Dx * UV.col(0);		
  Eigen::VectorXd Dxv = Dx * UV.col(1);		
  Eigen::VectorXd Dyu = Dy * UV.col(0);		
  Eigen::VectorXd Dyv = Dy * UV.col(1);

  double total_energy = 0;
  for(int i=0; i<F.rows();++i){
    Eigen::Matrix2d J;
		J << Dxu(i), Dyu(i), Dxv(i), Dyv(i);
    double temp = energy(J)*areas(i);
    total_energy += energy(J)*areas(i);
  }

  return total_energy/areas.sum();

}


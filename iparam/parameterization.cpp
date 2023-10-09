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

void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V) {
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
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.
	
	// size of C : #constraints x 2#V
	// C contains one 1 per row and it corresponds to one u or v coordinate
	
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

void computeConstraints(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, bool isFreeBoundary, Eigen::SparseMatrix<double> &C, Eigen::VectorXd &d) {
  Eigen::VectorXi fixed_UV_indices;
  Eigen::MatrixXd fixed_UV_positions;
  if (!isFreeBoundary) {
    // The boundary vertices should be fixed to positions on the unit disc. Find these position and
    // save them in the #V x 2 matrix fixed_UV_position.
    igl::boundary_loop(F, fixed_UV_indices);
    igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
  }
  else {
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
  }
	
	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions,V.rows(), C, d); 

}

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

void computeParameterization(DataGeo &data_mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &UV, Eigen::MatrixXd &new_UV,
    bool isFreeBoundary, bool igrad, const FlipType &ft, const EnergyType &et, int type) {

  auto flip_func = edgeorder_flip;
  double (*energy)(Eigen::Matrix2d);

  if(ft == FlipType::EDGEORDER) flip_func = edgeorder_flip;
  if(ft == FlipType::GREEDY) flip_func = greedy_flip;
  if(ft == FlipType::RANDOM) flip_func = random_flip;
  if(ft == FlipType::HEURISTIC) flip_func = heuristic_flip;

  if(et == EnergyType::DIRICHLET) energy = dirichlet;
  if(et == EnergyType::ASAP) energy = asap;
  if(et == EnergyType::ARAP) energy = arap;
  if(et == EnergyType::SYMMETRIC_DIRICHLET) energy = symmetric_dirichlet;

  if(igrad && UV.size()!=0) while(flip_func(data_mesh, UV, energy));
  
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
	Eigen::SparseMatrix<double> C;
  Eigen::VectorXd d;
  
  computeConstraints(V, F, isFreeBoundary, C, d);

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
    Eigen::SparseMatrix<double> L = compute_L_uniform(V, F);
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

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~

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
	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!
		
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
		
		// A = (DxADx + DyADy    DxADy - DyADx)
		//     (-DxADy + DyADx   DxADx + DyADy)
    Eigen::SparseMatrix<double> B1 = Dx.transpose()*areas.asDiagonal()*Dx + Dy.transpose()*areas.asDiagonal()*Dy; 
    Eigen::SparseMatrix<double> B2 = Dx.transpose()*areas.asDiagonal()*Dy - Dy.transpose()*areas.asDiagonal()*Dx; 
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

	if (type == '4') {
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices
		
		if(UV.size()==0) {
			computeParameterization(data_mesh, V, F, UV, new_UV, isFreeBoundary, igrad, ft, et, '3');
      UV = new_UV;
		}
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

  // Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail

	// build (A C^t; C 0)
  Eigen::SparseMatrix<double> Ct, temp1, temp2, res;
	Ct = C.transpose();
	igl::cat(2, A, Ct, temp1);
	C.conservativeResize(C.rows(), C.rows()+C.cols());
	igl::cat(1, temp1, C, res);

  Eigen::VectorXd rhs(2*V.rows()+C.rows());
	rhs << b,d;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(res);
	solver.factorize(res);
	Eigen::VectorXd x = solver.solve(rhs);

	new_UV.resize(V.rows(),2);
	new_UV.col(0) = x.segment(0,V.rows());
 	new_UV.col(1) = x.segment(V.rows(),V.rows());

}


#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
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

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

Eigen::MatrixXd edge_l;
bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

bool reset=false;
bool prevFreeBoundary=false;
VectorXi fixed_UV_indices;
MatrixXd fixed_UV_positions;

MatrixXd colors;
unsigned iterations=1;
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

static void computeGrad_intrinsic(SparseMatrix<double> & D1, SparseMatrix<double> & D2){
  SparseMatrix<double> G;
  MatrixXd L,L_out;
  MatrixXi F_out;
  igl::edge_lengths(V,F,L);
  igl::intrinsic_delaunay_triangulation(L,F,L_out,F_out);
  igl::grad_intrinsic(L_out, F_out, G);
  D1.resize(F.rows(),V.rows());
  D2.resize(F.rows(),V.rows());
  D1=G.block(0,0,F.rows(), V.rows());
  D2=G.block(F.rows(),0,F.rows(), V.rows());
  
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

void computeParameterization(int type)
{

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
		SparseMatrix<double> L;
		igl::cotmatrix(V,F,L);
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
		computeGrad_intrinsic(Dx,Dy);
		VectorXd areas;
		igl::doublearea(V,F,areas);
		// areas/=2; //since we always use areas and not its square expilicitly, scaling does not change anything
		
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
			computeParameterization('3');
		}
		VectorXd areas;
		igl::doublearea(V,F,areas);
		// areas/=2; //since we always use areas and not its square expilicitly, scaling does not change anything
		SparseMatrix<double> Dx, Dy;
		computeGrad_intrinsic(Dx,Dy);
		// to compute the SVD from the previous iteration
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
	// cout << solver.info() << endl;
	Eigen::VectorXd x = solver.solve(rhs);

	UV.resize(V.rows(),2);
	UV.col(0) = x.segment(0,V.rows());
 	UV.col(1) = x.segment(V.rows(),V.rows());

 	prevFreeBoundary = freeBoundary;
}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
	case '4':
		reset=true;
		if(key=='4'){
			unsigned its = iterations;
			while(its--) computeParameterization(key);
		}
		else computeParameterization(key);
		break;
	case '5':
			// Add your code for detecting and displaying flipped triangles in the
			// UV domain here
		
		break;
	
	case '6': // conformal
		{
		reset=false;
		cout << "Conformal distortion" << endl;
		SparseMatrix<double> Dx, Dy;
		computeGrad_intrinsic(Dx,Dy);
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
		computeGrad_intrinsic(Dx,Dy);
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
		computeGrad_intrinsic(Dx,Dy);
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
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();	
	return true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
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
    cout << "Usage ex4_bin <mesh.off/obj>" << endl;
    load_mesh("../data/cathead.obj");
  }
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }
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

			// TODO: Add more parameters to tweak here...
		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_init = callback_init;

  viewer.launch();
}

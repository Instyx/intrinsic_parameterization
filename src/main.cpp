#include <asm-generic/errno-base.h>
#include <cinttypes>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>


/*** insert any necessary libigl headers here ***/
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/slim.h>

#include <datageo.hpp>
#include "geometrycentral/surface/edge_length_geometry.h"
#include <fstream>
#include <parameterization.hpp>

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
int option_flip;

bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

bool reset=false;
bool prevFreeBoundary=false;
VectorXi fixed_UV_indices;
MatrixXd fixed_UV_positions;


bool igrad;

MatrixXd colors;
MatrixXd C1, P1;
MatrixXd C2, P2;
vector<Eigen::Vector3d> colored_points;
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

/*
void slim_parameterization(){
	if(UV.size()==0) {
      cout <<  " start UV " << endl;
			computeParameterization('2');
      cout << " total energy: " << compute_total_energy() << endl;
		}
  if(igrad && UV.size()!=0) while(flip_func());
	igl::SLIMData slimdata;

  slimdata.V = V;
  if(igrad)
    slimdata.F = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
  else
    slimdata.F = F;
  slimdata.V_o = UV;

  slimdata.v_num = V.rows();
  slimdata.f_num = F.rows(); 
  slimdata.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;

  slimdata.b = fixed_UV_indices;
  slimdata.bc = fixed_UV_positions;
  slimdata.soft_const_p = 0;

  slimdata.proximal_p = 0.0001;
  VectorXd areas;
	SparseMatrix<double> Dx, Dy;
  cout << "befre grad " << endl;
  if(igrad){
    computeGrad_intrinsic(Dx, Dy, areas);
  }
  else{
    computeSurfaceGradientMatrix(Dx,Dy);
    igl::doublearea(V,F,areas);
    areas/=2;
  }
  cout << "after grad" << endl;
  slimdata.Dx = Dx;
  slimdata.Dy = Dy;
  slimdata.M = areas;
  slimdata.mesh_area = slimdata.M.sum();
  slimdata.mesh_improvement_3d = false; // whether to use a jacobian derived from a real mesh or an abstract regular mesh (used for mesh improvement)
  cout << " before energy func " << endl;
  slimdata.energy = igl::slim::compute_energy(slimdata,slimdata.V_o) / slimdata.mesh_area;
  cout << " before solve " << endl; 
  MatrixXd newUV = igl::slim_solve(slimdata,1);
  UV = newUV;
}
*/

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
  EnergyType et;
  FlipType ft;
  if(option_en==0) et = EnergyType::DIRICHLET;
  if(option_en==1) et = EnergyType::SYMMETRIC_DIRICHLET;
  if(option_en==2) et = EnergyType::ASAP;
  if(option_en==3) et = EnergyType::ARAP;

  if(option_flip==0) ft = FlipType::EDGEORDER; 
  if(option_flip==1) ft = FlipType::GREEDY; 
  if(option_flip==2) ft = FlipType::RANDOM; 
  if(option_flip==3) ft = FlipType::HEURISTIC; 

	switch (key) {
  case '1':
	case '2':
	case '3':
	case '4':
    {
		MatrixXd new_UV;
    reset=true;
		if(key=='4'){
			unsigned its = iterations;
			while(its--){ 
        computeParameterization(data_mesh, V, F, UV, new_UV, freeBoundary, igrad, ft, et, key);
        UV = new_UV; 
      }
		}
    else{
      computeParameterization(data_mesh, V, F, UV, new_UV, freeBoundary, igrad, ft, et, key);
      UV = new_UV;
    }
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
  data_mesh.V=V;
  data_mesh.F=F;
  data_mesh.inputMesh.reset(new gcs::ManifoldSurfaceMesh(F));
  data_mesh.inputGeometry.reset(new gcs::VertexPositionGeometry(*data_mesh.inputMesh, V));
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->requireFaceAreas();
  cout << " Edges: " << data_mesh.intTri->intrinsicMesh->nEdges() << endl;
  cout << " Faces: " << data_mesh.intTri->intrinsicMesh->nFaces() << endl;
 
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
     return 0;
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
      ImGui::Checkbox("intrinsic grad", &igrad);
      ImGui::RadioButton("drichlet", &option_en, 0); 
      ImGui::RadioButton("symmetric drichlet", &option_en, 1); 
      ImGui::RadioButton("asap", &option_en, 2); 
      ImGui::RadioButton("arap", &option_en, 3);
      ImGui::RadioButton("edge order", &option_flip, 0); 
      ImGui::RadioButton("greedy", &option_flip, 1); 
      ImGui::RadioButton("random", &option_flip, 2); 
      ImGui::RadioButton("heuristic", &option_flip, 3);

		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_init = callback_init;

  viewer.launch();
}

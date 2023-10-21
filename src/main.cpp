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
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/slim.h>

#include <datageo.hpp>
#include "geometrycentral/surface/edge_length_geometry.h"
#include <fstream>
#include <parameterization.hpp>
#include <iglslim.hpp>
#include <intrinsicflip.hpp>
#include <dirent.h>
#include <igl/lscm.h>

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
Eigen::MatrixXd UV_o;

DataGeo data_mesh;

igl::SLIMData slimdata;

int option_en;
int option_flip;

bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
//igl::opengl::ViewerCore temp3D;
//igl::opengl::ViewerCore temp2D;

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
unsigned flip_granularity=1;

void Redraw()
{
  cout << " in redraw " << endl;
	viewer.data().clear();
	if (!showingUV)
	{
	  viewer.data().set_mesh(V, F);
    //viewer.data().compute_normals();
    viewer.core().align_camera_center(V,F);
		viewer.data().set_face_based(false);
    if(UV.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
    }
	}
	else
	{
	  viewer.data().set_mesh(UV, F);
    viewer.data().compute_normals();
    viewer.core().align_camera_center(UV,F);
		viewer.data().show_texture = false;
	}
}

bool load_mesh(string filename, bool istest)
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
 
  if(!istest){
    Redraw();
    viewer.core().align_camera_center(V);
    showingUV = false;
  }
  return true;
}

void reset_datageo(DataGeo &data_mesh){
  data_mesh.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  data_mesh.intTri->requireEdgeLengths();
  data_mesh.intTri->requireVertexIndices();
  data_mesh.intTri->requireFaceAreas();

}


void intrinsic(fstream &fout, bool free_boundary, unsigned iterations, char key, const EnergyType &et){
  cout << "in intrinsic" << endl;
  MatrixXd new_UV;
  computeParameterization(data_mesh, V, F, UV, new_UV, false, freeBoundary, key);
  UV = new_UV;
  fout << key << "," << free_boundary << "," << compute_total_energy(data_mesh, UV, et) << '\n';

  fout << "edge order,";
  for(unsigned i=0;i<iterations;++i){
    cout << " edge order " << i << endl;
    unsigned flips = 1;
    unsigned flip_iterations = 0;
    unsigned max_flip_iterations = 20;
    while(flip_iterations<max_flip_iterations && (flips=edgeorder_flip(data_mesh, UV, et))){
      cout << flips << endl;
      ++flip_iterations;
      fout << "(" <<flips << ";" << compute_total_energy(data_mesh, UV, et) << ")," ;
    }
    computeParameterization(data_mesh, V, F, UV, new_UV, true, free_boundary, key);
    fout << compute_total_energy(data_mesh, new_UV, et) << '\n';
  }
  reset_datageo(data_mesh);
  fout << "greedy,";
  for(unsigned i=0;i<iterations;++i){
    cout << " greedy " << i << endl;
    unsigned flips = 1;
    unsigned flip_iterations = 0;
    unsigned max_flip_iterations = 20;
    while(flip_iterations<max_flip_iterations && (flips=greedy_flip(data_mesh, UV, et))){
      cout << flips << endl;
      ++flip_iterations;
      fout << "(" <<flips << ";" << compute_total_energy(data_mesh, UV, et) << ")," ;
    }
    computeParameterization(data_mesh, V, F, UV, new_UV, true, free_boundary, key);
    fout << compute_total_energy(data_mesh, new_UV, et) << '\n';
  }
  reset_datageo(data_mesh);
  fout << "random,";
  for(unsigned i=0;i<iterations;++i){
    cout << " random " << i << endl;
    unsigned flips = 1;
    unsigned flip_iterations = 0;
    unsigned max_flip_iterations = 20;
    while(flip_iterations<max_flip_iterations && (flips=random_flip(data_mesh, UV, et))){
      cout << flips << endl;
      ++flip_iterations;
      fout << "(" <<flips << ";" << compute_total_energy(data_mesh, UV, et) << ")," ;
    }
    computeParameterization(data_mesh, V, F, UV, new_UV, true, free_boundary, key);
    fout << compute_total_energy(data_mesh, new_UV, et) << '\n';
  }
  reset_datageo(data_mesh);
  fout << "heuristic,";
  for(unsigned i=0;i<iterations;++i){
    cout << " heuristic " << i << endl;
    unsigned flips = 1;
    unsigned flip_iterations = 0;
    unsigned max_flip_iterations = 20;
    while(flip_iterations<max_flip_iterations && (flips=heuristic_flip(data_mesh, UV, et))){
      cout << flips << endl;
      ++flip_iterations;
      fout << "(" <<flips << ";" << compute_total_energy(data_mesh, UV, et) << ")," ;
    }
    computeParameterization(data_mesh, V, F, UV, new_UV, true, free_boundary, key);
    fout << compute_total_energy(data_mesh, new_UV, et) << '\n';
  }
  reset_datageo(data_mesh);
}

void test(){
  const char* folderPath = "../data";
  DIR* directory = opendir(folderPath);
  if (directory == NULL) {
    cerr << "Failed to open directory." << endl;
    return;
  }
  fstream fout; 
  // opens an existing csv file or creates a new file. 
  fout.open("results.csv", ios::out | ios::app); 

  // Read directory entries
  struct dirent* entry;
  
  while ((entry = readdir(directory)) != NULL) {
    // Skip "." and ".." entries
    if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
        continue;
    }
    string filePath = string(folderPath) + "/" + string(entry->d_name);
    
    fout << string(entry->d_name) << '\n';
    iterations = 2;

    load_mesh(filePath, true);
    intrinsic(fout, false, iterations, '1', EnergyType::DIRICHLET); 
    intrinsic(fout, false, iterations, '2', EnergyType::DIRICHLET); 
    intrinsic(fout, false, iterations, '3', EnergyType::ASAP);
    reset_constraints();
    //intrinsic(fout, true, iterations, '3', EnergyType::ASAP); 
    //reset_constraints();
  }
  
  fout.close();
  // Close the directory
  closedir(directory);
}

void intrinsicUV(const std::unique_ptr<gcs::IntrinsicTriangulation>& intTri, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2){
  vector<Eigen::Vector2d> points;
  vector<array<int,2> > edges;
  for(gcs::Edge e : intTri->intrinsicMesh->edges()){
    points.push_back(UV.row(e.firstVertex().getIndex()).transpose());
    points.push_back(UV.row(e.secondVertex().getIndex()).transpose());
    int t = points.size();
    std::array<int, 2> tmp{t-2, t-1};
    edges.push_back(tmp);

  }
  P1.resize(edges.size(),2);
  P2.resize(edges.size(),2);
  int i = 0;
  for(auto e : edges){
    P1.row(i) = points[e[0]].transpose();
    P2.row(i) = points[e[1]].transpose();
    i++;
  }

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
    cout << "FACEEE in b" << endl;
  }
  return result;
}

void intrinsicEdges(const std::unique_ptr<gcs::IntrinsicTriangulation>& intTri, const Eigen::MatrixXd& V, Eigen::MatrixXd& P1, Eigen::MatrixXd& P2) {
  vector<Eigen::Vector3d> points;
  vector<array<int,2> > edges;
  
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
  
  int i = 0;
  for(auto e : edges){
    P1.row(i) = points[e[0]].transpose();
    P2.row(i) = points[e[1]].transpose();
    i++;
  }
}


bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
  return false;
}

// TODO: datageo da zaten V ve F varmis, ayreten fonksiyonlara parametre olarak koyma iparamda

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
  EnergyType et;
  if(option_en==0) et = EnergyType::DIRICHLET;
  if(option_en==1) et = EnergyType::SYMMETRIC_DIRICHLET;
  if(option_en==2) et = EnergyType::ASAP;
  if(option_en==3) et = EnergyType::ARAP;

  auto flip_func = edgeorder_flip;
  if(option_flip==0) flip_func = edgeorder_flip; 
  if(option_flip==1) flip_func = greedy_flip; 
  if(option_flip==2) flip_func = random_flip;
  if(option_flip==3) flip_func = heuristic_flip; 

	switch (key) {
  case '1':
	case '2':
	case '3':
  {
		MatrixXd new_UV;
    reset=true;
    computeParameterization(data_mesh, V, F, UV, new_UV, freeBoundary, igrad, key);
    UV = new_UV;
    cout << "energy: " << compute_total_energy(data_mesh, UV, et) << endl;
    if(igrad){
      int flips = 1; 
      while(flips){
        flips = flip_func(data_mesh, UV, et);
        cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh, UV, et) <<  endl;
      }
    }
		break;
  }

	case '4':
  {
    MatrixXd new_UV;
    if(UV.size()==0) {
      computeParameterization(data_mesh, V, F, UV, new_UV, freeBoundary, igrad, '3'); 
      UV = new_UV;
      cout << "Initial energy: " << compute_total_energy(data_mesh, UV, et) << endl;
    }
    reset=true;
    unsigned its = iterations;
    while(its--){ 
      computeParameterization(data_mesh, V, F, UV, new_UV, freeBoundary, igrad, key);
      UV = new_UV;
      if(igrad && its%flip_granularity==0){
        int flips = 1; 
        while(flips){
          flips = flip_func(data_mesh, UV, et);
          cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh, UV, et) << endl;
        }
      }
      cout << "energy: " << compute_total_energy(data_mesh, UV, et) << endl;
    }

		break;
  }
  case '5': {
    if(UV.size()==0) {
      MatrixXd new_UV;
      computeParameterization(data_mesh, V, F, UV, new_UV, false, igrad, '2'); 
      UV = new_UV;
      cout << "Initial energy: " << compute_total_energy(data_mesh, UV, et) << endl;
    }
    unsigned its = iterations;
    while(its--){
      slim_parameterization(data_mesh, slimdata, UV, igrad, freeBoundary);
      cout << "energy: " << compute_total_energy(data_mesh, UV, et) << endl;
      cout << " energy in slim: " << slimdata.energy << endl;
      if(igrad && (its+1)%flip_granularity==0){
        int flips = 1; 
        while(flips){
          flips = flip_func(data_mesh, UV, et);
          cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh, UV, et) << endl;
        }
      }
    }
    break;
  }
  case '6': {
    if(UV.size()==0) {
      MatrixXd new_UV;
      computeParameterization(data_mesh, V, F, UV, new_UV, false, igrad, '2'); 
      UV = new_UV;
      cout << "Initial energy: " << compute_total_energy(data_mesh, UV, et) << endl;
    }
    if(!slimdata.has_pre_calc){
      cout << " aaa " << endl;
      Eigen::MatrixXd fixed_UV_positions;
      Eigen::VectorXi fixed_UV_indices;
      boundary(V,F,freeBoundary,fixed_UV_indices, fixed_UV_positions);
      igl::slim_precompute(V,F,UV,slimdata, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, fixed_UV_indices, fixed_UV_positions, 0);
    }
    cout << "it: " << iterations << endl;
    igl::slim_solve(slimdata, iterations);
    UV=slimdata.V_o;
    cout << slimdata.energy << endl;
    break;
  }
  case 'f': {
    if(UV.size()==0) {
      MatrixXd new_UV;
      computeParameterization(data_mesh, V, F, UV, new_UV, false, igrad, '2'); 
      UV = new_UV;
      cout << "Initial energy: " << compute_total_energy(data_mesh, UV, et) << endl;
      
    }
    int flips = 1; 
    while(flips){
      flips = flip_func(data_mesh, UV, et);
      cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh, UV, et) << endl;
    }
    break;
  }
  case 'l': 
  {
    igl::lscm(V,F,UV);
    reset = true;
    break;
  }
  case 's':
  {
    MatrixXd P1, P2;
    intrinsicUV(data_mesh.intTri, P1, P2);
    viewer.data().clear();
    viewer.data().set_mesh(UV, F);
    viewer.core().align_camera_center(UV,F);
    viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1,0,0));
    return true;
  }
  case 'c':
  {
    MatrixXd P1, P2;
    intrinsicEdges(data_mesh.intTri, V, P1, P2);
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.core().align_camera_center(V,F);
    viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1,0,0));
    return true;
  }
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
  case ' ': // space bar -  switches view between mesh and parameterization
    cout << " in space " << endl;
    if(showingUV)
    {
      cout << " sa " << endl;
      //temp2D = viewer.core();
      //viewer.core() = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        cout << " as " << endl;
        //temp3D = viewer.core();
        //viewer.core() = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
  cout << " before redraw " << endl;
	Redraw();	
  cout << " after redraw " << endl;
	return true;
}
/*
bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core();
	temp2D = viewer.core();
	temp2D.orthographic = true;

	return false;
}
*/
int main(int argc,char *argv[]) {
  if(argc != 2) {
    test();
    return 0;
  }
 
  else
  {
    // Read points and normals
    load_mesh(argv[1], false);

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
      ImGui::InputScalar("ARAP/SLIM iterations", ImGuiDataType_U32, &iterations, 0, 0);
      ImGui::InputScalar("Flip Remesh Granularity", ImGuiDataType_U32, &flip_granularity, 0, 0);
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
  //viewer.callback_init = callback_init;

  viewer.launch();
}

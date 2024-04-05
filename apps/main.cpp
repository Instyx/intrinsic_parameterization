#include <asm-generic/errno-base.h>
#include <cinttypes>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>


#include <igl/slim.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>

#include <datageo.hpp>
#include <fstream>
#include <parameterization.hpp>
#include <iglslim.hpp>
#include <intrinsicflip.hpp>
#include <intrinsicslim.hpp>
#include <test.hpp>
#include <dirent.h>
#include <igl/lscm.h>
#include <igl/arap.h>
#include <metrics.hpp>
#include <string>
#include <chrono>

using namespace std;
using namespace Eigen;

ofstream outfile("output.txt");

using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;
Eigen::MatrixXi F_o;
// UV coordinates, #V x2
Eigen::MatrixXd UV;
Eigen::MatrixXd UV_o;

DataGeo data_mesh;
DataGeo data_mesh_o;

igl::SLIMData slimdata;

int option_en;
int option_flip;

bool showingUV = false;
bool freeBoundary = false;
bool onlyDelaunay= false;
double TextureResolution = 10;
//igl::opengl::ViewerCore temp3D;
//igl::opengl::ViewerCore temp2D;

bool reset=false;
bool prevFreeBoundary=false;
VectorXi fixed_UV_indices;
MatrixXd fixed_UV_positions;


bool igrad;
bool intrinsic_edges;

MatrixXd colors;
MatrixXd C1, P1;
MatrixXd C2, P2;
vector<Eigen::Vector3d> colored_points;
unsigned iterations=1;
unsigned flip_granularity=1;

void Redraw()
{
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
  if(colors.size()!=0) viewer.data().set_colors(colors);
}

bool load_mesh(DataGeo &datageo, string filename, bool istest)
{
  igl::read_triangle_mesh(filename,V,F);
  datageo.V=V;
  datageo.F=F;
  datageo.inputMesh.reset(new gcs::ManifoldSurfaceMesh(F));
  datageo.inputGeometry.reset(new gcs::VertexPositionGeometry(*datageo.inputMesh, V));
  datageo.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*datageo.inputMesh, *datageo.inputGeometry));
  datageo.intTri->requireEdgeLengths();
  datageo.intTri->requireVertexIndices();
  datageo.intTri->requireFaceAreas();

  if(!istest){
    Redraw();
    viewer.core().align_camera_center(V);
    showingUV = false;
  }
  return true;
}

void reset_datageo(DataGeo &datageo){
  datageo.intTri.reset(new gcs::SignpostIntrinsicTriangulation(*data_mesh.inputMesh, *data_mesh.inputGeometry));
  datageo.intTri->requireEdgeLengths();
  datageo.intTri->requireVertexIndices();
  datageo.intTri->requireFaceAreas();

}

// to visualize intrinsic edges on the UV domain
void intrinsicUV(const std::unique_ptr<gcs::IntrinsicTriangulation>& intTri,Eigen::MatrixXd &UV, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2){
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
    cout << "FACE in b" << endl;
  }
  return result;
}

// to visualize intrinsic edges on the mesh
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
  // uniform laplacian 
  case '1':
  // harmonic - cotangent laplacian
	case '2':
  // LSCM
	case '3':
  {
		MatrixXd new_UV;
    reset=true;
    computeParameterization(data_mesh, V, F, UV, new_UV, freeBoundary, igrad, key);
    UV = new_UV;
    cout << "energy sadf: " << compute_total_energy(data_mesh, UV, et, true) << endl;
    // cout << "energy: " << compute_energy_ext(V, F, UV, et) << endl;
    if(F_o.size()!=0){
      computeParameterization(data_mesh, V, F_o, UV_o, new_UV, freeBoundary, igrad, key);
      UV_o = new_UV;
      cout << "energy F_new: " << compute_total_energy(data_mesh, UV_o, et, true) << endl;
      //cout << "energy F_new: " << compute_energy_ext(V, F_o, UV_o, et) << endl;
    }
		break;
  }
  // ARAP
	case '4':
  {
    MatrixXd new_UV;
    if(UV.size()==0) {
      computeParameterization(data_mesh, V, F, UV, new_UV, false , igrad, '2');
      UV = new_UV;
      //cout << "Initial energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
      cout << "Initial energy: " << compute_energy_ext(V, F, UV, et) << endl;
      reset_constraints();
      computeParameterization(data_mesh, V, F, UV, new_UV, freeBoundary, igrad, key);
      UV = new_UV;
      //cout << "energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
      cout << "energy: " << compute_energy_ext(V, F, UV, et) << endl;
    }
    //reset=true;
    unsigned its = iterations;
    for(unsigned i=0;i<its-1;++i){
      if(igrad && (i+1)%flip_granularity==0){
        int flips = 1;
        unsigned del = 0;
        while(flips){
          flips = flip_func(data_mesh, UV, del, et);
          cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
        }
      }
      computeParameterization(data_mesh, V, F, UV, new_UV, freeBoundary, igrad, key);
      UV = new_UV;
      //cout << "energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
      cout << "energy: " << compute_energy_ext(V, F, UV, et) << endl;
    }

		break;
  }
  case 'e':
  {
    MatrixXd new_UV;
    if(UV_o.size()==0) {
      reset_constraints();
      computeParameterization(data_mesh_o, V, F_o, UV_o, new_UV, false , igrad, '2');
      UV_o = new_UV;
      //cout << "Initial energy: " << compute_total_energy(data_mesh_o, UV_o, et, true) << endl;
      cout << "Initial energy: " << compute_energy_ext(V, F_o, UV_o, et) << endl;
      reset_constraints();
      computeParameterization(data_mesh_o, V, F_o, UV_o, new_UV, freeBoundary, igrad, '4');
      UV_o = new_UV;
      //cout << "energy: " << compute_total_energy(data_mesh_o, UV_o, et, true) << endl;
      cout << "energy: " << compute_energy_ext(V, F_o, UV_o, et) << endl;
    }
    //reset=true;
    unsigned its = iterations;
    for(unsigned i=0;i<its-1;++i){
      if(igrad && (i+1)%flip_granularity==0){
        int flips = 1;
        unsigned del;
        while(flips){
          flips = flip_func(data_mesh_o, UV_o, del, et);
          cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh_o, UV_o, et, true) << endl;
        }
      }
      computeParameterization(data_mesh_o, V, F_o, UV_o, new_UV, freeBoundary, igrad, '4');
      UV_o = new_UV;
      //cout << "energy: " << compute_total_energy(data_mesh_o, UV_o, et, true) << endl;
      cout << "energy: " << compute_energy_ext(V, F_o, UV_o, et) << endl;
    }

		break;
  }
  // SLIM
  case '5': 
  {
    if(UV.size()==0) {
      MatrixXd new_UV;
      computeParameterization(data_mesh, V, F, UV, new_UV, false, igrad, '2');
      UV = new_UV;
      //cout << "Initial energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
      cout << "Initial energy: " << compute_energy_ext(V, F, UV, et) << endl;
      slim_parameterization(data_mesh, slimdata, V, F, UV, igrad, freeBoundary);
      //cout << "energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
      cout << "energy: " << compute_energy_ext(V, F, UV, et) << endl;
    }
    unsigned its = iterations;
    for(unsigned i=0;i<its-1;++i){
      if(igrad && (i+1)%flip_granularity==0){
        int flips = 1;
        unsigned del;
        while(flips){
          flips = flip_func(data_mesh, UV, del, et);
          cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
        }
      }
      slim_parameterization(data_mesh, slimdata, V, F, UV, igrad, freeBoundary);
      //cout << "energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
      cout << "energy: " << compute_energy_ext(V, F, UV, et) << endl;
      cout << " energy in slim: " << slimdata.energy << endl;
    }
    break;
  }
  // for comparison
  case 'r': 
  {
    reset_datageo(data_mesh_o);
    slimdata.has_pre_calc = false;
    if(UV_o.size()==0) {
      reset_constraints();
      MatrixXd new_UV;
      computeParameterization(data_mesh_o, V, F_o, UV_o, new_UV, false, igrad, '2');
      UV_o = new_UV;
      //cout << "Initial energy: " << compute_total_energy(data_mesh_o, UV_o, et, true) << endl;
      cout << "Initial energy: " << compute_energy_ext(V, F_o, UV_o, et) << endl;
      slim_parameterization(data_mesh_o, slimdata, V, F_o, UV_o, igrad, freeBoundary);
      //cout << "energy: " << compute_total_energy(data_mesh_o, UV_o, et, true) << endl;
      cout << "energy: " << compute_energy_ext(V, F_o, UV_o, et) << endl;
    }
    unsigned its = iterations;
    for(unsigned i=0;i<its-1;++i){
      if(igrad && (i+1)%flip_granularity==0){
        int flips = 1;
        unsigned del;
        while(flips){
          flips = flip_func(data_mesh_o, UV_o, del, et);
          cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh_o, UV_o, et, true) << endl;
        }
      }
      slim_parameterization(data_mesh_o, slimdata, V, F_o, UV_o, igrad, freeBoundary);
      //cout << "energy: " << compute_total_energy(data_mesh_o, UV_o, et, true) << endl;
      cout << "energy: " << compute_energy_ext(V, F_o, UV_o, et) << endl;
      cout << " energy in slim: " << slimdata.energy << endl;
    }
    break;
  }
  case '6': {
    //test_ARAP();
    Eigen::MatrixXd UV_init = tutte(data_mesh, true);
    std::fstream foutt;
    intrinsicslim(data_mesh, UV_init, UV, 1000, 50, foutt);
    break;
  }

  case '7': {
    break;
  }
  case '8': {
    //compareIDTvsIPARAM(data_mesh, freeBoundary, EnergyType::ARAP, UV_o, UV);
    break;
  }
  case 'l':
  {
    //igl::lscm(V,F,UV);
    //reset = true;
   
    data_mesh.intTri->flipToDelaunay();
    MatrixXd new_UV;
    computeParameterization(data_mesh, V, F, UV_o, new_UV, false, true, '2');
    UV_o = new_UV;
    cout << "IDT: " << compute_total_energy(data_mesh, UV_o, EnergyType::DIRICHLET, true) << endl;
    reset_datageo(data_mesh);
    computeParameterization(data_mesh, V, F, UV, new_UV, false, false, '2');
    UV = new_UV;
    unsigned del;
    while(greedy_flip(data_mesh, UV, del, EnergyType::DIRICHLET));
    computeParameterization(data_mesh, V, F, UV, new_UV, false, true, '2');
    UV = new_UV;
    cout << "greedy: " << compute_total_energy(data_mesh, UV, EnergyType::DIRICHLET, true) << endl;
    
   // reset_datageo(data_mesh);
    return true;
    
  }
  // visualize intrinsic edges on the UV domain
  case 's':
  {
    viewer.data().clear();
    viewer.data().set_mesh(UV, F);
    viewer.core().align_camera_center(UV,F);
    if(intrinsic_edges) {
      intrinsicUV(data_mesh.intTri, UV, C1, C2);
      viewer.data().add_edges(C1, C2, Eigen::RowVector3d(1,0,0));
    }
    //reset_datageo(data_mesh);
    return true;
  }
  // visualize edges on the UV domain with UV_o for comparison
  case 'b':
  {
    viewer.data().clear();
    viewer.data().set_mesh(UV_o, F);
    viewer.core().align_camera_center(UV_o,F);
    if(intrinsic_edges) {
      intrinsicUV(data_mesh_o.intTri, UV_o, P1, P2);
      viewer.data().add_edges(P1, P2, Eigen::RowVector3d(0,0,1));
    }
    //reset_datageo(data_mesh);
    return true;

  }

  // visualize edges on the mesh and print metrics
  case 'c':
  {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(true);
    if(UV.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
      /*
      vector<double> res, res1, res2;
      compute_metrics(data_mesh, UV, res);
      cout << "extrinsic_flipped: " << res[0] << endl
          << "extrinsic_max_area_dist: " << res[1] << endl
          << "extrinsic_average_area_error: " << res[2] << endl
          << "extrinsic_max_angle_dist: " << res[3] << endl
          << "extrinsic_average_angle_error: " << res[4] << endl
          << "intrinsic_flipped: " << res[5] << endl
          << "intrinsic_max_area_dist: " << res[6] << endl
          << "intrinsic_average_area_error: " << res[7] << endl
          << "intrinsic_max_angle_dist: " << res[8] << endl
          << "intrinsic_average_angle_error: " << res[9] << endl;
      minmax_distortions(V,F, UV, res1);
      minmax_distortions_intri(data_mesh, UV, res2);
      cout << " conformal min/max: " << res1[0] << "   " << res1[1] << endl;
      cout << " intri conformal min/max: " << res2[0] << "   " << res2[1] << endl;
      cout << " isometric min/max: " << res1[2] << "   " << res1[3] << endl;
      cout << " intri isometric min/max: " << res2[2] << "   " << res2[3] << endl;
      cout << " authalic min/max: " << res1[4] << "   " << res1[5] << endl;
      cout << " intri authalic min/max: " << res2[4] << "   " << res2[5] << endl;
      cout << endl;
      */

    }
    viewer.core().align_camera_center(V,F);
    if(intrinsic_edges) {
      intrinsicEdges(data_mesh.intTri, V, C1, C2);
      viewer.data().add_edges(C1, C2, Eigen::RowVector3d(1,0,0));
    }
   
    //reset_datageo(data_mesh);
    return true;
  }
  // visualize intrinsic edges on the mesh with UV_o for comparison and print metrics
  case 'v':
  {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(true);
    if(UV_o.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV_o);
      viewer.data().show_texture = true;
      /*
      vector<double> res, res1, res2;
      compute_metrics(data_mesh, UV_o, res);
      cout << "extrinsic_flipped: " << res[0] << endl
          << "extrinsic_max_area_dist: " << res[1] << endl
          << "extrinsic_average_area_error: " << res[2] << endl
          << "extrinsic_max_angle_dist: " << res[3] << endl
          << "extrinsic_average_angle_error: " << res[4] << endl
          << "intrinsic_flipped: " << res[5] << endl
          << "intrinsic_max_area_dist: " << res[6] << endl
          << "intrinsic_average_area_error: " << res[7] << endl
          << "intrinsic_max_angle_dist: " << res[8] << endl
          << "intrinsic_average_angle_error: "  << res[9] << endl;
      minmax_distortions(V,F, UV_o, res1);
      minmax_distortions_intri(data_mesh, UV_o, res2);
      cout << " conformal min/max: " << res1[0] << "   " << res1[1] << endl;
      cout << " intri conformal min/max: " << res2[0] << "   " << res2[1] << endl;
      cout << " isometric min/max: " << res1[2] << "   " << res1[3] << endl;
      cout << " intri isometric min/max: " << res2[2] << "   " << res2[3] << endl;
      cout << " authalic min/max: " << res1[4] << "   " << res1[5] << endl;
      cout << " intri authalic min/max: " << res2[4] << "   " << res2[5] << endl;
      cout << endl;
    */
    }
    viewer.core().align_camera_center(V,F);
    if(intrinsic_edges) {
      intrinsicEdges(data_mesh_o.intTri, V, P1, P2);
      viewer.data().add_edges(P1, P2, Eigen::RowVector3d(0,1,0));
    }

    //reset_datageo(data_mesh);
    return true;
  }
  case 'h':
  {
    if(isDelaunayFlipBad(data_mesh, UV)) cout << "Delaunay increased energy, Baddo" << endl;
    else cout << "Gooddo" << endl;
    /*
    cout << " flipso: " << flipEdgesifCoplanar(data_mesh, onlyDelaunay) << endl;
    F_o = data_mesh.intTri->intrinsicMesh->getFaceVertexMatrix<int>();
    viewer.data().clear();
    viewer.data().set_mesh(V, F_o);
		viewer.data().set_face_based(true);
    //F=F_o;
    //
    */
    return true;
  }
  case 'j':
  {
    data_mesh_o = compareIDTvsGreedy(data_mesh);
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
      //temp2D = viewer.core();
      //viewer.core() = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        //temp3D = viewer.core();
        //viewer.core() = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();
	return true;
}


void print_usage(){
  cout << "  - tick 'free boundary' to compute parameterization boundary-free" << endl;
  cout << "  - tick 'igrad' to compute parameterization intrinsically" << endl;
  cout << "  - tick 'intrinsic edges' to visualize edges when using keys 's,b,c,v'. " << endl;
  cout << "  - energy type and flipping order can be selected from the menu on left " << endl;
  cout << "  - adjust 'ARAP/SLIM iterations' for the change the total iteration of the global/local approach when minimizing ARAP, Symmetric Dirichlet" << endl;
  cout << "  - adjust 'Flip Granularity' to change the flip granularity for intrinsic parameterization for ARAP, Symmetric Dirichlet" << endl; 
  cout << "  - to compute parameterization (result is saved in UV) :" << endl;
  cout << "    '1' for uniform Laplacian" << endl;
  cout << "    '2' for cotangent Laplacian (Dirichlet Energy)" << endl;
  cout << "    '3' for LSCM (ASAP)" << endl;
  cout << "    '4' for ARAP" << endl;
  cout << "    '5' for SLIM (Symmetric Dirichlet)" << endl;
  cout << "  - to compare different UV mappings when minimizing ARAP and Symetric Dirichlet : " << endl;
  cout << "    'e' for ARAP and result is saved in UV_o: " << endl;
  cout << "    'r' for SLIM and result is saved in UV_o: " << endl;
  cout << "  - 'f' for flipping algortihm with the selected order and energy" << endl;
  cout << "    when 'f' is used, the current UV mapping that is used for the flipping is saved as UV_o" << endl;
  cout << "  - 'space' for switching between UV domain and 3D mesh surface" << endl;
  cout << "  - 'l' for resetting intrinsic triangulation to input extrinsic triangulation" << endl;
  cout << "  - 's', 'b' for showing respectively UV, UV_o" << endl;
  cout << "  - 'c', 'v' for showing the textured 3D mesh and printing the metrics respectively for UV, UV_o  " << endl;
}

int main(int argc,char *argv[]) {
  // evaluation
  if(argc != 2) {
    test_all_withtextures();
    return 0;
  }

  // manual mesh selection with GUI
  else
  {
    // Read points and normals
    load_mesh(data_mesh, argv[1], false);
    load_mesh(data_mesh_o, argv[1], false);
    print_usage();
    cout << endl;
    cout << " Vertices: " << data_mesh.intTri->intrinsicMesh->nVertices() << endl;
    cout << " Edges: " << data_mesh.intTri->intrinsicMesh->nEdges() << endl;
    cout << " Faces: " << data_mesh.intTri->intrinsicMesh->nFaces() << endl;

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
		if (ImGui::CollapsingHeader("Parameterization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::Checkbox("Free boundary", &freeBoundary);
      ImGui::InputScalar("ARAP/SLIM iterations", ImGuiDataType_U32, &iterations, 0, 0);
      ImGui::InputScalar("Flip Remesh Granularity", ImGuiDataType_U32, &flip_granularity, 0, 0);
      ImGui::Checkbox("intrinsic grad", &igrad);
      ImGui::Checkbox("intrinsic edges", &intrinsic_edges);
      ImGui::Checkbox("only delaunay", &onlyDelaunay);
      ImGui::RadioButton("Dirichlet", &option_en, 0);
      ImGui::RadioButton("symmetric Dirichlet", &option_en, 1);
      ImGui::RadioButton("ASAP", &option_en, 2);
      ImGui::RadioButton("ARAP", &option_en, 3);
      ImGui::RadioButton("edge order", &option_flip, 0);
      ImGui::RadioButton("greedy", &option_flip, 1);
      ImGui::RadioButton("random", &option_flip, 2);
      ImGui::RadioButton("heuristic", &option_flip, 3);
      ImGui::RadioButton("delaunay", &option_flip, 4);

		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;

  viewer.launch();
}

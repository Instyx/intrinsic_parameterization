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
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>

#include <datageo.hpp>
#include "geometrycentral/surface/edge_length_geometry.h"
#include <fstream>
#include <parameterization.hpp>
#include <iglslim.hpp>
#include <intrinsicflip.hpp>
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
bool intrinsic_edges;

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
  if(colors.size()!=0) viewer.data().set_colors(colors);
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
  cout << " Vertices: " << data_mesh.intTri->intrinsicMesh->nVertices() << endl;
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


void flip_res(fstream &fout, bool free_boundary, char key, int start_iterations, string mesh_name, double start_energy, const EnergyType &et){

  auto start = std::chrono::high_resolution_clock::now();

  auto end = std::chrono::high_resolution_clock::now();

  unsigned flips = 1;
  unsigned i = 0;
  fout << mesh_name << "," << key << "," << free_boundary << "," << "start," << i << "," << start_iterations << "," << start_energy << "," << -1 << '\n';
  while(flips){
    start = chrono::high_resolution_clock::now();
    flips = edgeorder_flip(data_mesh, UV, et);
    end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    ++i;
    fout << mesh_name << "," << key << "," << free_boundary << "," << "edge order," << i << "," << flips << ","
      << compute_total_energy(data_mesh, UV, et, true) << "," << duration << '\n';
  }
  reset_datageo(data_mesh);
  i = 0;
  flips=1;
  while(flips){
    start = chrono::high_resolution_clock::now();
    flips = greedy_flip(data_mesh, UV, et);
    end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    ++i;
    fout << mesh_name << "," << key << "," << free_boundary << "," << "greedy," << i << "," << flips << "," 
      << compute_total_energy(data_mesh, UV, et, true) << "," << duration << '\n';
  }
  reset_datageo(data_mesh);
  i = 0;
  flips=1;
  while(flips){
    start = chrono::high_resolution_clock::now();
    flips = random_flip(data_mesh, UV, et);
    end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    ++i;
    fout << mesh_name << "," << key << "," << free_boundary << "," << "random," << i << "," << flips << "," 
      << compute_total_energy(data_mesh, UV, et, true) << "," << duration << '\n';
  }
  reset_datageo(data_mesh);
  i = 0;
  flips=1;
  while(flips){
    start = chrono::high_resolution_clock::now();
    flips = heuristic_flip(data_mesh, UV, et);
    end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    ++i;
    fout << mesh_name << "," << key << "," << free_boundary << "," << "heuristic," << i << "," << flips << ","
      << compute_total_energy(data_mesh, UV, et, true) << "," << duration << '\n';
  }
  reset_datageo(data_mesh);
}


void intrinsic(fstream &fout, bool free_boundary, unsigned iterations, char key, const EnergyType &et){
  cout << "in intrinsic" << endl;
  MatrixXd new_UV;
  computeParameterization(data_mesh, V, F, UV, new_UV, false, freeBoundary, key);
  UV = new_UV;
  fout << key << "," << free_boundary << "," << compute_total_energy(data_mesh, UV, et, true) << '\n';

  fout << "edge order,";
  for(unsigned i=0;i<iterations;++i){
    cout << " edge order " << i << endl;
    unsigned flips = 1;
    unsigned flip_iterations = 0;
    unsigned max_flip_iterations = 20;
    while(flip_iterations<max_flip_iterations && (flips=edgeorder_flip(data_mesh, UV, et))){
      cout << flips << endl;
      ++flip_iterations;
      fout << "(" <<flips << ";" << compute_total_energy(data_mesh, UV, et, true) << ")," ;
    }
    computeParameterization(data_mesh, V, F, UV, new_UV, true, free_boundary, key);
    fout << compute_total_energy(data_mesh, new_UV, et, true) << '\n';
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
      fout << "(" <<flips << ";" << compute_total_energy(data_mesh, UV, et, true) << ")," ;
    }
    computeParameterization(data_mesh, V, F, UV, new_UV, true, free_boundary, key);
    fout << compute_total_energy(data_mesh, new_UV, et, true) << '\n';
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
      fout << "(" <<flips << ";" << compute_total_energy(data_mesh, UV, et, true) << ")," ;
    }
    computeParameterization(data_mesh, V, F, UV, new_UV, true, free_boundary, key);
    fout << compute_total_energy(data_mesh, new_UV, et, true) << '\n';
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
      fout << "(" <<flips << ";" << compute_total_energy(data_mesh, UV, et, true) << ")," ;
    }
    computeParameterization(data_mesh, V, F, UV, new_UV, true, free_boundary, key);
    fout << compute_total_energy(data_mesh, new_UV, et, true) << '\n';
  }
  reset_datageo(data_mesh);
}


// test for intrinsic flipping
void test_intflip(){
  const char* folderPath = "../res_data";
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
  fout << "mesh name," << "type," << "free boundary," << "flipping order," << "flip iteration," << "flips / start iterations," << "energy," << "time"<< '\n';
  while ((entry = readdir(directory)) != NULL) {
    // Skip "." and ".." entries
    if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
        continue;
    }
    string filePath = string(folderPath) + "/" + string(entry->d_name);
    

    load_mesh(filePath, true);
    MatrixXd new_UV;
    // dirichlet 
    computeParameterization(data_mesh, V, F, UV, new_UV, false, false, '2');
    UV = new_UV;
    double start_energy = compute_total_energy(data_mesh, UV, EnergyType::DIRICHLET, true);
    flip_res(fout, false, '2', -1, string(entry->d_name), start_energy, EnergyType::DIRICHLET);

    // LSCM fixed
    computeParameterization(data_mesh, V, F, UV, new_UV, false, false, '3');
    UV = new_UV;
    start_energy = compute_total_energy(data_mesh, UV, EnergyType::ASAP, true);
    flip_res(fout, false, '3', -1, string(entry->d_name), start_energy, EnergyType::ASAP);

    // ARAP fixed
    MatrixXd start_UV = UV;
    computeParameterization(data_mesh, V, F, start_UV, new_UV, false, false, '4');
    UV = new_UV;
    start_energy = compute_total_energy(data_mesh, UV, EnergyType::ARAP, true);
    flip_res(fout, false, '4', 1, string(entry->d_name), start_energy, EnergyType::ARAP);

    int iterations = 2;
    while(iterations--){
      computeParameterization(data_mesh, V, F, UV, new_UV, false, false, '4');
      UV = new_UV;
    }
    start_energy = compute_total_energy(data_mesh, UV, EnergyType::ARAP, true);
    flip_res(fout, false, '4', 3, string(entry->d_name), start_energy, EnergyType::ARAP);
    
    iterations = 7;
    while(iterations--){
      computeParameterization(data_mesh, V, F, UV, new_UV, false, false, '4');
      UV = new_UV;
    }
    start_energy = compute_total_energy(data_mesh, UV, EnergyType::ARAP, true);
    flip_res(fout, false, '4', 10, string(entry->d_name), start_energy, EnergyType::ARAP);

    // ARAP free
    reset_constraints();
    computeParameterization(data_mesh, V, F, start_UV, new_UV, true, false, '4');
    UV = new_UV;
    start_energy = compute_total_energy(data_mesh, UV, EnergyType::ARAP, true);
    flip_res(fout, true, '4', 1, string(entry->d_name), start_energy, EnergyType::ARAP);

    iterations = 2;
    while(iterations--){
      computeParameterization(data_mesh, V, F, UV, new_UV, true, false, '4');
      UV = new_UV;
    }
    start_energy = compute_total_energy(data_mesh, UV, EnergyType::ARAP, true);
    flip_res(fout, true, '4', 3, string(entry->d_name), start_energy, EnergyType::ARAP);
    
    iterations = 7;
    while(iterations--){
      computeParameterization(data_mesh, V, F, UV, new_UV, true, false, '4');
      UV = new_UV;
    }
    start_energy = compute_total_energy(data_mesh, UV, EnergyType::ARAP, true);
    flip_res(fout, true, '4', 10, string(entry->d_name), start_energy, EnergyType::ARAP);


    // LSCM free
    reset_constraints();
    computeParameterization(data_mesh, V, F, UV, new_UV, true, false, '3');
    UV = new_UV;
    start_energy = compute_total_energy(data_mesh, UV, EnergyType::ASAP, true);
    flip_res(fout, true, '3', -1, string(entry->d_name), start_energy, EnergyType::ASAP);

    reset_constraints();
  }
  
  fout.close();
  // Close the directory
  closedir(directory);
}

void single_intri(fstream &fout, string mesh_name, bool free_boundary, int type, char key, const EnergyType &et){
    Eigen::MatrixXd new_UV;
    
    // greedy
    auto start = chrono::high_resolution_clock::now();
    computeParameterization(data_mesh, V, F, UV, new_UV, free_boundary, false, key);
    auto end = chrono::high_resolution_clock::now();
    UV = new_UV;
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::vector<double> res_ext, res_int, res;
    minmax_distortions(V, F, UV, res_ext);
    minmax_distortions_intri(data_mesh, UV, res_int);
    compute_metrics(data_mesh, UV, res);

    fout << mesh_name << "," << type << "," << free_boundary << ",greedy," << 0 << "," << false << "," << 
      compute_total_energy(data_mesh, UV, et, true) << "," <<  compute_total_energy(data_mesh, UV, et, false) << ","
      << res_ext[0] << "," << res_ext[1] << "," << res_ext[2] << "," << res_int[0] << "," <<  res_int[1] << "," << res_int[2] << ","
      << res_ext[3] << "," << res_ext[4] << "," << res_ext[5] << "," << res_int[3] << "," <<  res_int[4] << "," << res_int[5] << ","
      << res[0] << "," << res[1] << "," << res[2] << "," << res[3] << "," << res[4] << "," 
      << res[5] << "," << res[6] << "," << res[7] << "," << res[8] << "," << res[9] << ","
      << duration << '\n';
    start = chrono::high_resolution_clock::now();
    while(greedy_flip(data_mesh, UV, et));
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    minmax_distortions_intri(data_mesh, UV, res_int);
    compute_metrics(data_mesh, UV, res);
    fout << mesh_name << "," << type << "," << free_boundary << ",greedy," << 1 << "," << true << "," <<
      compute_total_energy(data_mesh, UV, et, true) << "," <<  compute_total_energy(data_mesh, UV, et, false) << ","
      << res_ext[0] << "," << res_ext[1] << "," << res_ext[2] << "," << res_int[0] << "," <<  res_int[1] << "," << res_int[2] << ","
      << res_ext[3] << "," << res_ext[4] << "," << res_ext[5] << "," << res_int[3] << "," <<  res_int[4] << "," << res_int[5] << ","
      << res[0] << "," << res[1] << "," << res[2] << "," << res[3] << "," << res[4] << "," 
      << res[5] << "," << res[6] << "," << res[7] << "," << res[8] << "," << res[9] << ","
      << duration << '\n';
    start = chrono::high_resolution_clock::now();
    computeParameterization(data_mesh, V, F, UV, new_UV, free_boundary, true, key);
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double after_igrad_energy = compute_total_energy(data_mesh, new_UV, et, true);
    minmax_distortions(V, F, new_UV, res_ext);
    minmax_distortions_intri(data_mesh, new_UV, res_int);
    compute_metrics(data_mesh, new_UV, res);
    fout << mesh_name << "," << type << "," << free_boundary << ",greedy," << 2 << "," << false << "," << 
      compute_total_energy(data_mesh, new_UV, et, true) << "," <<  compute_total_energy(data_mesh, new_UV, et, false) << ","
    << res_ext[0] << "," << res_ext[1] << "," << res_ext[2] << "," << res_int[0] << "," <<  res_int[1] << "," << res_int[2] << ","
      << res_ext[3] << "," << res_ext[4] << "," << res_ext[5] << "," << res_int[3] << "," <<  res_int[4] << "," << res_int[5] << ","
      << res[0] << "," << res[1] << "," << res[2] << "," << res[3] << "," << res[4] << "," 
      << res[5] << "," << res[6] << "," << res[7] << "," << res[8] << "," << res[9] << ","
      << duration << '\n';


    reset_datageo(data_mesh);
}

void single_intri_arap(fstream &fout, string mesh_name, bool free_boundary, int type, unsigned iterations, unsigned granularity, const EnergyType &et){
    Eigen::MatrixXd new_UV, start_UV;
    
    // greedy
    auto start = chrono::high_resolution_clock::now();
    if(free_boundary) reset_constraints();
    computeParameterization(data_mesh, V, F, UV, start_UV, false, false, '2');
    if(free_boundary) reset_constraints();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::vector<double> res_ext, res_int, res;
    minmax_distortions(V,F, start_UV, res_ext);
    minmax_distortions_intri(data_mesh, start_UV, res_int);
    compute_metrics(data_mesh, start_UV, res);
    fout << mesh_name << "," << type << "," << free_boundary << ",greedy," << 0 << "," << false << "," << 
      compute_total_energy(data_mesh, start_UV, et, true) << "," <<  compute_total_energy(data_mesh, start_UV, et, false) << ","
      << res_ext[0] << "," << res_ext[1] << "," << res_ext[2] << "," << res_int[0] << "," <<  res_int[1] << "," << res_int[2] << ","
      << res_ext[3] << "," << res_ext[4] << "," << res_ext[5] << "," << res_int[3] << "," <<  res_int[4] << "," << res_int[5] << ","
      << res[0] << "," << res[1] << "," << res[2] << "," << res[3] << "," << res[4] << "," 
      << res[5] << "," << res[6] << "," << res[7] << "," << res[8] << "," << res[9] << ","
      << duration << '\n';
    UV = start_UV;
    unsigned to_add=0;
  for(unsigned i=0; i<iterations; ++i){
      if(i!=0 && (i-1)%granularity==0){
        start = chrono::high_resolution_clock::now();
        while(greedy_flip(data_mesh, UV, et));
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        minmax_distortions_intri(data_mesh, start_UV, res_int);
        compute_metrics(data_mesh, start_UV, res);

        fout << mesh_name << "," << type << "," << free_boundary << ",greedy," << i+to_add+1 << "," << true << "," <<  
      compute_total_energy(data_mesh, UV, et, true) << "," <<  compute_total_energy(data_mesh, UV, et, false) << ","
          << res_ext[0] << "," << res_ext[1] << "," << res_ext[2] << "," << res_int[0] << "," <<  res_int[1] << "," << res_int[2] << ","
      << res_ext[3] << "," << res_ext[4] << "," << res_ext[5] << "," << res_int[3] << "," <<  res_int[4] << "," << res_int[5] << ","
          << res[0] << "," << res[1] << "," << res[2] << "," << res[3] << "," << res[4] << "," 
          << res[5] << "," << res[6] << "," << res[7] << "," << res[8] << "," << res[9] << ","
          << duration << '\n';
        ++to_add;
      }

    start = chrono::high_resolution_clock::now();
    computeParameterization(data_mesh, V, F, UV, new_UV, free_boundary, true, '4');
    UV = new_UV;
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    minmax_distortions(V,F, UV, res_ext);
    minmax_distortions_intri(data_mesh, UV, res_int);
    compute_metrics(data_mesh, UV, res);
    fout << mesh_name << "," << type << "," << free_boundary << ",greedy," << i+to_add+1 << "," << false << "," <<  
      compute_total_energy(data_mesh, UV, et, true) << "," <<  compute_total_energy(data_mesh, UV, et, false) << ","
    << res_ext[0] << "," << res_ext[1] << "," << res_ext[2] << "," << res_int[0] << "," <<  res_int[1] << "," << res_int[2] << ","
      << res_ext[3] << "," << res_ext[4] << "," << res_ext[5] << "," << res_int[3] << "," <<  res_int[4] << "," << res_int[5] << ","
      << res[0] << "," << res[1] << "," << res[2] << "," << res[3] << "," << res[4] << "," 
      << res[5] << "," << res[6] << "," << res[7] << "," << res[8] << "," << res[9] << ","
      << duration << '\n';


  } 
    
    reset_datageo(data_mesh);
}



void test_intri(){
  const char* folderPath = "../res_data";
  DIR* directory = opendir(folderPath);
  if (directory == NULL) {
    cerr << "Failed to open directory." << endl;
    return;
  }
  fstream fout; 
  // opens an existing csv file or creates a new file. 
  fout.open("results_intri.csv", ios::out | ios::app); 

  // Read directory entries
  struct dirent* entry;
  fout << "mesh_name," << "type," << "is_free_boundary," << "flipping_order," << "iteration," << "is_flip," << 
    "int_energy," << "ext_energy," 
    << "conformal_min," << "conformal_max," << "isometric_min," << "conformal_min_int," << "conformal_max_int," << "isometric_min_int," 
    << "isometric_max," << "auhalic_min," << "authalic_max," << "isometric_max_int," << "auhalic_min_int," << "authalic_max_int,"
    << "ext_flip_percentage," << "ext_max_area_dist," << "ext_avg_area_err," << "ext_max_angle_dist," << "ext_avg_angle_error,"
    << "int_flip_percentage," << "int_max_area_dist," << "int_avg_area_err," << "int_max_angle_dist," << "int_avg_angle_error,"
    << "time"<< '\n';
  while ((entry = readdir(directory)) != NULL) {
    // Skip "." and ".." entries
    if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
        continue;
    }
    string filePath = string(folderPath) + "/" + string(entry->d_name);
    string mesh_name = string(entry->d_name);

    load_mesh(filePath, true);
    // dirichlet
    single_intri(fout, mesh_name, false, 0, '2', EnergyType::DIRICHLET);
    single_intri(fout, mesh_name, false, 1, '3', EnergyType::ASAP);
    reset_constraints();
    single_intri(fout, mesh_name, true, 2, '3', EnergyType::ASAP);
    reset_constraints();
    single_intri_arap(fout, mesh_name, false, 3, 50, 1, EnergyType::ARAP);
    single_intri_arap(fout, mesh_name, false, 4, 50, 2, EnergyType::ARAP);
    single_intri_arap(fout, mesh_name, false, 5, 50, 6, EnergyType::ARAP);
    single_intri_arap(fout, mesh_name, false, 6, 50, 10, EnergyType::ARAP);
    single_intri_arap(fout, mesh_name, true, 7, 50, 1, EnergyType::ARAP);
    single_intri_arap(fout, mesh_name, true, 8, 50, 2, EnergyType::ARAP);
    single_intri_arap(fout, mesh_name, true, 9, 50, 6, EnergyType::ARAP);
    single_intri_arap(fout, mesh_name, true, 10, 50, 10, EnergyType::ARAP);
    reset_constraints();
  }
  
  fout.close();
  // Close the directory
  closedir(directory);


}


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
    cout << "energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
		break;
  }

	case '4':
  {
    MatrixXd new_UV;
    if(UV.size()==0) {
      computeParameterization(data_mesh, V, F, UV, new_UV, false , igrad, '2'); 
      UV = new_UV;
      UV_o = new_UV;
      cout << "Initial energy: " << compute_total_energy(data_mesh, UV_o, et, true) << endl;
      reset_constraints();
    }
    //reset=true;
    unsigned its = iterations;
    while(its--){ 
      computeParameterization(data_mesh, V, F, UV, new_UV, freeBoundary, igrad, key);
      UV = new_UV;
      
      cout << "energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
    }

		break;
  }
  case '5': {
    if(UV.size()==0) {
      MatrixXd new_UV;
      computeParameterization(data_mesh, V, F, UV, new_UV, false, igrad, '2'); 
      UV = new_UV;
      cout << "Initial energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
    }
    unsigned its = iterations;
    while(its--){
      slim_parameterization(data_mesh, slimdata, UV, igrad, freeBoundary);
      cout << "energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
      cout << " energy in slim: " << slimdata.energy << endl;
      /*
      if(igrad && (its+1)%flip_granularity==0){
        int flips = 1; 
        while(flips){
          flips = flip_func(data_mesh, UV, et);
          cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
        }
      }*/
    }
    break;
  }
  case '6': {
    if(UV.size()==0) {
      MatrixXd new_UV;
      computeParameterization(data_mesh, V, F, UV, new_UV, false, igrad, '2'); 
      UV = new_UV;
      cout << "Initial energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
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
  case '7':
    {
    Eigen::VectorXi fixed_UV_indices;
    Eigen::MatrixXd fixed_UV_positions;
    igl::boundary_loop(F, fixed_UV_indices);
    igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
    igl::ARAPData arapdata;
    igl::arap_precomputation(V,F,2,fixed_UV_indices, arapdata);
    arapdata.max_iter = 10;
    MatrixXd new_UV;
    igl::harmonic(V,F, fixed_UV_indices, fixed_UV_positions, 1 , new_UV);  
    igl::arap_solve(fixed_UV_positions, arapdata, new_UV);
    UV = new_UV;
    cout << " energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
    break;
    }
  case '8':
  case '9':{
    vector<double> dists;
    MatrixXd colors_all;
    minmax_distortions(V,F,UV,dists);
    if(key == '7'){
      cout << " conformal min/max: " << dists[0] << "   " << dists[1] << endl;
      colors = colors_all.block(0,0,F.rows(),3);
    }
    if(key == '8'){
      cout << " isometric min/max: " << dists[2] << "   " << dists[3] << endl;
      colors = colors_all.block(F.rows(),0,F.rows(),3);
    }
    if(key == '9'){
      cout << " authalic min/max: " << dists[4] << "   " << dists[5] << endl;
      colors = colors_all.block(2*F.rows(),0,F.rows(),3);
    }
    break;
  }
  case 'h':
  case 'j':
  case 'k':{
    vector<double> dists;
    MatrixXd colors_all;
    minmax_distortions_intri(data_mesh,UV,dists);
    if(key == 'h'){
      cout << " conformal min/max: " << dists[0] << "   " << dists[1] << endl;
    }
    if(key == 'j'){
      cout << " isometric min/max: " << dists[2] << "   " << dists[3] << endl;
    }
    if(key == 'k'){
      cout << " authalic min/max: " << dists[4] << "   " << dists[5] << endl;
    }
    break;
  }
  case 'f': {
    if(UV.size()==0) {
      MatrixXd new_UV;
      computeParameterization(data_mesh, V, F, UV, new_UV, false, igrad, '2'); 
      UV = new_UV;
      cout << "Initial energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
    }
    UV_o = UV;
    int flips = 1; 
    while(flips){
      flips = flip_func(data_mesh, UV, et);
      cout << "  total flips: " << flips << " ;  new energy: " << compute_total_energy(data_mesh, UV, et, true) << endl;
    }
    return true;
  }
  case 'l': 
  {
    //igl::lscm(V,F,UV);
    //reset = true;
    data_mesh.intTri->flipToDelaunay();
    return true;
  }
  case 's':
  {
    viewer.data().clear();
    viewer.data().set_mesh(UV, F);
    viewer.core().align_camera_center(UV,F);
    if(intrinsic_edges) {
      intrinsicUV(data_mesh.intTri, UV, C1, C2);
      viewer.data().add_edges(C1, C2, Eigen::RowVector3d(0,0,0));
    }
    reset_datageo(data_mesh);
    return true;
  }
  case 'b':
  {
    if(intrinsic_edges) {
      intrinsicUV(data_mesh.intTri, UV, P1, P2);
      viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1,0,0));
    }
    return true;
  }
  case 'c':
  {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    if(UV.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
    }
    viewer.core().align_camera_center(V,F);
    if(intrinsic_edges) {
      intrinsicEdges(data_mesh.intTri, V, C1, C2);
      viewer.data().add_edges(C1, C2, Eigen::RowVector3d(0,0,0));
    }
    reset_datageo(data_mesh);
    return true;
  }
  case 'v':
  {
    if(intrinsic_edges){
      intrinsicEdges(data_mesh.intTri, V, P1, P2);
      viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1,0,0));
    }
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
    test_intri();
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
      ImGui::Checkbox("intrinsic edges", &intrinsic_edges);
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

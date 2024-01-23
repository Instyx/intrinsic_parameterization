import polyscope as ps
import meshio
from plyfile import PlyData, PlyElement
import numpy as np
import igl

ps.init()

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


cblack = (0.0, 0.0, 0.0)
cA = (0.082, 0.463, 0.98)
cB = (0.957, 0.733, 0.043)
cC = (0.4, 0.4 , 0.44)


def parse_edgefile(path):
    lines = open(path).read().split("\n")
    nV, nOGE, nE = [int(i) for i in lines[0].split(" ")]
    V = []
    OGE = []
    E = []
    for v in range(1,nV+1):
        V += [[float(i) for i in lines[v].split(",")]]
    for e in range(nV+1,nV+nOGE+1):
        OGE += [[int(i) for i in lines[e].split(",")]]
    for e in range(nV+nOGE+2,nV+nOGE+nE+2):
        E += [[int(i) for i in lines[e].split(",")]]
    return np.array(V), np.array(OGE), np.array(E)


def create_intrinsic_mesh(path, mesh_name):
    mesh = PlyData.read(f"{path}/{mesh_name}.int.ply")
    eV, OGE, E = parse_edgefile(f"{path}/{mesh_name}.int.edges")
    colors = []
    V = np.array([[i["x"],i["y"],i["z"]] for i in mesh["vertex"].data])
    F = np.array([i[0] for i in mesh["face"].data], dtype=int)
    color = np.array([i[1] for i in mesh["face"].data])
    UV = np.array([[i["u"],i["v"]] for i in mesh["vertex"].data])
    return dotdict(mesh=mesh,V=V,F=F,UV=UV,colors=color, eV=eV, OGE=OGE, E=E)


def render_intrinsic_mesh(name, mesh):
    ps_mesh = ps.register_surface_mesh(name, mesh.V, mesh.F, color=cC)
    ps_mesh.add_scalar_quantity("intrisic triangles", mesh.colors, enabled=False, defined_on='faces')
    ps_mesh.add_parameterization_quantity("parametrization", mesh.UV, defined_on="vertices",
                                   coords_type="unit", viz_style="checker", checker_colors=(cblack, cA), checker_size=0.05, enabled=True)
    mask = np.ones(mesh.eV.shape[0], dtype=bool)
    mask[np.array(list(set(mesh.OGE.ravel())))] = False
    og = np.array(mesh.eV)
    og[mask] = 0,0,0
    network1 = ps.register_curve_network(f"{name}_OG", og, mesh.OGE, material="flat", color=(32/255,158/255,38/255))
    mask = np.ones(mesh.eV.shape[0], dtype=bool)
    mask[np.array(list(set(mesh.E.ravel())))] = False
    intr = np.array(mesh.eV)
    intr[mask] = 0,0,0
    network2 = ps.register_curve_network(f"{name}_Flipped", intr, mesh.E, material="flat", color=(156/255,22/255,22/255))
    return ps_mesh, network1, network2


#
#
# ext_mesh = meshio.read("output/spike/dirichlet/spike_ext.obj")
# intr_mesh = meshio.read("output/spike/dirichlet/spike_iparam.obj")
# V = ext_mesh.points
# F = ext_mesh.cells[0].data
# ps_ext_mesh = ps.register_surface_mesh("Ext Mesh", V, F)
# ps_ext_mesh.add_parameterization_quantity("Normal", ext_mesh.point_data["obj:vt"], defined_on="vertices",
#                                    coords_type="unit", viz_style="checker",checker_colors=(cblack, cB), checker_size=0.05, enabled=False)
# ps_ext_mesh.add_parameterization_quantity("Normal", intr_mesh.point_data["obj:vt"], defined_on="vertices",
#                                    coords_type="unit", viz_style="checker",checker_colors=(cblack, cC), checker_size=0.05, enabled=False)
# eV, OGE, E = parse_edgefile("output/spike/symdirichlet/spike.int.edges")
# mask = np.ones(eV.shape[0], dtype=bool)
# mask[np.array(list(set(OGE.ravel())))] = False
# og = np.array(eV)
# og[mask] = 0,0,0
# ps.register_curve_network("OGintrinsic", og, OGE)
# mask = np.ones(eV.shape[0], dtype=bool)
# mask[np.array(list(set(E.ravel())))] = False
# intr = np.array(eV)
# intr[mask] = 0,0,0
# ps.register_curve_network("intrinsic", intr, E)
#
#
#
#
# mesh = PlyData.read("output/spike/dirichlet/spike.int.ply")
#
# colors = []
# V = np.array([[i["x"],i["y"],i["z"]] for i in mesh["vertex"].data])
# F = np.array([i[0] for i in mesh["face"].data], dtype=int)
# color = np.array([i[1] for i in mesh["face"].data])
# UV = np.array([[i["u"],i["v"]] for i in mesh["vertex"].data])
#
# ps_mesh = ps.register_surface_mesh("my mesh", V, F)
# ps_mesh.add_scalar_quantity("rand vals", color, enabled=True,defined_on='faces')
# ps_mesh.add_parameterization_quantity("Intrinsic", UV, defined_on="vertices",
#                                coords_type="unit", viz_style="checker",checker_colors=(cblack, cA), checker_size=0.05, enabled=True)
#
# ps.show()
# TT, TTi = igl.triangle_triangle_adjacency(F)
#
# edges = []
# all = []
# for i,f in enumerate(F):
#     for k in range(3):
#         if f[k]<f[(k+1)%3]:
#             if color[i] != color[TT[i,k]]:
#                 edges.append((f[k], f[(k+1)%3]))
#             all.append((f[k], f[(k+1)%3]))
#
#
# ps.register_curve_network("intrinsic", V, np.array(list(set(all)-set(edges))))
# ps.register_curve_network("All", V, np.array(edges))
# ps.show()
# ### Register a point cloud
# # `my_points` is a Nx3 numpy array
# ps.register_point_cloud("my points", my_points)
#
# PlyData.read("build/filename.ply").elements[0]
#
# np.array(PlyData.read("build/filename.ply")["vertex"].data)

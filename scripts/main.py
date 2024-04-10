from glob import glob
import interactive_polyscope


%gui polyscope

con, har, sym = [set(open(i).read().split("\n")[:-1]) for i in sorted(glob("*.txt"))]

con.intersection(har)
har.intersection(sym)
sym.intersection(con)

# import os
# i = list(con)[0]
# for i in con:
#     os.system(f"cp ../cutted/{i} meshes/con/{i}")
# for i in har:
#     os.system(f"cp ../cutted/{i} meshes/har/{i}")
# for i in sym:
#     os.system(f"cp ../cutted/{i} meshes/sym/{i}")


import igl
import polyscope as ps
import ipykernel
import numpy as np
import pyiparam

def load_mesh(i):
    global V,F,B
    V, F = igl.read_triangle_mesh(f"meshes/con/{i}")
    B = igl.boundary_loop(F)
    E = np.array([np.arange(len(B))[:-1],np.arange(len(B))[1:]]).T
    # ps.register_surface_mesh("mesh", V, F)
    # ps.register_curve_network("boundary", V[B], "loop")
    DG = pyiparam.create_datastructure(V, F.astype("int32"))
    uv = pyiparam.harmonic(DG, False)
    ps.register_surface_mesh("harmonic",uv, F)
    uv = pyiparam.tutte(DG)
    ps.register_surface_mesh("tutte",uv, F)
    uv = pyiparam.conformal(DG, True, False)
    ps.register_surface_mesh("conformal",uv, F)

b = np.array([2, 1])

bnd = igl.boundary_loop(F)
b[0] = bnd[0]
b[1] = bnd[int(bnd.size / 2)]

bc = np.array([[0.0, 0.0], [1.0, 0.0]])

# LSCM parametrization
_, uv = igl.lscm(v, f, b, bc)

dir(pyiparam)

cnt = 0

i = list(con)[cnt]
load_mesh(i)
cnt+=1

ps.reset_camera_to_home_view()
ps.remove_all_structures()

ps.show()

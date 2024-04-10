import numpy as np
import pyiparam
import igl


V,F = igl.read_triangle_mesh("../meshes/thingi_data/40009.off")
F = F.astype("int32")
L = igl.cotmatrix(V,F)

VT, VTi = pyiparam.vt_adjacency(F, V)
TT = pyiparam.tt_adjacency(F, VT, VTi)
VV, VVi = pyiparam.vv_adjacency(F, V, TT)
B = pyiparam.bdy_loop(F, TT, VT, VTi)

pyiparam.harmonic(L, B)
3+3

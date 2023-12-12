# Exploring intrinsic triangulations for optimizing triangle mesh parameterizations
This repository contains the code used for the bachelor thesis "Exploring intrinsic triangulations for optimizing triangle mesh parameterizations" supervised by Prof. Dr. Marc Alexa and Prof. Dr. Olga Sorkine-Hornung. 

Some parts are adapted from my submission to the assigment 4 for the course Shape Modeling and Geometry Processing in 2022 Autumn Semester at ETH Zürich[^1] and from the code used for iARAP (ARAP Revisited: Discretizing the Elastic Energy using Intrinsic Voronoi Cells)[^2].

The implementation is in the folder iparam.

To clone and build:

```
git clone --recursive git@github.com:Instyx/intrinsic_parameterization.git
mkdir build
cd build 
cmake ..
make -j4
```
The executable `iparam` is built.

## Usage

### Menu Items

On the menu left, under the category "Parameterization":
- `Free boundary`[^3]
    * Unmarked: fixed boundary 
    * Marked  : free boundary
- `ARAP/SLIM iterations`: the number of iterations used for the global/local approach when minimized w.r.t. ARAP, symmetric Dirichlet energies
- `Flip granularity`: represents after how many iterations of the global/local approaches intrinsic flipping is performed when computing intrinsic parameterization
- `intrinsic grad`
    * Unmarked: Use extrinsic triangulation for mesh parameterization
    * Marked  : Use intrinsic triangulation for mesh parameterization
- `intrinsic edges`
    * Unmarked: Intrinsic edges are not displayed when keys `s,b,c,v` are pressed
    * Marked  : Intrinsic edges are displayed as red edges when keys `s,b,c,v` are pressed
- select one of the energies Dirichlet, symmetric Dirichlet, ASAP, ARAP for parameterization and computing total energy
- select one of the flipping orders Edge Order, Greedy, Random, Heuristic for the intrinsic flipping algorithm

### Keys

- `1`, `2`, `3`, `4`, `5` are used for parameterization respectively for minimizing spring energy (uniform laplacian), Dirichlet energy (cotangent Laplacian), ASAP energy (LSCM), ARAP energy (global/local approach), symmetric Dirichlet energy (SLIM) and the result is saved in the matrix `UV`.[^4] 
- For ARAP and symmetric Dirichlet energies, respectively the keys `e`, `r` compute the parameterization but save the result in a different matrix `UV_o`
- `f` for intrinsic flipping algorithm. The UV mapping that is used for the flipping algorithm is saved in `UV_o`
- `space` for switching between the UV domain and 3D mesh surface.
- `s`, `b` for showing the UV domain using respectively `UV`, `UV_o` and extrinsic geometry. If "intrinsic edges" is marked, then the intrinsic edges are also visualized in th UV domain with different color. 
- `c`, `v` for showing the textured 3D mesh and printing the metrics respectively for `UV`, `UV_o`. If "intrinsic edges" is marked, then the intrinsic edges are also visualized on the mesh surface with different color. 
- `l` for resetting intrinsic triangulation to input extrinsic triangulation

### Example


 `./iparam ../res_data/Octo_cut2.obj` to run the code with the mesh `Octo_cut2.obj`.


- Select `ASAP` from the menu and press `3`. 
- Then the parameterization is computed, the textured mesh is shown and the total energy normalized with respect to the mesh area is printed out.
- Press `f` to apply intrinsic flipping algorithm.
- Mark 'intrinsic grad' from the menu left and press `3` again to compute the intrinsic parametrization using the intrinsic triangulation resulting from the intrinsic flipping algorithm.
- Press `s`, `b` to compare the parametrization with intrinsic parameterization.

**OR**

- Select `ARAP` from the menu, adjust the `ARAP\SLIM iterations` and `Flip granularity`.
- Mark `intrinsic grad`.
- Press `4` to compute intrinsic parameterization.
- Mark `intrinsic edges` and press `s` to visualize intrinsic edges on top of the extrinsic geometry. Unmark `Wireframe` from the menu left, if you only want to see the intrinsic edges
- Adjust the menu and press `e` for computing another ARAP parameterization. The result is saved in another matrix and data structure for intrinsic triangulation.
- To compare the 2 different methods, press `s`, `b`, `c`, `v`.

[^1]: https://igl.ethz.ch/teaching/shape-modeling/sm2022/
[^2]: https://cybertron.cg.tu-berlin.de/projects/iARAP/
[^3]: SLIM always uses free boundary, and the initial UV parameterizations for minimizing ARAP and symmetric Dirichlet are computed with fixed boundary.
[^4]: The keys `4`, `5` applies the intrinsic flipping algorithm when the 'intrinsic grad' is marked. Other keys require using `f` to update the intrinsic triangulation. 

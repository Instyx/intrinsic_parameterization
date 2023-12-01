# Exploring intrinsic triangulations for optimizing triangle mesh parameterizations

This is the code that is used for the bachelor thesis "Exploring intrinsic triangulations for optimizing triangle mesh parameterizations" supervised by Prof. Dr. Marc Alexa and Prof. Dr. Olga Sorkine-Hornung. 

Some parts are adapted from my submission to the assigment 4 for the course Shape Modeling and Geometry Processing in 2022 Autumn Semester at ETH Zürich and from the code that is used for iARAP (ARAP Revisited: Discretizing the Elastic Energy using Intrinsic Voronoi Cells).

The implementation is in the folder iparam.

To clone and build:

```
git clone --recursive git@github.com:Instyx/intrinsic_parameterization.git
mkdir build
cd build 
cmake ..
make -j6
```

## Usage

On the menu left, under the category "Parameterization":
- Free boundary
    * Unmarked: fixed boundary 
    * Marked  : free boundary
- ARAP/SLIM iterations: the number of iterations used for the global/local approach when minimized w.r.t. ARAP, symmetric Dirichlet energies
- Flip granularity: represents after how many iterations of the global/local approaches intrinsic flipping is performed when computing intrinsic parameterization
- intrinsic grad
    * Unmarked: Use extrinsic triangulation for mesh parameterization
    * Marked  : Use intrinsic triangulation for mesh parameterization
- intrinsic edges
    * Unmarked: Intrinsic edges are not displayed when keys `s,b,c,v` are pressed
    * Marked  : Intrinsic edges are displayed as red edges when keys `s,b,c,v` are pressed




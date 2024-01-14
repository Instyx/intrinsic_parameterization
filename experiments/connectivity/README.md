Tolerance tol is used for checking co-planarity comparison. Delaunay means all possible co-planar Delaunay flips are performed.
Orig stands for original, mod stands for modified (after the flips).

You can click on the images to view full size.

| Method      | Co-planar Flips         | Energy Orig | Energy Mod | Texture Orig                                                              | Texture Mod                                                              |
| :---------: | :------------:          | :---------: | :--------: | :----------:                                                              | :---------:                                                              |
| ARAP, free  | 50206, tol=1e-8         | 0.28221     | 0.28228    | <img align="center" src="./superman_arap1_orig_tol8.png" width="300">     | <img align="center" src="./superman_arap1_mod_tol8.png" width="300">     |
| ASAP        | 38766, tol=1e-9         | 5.13079     | 5.1197     | <img align="center" src="./superman_asap1_orig_tol9.png" width="300">     | <img align="center" src="./superman_asap1_mod_tol9.png" width="300">     |
| ASAP        | 13209, tol=1e-8         | 4.62307     | 4.62755    | <img align="center" src="./gargoyle_asap1_orig_tol8.png" width="300">     | <img align="center" src="./gargoyle_asap1_mod_tol8.png" width="300">     |
| Dirichlet   | 22575,  tol=1e-8        | 6.30999     | 6.31236    | <img align="center" src="./vaselion_dirichlet_orig_tol8.png" width="300"> | <img align="center" src="./vaselion_dirichlet_mod_tol8.png" width="300"> |
| ARAP, fixed | 22575,  tol=1e-8        | 2.90685     | 2.9089     | <img align="center" src="./vaselion_arap0_orig_tol8.png" width="300">     | <img align="center" src="./vaselion_arap0_mod_tol8.png" width="300">     |
| ARAP, free  | 22575,  tol=1e-8        | 0.530081    | 0.530254   | <img align="center" src="./vaselion_arap1_orig_tol8.png" width="300">     | <img align="center" src="./vaselion_arap1_mod_tol8.png" width="300">     |
| Dirichlet   | 464, Delaunay, tol=1e-9 | 0.00012     | 0.00012    | <img align="center" src="./stripe_dirichlet_orig_tol9_d.png" width="300"> | <img align="center" src="./stripe_dirichlet_mod_tol9_d.png" width="300"> |
| ARAP, fixed | 2166, tol=1e-8          | 1.02528     | 1.02520    | <img align="center" src="./stripe_arap0_orig_tol8.png" width="300">       | <img align="center" src="./stripe_arap0_mod_tol8.png" width="300">       |
| Sym. Dirichlet  | 2166, tol=1e-8            | 5.24653| 5.24656    | <img align="center" src="./stripe_symdir_orig_tol8.png" width="300"> | <img align="center" src="./stripe_symdir_mod_tol8.png" width="300"> |

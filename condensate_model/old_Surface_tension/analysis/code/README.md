### Code to analyze constant concentration simulations of the Nup62 complex

- write_final_box.py - Script to write a .gro file (final.gro) of the last frame of a simulation, centered at the largest cluster.
- cluster_FG_box_mindist.py - Script to calculate density and cluster size from constant concentration simulations.
- surface_tension.py - Script to surface tension from constant concentration simulations. (uses Eq. 27 and Eq. 28 of Benayad, Z et al. Simulation of FUS Protein Condensates with an Adapted Coarse-Grained Model, J. Chem. Theory Comput., 2021)
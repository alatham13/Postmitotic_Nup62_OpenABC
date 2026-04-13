### This directory contains the inputs for running slab simulations of the Nup62-Nup58-Nup54 complex in HXD. HXD is modeled by modifying the interaction strength, with epsilon=0.9 corresponding to 1.5% HXD and epsilon=0.7 corresponding to 5% HXD.

### These simulations are done in 4 steps (Note, simulations can take up to 11 days to run on one A40 GPU):

1. To NPT equilibration simulations:
- python run_npt.py

2. To NPT equilibration simulations:
- python run_eq1.py

3. To start production simulations:
- python run_prod.py

4. To finish production simulations (repeat until simulation is finished):
- python run_restart.py

### To analyze the simulations:

1. calculate the density profile of the protein as a function of Z-axis (where XXX is the name of the protein):
- python ../../../../Slab_simulations/analysis/code/cluster_FG_mindist.py XXX_system.xml

2. write the final configuration to a .gro file (where XXX is the name of the protein):
- python ../../../../Slab_simulations/analysis/code/write_final.py XXX_system.xml

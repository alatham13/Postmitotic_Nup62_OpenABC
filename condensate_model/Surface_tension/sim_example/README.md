### This directory contains the inputs for running constant concentration simulations of the Nup62 complexes.
Final results include the results of 15 simulations with different random seeds.

To start simulations:
- python run_prod.py

To continue unfinished simulations:
- python run_restart.py

To analyze simulations (where XXX is the name of the system):
- python ../analysis/code/write_final_box.py XXX.xml
- python ../analysis/code/cluster_FG_box_mindist.py XXX.xml
- python ../analysis/code/surface_tension.py XXX.xml

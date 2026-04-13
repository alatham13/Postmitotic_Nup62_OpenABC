# Function to calculate radius of gyration
# Written by Andrew Latham
# Note the inputs at top.

import sys
import os
import math
import numpy
import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
import xml.etree.ElementTree as ET

# Inputs
# GRO file used to run the simulation
GRO='npt.gro'
# XTC file used to run the simulation
XTC=sys.argv[1]
# number of equilibration steps
eq=-1

# exclude=1 creates exclusion for the same chain
def calc_rdf(GRO_file,xtc_file,start):
    u1 = mda.Universe(GRO_file, xtc_file)

    # length of trajectory
    timesteps=len(u1.trajectory)
    print('Number of timesteps: '+str(timesteps))
    # select the entire protein
    sel_all = u1.select_atoms("protein")


    Rg=[]
    for ts in u1.trajectory:
        print('frame: ' + str(ts.frame))
        if ts.frame < start:
            pass
        else:
            # Select the entire chain and calculate the radius of gyration
            Rg_temp=sel_all.radius_of_gyration(wrap=False)
            Rg.append(Rg_temp)

    Rg=numpy.asarray(Rg)*0.1 # convert Ang to nm

    return Rg



Rg1=calc_rdf(GRO,XTC,eq)

numpy.savetxt('Rg_v1.txt',Rg1)




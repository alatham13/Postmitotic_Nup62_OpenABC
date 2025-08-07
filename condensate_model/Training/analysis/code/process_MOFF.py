# Function to calculate radius of gyration and distance between 2nd and 2nd-last residue (to compare to FRET)
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
# pdb file used to run the simulation
PDB=sys.argv[1]
# xml file used to define the system
SYSTEM=sys.argv[2]
# dcd file used to run the simulation
DCD='output.dcd'
# number of equilibration steps
eq=-1

# function to find the mass of particles from the system XML file
def find_mass(xml_file):
    # read in system XML file
    tree=ET.parse(xml_file)
    root=tree.getroot()
    # list for masses of each particle
    mass_list=[]
    # go over all branches of the XML
    for branch in root:
        # the branch that stores mass is named 'Particles'
        if branch.tag=='Particles':
            for leaf in branch:
                # append these masses to the mass_list
                mass_list.append(float(leaf.attrib['mass']))
    return mass_list

# exclude=1 creates exclusion for the same chain
def calc_rdf(pdb_file,mass_xml,dcd_file,start):
    u1 = mda.Universe(pdb_file, dcd_file)

    # length of trajectory
    timesteps=len(u1.trajectory)
    print('Number of timesteps: '+str(timesteps))
    # select the entire chain
    sel1_str = 'all'
    sel_all = u1.select_atoms(sel1_str)
    # select the start and end of the chain
    N_atoms=len(sel_all)
    begin_str='index 1'
    end_str='index '+str(N_atoms-2)
    begin=u1.select_atoms(begin_str)
    end=u1.select_atoms(end_str)

    # set mass of atoms
    atom_masses = find_mass(mass_xml)
    for i in range(N_atoms):
        u1.atoms[i].mass = atom_masses[i]


    Rg=[]
    FRET=[]
    for ts in u1.trajectory:
        print('frame: ' + str(ts.frame))
        if ts.frame < start:
            pass
        else:
            # Select the entire chain and calculate the radius of gyration
            Rg_temp=sel_all.radius_of_gyration(wrap=False)
            Rg.append(Rg_temp)
            # Select the start and end of the chain to calculate the FRET distance
            FRET_temp=dist(begin,end)
            FRET.append(FRET_temp[2])

    Rg=numpy.asarray(Rg)*0.1 # convert Ang to nm
    FRET=numpy.asarray(FRET)*0.1

    return Rg,FRET



Rg1,FRET1=calc_rdf(PDB,SYSTEM,DCD,eq)

numpy.savetxt('Rg_v2.txt',Rg1)
numpy.savetxt('FRET_v2.txt',FRET1)




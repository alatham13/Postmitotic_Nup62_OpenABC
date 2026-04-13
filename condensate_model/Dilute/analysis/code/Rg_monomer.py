# Script to calculate polymer properties from monomeric protein simulations.
# Calculates:
# the radius of gyration, in nm
# Written by Andrew Latham
# Note the inputs at top. These need to be customized for the system.

# Output files:
# Rg.txt, Ree.txt, asphericity.txt, scaling_exponent.txt, A_0.txt
# Each is written as a matrix where rows correspond to each time frame and columns correspond to protein copy numbers

import sys
import os
import math
import numpy
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt


# Inputs
# dcd file used to run the simulation
DCD='prod.dcd'
# pdb file used to set the topology
PDB='start.pdb'
# xml file used to define the system
SYSTEM=sys.argv[1]
# number of protein chains in the simulation box
nchain=1
# number of equilibration frames
eq=50000

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



def polymer_ana(dcd_file,pdb_file,mass_xml,nchain,start):
    u1 = mda.Universe(pdb_file,dcd_file)

    # caclulate number of atoms per protein chain and number of timesteps
    protien=u1.select_atoms("all")
    print(protien.atoms)
    N1=int(len(protien)/nchain)
    print('Number of protein atoms per chain: '+str(N1))
    timesteps=len(u1.trajectory)-(start)
    print('Number of timesteps: '+str(timesteps))
    N = len(u1.atoms)

    # set mass of atoms
    atom_masses=find_mass(mass_xml)
    for i in range(N):
        u1.atoms[i].mass=atom_masses[i]

    # counter for number of timesteps, initial empty pas
    count=0
    Rg = numpy.zeros((timesteps, nchain))

    for ts in u1.trajectory:
        print('frame: ' + str(ts.frame))
        if ts.frame < start:
            pass
        else:
            # calculate com of each protein chain
            for chain in range(nchain):
                # select atoms in each chain
                index1 = chain * N1
                index2 = ((chain + 1) * N1)
                atoms1 = u1.atoms[index1:index2]
                Rg[count,chain]=atoms1.radius_of_gyration(wrap=False)

            count = count + 1

    return Rg


gyr=polymer_ana(DCD,PDB,SYSTEM,nchain,eq)

numpy.savetxt('Rg.txt',gyr/10,header='Rg (nm)')


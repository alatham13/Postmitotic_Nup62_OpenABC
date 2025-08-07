# Function to calculate density and cluster size from constant concentration simulations
# Written by Andrew Latham
# Note the inputs at top. These need to be customized for the system.
# Here, we assume 32 molecules, with 9999 equilibration frames (i.e. 10000 total frames), and a minimum distance cutoff of 10 Ang

# Output files:
# cluster_size_mindist.txt - size of the largest 4 clusters, in number of molecules
# dist_from_center.txt - histogram of radial protein density, centered at the largest cluster

import sys
import os
import math
import numpy
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import xml.etree.ElementTree as ET

# Inputs
# dcd file used to run the simulation
DCD='prod.dcd'
# pdb file used to set the topology
PDB='start.pdb'
# xml file used to define the system
SYSTEM=sys.argv[1]
# number of protein chains in the simulation box
nchain=32
# number of equilibration frames
eq=5000
# distance cutoff
cut=10




# function that runs depth first search algorithm on a graph
def dfs_iterative(graph, start):
    stack, path = [start], []
    while stack:
        vertex = stack.pop()
        if vertex in path:
            continue
        path.append(vertex)
        if vertex in graph:
            for neighbor in graph[vertex]:
                stack.append(neighbor)
    return path


def mat_to_dict(mat):
    dict1={}
    for i in range(0,len(mat)):
        for j in range(0,len(mat[i])):
            if i in dict1 and mat[i][j]==1:
                dict1[i].append(j)
            elif i not in dict1 and mat[i][j]==1:
                dict1[i] = [j]
    return dict1

def wrap_coord(pos,side):
    while pos<0:
        pos=pos+side
    while pos>side/2:
        pos=pos-side
    return pos

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

# function to find the X/Y box edges from the topology XML file
def find_box(xml_file,find_z=False,default_Z=5000):
    # read in system XML file
    tree=ET.parse(xml_file)
    root=tree.getroot()
    box=numpy.zeros(6)
    # go over all branches of the XML
    for branch in root:
        # the branch that stores box length is named 'PeriodicBoxVectors'
        if branch.tag=='PeriodicBoxVectors':
            # multiply by 10 to convert nm to Ang.
            box[0]=float(branch[0].attrib['x'])*10
            box[1] = float(branch[1].attrib['y'])*10
            # Z has not yet been resized. Save only if asked
            if find_z:
                box[2] = float(branch[2].attrib['z'])*10
            else:
                box[2]=default_Z
    # set angles for box edges. This assumes a rectangular box
    box[3]=90
    box[4]=90
    box[5]=90
    return box


def contact_mat(dcd_file,pdb_file,mass_xml,nchain,start,cutoff):
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

    # read in box dimensions
    box=find_box(mass_xml,find_z=True)

    # counter for number of timesteps, initial empty pas
    count=0
    path_tot_timestep = []
    cluster_size = numpy.zeros((timesteps, 4))
    X_pos = numpy.zeros((N,timesteps))
    Y_pos = numpy.zeros((N,timesteps))
    Z_pos = numpy.zeros((N,timesteps))
    overall_pos = numpy.zeros((N,timesteps))


    for ts in u1.trajectory:
        print('frame: ' + str(ts.frame))
        if ts.frame < start:
            pass
        else:
            com = numpy.zeros((nchain, 3))
            print(box)
            mass = numpy.zeros((nchain, 1))
            # convert to adjacency matrix
            mat = numpy.zeros((nchain, nchain))
            # calculate com of each protein chain
            for i in range(0, nchain-1):
                # select atoms in each chain
                index1 = i * N1
                index2 = ((i + 1) * N1)
                atoms1 = u1.atoms[index1:index2]
                com[i, :] = atoms1.center_of_mass()
                mass[i] = atoms1.total_mass()
                for j in range(i+1,nchain):
                    # find all pairs of chain
                    index3 = j * N1
                    index4 = ((j + 1) * N1)
                    atoms2 = u1.atoms[index3:index4]
                    # caclulate total distance matrix
                    dist2 = distance_array(atoms1.positions, atoms2.positions, box)
                    # use minimum of distance matrix to calculate contact matrix
                    dist = numpy.amin(dist2[:])
                    if dist < cutoff:
                        mat[i, j] = mat[i, j] + 1
                        mat[j, i] = mat[i, j]
            # need com/mass for all chains. Add for the last chain here
            chain_index=nchain-1
            index1 = chain_index * N1
            index2 = ((chain_index + 1) * N1)
            atoms1 = u1.atoms[index1:index2]
            com[chain_index, :] = atoms1.center_of_mass()
            mass[chain_index] = atoms1.total_mass()

            # Calculate dfs at each timestep from matrix
            graph2 = mat_to_dict(mat)
            visited = []
            path_tot = []
            for i in range(0, len(mat)):
                flag = 0
                for k in range(0, len(visited)):
                    if visited[k] == i:
                        flag = 1

                if flag == 0:
                    path1 = []
                    path1 = dfs_iterative(graph2, i)
                    for j in range(0, len(path1)):
                        visited.append(path1[j])
                    path_tot.append(path1)
                path_tot_timestep.append(path_tot)
                #print(path_tot_timestep)

            # find atoms in largest cluster
            n = len(path_tot)
            max = 0
            l2 = -1
            index = -1
            #print(path_tot)
            l_list=[]
            for i in range(0, n):
                l = len(path_tot[i])
                l_list.append(l)
                if l > l2:
                    # reset old variables
                    l2 = l
                    index = i
            atoms_in_cluser = path_tot[index]
            mass_tot = 0
            for i in range(0, len(atoms_in_cluser)):
                index = atoms_in_cluser[i]
                mass_tot = mass_tot + mass[index]
            # X ----------------------------------------------------------------------------
            comX2 = (mass[:, 0] / mass_tot) * numpy.cos((com[:, 0] / box[0]) * 2 * numpy.pi)
            comX3 = (mass[:, 0] / mass_tot) * numpy.sin((com[:, 0] / box[0]) * 2 * numpy.pi)
            # Y ----------------------------------------------------------------------------
            comY2 = (mass[:, 0] / mass_tot) * numpy.cos((com[:, 1] / box[1]) * 2 * numpy.pi)
            comY3 = (mass[:, 0] / mass_tot) * numpy.sin((com[:, 1] / box[1]) * 2 * numpy.pi)
            # Z ----------------------------------------------------------------------------
            comZ2 = (mass[:, 0] / mass_tot) * numpy.cos((com[:, 2] / box[2]) * 2 * numpy.pi)
            comZ3 = (mass[:, 0] / mass_tot) * numpy.sin((com[:, 2] / box[2]) * 2 * numpy.pi)


            # Calculate COM of cluster. Use angles to avoid issues with periodicity
            X1 = 0
            X2 = 0
            Y1 = 0
            Y2 = 0
            Z1 = 0
            Z2 = 0
            # l_list: list of cluster sizes
            l_list=numpy.asarray(l_list)
            l_list[::-1].sort()
            l_list.resize((4))
            cluster_size[count,:]=l_list
            print(cluster_size[count,:])
            for i in range(0, l2):
                index = atoms_in_cluser[i]
                # X ----------------------------------------------------------------------------
                X1 = X1 + comX2[index]
                X2 = X2 + comX3[index]
                # Y ----------------------------------------------------------------------------
                Y1 = Y1 + comY2[index]
                Y2 = Y2 + comY3[index]
                # Z ----------------------------------------------------------------------------
                Z1 = Z1 + comZ2[index]
                Z2 = Z2 + comZ3[index]
            # X ----------------------------------------------------------------------------
            X1 = X1 / l2
            X2 = X2 / l2
            thetaX = numpy.arctan2(-1 * X2, -1 * X1) + numpy.pi
            X_com = (box[0] / (2 * numpy.pi)) * thetaX
            # Y ----------------------------------------------------------------------------
            Y1 = Y1 / l2
            Y2 = Y2 / l2
            thetaY = numpy.arctan2(-1 * Y2, -1 * Y1) + numpy.pi
            Y_com = (box[1] / (2 * numpy.pi)) * thetaY
            # Z ----------------------------------------------------------------------------
            Z1 = Z1 / l2
            Z2 = Z2 / l2
            thetaZ = numpy.arctan2(-1 * Z2, -1 * Z1) + numpy.pi
            Z_com = (box[2] / (2 * numpy.pi)) * thetaZ


            #calculate wrapped list of atomistic distances to the largest cluster
            for i in range(0,N):
                atom = u1.atoms[i]
                # X ----------------------------------------------------------------------------
                X = atom.position[0]
                X_diff = X - X_com
                X_final = wrap_coord(X_diff, box[0])
                X_pos[i, count] = X_final
                # Z ----------------------------------------------------------------------------
                Y = atom.position[1]
                Y_diff = Y - Y_com
                Y_final = wrap_coord(Y_diff, box[1])
                Y_pos[i, count] = Y_final
                # Z ----------------------------------------------------------------------------
                Z = atom.position[2]
                Z_diff=Z-Z_com
                Z_final = wrap_coord(Z_diff, box[2])
                Z_pos[i,count]=Z_final
                # overall
                overall_pos[i,count]=numpy.sqrt(X_final**2+Y_final**2+Z_final**2)

            count = count + 1

    # Calculate density in X, Y, and Z separately, for debugging
    """# calculate the protein density as a function of X-position
    nbins = 100
    maxX = box[0] / 2
    minX = -1 * box[0] / 2
    dX = (maxX - minX) / nbins
    Xhist = numpy.zeros((nbins, 2))
    # set X-axis
    for i in range(0, nbins):
        Xhist[i, 0] = minX + dX * i + dX / 2
    Xhist[0, 0] = Xhist[0, 0] - dX / 2
    Xhist[nbins - 1, 0] = Xhist[nbins - 1, 0] + dX / 2
    # calculate Xhistogram of atomic masses
    for i in range(0, N):
        for j in range(0, count):
            bin = int((X_pos[i, j] - minX) / dX)
            # atoms in protein
            Xhist[bin, 1] = Xhist[bin, 1] + u1.atoms.masses[i]

    # calculate the protein density as a function of Y-position
    nbins = 100
    maxY = box[1] / 2
    minY = -1 * box[1] / 2
    dY = (maxY - minY) / nbins
    Yhist = numpy.zeros((nbins, 2))
    # set Y-aYis
    for i in range(0, nbins):
        Yhist[i, 0] = minY + dY * i + dY / 2
    Yhist[0, 0] = Yhist[0, 0] - dY / 2
    Yhist[nbins - 1, 0] = Yhist[nbins - 1, 0] + dY / 2
    # calculate Yhistogram of atomic masses
    for i in range(0, N):
        for j in range(0, count):
            bin = int((Y_pos[i, j] - minY) / dY)
            # atoms in protein
            Yhist[bin, 1] = Yhist[bin, 1] + u1.atoms.masses[i]

    # calculate the protein density as a function of Z-position
    nbins = 100
    maxZ = box[2] / 2
    minZ = -1 * box[2] / 2
    dZ = (maxZ - minZ) / nbins
    Zhist = numpy.zeros((nbins, 2))
    # set Z-axis
    for i in range(0, nbins):
        Zhist[i, 0] = minZ + dZ * i + dZ / 2
    Zhist[0, 0] = Zhist[0, 0] - dZ / 2
    Zhist[nbins - 1, 0] = Zhist[nbins - 1, 0] + dZ / 2
    # calculate Zhistogram of atomic masses
    for i in range(0, N):
        for j in range(0, count):
            bin = int((Z_pos[i, j] - minZ) / dZ)
            # atoms in protein
            Zhist[bin, 1] = Zhist[bin, 1] + u1.atoms.masses[i]"""

    # calculate the protein density as a function of Overall-distance
    nbins = 100
    maxO = numpy.sqrt(box[0] ** 2 + box[1] ** 2 + box[2] ** 2) / 2
    minO = 0
    dO = (maxO - minO) / nbins
    Ohist = numpy.zeros((nbins, 2))
    # set O-axis
    for i in range(0, nbins):
        Ohist[i, 0] = minO + dO * i + dO / 2
    Ohist[0, 0] = Ohist[0, 0] - dO / 2
    Ohist[nbins - 1, 0] = Ohist[nbins - 1, 0] + dO / 2
    # calculate Ohistogram of atomic masses
    for i in range(0, N):
        for j in range(0, count):
            bin = int((overall_pos[i, j] - minO) / dO)
            # atoms in protein
            Ohist[bin, 1] = Ohist[bin, 1] + u1.atoms.masses[i]

    # Header for output file
    key_string = '\tr\t\t\tProtein'


    #return cluster_size,Xhist,Yhist,Zhist,Ohist,key_string
    return cluster_size, Ohist, key_string



cluster_size,dist,key_string=contact_mat(DCD,PDB,SYSTEM,nchain,eq,cut)

numpy.savetxt('cluster_size_mindist.txt',cluster_size)
"""numpy.savetxt('proteinX.txt',X,header=key_string)
numpy.savetxt('proteinY.txt',Y,header=key_string)
numpy.savetxt('proteinZ.txt',Z,header=key_string)"""
numpy.savetxt('dist_from_center.txt',dist,header=key_string)




# run simulation of a protein.
import numpy as np
import pandas as pd
import os
import sys
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

sys.path.insert(0, "~/Postmitotic_Nup62_OpenABC")
from openabc2.forcefields.parsers import MOFFParser, MRGdsDNAParser
from openabc2.forcefields import MOFFMRGModel
from openabc2.utils.insert import insert_molecules

tot_t=1000000000
dt=100000

# set simulation platform
platform_name = 'CUDA'

print('Successful import')

# set the maximum allowed angle
max_angle=130*np.pi/180

# start from predicted PDB
FG_parser = MOFFParser.from_atomistic_pdb('Nup62_Nup58_Nup54_AF.pdb', 'Nup62_Nup58_Nup54_CA.pdb')

# Nup62_Nup58_Nup54 complex has an ordered domain that includes all 3 proteins, as well as 2 seperate helices, pass these as pairs
old_native_pairs = FG_parser.native_pairs.copy()
new_native_pairs = pd.DataFrame(columns=old_native_pairs.columns)
OD1=np.arange(110,507)
OD2=np.arange(752,932)
OD3=np.arange(1435,1627)
OD_tot=np.concatenate((OD1,OD2,OD3))
print(OD_tot)
helix1=np.arange(743,751)
helix2=np.arange(937,951)
helix3=np.arange(967,979)
for i, row in old_native_pairs.iterrows():
    a1, a2 = int(row['a1']), int(row['a2'])
    if a1 > a2:
        a1, a2 = a2, a1
    flag1=( (a1 in helix1) and (a2 in helix1) )
    flag2=( (a1 in helix2) and (a2 in helix2) )
    flag3=( (a1 in helix3) and (a2 in helix3) )
    flag4=( (a1 in OD_tot) and (a2 in OD_tot) )
    if flag1 or flag2 or flag3 or flag4:
        new_native_pairs.loc[len(new_native_pairs.index)] = row
FG_parser.native_pairs = new_native_pairs
FG_parser.parse_exclusions() # update exclusions based on the new native pairs

# check interaction parameters
print(FG_parser.native_pairs)

# loop to remove large angles
# save the original angles to CSV. For debugging
FG_parser.protein_angles.to_csv('old_angles.csv')
# copy the angle potential
old_angles=FG_parser.protein_angles.copy()
new_angles = pd.DataFrame(columns=old_angles.columns)
count_bad_angle=0
# loop over all angles
for i, row in old_angles.iterrows():
    # round down if greater than max_angle
    if max_angle<row['theta0']:
        #print('Warning!!! Protein angle value is greater than large angle: '+str(row['theta0']*180/np.pi))
        #print(row)
        count_bad_angle+=1
        row['theta0']=max_angle
        new_angles.loc[len(new_angles.index)] = row
    # Leave alone otherwise
    else:
        new_angles.loc[len(new_angles.index)] = row
        print(row)
# report warning based on the number of large angles
print('Warning!!! Found '+str(count_bad_angle)+' angles greater than '+str(max_angle*180/np.pi)+' degrees.')
print('This may effect the numeric stability of your simulation. Consider rounding down large angles.')
# replace angles with new angles
FG_parser.protein_angles=new_angles
FG_parser.protein_angles.to_csv('new_angles.csv')

# setup initial topology
n_mol = 32
box_a=81.00
box_b=81.00
box_c=81.00
insert_molecules('Nup62_Nup58_Nup54_CA.pdb', 'start.pdb', n_mol, box=[box_a, box_b, box_c])


# set up protein model
condensate = MOFFMRGModel()
# add protein to model
for i in range(n_mol):
    # append multiple FG parser instances
    condensate.append_mol(FG_parser)
# setup initial topology
top = app.PDBFile('start.pdb').getTopology()
condensate.create_system(top, box_a=box_a, box_b=box_b, box_c=box_c, remove_cmmotion=True)
# setup variables necessary for simulations
salt_concentration = 162*unit.millimolar
temperature = 300*unit.kelvin
# add force terms to simulation
condensate.add_protein_bonds(force_group=1)
condensate.add_protein_angles(force_group=2)
condensate.add_protein_dihedrals(force_group=3)
condensate.add_native_pairs(force_group=4)
condensate.add_contacts(force_group=5,scale=1.05)
condensate.add_elec_switch(salt_concentration, temperature, force_group=6)
# setup integrator
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
# setup initial coordinates
init_coord = app.PDBFile('start.pdb').getPositions()
# setup simulation
condensate.set_simulation(integrator, platform_name, init_coord=init_coord)
# save system to xml
system=condensate.simulation.context.getSystem()
with open('Nup62_Nup58_Nup54.xml', 'w') as output:
    output.write(mm.XmlSerializer.serialize(system))
# perform energy minimization
condensate.simulation.minimizeEnergy()
output_interval = dt
output_dcd = 'prod.dcd'
# Run short simulation
condensate.add_reporters(output_interval, output_dcd)
condensate.simulation.reporters.append(app.checkpointreporter.CheckpointReporter('prod.cpt',output_interval))
condensate.simulation.reporters.append(app.StateDataReporter('prod.csv', output_interval, step=True, potentialEnergy=True, temperature=True, volume=True, density=True, speed=True, time=True))
condensate.simulation.context.setVelocitiesToTemperature(temperature)
condensate.simulation.step(tot_t)

# print final box vectors
state = condensate.simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
with open('prod.xml', 'w') as f:
    f.write(mm.XmlSerializer.serialize(state))
# save final system
condensate.save_system('Nup62_Nup58_Nup54_system.xml')

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

tot_t=10000000
dt=int(tot_t/100)

# set simulation platform
platform_name = 'CUDA'

print('Successful import')

# set the maximum allowed angle
max_angle=130*np.pi/180

# start from predicted PDB
FG_parser = MOFFParser.from_atomistic_pdb('Nup58_CT_AF.pdb', 'Nup58_CT_CA.pdb')

# Nup58 has an OD and 2 helicies - pass these as native pairs
old_native_pairs = FG_parser.native_pairs.copy()
new_native_pairs = pd.DataFrame(columns=old_native_pairs.columns)
helix1=np.arange(243,425)-418 # 244-425
helix2=np.arange(430,444)-418 # 431-444
helix3=np.arange(460,471)-418 # 461-471
for i, row in old_native_pairs.iterrows():
    a1, a2 = int(row['a1']), int(row['a2'])
    if a1 > a2:
        a1, a2 = a2, a1
    flag1=( (a1 in helix1) and (a2 in helix1) )
    flag2=( (a1 in helix2) and (a2 in helix2) )
    flag3=( (a1 in helix3) and (a2 in helix3) )
    if flag1 or flag2 or flag3:
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
n_mol = 100
box_a=100
box_b=100
box_c=100
insert_molecules('Nup58_CT_CA.pdb', 'start.pdb', n_mol, box=[box_a, box_b, box_c])


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
temperature = 150*unit.kelvin
# add force terms to simulation
condensate.add_protein_bonds(force_group=1)
condensate.add_protein_angles(force_group=2)
condensate.add_protein_dihedrals(force_group=3)
condensate.add_native_pairs(force_group=4)
condensate.add_contacts(force_group=5,scale=1.05)
condensate.add_elec_switch(salt_concentration, 300*unit.kelvin, force_group=6)
pressure = 1*unit.bar
condensate.system.addForce(mm.MonteCarloBarostat(pressure, temperature))
# setup integrator
friction_coeff = 1/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
# setup initial coordinates
init_coord = app.PDBFile('start.pdb').getPositions()
# setup simulation
condensate.set_simulation(integrator, platform_name, init_coord=init_coord)
# perform energy minimization
condensate.simulation.minimizeEnergy()
output_interval = dt
output_dcd = 'npt.dcd'
# Run short simulation
condensate.add_reporters(output_interval, output_dcd)
condensate.simulation.reporters.append(app.checkpointreporter.CheckpointReporter('npt.cpt',output_interval))
condensate.simulation.reporters.append(app.StateDataReporter('npt.csv', output_interval, step=True, potentialEnergy=True, temperature=True, volume=True, density=True, speed=True, time=True))
condensate.simulation.context.setVelocitiesToTemperature(temperature)
condensate.simulation.step(tot_t)

# print final box vectors
state = condensate.simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
with open('npt.xml', 'w') as f:
    f.write(mm.XmlSerializer.serialize(state))
# save final system
condensate.save_system('Nup58_CT_system.xml')

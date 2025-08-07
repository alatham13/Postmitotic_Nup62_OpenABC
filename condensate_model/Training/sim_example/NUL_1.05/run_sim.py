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

# set simulation platform
platform_name = 'CUDA'

print('Successful import')

# start from predicted PDB
Nup_parser = MOFFParser.from_atomistic_pdb('NUL_AF.pdb', 'NUL_CA.pdb')

# ProTa has no ordered domains- pass empty native pairs
old_native_pairs = Nup_parser.native_pairs.copy()
new_native_pairs = pd.DataFrame(columns=old_native_pairs.columns)
Nup_parser.native_pairs = new_native_pairs
Nup_parser.parse_exclusions() # update exclusions based on the new native pairs

# loop to remove large angles
max_angle=130*np.pi/180
# copy the angle potential
old_angles=Nup_parser.protein_angles.copy()
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
print('These angles have been rounded down to '+str(count_bad_angle)+' degrees to ensure numeric stability.')
# replace angles with new angles
Nup_parser.protein_angles=new_angles

# check interaction parameters
print(Nup_parser.native_pairs)


# setup initial topology
#box_a, box_b, box_c = 200, 200, 200
# set up protein model
protein = MOFFMRGModel()
# add protein to model
protein.append_mol(Nup_parser)
# setup initial topology
top = app.PDBFile('NUL_CA.pdb').getTopology()
protein.create_system(top, remove_cmmotion=True)
# setup variables necessary for simulations
salt_concentration = 162*unit.millimolar
temperature = 300*unit.kelvin
# add force terms to simulation
protein.add_protein_bonds(force_group=1)
protein.add_protein_angles(force_group=2)
protein.add_protein_dihedrals(force_group=3)
protein.add_native_pairs(force_group=4)
protein.add_contacts(force_group=5,scale=1.05)
protein.add_elec_switch(salt_concentration, temperature, force_group=6)
# save initial system
protein.save_system('NUL_system.xml')
# setup integrator
friction_coeff = 1/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
# setup initial coordinates
init_coord = app.PDBFile('NUL_CA.pdb').getPositions()
# setup simulation
protein.set_simulation(integrator, platform_name, init_coord=init_coord)
# perform energy minimization
protein.simulation.minimizeEnergy()
output_interval = 10000
output_dcd = 'output.dcd'
# Run short simulation
protein.add_reporters(output_interval, output_dcd)
protein.simulation.reporters.append(app.checkpointreporter.CheckpointReporter('checkpoint.cpt',output_interval))
protein.simulation.reporters.append(app.StateDataReporter('data.csv', output_interval, step=True, potentialEnergy=True, temperature=True, volume=True, density=True, speed=True, time=True))
protein.simulation.context.setVelocitiesToTemperature(temperature)
protein.simulation.step(200000000)


# run simulation of a protein.
import numpy as np
import pandas as pd
import os
import sys
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

sys.path.insert(0, "~/Postmitotic_Nup62_OpenABC")
# from openabc2.forcefields.parsers import MOFFParser, MRGdsDNAParser
from openabc2.forcefields import MOFFMRGModel
from openabc2.utils.insert import insert_molecules

# start revised MOFFParser. This is necessary due to a bug that arises in the original code when no native pairs are found in the PDB -------
import mdtraj
from openabc2.utils import helper_functions
from openabc2.utils.shadow_map import find_ca_pairs_from_atomistic_pdb

_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                'SER', 'THR', 'TRP', 'TYR', 'VAL']

_moff_amino_acid_mass_dict = dict(ALA=71.08, ARG=156.20, ASN=114.10, ASP=115.10, CYS=103.10,
                                  GLN=128.10, GLU=129.10, GLY=57.05, HIS=137.10, ILE=113.20,
                                  LEU=113.20, LYS=128.20, MET=131.20, PHE=147.20, PRO=97.12,
                                  SER=87.08, THR=101.10, TRP=186.20, TYR=163.20, VAL=99.07)

_moff_amino_acid_charge_dict = dict(ALA=0.0, ARG=1.0, ASN=0.0, ASP=-1.0, CYS=0.0,
                                    GLN=0.0, GLU=-1.0, GLY=0.0, HIS=0.25, ILE=0.0,
                                    LEU=0.0, LYS=1.0, MET=0.0, PHE=0.0, PRO=0.0,
                                    SER=0.0, THR=0.0, TRP=0.0, TYR=0.0, VAL=0.0)

_kcal_to_kj = 4.184


class MOFFParser(object):
    """
    MOFF protein parser.
    """

    def __init__(self, atomistic_pdb, ca_pdb, default_parse=True):
        """
        Initialize a protein with MOFF model.

        Parameters
        ----------
        atomistic_pdb : str
            Path for the atomistic pdb file.

        ca_pdb : str
            path for the CA pdb file.

        default_parse : bool
            Whether to parse with default settings.

        """
        self.atomistic_pdb = atomistic_pdb
        self.pdb = ca_pdb
        self.atoms = helper_functions.parse_pdb(ca_pdb)
        # check if all the atoms are protein CA atoms
        assert ((self.atoms['resname'].isin(_amino_acids)).all() and self.atoms['name'].eq('CA').all())
        if default_parse:
            print('Parse molecule with default settings.')
            self.parse_mol()

    @classmethod
    def from_atomistic_pdb(cls, atomistic_pdb, ca_pdb, write_TER=False, default_parse=True):
        """
        Initialize an HPS model protein from atomistic pdb.

        Parameters
        ----------
        atomistic_pdb : str
            Path for the atomistic pdb file.

        ca_pdb : str
            path for the CA pdb file.

        write_TER : bool
            Whether to write TER between two chains.

        default_parse : bool
            Whether to parse with default settings.

        """
        helper_functions.atomistic_pdb_to_ca_pdb(atomistic_pdb, ca_pdb, write_TER)
        return cls(atomistic_pdb, ca_pdb, default_parse)

    def parse_exclusions(self, exclude12=True, exclude13=True, exclude14=True, exclude_native_pairs=True):
        """
        Parse nonbonded exclusions based on bonds, angles, dihedrals, and native pairs.

        Parameters
        ----------
        exclude12 : bool
            Whether to exclude nonbonded interactions between 1-2 atom pairs.

        exclude13 : bool
            Whether to exclude nonbonded interactions between 1-3 atom pairs.

        exclude14 : bool
            Whether to exclude nonbonded interactions between 1-4 atom pairs.

        exclude_native_pairs : bool
            Whether to exclude nonbonded interactions between native pairs.

        """
        exclusions = []
        if exclude12 and hasattr(self, 'protein_bonds'):
            for i, row in self.protein_bonds.iterrows():
                exclusions.append((int(row['a1']), int(row['a2'])))
        if exclude13 and hasattr(self, 'protein_angles'):
            for i, row in self.protein_angles.iterrows():
                exclusions.append((int(row['a1']), int(row['a3'])))
        if exclude14 and hasattr(self, 'protein_dihedrals'):
            for i, row in self.protein_dihedrals.iterrows():
                exclusions.append((int(row['a1']), int(row['a4'])))
        if exclude_native_pairs and hasattr(self, 'native_pairs'):
            for i, row in self.native_pairs.iterrows():
                exclusions.append((int(row['a1']), int(row['a2'])))
        exclusions = np.array(sorted(exclusions))
        self.exclusions = pd.DataFrame(exclusions, columns=['a1', 'a2']).drop_duplicates(ignore_index=True)

    def parse_mol(self, get_native_pairs=False, frame=0, radius=0.1, bonded_radius=0.05, cutoff=0.6, box=None,
                  use_pbc=False, exclude12=True, exclude13=True, exclude14=True, exclude_native_pairs=True,
                  mass_dict=_moff_amino_acid_mass_dict, charge_dict=_moff_amino_acid_charge_dict):
        """
        Parse molecule.

        Parameters
        ----------
        get_native_pairs : bool
            Whether to get native pairs from atomistic pdb with shadow algorithm.

        frame : int
            Frame index to compute native structure parameters for both CG and atomistic models.

        radius : float or int
            Shadow algorithm radius.

        bonded_radius : float or int
            Shadow algorithm bonded radius.

        cutoff : float or int
            Shadow algorithm cutoff.

        box : None or np.ndarray
        Specify box shape.
        Note `use_pbc` = False is only compatible with orthogonal box.
        If `use_pbc` = False, then set `box` as None.
        If `use_pbc` = True, then `box` = np.array([lx, ly, lz, alpha1, alpha2, alpha3]).

        use_pbc : bool
            Whether to use periodic boundary condition (PBC) for searching native pairs.

        exclude12 : bool
            Whether to exclude nonbonded interactions between 1-2 atom pairs.

        exclude13 : bool
            Whether to exclude nonbonded interactions between 1-3 atom pairs.

        exclude14 : bool
            Whether to exclude nonbonded interactions between 1-4 atom pairs.

        exclude_native_pairs : bool
            Whether to exclude nonbonded interactions between native pairs.

        mass_dict : dict
            Mass dictionary.

        charge_dict : dict
            Charge dictionary.

        """
        # set bonds, angles, and dihedrals
        bonds, angles, dihedrals = [], [], []
        n_atoms = len(self.atoms.index)
        for atom1 in range(n_atoms):
            chain1 = self.atoms.loc[atom1, 'chainID']
            if atom1 < n_atoms - 1:
                atom2 = atom1 + 1
                chain2 = self.atoms.loc[atom2, 'chainID']
                if chain1 == chain2:
                    bonds.append([atom1, atom2])
            if atom1 < n_atoms - 2:
                atom3 = atom1 + 2
                chain3 = self.atoms.loc[atom3, 'chainID']
                if (chain1 == chain2) and (chain1 == chain3):
                    angles.append([atom1, atom2, atom3])
            if atom1 < n_atoms - 3:
                atom4 = atom1 + 3
                chain4 = self.atoms.loc[atom4, 'chainID']
                if (chain1 == chain2) and (chain1 == chain3) and (chain1 == chain4):
                    dihedrals.append([atom1, atom2, atom3, atom4])
        bonds, angles, dihedrals = np.array(bonds), np.array(angles), np.array(dihedrals)
        traj = mdtraj.load_pdb(self.pdb)
        self.protein_bonds = pd.DataFrame(bonds, columns=['a1', 'a2'])
        self.protein_bonds.loc[:, 'r0'] = 0.38
        self.protein_bonds.loc[:, 'k_bond'] = 1000
        self.protein_angles = pd.DataFrame(angles, columns=['a1', 'a2', 'a3'])
        self.protein_angles['theta0'] = mdtraj.compute_angles(traj, angles, periodic=use_pbc)[frame]
        self.protein_angles.loc[:, 'k_angle'] = 120
        self.protein_dihedrals = pd.DataFrame(columns=['a1', 'a2', 'a3', 'a4', 'periodicity', 'phi0', 'k_dihedral'])
        phi = mdtraj.compute_dihedrals(traj, dihedrals, periodic=use_pbc)[frame]
        for i in range(dihedrals.shape[0]):
            row = dihedrals[i].tolist() + [1, phi[i] + np.pi, 3.0]
            self.protein_dihedrals.loc[len(self.protein_dihedrals.index)] = row
            row = dihedrals[i].tolist() + [3, 3 * (phi[i] + np.pi), 1.5]
            self.protein_dihedrals.loc[len(self.protein_dihedrals.index)] = row
        # set native pairs
        if get_native_pairs:
            print('Get native pairs with shadow algorithm.')
            self.native_pairs = find_ca_pairs_from_atomistic_pdb(self.atomistic_pdb, frame, radius, bonded_radius,
                                                                 cutoff, box, use_pbc)
            if len(self.native_pairs) > 0:
                self.native_pairs.loc[:, 'epsilon'] = 3.0
        else:
            self.native_pairs = pd.DataFrame(columns=['a1', 'a2', 'mu', 'epsilon'])
        # set exclusions
        self.parse_exclusions(exclude12, exclude13, exclude14, exclude_native_pairs)
        # set mass and charge
        for i, row in self.atoms.iterrows():
            self.atoms.loc[i, 'mass'] = mass_dict[row['resname']]
            self.atoms.loc[i, 'charge'] = charge_dict[row['resname']]

# Done -----------------------------------------------------------------------------------------------------------------

tot_t=10000000
dt=int(tot_t/100)

# set simulation platform
platform_name = 'CUDA'

print('Successful import')

# set the maximum allowed angle
max_angle=130*np.pi/180

# start from predicted PDB
FG_parser = MOFFParser.from_atomistic_pdb('Nup58_NT_AF.pdb', 'Nup58_NT_CA.pdb')

# Nup58_NT has no native pairs
old_native_pairs = FG_parser.native_pairs.copy()
new_native_pairs = pd.DataFrame(columns=old_native_pairs.columns)
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
insert_molecules('Nup58_NT_CA.pdb', 'start.pdb', n_mol, box=[box_a, box_b, box_c])


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
condensate.save_system('Nup58_NT_system.xml')

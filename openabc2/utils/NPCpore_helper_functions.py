import numpy as np
import pandas as pd

try:
    import openmm.unit as unit
except ImportError:
    import simtk.unit as unit
import sys
import os
from openabc2.lib import (_amino_acids, _amino_acid_1_letter_to_3_letters_dict,
                         _dna_nucleotides, _dna_WC_pair_dict,
                         _rna_nucleotides)

"""
Some code is adapted from Open3SPN2. 

Open3SPN2 and OpenAWSEM paper: 
Lu, Wei, et al. "OpenAWSEM with Open3SPN2: A fast, flexible, and accessible framework for large-scale coarse-grained biomolecular simulations." PLoS computational biology 17.2 (2021): e1008308.
"""

_nucleotides = _dna_nucleotides + _rna_nucleotides


def parse_pdb_long(pdb_file):
    """
    Load pdb file as pandas dataframe.

    Parameters
    ----------
    pdb_file : str
        Path for the pdb file.

    Returns
    -------
    pdb_atoms : pd.DataFrame
        A pandas dataframe includes atom information.

    """

    def pdb_line(line):
        return dict(recname=str(line[0:6]).strip(),
                    serial=str(line[6:11]),
                    name=str(line[12:16]).strip(),
                    altLoc=str(line[16:17]),
                    resname=str(line[17:20]).strip(),
                    chainID=str(line[20:22]),
                    resSeq=int(line[22:26]),
                    iCode=str(line[26:27]),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=0.0 if line[54:60].strip() == '' else float(line[54:60]),
                    tempFactor=0.0 if line[60:66].strip() == '' else float(line[60:66]),
                    element=str(line[76:78].strip()),
                    charge=str(line[78:80].strip()))

    with open(pdb_file, 'r') as pdb:
        lines = []
        for line in pdb:
            if (len(line) > 6) and (line[:6] in ['ATOM  ', 'HETATM']):
                lines += [pdb_line(line)]
    pdb_atoms = pd.DataFrame(lines)
    pdb_atoms = pdb_atoms[['recname', 'serial', 'name', 'altLoc',
                           'resname', 'chainID', 'resSeq', 'iCode',
                           'x', 'y', 'z', 'occupancy', 'tempFactor',
                           'element', 'charge']]
    return pdb_atoms

def write_pdb_long(pdb_atoms, pdb_file, write_TER=False):
    """
    Write pandas dataframe to pdb file.

    Parameters
    ----------
    pdb_atoms : pd.DataFrame
        A pandas dataframe includes atom information.

    pdb_file : str
        Output path for the pdb file.

    write_TER : bool
        Whether to write TER between two chains.

    """
    chainID = None
    with open(pdb_file, 'w') as pdb:
        for i, atom in pdb_atoms.iterrows():
            if chainID is not None:
                if write_TER and (atom['chainID'] != chainID):
                    pdb.write('TER\n')
            chainID = atom['chainID']
            #print(atom.recname)
            #print(atom.serial)
            pdb_line = f'{atom.recname:<6}{int(atom.serial):>6}{atom["name"]:^4}{atom.altLoc:1}' + \
                       f'{atom.resname:<3} {atom.chainID[1]:1}{int(atom.resSeq):>4}{atom.iCode:1}   ' + \
                       f'{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}' + \
                       f'{atom.occupancy:>6.2f}{atom.tempFactor:>6.2f}' + f' {atom.chainID:2}' + ' ' * 7 + \
                       f'{atom.element:>2}{atom.charge:>2}'
            assert len(
                pdb_line) == 80, f'An item in the atom table is longer than expected ({len(pdb_line)})\n{pdb_line}'
            pdb.write(pdb_line + '\n')
        pdb.write('END\n')


def atomistic_pdb_to_ca_pdb(atomistic_pdb, ca_pdb, write_TER=False):
    """
    Convert atomistic pdb to protein CA pdb.

    Parameters
    ----------
    atomistic_pdb : str
        Path for the atomistic pdb file.

    ca_pdb : str
        Output path for the CA pdb file.

    write_TER : bool
        Whether to write TER between two chains.

    """
    atomistic_pdb_atoms = parse_pdb_long(atomistic_pdb)
    ca_pdb_atoms = pd.DataFrame(columns=atomistic_pdb_atoms.columns)
    for i, row in atomistic_pdb_atoms.iterrows():
        if (row['resname'] in _amino_acids) and (row['name'] == 'CA'):
            ca_pdb_atoms.loc[len(ca_pdb_atoms.index)] = row
    print('Read '+str(len(ca_pdb_atoms))+' atoms.')
    ca_pdb_atoms['serial'] = list(range(1, len(ca_pdb_atoms.index) + 1))
    ca_pdb_atoms.loc[:, 'charge'] = ''  # remove charge
    write_pdb_long(ca_pdb_atoms, ca_pdb, write_TER)



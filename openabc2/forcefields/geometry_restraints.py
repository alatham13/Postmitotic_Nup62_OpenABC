import numpy as np
import pandas as pd
import os
import sys
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

def spherical_boundary_term(CGModel,k,x0,y0,z0):
    membrane = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    membrane.addGlobalParameter('x0', x0*unit.angstrom)
    membrane.addGlobalParameter('y0', y0*unit.angstrom)
    membrane.addGlobalParameter('z0', z0*unit.angstrom)
    membrane.addGlobalParameter('k', k*(unit.kilojoules_per_mole/(unit.angstrom*unit.angstrom)))
    CGModel.system.addForce(membrane)
    for i in range(CGModel.system.getNumParticles()):
        membrane.addParticle(i, [])

# Use step for if statements
def flat_membrane(CGModel,k,h):
    membrane = mm.CustomExternalForce("k * ( step(z) * step( (h/2)-z) * (z-(h/2))^2) + k * ( step(-z) * step((h/2)+z) * (z+(h/2))^2)")
    membrane.addGlobalParameter('h', h*unit.angstrom)
    membrane.addGlobalParameter('k', k*(unit.kilojoules_per_mole/(unit.angstrom*unit.angstrom)))
    CGModel.system.addForce(membrane)
    for i in range(CGModel.system.getNumParticles()):
        membrane.addParticle(i, [])
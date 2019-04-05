#!/usr/bin/env python
from collections import defaultdict

import sys

# OpenMM Imports
import simtk.openmm as mm
from simtk.openmm import KcalPerKJ
import simtk.openmm.app as app 
import parmed as pmd
# ParmEd Imports
from parmed import load_file, unit as u
from parmed.charmm import CharmmParameterSet
from parmed.openmm import StateDataReporter, energy_decomposition_system
from sys import setrecursionlimit

setrecursionlimit(1500)

pdb = app.PDBFile('QUBE_pro.pdb')
forcefield = app.ForceField('QUBE_pro.xml')


system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,nonbondedCutoff=500.0*u.angstroms, switchDistance=496.0*u.angstroms)

#serializing the system
from simtk.openmm import XmlSerializer
serialized_system = XmlSerializer.serialize(system)
outfile = open('2.xml','w')
outfile.write(serialized_system)
outfile.close()

app.PDBFile.writeFile(pdb.topology, pdb.positions, open('output.pdb', 'w+'))
# Create the integrator to do Langevin dynamics
integrator = mm.LangevinIntegrator(
                        300*u.kelvin,       # Temperature of heat bath
                        1.0/u.picoseconds,  # Friction coefficient
                        2.0*u.femtoseconds, # Time step
)





platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
print ('energy from openmm library')
print (simulation.context.getState(getEnergy=True).getPotentialEnergy()) 
#platform.setPropertyValue(simulation.context, property='Precision', value='double')
# Minimize the energy
#print('Minimizing energy')

#Minimize(simulation,iters=0)

struct=pmd.load_file('QUBE_pro.pdb')
ecomps=(pmd.openmm.energy_decomposition_system(struct,system))
tot_ene=0.0
for i in range(0,len(ecomps)):
        tot_ene+=ecomps[i][1]
        print(ecomps[i][0],ecomps[i][1])
print('Total-energy %6.6f'%tot_ene)

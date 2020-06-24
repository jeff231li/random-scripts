import os 
import pandas as pd
import matplotlib.pyplot as plt
import simtk.unit as unit
from MDAnalysis.coordinates.DCD import DCDFile
from simtk.openmm.app import PDBFile
from simtk.openmm import XmlSerializer
from simtk.openmm.app import AmberPrmtopFile, Topology
from simtk import openmm

# Load system.xml
with open("system.xml", "r") as f:
    system = XmlSerializer.deserialize(f.read())
prmtop = AmberPrmtopFile('MGLab8-MCH.prmtop')

dt   = 2.0 * unit.femtoseconds
Temp = 298.15 * unit.kelvin
Pres = 1.01325 * unit.bar
integrator = openmm.LangevinIntegrator(Temp, 1.0 / unit.picoseconds, dt)
barostat   = openmm.MonteCarloBarostat(Pres, Temp, 100)
system.addForce(barostat)

# Create Simulation Object
simulation = openmm.app.Simulation(
    prmtop.topology,
    system,
    integrator,
    openmm.Platform.getPlatformByName('CUDA')
)

# Extract Collective variables from trajectory
colvars_list = []
with DCDFile("MD-Equil.dcd") as dcd:
    header = dcd.header
    for frame in dcd.readframes()[0]:
        simulation.context.setPositions(frame * unit.angstrom)
        colvars_list.append(list(system.getForce(28).getCollectiveVariableValues(simulation.context)))

# Write Collective variables to file
with open('COLVAR.txt', 'w') as file:
    for colvars in colvars_list:
        file.writelines(' '.join(f"{value*10:10.5f}" for value in colvars))
        file.writelines('\n')


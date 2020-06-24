#!/usr/bin/env python
import numpy as np
import parmed as pmd
from openmmtools.alchemy import *
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from time import time

# Initialize
temperature = 298.15
pressure    = 1.01325
timestep    = 2.0
equil_steps = 500000

platform = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0',
              'CudaPrecision'  : 'mixed'}

# Load System
prmtop = AmberPrmtopFile('MGLab8-MCH.prmtop')
inpcrd = AmberInpcrdFile('MGLab8-MCH.rst7')
with open("system.xml", "r") as file:
    system = XmlSerializer.deserialize(file.read())

# Thermostat and Barostat
integrator = LangevinIntegrator(temperature*kelvin, 1.0/picosecond, timestep*femtoseconds)
barostat = MonteCarloBarostat(pressure*bar, temperature*kelvin)
system.addForce(barostat)

# Simulation Object
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)

# Reporters
simulation.reporters.append(
    StateDataReporter(
        f'openmm-statistics.csv',
        5000,
        step=True,
        kineticEnergy=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
        separator=",",
    )
)
simulation.reporters.append(DCDReporter(f'MD-Equil.dcd', 5000))

# Minimize
print("Minimizing Energy")
simulation.minimizeEnergy(tolerance=0.01*kilocalories_per_mole, maxIterations=5000)
simulation.context.setVelocitiesToTemperature(temperature*kelvin)

# Equilibration
print("Equilibration")
simulation.step(equil_steps)


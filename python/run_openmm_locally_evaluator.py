#!/usr/bin/env python
# coding: utf-8

import os

import nglview as nv
import pytraj as pt

from progressbar import ProgressBar, Bar, Percentage, Timer, ETA, progressbar

import simtk.unit as unit
from simtk.openmm.app import PDBFile, DCDReporter, StateDataReporter, Simulation
from simtk.openmm import XmlSerializer, LangevinIntegrator, MonteCarloBarostat, Platform

# Initialization
with open("restrained.xml", "r") as f:
    system = XmlSerializer.deserialize(f.read())
coords = PDBFile(os.path.join('simulations', 'npt_production', 'input.pdb'))

dt   = 4.0 * unit.femtoseconds
Temp = 298.15 * unit.kelvin
Pres = 1.01325 * unit.bar
out_freq = 2500
print_freq = 50000
MD_steps = 12500000
bar_freq = int(MD_steps/print_freq)

integrator = LangevinIntegrator(Temp, 1.0 / unit.picoseconds, dt)
barostat   = MonteCarloBarostat(Pres, Temp, 100)
system.addForce(barostat)

# Simulation Object
simulation = Simulation(
    coords.topology, 
    system, 
    integrator,
    Platform.getPlatformByName('CUDA'),
    {'CudaDeviceIndex': '0', 'CudaPrecision': 'single'}
)
simulation.context.setPositions(coords.positions)

simulation.reporters.append(DCDReporter(os.path.join('simulations', 'npt_production', 'test.dcd'), out_freq, False))
simulation.reporters.append(StateDataReporter(
    os.path.join('simulations', 'npt_production', 'test.log'),
    out_freq,
    step=True,
    kineticEnergy=True,
    potentialEnergy=True,
    totalEnergy=True,
    temperature=True,
    totalSteps=MD_steps,
    progress=True,
    remainingTime=True,
    speed=True,
    separator=",",
))

# Run OpenMM
widgets = [' [', Timer(), '] ', Bar(), ' (', ETA(), ')']
bar = ProgressBar(max=bar_freq, widgets=widgets).start()
restr_potential = []

try:
    for i in progressbar(range(bar_freq), widgets=widgets):
        simulation.step(print_freq)

        static_restraint_energy = simulation.context.getState(getEnergy=True, groups={10})
        static_restraint_energy = static_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

        conformational_restraint_energy = simulation.context.getState(getEnergy=True, groups={11})
        conformational_restraint_energy = conformational_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

        guest_restraint_energy = simulation.context.getState(getEnergy=True, groups={12})
        guest_restraint_energy = guest_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

        symmetry_restraint_energy = simulation.context.getState(getEnergy=True, groups={13})
        symmetry_restraint_energy = symmetry_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

        wall_restraint_energy = simulation.context.getState(getEnergy=True, groups={14})
        wall_restraint_energy = wall_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

        restr_potential.append([
            static_restraint_energy,
            conformational_restraint_energy,
            guest_restraint_energy,
            symmetry_restraint_energy,
            wall_restraint_energy,
        ])

except:
    print('\nExtracting coordinates')
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(os.path.join('simulations', 'npt_production', 'test.pdb'), 'w'))
    
    print('Extracting info on restraints')
    static_restraint_energy = simulation.context.getState(getEnergy=True, groups={10})
    static_restraint_energy = static_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

    conformational_restraint_energy = simulation.context.getState(getEnergy=True, groups={11})
    conformational_restraint_energy = conformational_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

    guest_restraint_energy = simulation.context.getState(getEnergy=True, groups={12})
    guest_restraint_energy = guest_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

    symmetry_restraint_energy = simulation.context.getState(getEnergy=True, groups={13})
    symmetry_restraint_energy = symmetry_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

    wall_restraint_energy = simulation.context.getState(getEnergy=True, groups={14})
    wall_restraint_energy = wall_restraint_energy.getPotentialEnergy() / unit.kilocalorie_per_mole

    restr_potential.append([
        static_restraint_energy,
        conformational_restraint_energy,
        guest_restraint_energy,
        symmetry_restraint_energy,
        wall_restraint_energy,
    ])
    print(f"{restr_potential[-1]}")

# Write out restraint potential
f = open(os.path.join('simulations', 'npt_production', 'restraint_potential.csv'), "w")
f.writelines(f'"frame","static","conformation","guest","symmetry","wall"\n')
i = 0
for restr in restr_potential:
    f.writelines(f"{i},{restr[0]:.5f},{restr[1]:.5f},{restr[2]:.5f},{restr[3]:.5f},{restr[4]:.5f}\n")
    i += 1
f.close()


#!/usr/bin/env python
# coding: utf-8

import logging
import time
from importlib import reload

import numpy as np
import simtk.openmm as openmm
import simtk.openmm.app as app
import simtk.unit as openmm_unit
from MDAnalysis import Universe
from mpi4py import MPI
from scipy.integrate import cumtrapz
from tqdm import tqdm

reload(logging)

logger = logging.getLogger()
logging.basicConfig(
    filename="analysis-mpi.log",
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.INFO,
)

# MPI setting
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def add_wca_repulsive(
    system,
    particle_selection,
    environment_selection,
    particle_radius=4.0 * openmm_unit.angstrom,
    wca_sigma=4.0 * openmm_unit.angstrom,
    wca_epsilon=1.0 * openmm_unit.kilocalorie_per_mole,
    coeff_r=50,
    coeff_a=49,
    force_group=10,
):
    nonbonded = [
        force
        for force in system.getForces()
        if isinstance(force, openmm.NonbondedForce)
    ][0]

    wca_repulsive = openmm.CustomNonbondedForce(
        "U_repulsive;"
        "U_repulsive = step(R_particle - r) * (U_Mie + epsilon_wall);"
        "U_Mie = prefactor * epsilon_wall * ((sigma_wall/r_prime)^coeff_r - (sigma_wall/r_prime)^coeff_a);"
        "prefactor = coeff_r/(coeff_r-coeff_a)*(coeff_r/coeff_a)^(coeff_r/(coeff_r-coeff_a));"
        "r_prime = r - (R_particle - R_min);"
        "R_min = sigma_wall * (coeff_r/coeff_a)^(1/(coeff_r-coeff_a));"
    )
    wca_repulsive.addGlobalParameter("R_particle", particle_radius)
    wca_repulsive.addGlobalParameter("coeff_r", coeff_r)
    wca_repulsive.addGlobalParameter("coeff_a", coeff_a)
    wca_repulsive.addGlobalParameter("sigma_wall", wca_sigma)
    wca_repulsive.addGlobalParameter("epsilon_wall", wca_epsilon)
    wca_repulsive.addPerParticleParameter("sigma")
    wca_repulsive.addPerParticleParameter("epsilon")
    wca_repulsive.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    wca_repulsive.setCutoffDistance(nonbonded.getCutoffDistance())
    wca_repulsive.setUseLongRangeCorrection(False)
    wca_repulsive.setForceGroup(force_group)

    wca_repulsive.addInteractionGroup(particle_selection, environment_selection)

    # Set LJ parameters
    for atom in range(nonbonded.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(atom)
        if atom in particle_selection:
            wca_repulsive.addParticle([wca_sigma, wca_epsilon])
        else:
            wca_repulsive.addParticle([sigma, epsilon])

    # Set LJ parameters to zero
    for atom in particle_selection:
        [charge, sigma, epsilon] = nonbonded.getParticleParameters(atom)
        nonbonded.setParticleParameters(atom, 0.0, wca_sigma, 0.0)

    # Transfer Exclusion
    for exception_index in range(nonbonded.getNumExceptions()):
        iatom, jatom, chargeprod, sigma, epsilon = nonbonded.getExceptionParameters(
            exception_index
        )
        wca_repulsive.addExclusion(iatom, jatom)

    # Add Force to System
    system.addForce(wca_repulsive)


if rank == 0:
    logger.info("Initializing calculations - loading trajectories etc...")

particle_radius = 4.5 * openmm_unit.angstrom
wca_sigma = 3.0 * openmm_unit.angstrom
wca_epsilon = 0.1 * openmm_unit.kilocalorie_per_mole
host_resname = "CB7"
guest_resname = "DM1"
solvent_resname = "HOH"
base_name = "restrained"

# Create Universe
window_list = [f"p{i:03}" for i in range(41)]
universe = Universe(
    f"p000/{base_name}.pdb",
    [f"{window}/production.dcd" for window in window_list],
)
all_atoms = universe.select_atoms("all")
guest = universe.select_atoms(f"resname {guest_resname}")
host = universe.select_atoms(f"resname {host_resname} and element C")

# Load in PDB and XML file
pdbfile = app.PDBFile(f"p000/{base_name}.pdb")
with open(f"p000/{base_name}.xml", "r") as f:
    system = openmm.XmlSerializer.deserialize(f.read())

# Atom indices
dummy_particle = [
    atom.index
    for atom in pdbfile.topology.atoms()
    if atom.residue.name == guest_resname
]
host_indices = [
    atom.index for atom in pdbfile.topology.atoms() if atom.residue.name == host_resname
]
solvent_indices = [
    atom.index
    for atom in pdbfile.topology.atoms()
    if atom.residue.name == solvent_resname
]

# Add in WCA repulsive core
add_wca_repulsive(
    system,
    dummy_particle,
    host_indices,
    particle_radius=particle_radius,
    wca_sigma=wca_sigma,
    wca_epsilon=wca_epsilon,
    coeff_r=50,
    coeff_a=49,
    force_group=4,
)
add_wca_repulsive(
    system,
    dummy_particle,
    solvent_indices,
    particle_radius=particle_radius,
    wca_sigma=wca_sigma,
    wca_epsilon=wca_epsilon,
    coeff_r=50,
    coeff_a=49,
    force_group=7,
)

# Create Integrator and Simulation object again
thermostat = openmm.LangevinIntegrator(
    298.15 * openmm_unit.kelvin,
    1.0 / openmm_unit.picosecond,
    2.0 * openmm_unit.femtosecond,
)
simulation = app.Simulation(
    pdbfile.topology,
    system,
    thermostat,
    openmm.Platform.getPlatformByName("CPU"),
)

# Specify array sizes for MPI nodes
n_frames = universe.trajectory.n_frames
num_per_rank = int(n_frames / size)
lower_bound = rank * num_per_rank
upper_bound = (rank + 1) * num_per_rank
report_freq = int(num_per_rank / 100)

# Chunk arrays
z_axis = np.array([0, 0, 1])
r = np.linspace(0, 15, 61)

chunk_n_bins = np.zeros(len(r))
chunk_mean_force_host = np.zeros(len(r))
chunk_mean_force_all_solvent = np.zeros(len(r))

if rank == 0:
    n_bins = np.zeros(len(r))
    mean_force_host = np.zeros(len(r))
    mean_force_all_solvent = np.zeros(len(r))
else:
    n_bins = None
    mean_force_host = None
    mean_force_all_solvent = None


comm.Barrier()


if rank == 0:
    start = time.time()
    logger.info("Analyzing trajectories...")

for i, frame in enumerate(tqdm(universe.trajectory[lower_bound:upper_bound])):
    # Report Completion
    if rank == 0 and i % report_freq == 0:
        logger.info(f"Completed: {i/num_per_rank*100:5.2f}%")

    # Set Box vectors and Coordinates
    box = universe.dimensions
    box_vectors = openmm_unit.Quantity(
        value=[
            openmm.Vec3(box[0], 0, 0),
            openmm.Vec3(0, box[1], 0),
            openmm.Vec3(0, 0, box[2]),
        ],
        unit=openmm_unit.angstrom,
    )
    simulation.context.setPositions(all_atoms.positions * openmm_unit.angstrom)

    # Find location in bins
    vector = np.mean(guest.positions, axis=0) - np.mean(host.positions, axis=0)
    i_loc = np.where(
        np.min(abs(np.dot(vector, z_axis) - r)) == abs(np.dot(vector, z_axis) - r)
    )[0][0]
    chunk_n_bins[i_loc] += 1

    # Calculate Forces
    force_by_host = simulation.context.getState(getForces=True, groups={4}).getForces()
    force_by_all_solvent = simulation.context.getState(
        getForces=True, groups={7}
    ).getForces()

    # Vector projection of forces to reaction coordinate
    chunk_mean_force_host[i_loc] += np.dot(
        np.array(
            force_by_host[dummy_particle[0]].value_in_unit(
                openmm_unit.kilocalorie_per_mole / openmm_unit.angstrom
            )
        ),
        z_axis,
    )
    chunk_mean_force_all_solvent[i_loc] += np.dot(
        np.array(
            force_by_all_solvent[dummy_particle[0]].value_in_unit(
                openmm_unit.kilocalorie_per_mole / openmm_unit.angstrom
            )
        ),
        z_axis,
    )

comm.Barrier()


# Collect data - sum arrays from all nodes
comm.Reduce([chunk_n_bins, MPI.DOUBLE], [n_bins, MPI.DOUBLE], op=MPI.SUM, root=0)
comm.Reduce(
    [chunk_mean_force_host, MPI.DOUBLE],
    [mean_force_host, MPI.DOUBLE],
    op=MPI.SUM,
    root=0,
)
comm.Reduce(
    [chunk_mean_force_all_solvent, MPI.DOUBLE],
    [mean_force_all_solvent, MPI.DOUBLE],
    op=MPI.SUM,
    root=0,
)

comm.Barrier()


# Post-process and print results
if rank == 0:
    logger.info("Post-processing...")
    end = time.time()
    logger.info(f"Total time taken: {(end-start)} seconds")

    n_bins[n_bins == 0.0] = 1.0

    # Integrate Forces
    pmf_host = cumtrapz(-1 * mean_force_host / n_bins, x=r, initial=0.0)
    pmf_all_solvent = cumtrapz(-1 * mean_force_all_solvent / n_bins, x=r, initial=0.0)

    np.savetxt(
        "mean_force-split-mpi.txt",
        np.c_[
            r,
            pmf_host,
            pmf_all_solvent,
            pmf_host + pmf_all_solvent,
        ],
        fmt="%10.5f",
    )

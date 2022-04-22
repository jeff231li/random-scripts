import json
import logging
import sys
from copy import deepcopy
from importlib import reload

import numpy as np
import simtk.unit as unit
from MDAnalysis import Universe
from mpi4py import MPI
from paprika.io import NumpyEncoder
from simtk.openmm import (
    HarmonicAngleForce,
    HarmonicBondForce,
    NonbondedForce,
    PeriodicTorsionForce,
    Platform,
    Vec3,
    VerletIntegrator,
    XmlSerializer,
)
from simtk.openmm.app import PDBFile, Simulation
from tqdm import tqdm

reload(logging)

logger = logging.getLogger()
logging.basicConfig(
    filename="analysis.log",
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.INFO,
)

# MPI Settings
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    logger.info("Starting full analysis of Lennard-Jones interactions.")


def create_simulation(
    system,
    coords,
    dt=2.0 * unit.femtosecond,
):
    integrator = VerletIntegrator(dt)

    simulation = Simulation(
        coords.topology,
        system,
        integrator,
        Platform.getPlatformByName("CPU"),
        {"Threads": "1"},
    )

    return simulation


def get_energy(simulation, xyz, box, force_group=0):
    box_vectors = unit.Quantity(
        value=[Vec3(box[0], 0, 0), Vec3(0, box[1], 0), Vec3(0, 0, box[2])],
        unit=unit.angstrom,
    )

    simulation.context.setPositions(xyz * unit.angstrom)
    simulation.context.setPeriodicBoxVectors(*box_vectors)

    energy = simulation.context.getState(getEnergy=True, groups={force_group})
    energy = energy.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

    return energy


def turn_parm_off(
    nb_force, atoms_list, parm_type="all", parm_off=True, exception_off=True
):
    if parm_off:
        for atom in range(nb_force.getNumParticles()):
            if atom in atoms_list:
                charge, sigma, epsilon = nb_force.getParticleParameters(atom)
                if parm_type == "all":
                    nb_force.setParticleParameters(atom, 0.0, sigma, 0.0)
                elif parm_type == "elec":
                    nb_force.setParticleParameters(atom, 0.0, sigma, epsilon)
                elif parm_type == "vdw":
                    nb_force.setParticleParameters(atom, charge, sigma, 0.0)

    if exception_off:
        for exception in range(nb_force.getNumExceptions()):
            atomi, atomj, charge, sigma, epsilon = nb_force.getExceptionParameters(
                exception
            )
            if atomi in atoms_list or atomj in atoms_list:
                if parm_type == "all":
                    nb_force.setExceptionParameters(
                        exception, atomi, atomj, 0.0, sigma, 0.0
                    )
                elif parm_type == "elec":
                    nb_force.setExceptionParameters(
                        exception, atomi, atomj, 0.0, sigma, epsilon
                    )
                elif parm_type == "vdw":
                    nb_force.setExceptionParameters(
                        exception, atomi, atomj, charge, sigma, 0.0
                    )


# Load Coordinates and XML file
if rank == 0:
    logger.info("Loading XML and PDB files...")

window = sys.argv[1]
pdbfile = PDBFile(f"windows/{window}/restrained.pdb")
with open(f"windows/{window}/restrained.xml", "r") as f:
    system = XmlSerializer.deserialize(f.read())

# Atom indices
if rank == 0:
    logger.info("Defining atom indices of molecules in system...")

guest_resname = "AMT"
host_resname = "CB7"
guest_mol = [
    atom.index
    for atom in pdbfile.topology.atoms()
    if atom.residue.name == guest_resname
]
host_mol = [
    atom.index for atom in pdbfile.topology.atoms() if atom.residue.name == host_resname
]
solvent = [
    atom.index
    for atom in pdbfile.topology.atoms()
    if atom.residue.name not in [guest_resname, host_resname]
]
host_guest_mol = guest_mol + host_mol
guest_solvent_mol = guest_mol + solvent
host_solvent_mol = host_mol + solvent
all_mol = guest_mol + host_mol + solvent

if rank == 0:
    logger.info("Splitting LJ forces to different components...")

# Valence terms
harmonic_bond = [
    force for force in system.getForces() if isinstance(force, HarmonicBondForce)
][0]
harmonic_angle = [
    force for force in system.getForces() if isinstance(force, HarmonicAngleForce)
][0]
periodic_torsion = [
    force for force in system.getForces() if isinstance(force, PeriodicTorsionForce)
][0]
harmonic_bond.setForceGroup(5)
harmonic_angle.setForceGroup(6)
periodic_torsion.setForceGroup(7)

# Nonbonded force
nonbonded_Elec = [
    force for force in system.getForces() if isinstance(force, NonbondedForce)
][0]
nonbonded_Elec14 = deepcopy(nonbonded_Elec)
nonbonded_LJ = deepcopy(nonbonded_Elec)
nonbonded_LJ14 = deepcopy(nonbonded_Elec)

# Elec intermol - remove vdW and Elec 1-4 interactions
nonbonded_Elec.setForceGroup(1)
turn_parm_off(
    nonbonded_Elec, all_mol, parm_type="elec", parm_off=False, exception_off=True
)
turn_parm_off(
    nonbonded_Elec, all_mol, parm_type="vdw", parm_off=True, exception_off=True
)

# Elec intramol - remove vdW and Elec intermol interactions
nonbonded_Elec14.setForceGroup(2)
turn_parm_off(
    nonbonded_Elec14, all_mol, parm_type="elec", parm_off=True, exception_off=False
)
turn_parm_off(
    nonbonded_Elec14, all_mol, parm_type="vdw", parm_off=True, exception_off=True
)
system.addForce(nonbonded_Elec14)

# LJ intermol - remove vdW 1-4 interactions and Elec interactions
nonbonded_LJ.setForceGroup(3)
turn_parm_off(
    nonbonded_LJ, all_mol, parm_type="elec", parm_off=False, exception_off=False
)
turn_parm_off(
    nonbonded_LJ, all_mol, parm_type="vdw", parm_off=True, exception_off=False
)
system.addForce(nonbonded_LJ)

# LJ intramol - remove vdW and Elec 1-4 interactions
nonbonded_LJ14.setForceGroup(4)
turn_parm_off(
    nonbonded_LJ14, all_mol, parm_type="elec", parm_off=True, exception_off=True
)
turn_parm_off(
    nonbonded_LJ14, all_mol, parm_type="vdw", parm_off=True, exception_off=False
)
system.addForce(nonbonded_LJ14)

# Simulation Object
simulation = create_simulation(system, pdbfile)

# Load 1 microsecond long trajectory
if rank == 0:
    logger.info("Loading Trajectories")

traj_list = [f"windows/{window}/production-v{i+1}.dcd" for i in range(200)]
if rank == 0:
    logger.info(f"{traj_list}")
universe = Universe(f"windows/{window}/restrained.pdb", traj_list)

n_frames = universe.trajectory.n_frames

# Specify array sizes for MPI nodes
if rank == 0:
    logger.info(f"Total number of frames {n_frames}")

num_per_rank = int(n_frames / size)
lower_bound = rank * num_per_rank
upper_bound = (rank + 1) * num_per_rank
report_freq = int(num_per_rank / 100)

potential = {
    "bond": np.zeros(n_frames),
    "angle": np.zeros(n_frames),
    "dihedral": np.zeros(n_frames),
    "Elec14": np.zeros(n_frames),
    "LJ14": np.zeros(n_frames),
    "Elec": np.zeros(n_frames),
    "LJ": np.zeros(n_frames),
}
potential_chunk = {
    "bond": np.zeros(num_per_rank),
    "angle": np.zeros(num_per_rank),
    "dihedral": np.zeros(num_per_rank),
    "Elec14": np.zeros(num_per_rank),
    "LJ14": np.zeros(num_per_rank),
    "Elec": np.zeros(num_per_rank),
    "LJ": np.zeros(num_per_rank),
}

# Scatter array to child node
for decomp in potential.keys():
    comm.Scatter(
        [potential[decomp], num_per_rank, MPI.DOUBLE],
        [potential_chunk[decomp], num_per_rank, MPI.DOUBLE],
        root=0,
    )

# Calculate Potential Energy
if rank == 0:
    logger.info("Calculating Potential Energy...")

comm.Barrier()

for i, frame in enumerate(tqdm(universe.trajectory[lower_bound:upper_bound])):
    if rank == 0 and i % report_freq == 0:
        logger.info(f"Completed: {i/num_per_rank*100:5.2f}%")

    xyz = frame.positions
    box = frame.dimensions

    # Nonbonded interactions
    potential_chunk["Elec"][i] = get_energy(simulation, xyz, box, force_group=1)
    potential_chunk["Elec14"][i] = get_energy(simulation, xyz, box, force_group=2)
    potential_chunk["LJ"][i] = get_energy(simulation, xyz, box, force_group=3)
    potential_chunk["LJ14"][i] = get_energy(simulation, xyz, box, force_group=4)

    # Valence interactions
    potential_chunk["bond"][i] = get_energy(simulation, xyz, box, force_group=5)
    potential_chunk["angle"][i] = get_energy(simulation, xyz, box, force_group=6)
    potential_chunk["dihedral"][i] = get_energy(simulation, xyz, box, force_group=7)

comm.Barrier()

if rank == 0:
    logger.info("Completed: 100.00%")

# Combine data back from child node
for decomp in potential.keys():
    comm.Allgather(
        [potential_chunk[decomp], MPI.DOUBLE], [potential[decomp], MPI.DOUBLE]
    )

# Save results to JSON
if rank == 0:
    logger.info("Saving results to JSON file...")
    with open(f"enthalpy_decomposition-{window}.json", "w") as f:
        dumped = json.dumps(potential, cls=NumpyEncoder)
        f.write(dumped)

    logger.info("Analysis completed.")

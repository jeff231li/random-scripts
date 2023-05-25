import json
import logging
from copy import deepcopy
from importlib import reload

import numpy as np
import simtk.unit as unit
from MDAnalysis import Universe
from mpi4py import MPI
from paprika.io import NumpyEncoder
from simtk.openmm import (
    CustomNonbondedForce,
    NonbondedForce,
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


def add_wca_repulsive(system, nb_force, force_group):
    wca_repulsive = CustomNonbondedForce(
        "step(r_cut - r) * U_sterics;"
        "U_sterics=4*epsilon*x*(x-1.0) + epsilon;"
        "x=(sigma/r)^6;"
        "r_cut = sigma*2^(1/6);"
        "sigma=0.5*(sigma1+sigma2);"
        "epsilon=sqrt(epsilon1*epsilon2);"
    )
    wca_repulsive.addPerParticleParameter("sigma")
    wca_repulsive.addPerParticleParameter("epsilon")
    wca_repulsive.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    wca_repulsive.setCutoffDistance(nb_force.getCutoffDistance())
    wca_repulsive.setUseLongRangeCorrection(nb_force.getUseDispersionCorrection())
    wca_repulsive.setForceGroup(force_group)

    # Set LJ parameters
    for atom in range(nb_force.getNumParticles()):
        charge, sigma, epsilon = nb_force.getParticleParameters(atom)
        wca_repulsive.addParticle([sigma, epsilon])

    # Transfer Exclusion
    for exception_index in range(nb_force.getNumExceptions()):
        iatom, jatom, chargeprod, sigma, epsilon = nb_force.getExceptionParameters(
            exception_index
        )
        wca_repulsive.addExclusion(iatom, jatom)

    system.addForce(wca_repulsive)


def add_wca_dispersive(system, nb_force, force_group):
    wca_dispersive = CustomNonbondedForce(
        "(step(r_cut - r) * -1*epsilon) + (step(r - r_cut) * U_sterics);"
        "U_sterics=4*epsilon*x*(x-1.0);"
        "x=(sigma/r)^6;"
        "r_cut = sigma*2^(1/6);"
        "sigma=0.5*(sigma1+sigma2);"
        "epsilon=sqrt(epsilon1*epsilon2);"
    )
    wca_dispersive.addPerParticleParameter("sigma")
    wca_dispersive.addPerParticleParameter("epsilon")
    wca_dispersive.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    wca_dispersive.setCutoffDistance(nb_force.getCutoffDistance())
    wca_dispersive.setUseLongRangeCorrection(nb_force.getUseDispersionCorrection())
    wca_dispersive.setForceGroup(force_group)

    # Set LJ parameters
    for atom in range(nb_force.getNumParticles()):
        charge, sigma, epsilon = nb_force.getParticleParameters(atom)
        wca_dispersive.addParticle([sigma, epsilon])

    # Transfer Exclusion
    for exception_index in range(nb_force.getNumExceptions()):
        iatom, jatom, chargeprod, sigma, epsilon = nb_force.getExceptionParameters(
            exception_index
        )
        wca_dispersive.addExclusion(iatom, jatom)

    system.addForce(wca_dispersive)


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
    energy = energy.getPotentialEnergy() / unit.kilocalories_per_mole

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

pdbfile = PDBFile("simulations/restrained.pdb")
with open("simulations/restrained.xml", "r") as f:
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

# Nonbonded force
nb_force = [force for force in system.getForces() if isinstance(force, NonbondedForce)][
    0
]

# Remove electrostatics and all 1-4 interactions
turn_parm_off(nb_force, all_mol, parm_type="elec")
turn_parm_off(nb_force, all_mol, parm_type="vdw", parm_off=False, exception_off=True)

# Guest LJ
guest_LJ = deepcopy(nb_force)
turn_parm_off(guest_LJ, host_solvent_mol)

# Host LJ
host_LJ = deepcopy(nb_force)
turn_parm_off(host_LJ, guest_solvent_mol)

# Solvent LJ
sol_LJ = deepcopy(nb_force)
turn_parm_off(sol_LJ, host_guest_mol)

# Guest-sol LJ
guest_sol_LJ = deepcopy(nb_force)
turn_parm_off(guest_sol_LJ, host_mol)

# Host-sol LJ
host_sol_LJ = deepcopy(nb_force)
turn_parm_off(host_sol_LJ, guest_mol)

# Host-Guest LJ
host_guest_LJ = deepcopy(nb_force)
turn_parm_off(host_guest_LJ, solvent)

# 01) System LJ-rep
add_wca_repulsive(system, nb_force, force_group=1)

# 02) System LJ-dis
add_wca_dispersive(system, nb_force, force_group=2)

# 03) Guest LJ-rep
add_wca_repulsive(system, guest_LJ, force_group=3)

# 04) Guest LJ-dis
add_wca_dispersive(system, guest_LJ, force_group=4)

# 05) Host LJ-rep
add_wca_repulsive(system, host_LJ, force_group=5)

# 06) Host LJ-dis
add_wca_dispersive(system, host_LJ, force_group=6)

# 07) Solvent LJ-rep
add_wca_repulsive(system, sol_LJ, force_group=7)

# 08) Solvent LJ-dis
add_wca_dispersive(system, sol_LJ, force_group=8)

# 09) Guest-sol LJ-rep
add_wca_repulsive(system, guest_sol_LJ, force_group=9)

# 10) Guest-sol LJ-dis
add_wca_dispersive(system, guest_sol_LJ, force_group=10)

# 11) Host-sol LJ-rep
add_wca_repulsive(system, host_sol_LJ, force_group=11)

# 12) Host-sol LJ-dis
add_wca_dispersive(system, host_sol_LJ, force_group=12)

# 13) Host-Guest LJ-rep
add_wca_repulsive(system, host_guest_LJ, force_group=13)

# 14) Host-Guest LJ-dis
add_wca_dispersive(system, host_guest_LJ, force_group=14)

# Simulation Object
simulation = create_simulation(system, pdbfile)

# Load 1 microsecond long trajectory
if rank == 0:
    logger.info("Loading Trajectories")
traj_list = [f"simulations/production-{i+1:03}.dcd" for i in range(100)]
universe = Universe("simulations/restrained.pdb", traj_list)

n_frames = universe.trajectory.n_frames

# Specify array sizes for MPI nodes
if rank == 0:
    logger.info(f"Total number of frames {n_frames}")

num_per_rank = int(n_frames / size)
lower_bound = rank * num_per_rank
upper_bound = (rank + 1) * num_per_rank
report_freq = int(num_per_rank / 100)

potential = {
    "total_rep": np.zeros(n_frames),
    "total_dis": np.zeros(n_frames),
    "guest_rep": np.zeros(n_frames),
    "guest_dis": np.zeros(n_frames),
    "host_rep": np.zeros(n_frames),
    "host_dis": np.zeros(n_frames),
    "sol_rep": np.zeros(n_frames),
    "sol_dis": np.zeros(n_frames),
    "guest_sol_rep": np.zeros(n_frames),
    "guest_sol_dis": np.zeros(n_frames),
    "host_sol_rep": np.zeros(n_frames),
    "host_sol_dis": np.zeros(n_frames),
    "host_guest_rep": np.zeros(n_frames),
    "host_guest_dis": np.zeros(n_frames),
}
potential_chunk = {
    "total_rep": np.zeros(num_per_rank),
    "total_dis": np.zeros(num_per_rank),
    "guest_rep": np.zeros(num_per_rank),
    "guest_dis": np.zeros(num_per_rank),
    "host_rep": np.zeros(num_per_rank),
    "host_dis": np.zeros(num_per_rank),
    "sol_rep": np.zeros(num_per_rank),
    "sol_dis": np.zeros(num_per_rank),
    "guest_sol_rep": np.zeros(num_per_rank),
    "guest_sol_dis": np.zeros(num_per_rank),
    "host_sol_rep": np.zeros(num_per_rank),
    "host_sol_dis": np.zeros(num_per_rank),
    "host_guest_rep": np.zeros(num_per_rank),
    "host_guest_dis": np.zeros(num_per_rank),
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

    potential_chunk["total_rep"][i] = get_energy(simulation, xyz, box, force_group=1)
    potential_chunk["total_dis"][i] = get_energy(simulation, xyz, box, force_group=2)

    potential_chunk["guest_rep"][i] = get_energy(simulation, xyz, box, force_group=3)
    potential_chunk["guest_dis"][i] = get_energy(simulation, xyz, box, force_group=4)

    potential_chunk["host_rep"][i] = get_energy(simulation, xyz, box, force_group=5)
    potential_chunk["host_dis"][i] = get_energy(simulation, xyz, box, force_group=6)

    potential_chunk["sol_rep"][i] = get_energy(simulation, xyz, box, force_group=7)
    potential_chunk["sol_dis"][i] = get_energy(simulation, xyz, box, force_group=8)

    potential_chunk["guest_sol_rep"][i] = get_energy(
        simulation, xyz, box, force_group=9
    )
    potential_chunk["guest_sol_dis"][i] = get_energy(
        simulation, xyz, box, force_group=10
    )

    potential_chunk["host_sol_rep"][i] = get_energy(
        simulation, xyz, box, force_group=11
    )
    potential_chunk["host_sol_dis"][i] = get_energy(
        simulation, xyz, box, force_group=12
    )

    potential_chunk["host_guest_rep"][i] = get_energy(
        simulation, xyz, box, force_group=13
    )
    potential_chunk["host_guest_dis"][i] = get_energy(
        simulation, xyz, box, force_group=14
    )

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
    with open("enthalpy_WCA.json", "w") as f:
        dumped = json.dumps(potential, cls=NumpyEncoder)
        f.write(dumped)

    logger.info("Analysis completed.")

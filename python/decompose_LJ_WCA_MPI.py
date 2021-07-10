import json
from copy import deepcopy

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

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def add_wca_repulsive(system):
    nonbonded = [
        force for force in system.getForces() if isinstance(force, NonbondedForce)
    ][0]

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
    wca_repulsive.setCutoffDistance(nonbonded.getCutoffDistance())
    wca_repulsive.setUseLongRangeCorrection(nonbonded.getUseDispersionCorrection())

    # Set LJ parameters
    for atom in range(nonbonded.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(atom)
        wca_repulsive.addParticle([sigma, epsilon])

    # Transfer Exclusion
    for exception_index in range(nonbonded.getNumExceptions()):
        iatom, jatom, chargeprod, sigma, epsilon = nonbonded.getExceptionParameters(
            exception_index
        )
        wca_repulsive.addExclusion(iatom, jatom)

    system.addForce(wca_repulsive)


def add_wca_dispersive(system):
    nonbonded = [
        force for force in system.getForces() if isinstance(force, NonbondedForce)
    ][0]

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
    wca_dispersive.setCutoffDistance(nonbonded.getCutoffDistance())
    wca_dispersive.setUseLongRangeCorrection(nonbonded.getUseDispersionCorrection())

    # Set LJ parameters
    for atom in range(nonbonded.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(atom)
        wca_dispersive.addParticle([sigma, epsilon])

    # Transfer Exclusion
    for exception_index in range(nonbonded.getNumExceptions()):
        iatom, jatom, chargeprod, sigma, epsilon = nonbonded.getExceptionParameters(
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
    )

    return simulation


def get_potential_energy(simulation, xyz, box):
    box_vectors = unit.Quantity(
        value=[Vec3(box[0], 0, 0), Vec3(0, box[1], 0), Vec3(0, 0, box[2])],
        unit=unit.angstrom,
    )

    simulation.context.setPositions(xyz * unit.angstrom)
    simulation.context.setPeriodicBoxVectors(*box_vectors)

    energy = simulation.context.getState(getEnergy=True, groups={0})
    energy = energy.getPotentialEnergy() / unit.kilocalories_per_mole

    return energy


def turn_parameters_off(
    system, atoms_list, parm_type="all", parm_off=True, exception_off=True
):
    nb_force = [
        force for force in system.getForces() if isinstance(force, NonbondedForce)
    ][0]

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
pdbfile = PDBFile("simulations/restrained.pdb")
with open("simulations/restrained.xml", "r") as f:
    system = XmlSerializer.deserialize(f.read())

# Atom indices
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

# Remove all forces except `NonbondedForce`
system.removeForce(4)
system.removeForce(0)
system.removeForce(0)
system.removeForce(0)

# 01) System LJ
system_LJ = deepcopy(system)
turn_parameters_off(system_LJ, all_mol, parm_type="elec")
turn_parameters_off(system_LJ, all_mol, parm_type="vdw", parm_off=False, exception_off=True)
sim_total_LJ = create_simulation(system_LJ, pdbfile)

# 02) System LJ14
system_LJ14 = deepcopy(system)
turn_parameters_off(system_LJ14, all_mol, parm_type="elec")
turn_parameters_off(system_LJ14, all_mol, parm_type="vdw", parm_off=True, exception_off=False)
sim_total_LJ14 = create_simulation(system_LJ14, pdbfile)

# 03) System LJ-rep
system_rep = deepcopy(system_LJ)
add_wca_repulsive(system_rep)
turn_parameters_off(system_rep, all_mol)
sim_total_LJrep = create_simulation(system_rep, pdbfile)

# 04) System LJ-dis
system_dis = deepcopy(system_LJ)
add_wca_dispersive(system_dis)
turn_parameters_off(system_dis, all_mol)
sim_total_LJdis = create_simulation(system_dis, pdbfile)

# 05) Guest LJ
system_guest_LJ = deepcopy(system_LJ)
turn_parameters_off(system_guest_LJ, host_solvent_mol)
sim_guest_LJ = create_simulation(system_guest_LJ, pdbfile)

# 06) Guest LJ14
system_guest_LJ14 = deepcopy(system)
turn_parameters_off(system_guest_LJ14, all_mol, parm_type="elec")
turn_parameters_off(system_guest_LJ14, host_solvent_mol)
turn_parameters_off(
    system_guest_LJ14, guest_mol, parm_type="vdw", parm_off=True, exception_off=False
)
sim_guest_LJ14 = create_simulation(system_guest_LJ14, pdbfile)

# 07) Guest LJ-rep
system_guest_LJrep = deepcopy(system_guest_LJ)
add_wca_repulsive(system_guest_LJrep)
turn_parameters_off(system_guest_LJrep, all_mol)
sim_guest_LJrep = create_simulation(system_guest_LJrep, pdbfile)

# 08) Guest LJ-dis
system_guest_LJdis = deepcopy(system_guest_LJ)
add_wca_dispersive(system_guest_LJdis)
turn_parameters_off(system_guest_LJdis, all_mol)
sim_guest_LJdis = create_simulation(system_guest_LJdis, pdbfile)

# 09) Host LJ
system_host_LJ = deepcopy(system_LJ)
turn_parameters_off(system_host_LJ, guest_solvent_mol)
sim_host_LJ = create_simulation(system_host_LJ, pdbfile)

# 10) Host LJ14
system_host_LJ14 = deepcopy(system)
turn_parameters_off(system_host_LJ14, all_mol, parm_type="elec")
turn_parameters_off(system_host_LJ14, guest_solvent_mol)
turn_parameters_off(
    system_host_LJ14, host_mol, parm_type="vdw", parm_off=True, exception_off=False
)
sim_host_LJ14 = create_simulation(system_host_LJ14, pdbfile)

# 11) Host LJ-rep
system_host_LJrep = deepcopy(system_host_LJ)
add_wca_repulsive(system_host_LJrep)
turn_parameters_off(system_host_LJrep, all_mol)
sim_host_LJrep = create_simulation(system_host_LJrep, pdbfile)

# 12) Host LJ-dis
system_host_LJdis = deepcopy(system_host_LJ)
add_wca_dispersive(system_host_LJdis)
turn_parameters_off(system_host_LJdis, all_mol)
sim_host_LJdis = create_simulation(system_host_LJdis, pdbfile)

# 13) Solvent LJ
system_solvent_LJ = deepcopy(system_LJ)
turn_parameters_off(system_solvent_LJ, host_guest_mol)
sim_solvent_LJ = create_simulation(system_solvent_LJ, pdbfile)

# 14) Solvent LJ-rep
system_solvent_LJrep = deepcopy(system_solvent_LJ)
add_wca_repulsive(system_solvent_LJrep)
turn_parameters_off(system_solvent_LJrep, all_mol)
sim_solvent_LJrep = create_simulation(system_solvent_LJrep, pdbfile)

# 15) Solvent LJ-dis
system_solvent_LJdis = deepcopy(system_solvent_LJ)
add_wca_dispersive(system_solvent_LJdis)
turn_parameters_off(system_solvent_LJdis, all_mol)
sim_solvent_LJdis = create_simulation(system_solvent_LJdis, pdbfile)

# Load trajectory
traj_list = [f"simulations/production-{i+1:03}.dcd" for i in range(1)]
universe = Universe("simulations/restrained.pdb", traj_list)

n_frames = universe.trajectory.n_frames

# Configure MPI sizes for nodes
if rank == 0:
    print(f"Total number of frames {n_frames}")

num_per_rank = int(n_frames / size)
lower_bound = rank * num_per_rank
upper_bound = (rank + 1) * num_per_rank

potential = {
    "total": np.zeros(n_frames),
    "total_14": np.zeros(n_frames),
    "total_rep": np.zeros(n_frames),
    "total_dis": np.zeros(n_frames),
    "guest_total": np.zeros(n_frames),
    "guest_14": np.zeros(n_frames),
    "guest_rep": np.zeros(n_frames),
    "guest_dis": np.zeros(n_frames),
    "host_total": np.zeros(n_frames),
    "host_14": np.zeros(n_frames),
    "host_rep": np.zeros(n_frames),
    "host_dis": np.zeros(n_frames),
    "sol_total": np.zeros(n_frames),
    "sol_rep": np.zeros(n_frames),
    "sol_dis": np.zeros(n_frames),
}
potential_chunk = {
    "total": np.zeros(num_per_rank),
    "total_14": np.zeros(num_per_rank),
    "total_rep": np.zeros(num_per_rank),
    "total_dis": np.zeros(num_per_rank),
    "guest_total": np.zeros(num_per_rank),
    "guest_14": np.zeros(num_per_rank),
    "guest_rep": np.zeros(num_per_rank),
    "guest_dis": np.zeros(num_per_rank),
    "host_total": np.zeros(num_per_rank),
    "host_14": np.zeros(num_per_rank),
    "host_rep": np.zeros(num_per_rank),
    "host_dis": np.zeros(num_per_rank),
    "sol_total": np.zeros(num_per_rank),
    "sol_rep": np.zeros(num_per_rank),
    "sol_dis": np.zeros(num_per_rank),
}

# Scatter array to child node
for decomp in potential.keys():
    comm.Scatter(
        [potential[decomp], num_per_rank, MPI.DOUBLE],
        [potential_chunk[decomp], num_per_rank, MPI.DOUBLE],
        root=0,
    )

# Calculate Potential Energy
comm.Barrier()

for i, frame in enumerate(tqdm(universe.trajectory[lower_bound:upper_bound])):
    xyz = frame.positions
    box = frame.dimensions

    potential_chunk["total"][i]       = get_potential_energy(sim_total_LJ,      xyz, box)
    potential_chunk["total_14"][i]    = get_potential_energy(sim_total_LJ14,    xyz, box)
    potential_chunk["total_rep"][i]   = get_potential_energy(sim_total_LJrep,   xyz, box)
    potential_chunk["total_dis"][i]   = get_potential_energy(sim_total_LJdis,   xyz, box)

    potential_chunk["guest_total"][i] = get_potential_energy(sim_guest_LJ,      xyz, box)
    potential_chunk["guest_14"][i]    = get_potential_energy(sim_guest_LJ14,    xyz, box)
    potential_chunk["guest_rep"][i]   = get_potential_energy(sim_guest_LJrep,   xyz, box)
    potential_chunk["guest_dis"][i]   = get_potential_energy(sim_guest_LJdis,   xyz, box)

    potential_chunk["host_total"][i]  = get_potential_energy(sim_host_LJ,       xyz, box)
    potential_chunk["host_14"][i]     = get_potential_energy(sim_host_LJ14,     xyz, box)
    potential_chunk["host_rep"][i]    = get_potential_energy(sim_host_LJrep,    xyz, box)
    potential_chunk["host_dis"][i]    = get_potential_energy(sim_host_LJdis,    xyz, box)

    potential_chunk["sol_total"][i]   = get_potential_energy(sim_solvent_LJ,    xyz, box)
    potential_chunk["sol_rep"][i]     = get_potential_energy(sim_solvent_LJrep, xyz, box)
    potential_chunk["sol_dis"][i]     = get_potential_energy(sim_solvent_LJdis, xyz, box)


comm.Barrier()

# Combine data back from child node
for decomp in potential.keys():
    comm.Allgather(
        [potential_chunk[decomp], MPI.DOUBLE], [potential[decomp], MPI.DOUBLE]
    )

# Save results to JSON
if rank == 0:
    with open("enthalpy_WCA.json", "w") as f:
        dumped = json.dumps(potential, cls=NumpyEncoder)
        f.write(dumped)

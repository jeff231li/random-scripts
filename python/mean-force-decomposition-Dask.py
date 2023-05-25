#!/usr/bin/env python
# coding: utf-8

import time
import logging
from importlib import reload

import dask
import numpy as np
import simtk.openmm as openmm
import simtk.openmm.app as app
import simtk.unit as openmm_unit
from dask.distributed import Client
from MDAnalysis import Universe
from scipy.integrate import cumtrapz

reload(logging)

logger = logging.getLogger()
logging.basicConfig(
    filename="analysis-dask.log",
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    level=logging.INFO,
)


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


@dask.delayed
def analyze_block(
    blockslice,
    n_frames_per_block,
    universe,
    pdb_file_path,
    system_xml,
):
    # Load PDB and Load XML file
    pdbfile = app.PDBFile(pdb_file_path)
    with open(system_xml, "r") as f:
        system = openmm.XmlSerializer.deserialize(f.read())

    # Atom indices
    host_resname = "CB7"
    guest_resname = "DM1"
    solvent_resname = "HOH"

    all_atoms = universe.select_atoms("all")
    guest = universe.select_atoms(f"resname {guest_resname}")
    host = universe.select_atoms(f"resname {host_resname} and element C")

    dummy_particle = [
        atom.index
        for atom in pdbfile.topology.atoms()
        if atom.residue.name == guest_resname
    ]
    host_indices = [
        atom.index
        for atom in pdbfile.topology.atoms()
        if atom.residue.name == host_resname
    ]
    solvent_indices = [
        atom.index
        for atom in pdbfile.topology.atoms()
        if atom.residue.name == solvent_resname
    ]

    r = np.linspace(0, 15, 61)
    z_axis = np.array([0, 0, 1])
    particle_radius = 4.5 * openmm_unit.angstrom
    wca_sigma = 3.0 * openmm_unit.angstrom
    wca_epsilon = 0.1 * openmm_unit.kilocalorie_per_mole

    # Add WCA repulsive core
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

    # Create Bins
    n_bins = np.zeros(len(r))
    mean_force_host = np.zeros(len(r))
    mean_force_all_solvent = np.zeros(len(r))

    # Loop through trajectory
    for iframe, ts in enumerate(universe.trajectory[blockslice.start : blockslice.stop]):
        logger.info(f"Completed: {iframe/n_frames_per_block*100:5.2f}%")
        
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
        simulation.context.setPeriodicBoxVectors(*box_vectors)
        simulation.context.setPositions(all_atoms.positions * openmm_unit.angstrom)

        # Find location in bins
        vector = np.mean(guest.positions, axis=0) - np.mean(host.positions, axis=0)
        i_loc = np.where(
            np.min(abs(np.dot(vector, z_axis) - r)) == abs(np.dot(vector, z_axis) - r)
        )[0][0]
        n_bins[i_loc] += 1

        # Calculate Forces
        force_by_host = simulation.context.getState(
            getForces=True, groups={4}
        ).getForces()
        force_by_all_solvent = simulation.context.getState(
            getForces=True, groups={7}
        ).getForces()

        # Vector projection of forces to reaction coordinate
        mean_force_host[i_loc] += np.dot(
            np.array(
                force_by_host[dummy_particle[0]].value_in_unit(
                    openmm_unit.kilocalorie_per_mole / openmm_unit.angstrom
                )
            ),
            z_axis,
        )
        mean_force_all_solvent[i_loc] += np.dot(
            np.array(
                force_by_all_solvent[dummy_particle[0]].value_in_unit(
                    openmm_unit.kilocalorie_per_mole / openmm_unit.angstrom
                )
            ),
            z_axis,
        )

    return n_bins, mean_force_host, mean_force_all_solvent


def main(n_workers):
    # Create Universe
    window_list = [f"p{i:03}" for i in range(41)]
    universe = Universe(
        "p000/restrained.pdb",
        [f"{window}/production.dcd" for window in window_list],
    )

    # Specify array sizes for Dask nodes
    n_frames = universe.trajectory.n_frames
    n_blocks = n_workers
    n_frames_per_block = n_frames // n_blocks
    blocks = [
        range(i * n_frames_per_block, (i + 1) * n_frames_per_block)
        for i in range(n_blocks - 1)
    ]
    blocks.append(range((n_blocks - 1) * n_frames_per_block, n_frames))
    blocks

    # Setup Dask jobs
    jobs = []
    for bs in blocks:
        jobs.append(
            analyze_block(
                bs,
                n_frames_per_block,
                universe,
                "p000/restrained.pdb",
                "p000/restrained.xml",
            )
        )
    jobs = dask.delayed(jobs)

    # Do Analysis
    start = time.time()
    results = jobs.compute()
    end = time.time()

    logger.info(f"Total time taken: {(end-start)} seconds")

    # Gather Data from all processes
    for i in range(n_blocks):
        if i == 0:
            n_bins = results[0][0]
            mean_force_host = results[0][1]
            mean_force_all_solvent = results[0][2]
        else:
            n_bins += results[i][0]
            mean_force_host += results[i][1]
            mean_force_all_solvent += results[i][2]

    # Integrate Forces
    n_bins[n_bins == 0.0] = 1.0
    r = np.linspace(0, 15, 61)
    pmf_host = cumtrapz(-1 * mean_force_host / n_bins, x=r, initial=0.0)
    pmf_all_solvent = cumtrapz(-1 * mean_force_all_solvent / n_bins, x=r, initial=0.0)

    # Save to file
    np.savetxt(
        "mean_force-split-dask.txt",
        np.c_[
            r,
            pmf_host,
            pmf_all_solvent,
            pmf_host + pmf_all_solvent,
        ],
        fmt="%10.5f",
    )


if __name__ == "__main__":

    dask.config.set(scheduler="processes")

    n_workers = 20
    with Client(n_workers=n_workers, threads_per_worker=1) as client:
        main(n_workers=n_workers)

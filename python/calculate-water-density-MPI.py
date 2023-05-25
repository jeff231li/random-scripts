import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis import Universe
from MDAnalysis.analysis.align import AlignTraj
from mpi4py import MPI
from tqdm import tqdm

# MPI setting
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Specify Volume (Cylinder)
dr = 0.5
r_bins = np.arange(-10.0, 10.0 + dr, dr)
radius = 5.0
volume = np.pi * radius**2 * dr

# Allocate arrays
chunk_molecule_counts = np.zeros(len(r_bins))
molecule_counts = np.zeros(len(r_bins))

# Load trajectory and define atom selection
u = Universe("restrained.pdb", "production.dcd")
n_frames = u.trajectory.n_frames

ref = Universe("restrained.pdb")

cylinder = [
    u.select_atoms(
        f"resname HOH and byres cyzone {radius} {r_bin+dr} {r_bin} resname CB7 and element C N",
        updating=True,
    )
    for r_bin in r_bins
]

comm.Barrier()

# Align the system to CB7 as lab frame of reference
print("Aligning trajectory to reference frame ...")
aligner = AlignTraj(u, ref, select="resname CB7 and element C N", in_memory=True)
aligner.run()

comm.Barrier()

# Set upper and lower bound
num_per_rank = int(n_frames / size)
lower_bound = rank * num_per_rank
upper_bound = (rank + 1) * num_per_rank
report_freq = int(num_per_rank / 100)

comm.Barrier()

# Calculate the number of water molecules and density
if rank == 0:
    print("Calculating density profile ...")

for j, ts in enumerate(tqdm(u.trajectory[lower_bound:upper_bound])):
    for i in range(len(r_bins)):
        chunk_molecule_counts[i] += cylinder[i].n_atoms / 3.0

comm.Barrier()

# Sum up arrays from all processes
comm.Reduce(
    [chunk_molecule_counts, MPI.DOUBLE],
    [molecule_counts, MPI.DOUBLE],
    op=MPI.SUM,
    root=0,
)

comm.Barrier()

# Save results to file
if rank == 0:
    molecule_counts /= n_frames
    density = molecule_counts / volume

    np.savetxt(
        "water_density.txt", np.c_[r_bins, molecule_counts, density], fmt="%10.5f"
    )

    plt.figure(figsize=(8, 10))
    plt.subplot(211)
    plt.plot(r_bins, molecule_counts)
    plt.xlabel("z (Ang)")
    plt.ylabel("<N_water>")
    plt.subplot(212)
    plt.plot(r_bins, density)
    plt.xlabel("z (Ang)")
    plt.ylabel("water density (molecule/A^3)")
    plt.show()

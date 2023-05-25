from simtk.openmm.app import Modeller, PDBFile
from openff.toolkit.topology import Molecule

molecule = Molecule.from_file("TEMOA.mol2")
new_model = Modeller(
    molecule.to_topology().to_openmm(),
    molecule.conformers[0],
)

residues = [residue for residue in new_model.topology.residues()]
residues[0].id = 'OAM'
residues[0].index = 1

chains   = [chain for chain in new_model.topology.chains()]
chains[0].id = ' '
chains[0].index = 1

with open("TEMOA-conect.pdb", "w") as f:
    PDBFile.writeFile(
        new_model.topology,
        new_model.positions,
        f,
        keepIds=True,
    )

import sys
import parmed as pmd
from paprika.io as import load_restraints
from paprika.analysis import fe_calc
from paprika.restraints.utils import extract_guest_restraints

# Filenames etc.
structure_file = 'vac-dum'
guest_resname = 'AMT'
restraint_file = 'restraints.json'

# Extract restraints
structure = pmd.load_file(f'{structure_file}.prmtop', f'{structure_file}.rst7')
restraints = load_restraints(filepath=restraint_file)
guest_restraints = extract_guest_restraints(structure, guest_resname, restraints)

# Calculate ref state work
FE = fe_calc()
FE.compute_ref_state_work(guest_restraints, calc='ddm')
print(f"Ref state work: {FE.results['ref_state_work']:.2f} kcal/mol")

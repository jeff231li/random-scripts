# Load PBCTools plugin
package require pbctools

# Load (GROMACS) coordinates file
mol new xtc.gro

# Get PBC box vectors
set pbc_box [lindex [pbc get] 0]

# Extract (x,y,z) from pbc_box
set x [lindex $pbc_box 0]
set y [lindex $pbc_box 1]
set z [lindex $pbc_box 2]

# Get translational vector = half the box dimension
set vector "[expr $x/2] [expr $y/2] [expr $z/2]"

# Translate whole system by vector
set all [atomselect top "all"]
$all moveby $vector

# Re-wrap all particles in the system according to the PBC boundaries
pbc wrap

# Save shifted coordinates to a PDB file
$all writepdb xtc.shifted.pdb

# Save shifted coordinates to a GRO file
animate write gro xtc.shifted.gro beg 0 end -1 top

exit

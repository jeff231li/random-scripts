#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################
# Amber topology
amber              on
parmfile           $ambertopology
ambercoor          $ambercoor
readexclusions     yes
scnb               2.0

# Restarting from previous run option(s)
if {$restart == 0} {
  # Periodic Boundary Conditions
  # NOTE: Do not set the periodic cell basis if you have also 
  # specified an .xsc restart file!
  cellBasisVector1   30.0   0.0   0.0
  cellBasisVector2    0.0  30.0   0.0
  cellBasisVector3    0.0   0.0  30.0
  cellOrigin          0.0   0.0   0.0 

  # NOTE: Do not set the initial velocity temperature if you 
  # have also specified a .vel restart file!
  temperature       $temperature
  firsttimestep     0

} else {
  # Restart file(s)
  binCoordinates    $restartfile.restart.coor
  binVelocities     $restartfile.restart.vel   ;# remove the "temperature" entry if you use this!
  extendedSystem    $restartfile.restart.xsc

  # Initial temperature
  firsttimestep     $firsttimestep
  #temperature       $temperature
}

# Force-Field Parameters
exclude             scaled1-4
1-4scaling          0.833333
cutoff              9.0
switching           off
pairlistdist        11.0

# Integrator Parameters
timestep            $Dtime
rigidBonds          all
nonbondedFreq       1
fullElectFrequency  2
stepspercycle       $NStepsCycle

# PME (for full-system periodic electrostatics)
PME                 yes
PMEGridSpacing      1.0
margin              1.0

wrapWater           on
wrapAll             on

# Constant Temperature Control
langevin            on
langevinDamping     1.0
langevinTemp        $temperature
langevinHydrogen    off

# Constant Pressure Control (variable volume)
if {$npt == 1} {
  StrainRate            0.0 0.0 0.0
  useGroupPressure      yes     ;# needed for 2fs steps
  useFlexibleCell       yes     ;# no for water box, yes for membrane
  useConstantArea       yes     ;# no for water box, yes for membrane

  langevinPiston        on
  langevinPistonTarget  1.013250 ;#  in bar -> 1 atm
  langevinPistonPeriod  200.
  langevinPistonDecay   20.
  langevinPistonTemp    $temperature
}

# Output Options
outputName          $outputname

restartfreq         $out_step
dcdfreq             $out_step
xstFreq             $out_step
outputEnergies      $out_step
outputPressure      $out_step

#############################################################
## EXTRA PARAMETERS                                        ##
#############################################################
# Fixed Atoms (set PDB O-column to 1)
if {$fix == 1} {
  fixedAtoms         on
  fixedAtomsFile     $ref_fix
  fixedAtomsCol      O
  fixedAtomsForces   on
}

# SMD for Phosphate head group
if {$smd == 1} {
  SMD                on
  SMDFile            $ref_smd
  SMDDir             0 0 1
  SMDk               0.1
  SMDVel             0
  SMDOutputFreq      $out_step
}

# Harmonic Restraints
if {$cons == 1} {
  constraints        on
  consexp            2
  consref            $ref_umb
  conskfile          $ref_umb
  conskcol           B
  selectConstraints  on
  selectConstrX      on
  selectConstrY      on
  selectConstrZ      on
}

# Extra bonds
if {$ext == 1} {
  extraBonds         on
  extraBondsFile     $ref_ext
}

# Collective variables
if {$Col == 1} {
  colvars            on
  colvarsConfig      colvars.tcl
}


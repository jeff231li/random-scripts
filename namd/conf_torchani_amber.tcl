#############################################################
## USER OPTIONS                                            ##
#############################################################
set AmberPrmtopFile   solvate.prmtop
set AmberInpcrdFile   solvate.rst7

set outputName        out_eq01
set restartName       out_prev

set npt               1
set qmmm              1
set Min               1
set Heat              1
set MD                1
set restart           0

set dTime             2.0
set temperature       298.15
set NStepsCycle       20

set outSteps          500
set minSteps          5000   
set heatSteps         500    ;# Steps per T(i)
set equilSteps        50000  

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################
# Amber topology
amber              on
parmfile           $AmberPrmtopFile
ambercoor          $AmberInpcrdFile
readexclusions     yes
scnb               2.0

# Restarting from previous run option(s)
if {$restart == 0} {
  # Periodic Boundary Conditions
  # NOTE: Do not set the periodic cell basis if you have also 
  # specified an .xsc restart file!
  cellBasisVector1   35.0   0.0   0.0
  cellBasisVector2    0.0  35.0   0.0
  cellBasisVector3    0.0   0.0  35.0
  cellOrigin          0.0   0.0   0.0 

  # NOTE: Do not set the initial velocity temperature if you 
  # have also specified a .vel restart file!
  temperature       $temperature

} else {
  # Restart file(s)
  binCoordinates    $restartName.restart.coor
  binVelocities     $restartName.restart.vel   ;# remove the "temperature" entry if you use this!
  extendedSystem    $restartName.restart.xsc

  # Initial temperature
  #temperature       $temperature
}

# Force-Field Parameters
exclude             scaled1-4
1-4scaling          0.833333
cutoff              9.0
switching           off

# Integrator Parameters
timestep            $dTime
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
  useFlexibleCell       no      ;# no for water box, yes for membrane
  useConstantArea       no      ;# no for water box, yes for membrane

  langevinPiston        on
  langevinPistonTarget  1.013250 ;#  in bar -> 1 atm
  langevinPistonPeriod  200.
  langevinPistonDecay   20.
  langevinPistonTemp    $temperature
}

# Output Options
outputName          $outputName

restartfreq         $outSteps
dcdfreq             $outSteps
xstFreq             $outSteps
outputEnergies      $outSteps
outputPressure      $outSteps

#############################################################
## NNP/MM PROTOCOLS                                        ##
#############################################################
# "qmmm.pdb" should contain the molecule of interest tagged
# Make sure "client.py" and "server.py" is available in PATH
if {$qmmm == 1} {
  set rundir   $env(RUNDIR)
  qmforces     on
  qmParamPDB   qmmm.pdb
  qmSoftware   custom
  qmexecpath   client.py
  qmBaseDir    $rundir
  QMColumn     occ
  qmChargeMode none
  qmElecEmbed  off
}

#############################################################
## MD PROTOCOLS                                            ##
#############################################################
if {$Min == 1} {
  minimize   $minSteps
  reinitvels $temperature
}

if {$Heat == 1} {
  for {set tempi 50} {$tempi <= $temperature} {incr tempi 50} {
    langevinTemp $tempi
    reinitvels   $tempi
    run          $heatSteps
  }
  langevinTemp $temperature
  reinitvels   $temperature
}

if {$MD == 1} {
  run $equilSteps
}


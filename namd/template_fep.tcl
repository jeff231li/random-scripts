#############################################################
## ALCHEMICAL TRANSFORMATION                               ##
#############################################################
if {$Alch == 1} {
  alch                  on
  alchType              FEP
  alchFile              $ref_fep
  alchCol               B
  alchOutFile           $outputname.fepout
  alchOutFreq           50
  alchEquilSteps        $EquilSteps

  # Turn off Electrostatics Only
  alchVdwLambdaEnd      0.0
  alchElecLambdaStart   0.0
  alchVdwShiftCoeff     0.0
  alchDecouple          off

  # Turn off van der Waals Only (charges must be set to zero in topology)
  alchVdwLambdaEnd      1.0
  alchElecLambdaStart   0.0
  alchVdwShiftCoeff     5.0
  alchDecouple          off

  # Thermodynamics Integration
  set win 1
  foreach lambda $inp_lambda {
    set l1 $lambda
    set l2 [expr $lambda + $dlambda]

    print [format "Running FEP window %3s: Lambda1 %-6s Lambda2 %-6s \[dLambda %-6s\]"\
        $win $l1 $l2 [expr abs($l2 - $l1)]]

    alchLambda  $l1
    alchLambda2 $l2
    run         [expr $EquilSteps + $ProdSteps]
  }
}

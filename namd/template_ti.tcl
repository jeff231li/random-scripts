#############################################################
## ALCHEMICAL TRANSFORMATION                               ##
#############################################################
if {$Alch == 1} {
  alch                  on
  alchType              TI
  alchFile              $ref_fep
  alchCol               B
  alchOutFile           $outputname.tiout
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
  print [format "Running TI with lambda %-6s " $l1]
  alchLambda  $l1
  run         [expr $EquilSteps + $ProdSteps]
}


set f [open "angles.dat" w]

foreach angle [label list Angles] {
  # Angle atom indices
  set i1 [lindex [split [lindex $angle 0]] 1]
  set i2 [lindex [split [lindex $angle 1]] 1]
  set i3 [lindex [split [lindex $angle 2]] 1]
  
  # Angle atom names
  set a1 [[atomselect top "index $i1"] get name]
  set a2 [[atomselect top "index $i2"] get name]
  set a3 [[atomselect top "index $i3"] get name]
  
  # Angle atom resid
  set r1 [[atomselect top "index $i1"] get resid]
  set r2 [[atomselect top "index $i2"] get resid]
  set r3 [[atomselect top "index $i3"] get resid]
  
  # Angle
  set ang [lindex $angle 3]
  
  # Print
  puts $f ":${r1}@${a1} :${r2}@${a2} :${r3}@${a3} $ang"
}

close $f
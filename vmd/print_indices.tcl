set indices ""
foreach atom [label list Atoms] {
    set index [lindex [split [lindex $atom 0]] 1]
    set indices "${indices} ${index}"
}
puts $indices

#!/usr/bin/tclsh
# load.tcl

# display
display reposition 600 1100
display resize 500 500

# menus
menu graphics on
menu graphics move 1413 439
menu main move 1276 95

# molecule
set psfname Polymer.psf
set xyzname lB32.config
mol load psf ${psfname} xyz ${xyzname}
puts "PSF files: ${psfname}"
puts "XYZ files: ${xyzname}"

# repsentations
mol modselect 0 0 type 0
mol modcolor 0 0 Type
mol color Type
for {set x 1} {$x <= 3} {incr x} {
mol selection type $x
mol addrep 0
}
mol showrep 0 0 off
mol showrep 0 1 off

# style
mol modstyle 2 0 {CPK 0.1 0.3 10}
mol modstyle 3 0 {CPK 0.1 0.3 10}

# background
# color Display Background white
display backgroundgradient on
axes location off

# box
pbc set {4.0 4.0 4.0}
pbc box -color white -width 2 -material BrushedMetal


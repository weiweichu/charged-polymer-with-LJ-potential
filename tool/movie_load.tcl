#!/usr/bin/tclsh
# load.tcl

# data file
set psfname ../N96.psf
set xyzname config.merge

# display
display reposition 20 820
display resize 680 360
display rendermode GLSL
display depthcue on
display projection orthographic

# background and axes
#display backgroundgradient on
color Display Background gray
axes location off

# graphics window
menu graphics on
menu graphics move 800 350
menu main move 800 50

# molecule
mol load psf ${psfname} xyz ${xyzname}
puts "PSF files: ${psfname}"
puts "XYZ files: ${xyzname}"

# create representations for polycations, polyanions & microions
mol modselect 0 0 type 0
mol modcolor 0 0 Charge
mol modmaterial 0 0 Translucent

mol selection type 1
mol addrep 0
mol modcolor 1 0 Charge
mol modmaterial 1 0 Translucent

mol selection type 2
mol addrep 0
mol modcolor 2 0 Charge
mol modmaterial 2 0 Transparent
#mol showrep 0 2 off

# style
mol modstyle 0 0 {CPK 0.3 0.1 10}
mol modstyle 1 0 {CPK 0.3 0.1 10}
mol modstyle 2 0 {CPK 0.1 0.3 10}

# bounding box
set bx 2.0
set by 2.0
set bz 8.0

graphics 0 color silver
graphics 0 materials on
graphics 0 material BrushedMetal

graphics 0 line "0.0 0.0 0.0" "0.0 0.0 $bz" width 2
graphics 0 line "0.0 0.0 $bz" "$bx 0.0 $bz" width 2
graphics 0 line "$bx 0.0 $bz" "$bx 0.0 0.0" width 2
graphics 0 line "$bx 0.0 0.0" "0.0 0.0 0.0" width 2

graphics 0 line "0.0 $by 0.0" "0.0 $by $bz" width 2
graphics 0 line "0.0 $by $bz" "$bx $by $bz" width 2
graphics 0 line "$bx $by $bz" "$bx $by 0.0" width 2
graphics 0 line "$bx $by 0.0" "0.0 $by 0.0" width 2

graphics 0 line "0.0 0.0 0.0" "0.0 $by 0.0" width 2
graphics 0 line "0.0 0.0 $bz" "0.0 $by $bz" width 2
graphics 0 line "$bx 0.0 $bz" "$bx $by $bz" width 2
graphics 0 line "$bx 0.0 0.0" "$bx $by 0.0" width 2

# box
#pbc set {2.0 2.0 8.0}
#pbc box -color white -width 2 -material BrushedMetal

#!/usr/bin/env python
import sys, glob
nwave = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0

f = open('f.gp' if len(sys.argv) < 3 else sys.argv[2], 'w')

f.write('set term x11 enhanced font "Hevitica, 20"\n\n')
f.write('set xlabel "psi" font "Hevitica, 20"\n')
f.write('set ylabel "free energy" font "Hevitica, 20"\n\n')

f.write('w = 2\n')
f.write('plot [0:][] 0 ti "",\\\n')

wid = [ max ( ([ int( x [ x.index('.') + 1 : ] ) for x in glob.glob(str(i) + '/doshistory/weight.*') ] + [-1]) ) for i in xrange(8) ]
for i in xrange(8):
   f.write('   "' + str(i) + '/doshistory/weight.' + str(wid[i]) + '" u 1:($2 + ' + str(nwave) + '*log($1) + 6.859) w l lw w ti "' + str(i) + ('",\\\n' if i < 7 else '"\n'))

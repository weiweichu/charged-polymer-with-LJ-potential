#!/usr/bin/env python
from EnergyData import *

ind = [1, 2, 3, 5] if len(sys.argv) < 2 else map(int, sys.argv[1:])
efile = EnergyData()
sys.stdout.write('{0:<6s}'.format('chiN'))
efile.outlabel(ind)
sys.stdout.write('-'*(5 + len(ind)*11) + '\n')
for i in xrange(8):
   efile.parse(str(i) + '/status.dat')
   sys.stdout.write('{0:<6.2f}'.format(11.0 + i))
   efile.seloutput(ind)

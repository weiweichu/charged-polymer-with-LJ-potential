#!/usr/bin/env python
from subprocess import Popen, PIPE
import os

nBeads     = 4096
nBlockBead = [16, 17]
nChainUnit = 8

for i in xrange(1, 8+1):

   f = open('param.' + str(i), 'w')
   f.write('----- System parameters -----\n')
   f.write('2.000000  2.000000   8.00000        boxL\n')
   f.write('6    6    24                        nGrid\n')
   f.write('3                                   nType\n')
   f.write('{0:<5d}{1:5d}                          nBlock\n'.format(nBlockBead[0], nBlockBead[1]))
   f.write('{0:<5d}                               nPolymer\n'.format(nChainUnit*i))
   f.write('{0:<7d}                             nNeutralSolvent\n'.format(nBeads - nChainUnit*i*sum(nBlockBead)))
   f.write('7                                   chargeCode\n')
   f.close()

   cmd = ['/home/qin/psfgen.py', 'param.' + str(i), 'n' + str(nChainUnit*i) + '.psf']
   p = Popen(cmd, stdout=PIPE, stderr=PIPE)
   print reduce(lambda x, y: x + ' ' + y, cmd)
   print 'n = {0:4d}   PID = {1:<8d}'.format(nChainUnit*i, p.pid)

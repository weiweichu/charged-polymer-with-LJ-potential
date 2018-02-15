#!/usr/bin/env python
import os, sys
from subprocess import call

if len(sys.argv) < 2:
   sys.exit('dir name not provided')
dtdir = sys.argv[1]
if dtdir[-1] != '/': dtdir += '/'

if not os.path.isdir('tga'):
   os.mkdir('tga')

nChainUnit = 8
nBeadPerChain = 33
nTotal = 4096*2

for n in xrange(8):
   fload = open('load.tcl', 'w')
   for line in open(os.path.dirname(os.path.realpath(__file__)) + '/load_template', 'r'):
      line = line.replace('PSFID', 'n' + str((n+1)*nChainUnit))
      line = line.replace('DIR', dtdir + str(n))
      line = line.replace('FILE', 'tga/' + str(n))
      fload.write(line)
   fload.close()
   call(['vmd', '-e', 'load.tcl'])

for n in xrange(8):
   print n, '{0:4d} {1:6.3f} {2:6.3f}'.format( (n+1)*nChainUnit,\
             float(nBeadPerChain*(n+1)*nChainUnit) / float(nTotal),\
             float(nBeadPerChain*(nChainUnit*(n+1) + nChainUnit/2 ) ) / float(nTotal) )

call('montage -label %f -frame 5 -tile 2x4 -geometry 100% tga/*.tga all.tga'.split(' '))
os.remove('load.tcl')

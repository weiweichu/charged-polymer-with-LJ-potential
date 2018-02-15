#!/usr/bin/env python
'''
PURPOSE
   Parse the files param & mu.out and find salt concentration and chemical potential.
'''
import sys
for fld in [ str(i) + '/' for i in xrange(6 if len(sys.argv) < 2 else int(sys.argv[1])) ]:
   nSalts, nSolvents = [ int(x.strip().split()[0]) for x in open(fld + 'param', 'r').readlines()[5:7] ]
   mu = float( open(fld + 'mu.out', 'r').readlines()[0].strip().split()[1] )
   print '{0:<7.4f} {1:<8.4f}'.format(float(nSalts)/float(nSalts + nSolvents), mu)

#!/usr/bin/env python
import os
from subprocess import call

dt = [[0.0]*6 for x in xrange(16)]
counter = 0
for i in xrange(2, 5):
   call('glean16.py run' + str(i) + ' | tail -16 > psi.' + str(i), shell=True)
   fdt = open('psi.' + str(i), 'r').readlines()
   for j in xrange(16):
      line = fdt[j].strip().split()
      for k in xrange(6):
         dt[j][k] += float(line[k])
   os.remove('psi.' + str(i))
   counter += 1

for j in xrange(16):
   for k in xrange(6):
      dt[j][k] /= float(counter)
   print reduce(lambda x, y: x + ' ' + y, map(lambda x: '{0:>10.4f}'.format(x), dt[j]))

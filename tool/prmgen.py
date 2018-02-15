#!/usr/bin/env python
import sys, os, shutil
from ParamFile import *

# parameter file
prmfile = ParamFile()
for i in xrange(8):
   if not os.path.isdir(str(i)):
      os.mkdir(str(i))
   prmfile.nGrid[2] = 9
   prmfile.resetSystem()
   prmfile.chiN[0][0] = 13.0 + float(i)
   prmfile.write(open(str(i) + '/param', 'w'))

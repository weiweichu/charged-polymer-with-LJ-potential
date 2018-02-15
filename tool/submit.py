#!/usr/bin/env python
import os, sys, random
from subprocess import call

# check the number of argument
if len(sys.argv) < 3:
   sys.exit('source or destination folder not provided')

if not os.path.isdir(sys.argv[1]):
   sys.exit('source folder non-existing')
src = sys.argv[1]
if src[-1] != '/': src += '/'

if os.path.isdir(sys.argv[2]):
   sys.exit('destination folder exist')
dst = sys.argv[2]
if dst[-1] != '/': dst += '/'
os.mkdir(dst)

# copy and modify job script
f = open(dst + 'job', 'w')
for line in open(src + 'job', 'r'):
   line = line.strip()
   i = line.find('swp')
   if i >= 0 and i + 3 < len(line):
      line = line.replace(line[i:], 'swp' + str(int(line[i+3:])+1))
   f.write(line + '\n')
f.close()

# copy and modify param script
f = open(dst + 'prmgen.py', 'w')
for line in open(src + 'prmgen.py', 'r'):
   if line.find('configFile') >= 0:
      i1 = line.find('run')
      i2 = line.find('/', i1)
      line = line.replace(line[i1:i2], src[:-1])

   if line.find('.seed') >= 0:
      i1 = line.find('-')
      i2 = line.find(' ', i1)
      line = line.replace(line[i1+1:i2], str(random.randint(10000000, 99999999)))

   f.write(line)
f.close()
os.chmod(dst + 'prmgen.py', 0755)

# Submit job
os.chdir(dst)
call('./prmgen.py')
call(['sbatch', 'job'])

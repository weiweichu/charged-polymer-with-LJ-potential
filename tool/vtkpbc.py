#!/usr/bin/env python
import sys
if len(sys.argv) < 6:
   sys.exit('Insufficient number of arguments.')
nshift = [ int(x) for x in sys.argv[3:6] ]

# read the source data
dim = [1]*3
data = []
ndata = 0
spacing = ''
scalar = ''
i = 0
for line in open(sys.argv[1], 'r'):
   i += 1
   if i == 5:
      line = line.strip().split()[1:]
      for j in xrange(3): dim[j] = int(line[j])
      ndata = reduce(lambda x, y: x*y, dim)
   elif i == 7:
      spacing = line
   elif i == 10:
      scalar = line
   elif i > 11 and i <= 11 + ndata:
      data.append(line)

# calculate the compact
def wrap(ix, iy, iz):
   i1 = ix % dim[0]
   i2 = iy % dim[1]
   i3 = iz % dim[2]
   return (i3 + dim[2]*(i2 + dim[1]*i1))

# duplicate the source data
f = open(sys.argv[2], 'w')
f.write('# vtk DataFile Version 2.0\n')
f.write('GridMC configuration\n')
f.write('ASCII\n')
f.write('DATASET STRUCTURED_POINTS\n')
fmt = reduce(lambda x, y: x + ' ' + y,
      map(lambda x: '{0:<5d}'.format(x).strip(), [ dim[x]*nshift[x] for x in xrange(3) ]) )
f.write('DIMENSIONS ' + fmt + '\n')
f.write('ORIGIN 0 0 0\n')
f.write(spacing)
f.write('\n')
f.write('POINT_DATA {0:<16d}'.format(ndata*reduce(lambda x, y: x*y, nshift)).strip() + '\n')
f.write(scalar)
f.write('LOOKUP_TABLE default\n')

for ix in xrange(dim[0]*nshift[0]):
   for iy in xrange(dim[1]*nshift[1]):
      for iz in xrange(dim[2]*nshift[2]):
         f.write(data[wrap(ix, iy, iz)])
f.close()


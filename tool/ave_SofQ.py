#!/usr/bin/env python
'''
PURPOSE
   Read and average structure factor data.
'''

# read folder run1
header = reduce(lambda x, y: x + y, [ [ fd for i in xrange(8) ] for fd in ['c00_07', 'c08_15'] ] )
folder = reduce(lambda x, y: x + y, [ [ str(i) for i in xrange(8) ] for j in xrange(2) ] )
SofQ = [ [] for i in xrange(16) ]
for i in xrange(16):
   for line in open(header[i] + '/run1/' + folder[i] + '/Sq.dat', 'r'):
      line = line.strip().split()
      SofQ[i].append( [ float(line[1]), float(line[2]) / 3.0 ] )
weight = 1.0 / 3.0

# read folders run2, run3, ...
for j in xrange(2, 4):
   for i in xrange(16):
      iLine = -1
      for line in open(header[i] + '/run' + str(j) + '/' + folder[i] + '/Sq.dat', 'r'):
         iLine += 1
         SofQ[i][iLine][1] += float(line.strip().split()[2])
   weight += 1.0

# normalize and output data
for i in xrange(len(SofQ)):
   for j in xrange(len(SofQ[i])):
      SofQ[i][j][1] /= weight

   f = open('SofQ_average/chiN_' + str(i) + '.dat', 'w')
   for j in xrange(len(SofQ[i])):
      f.write('{0:<12.4f} {1:<12.4f}\n'.format(SofQ[i][j][0], -SofQ[i][j][1]))
   f.close()


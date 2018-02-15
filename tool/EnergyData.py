#!/usr/bin/env python
'''
PURPOSE
   Parse an energy data file.
'''
import sys

class EnergyData:

   def __init__(self):
      '''
      PURPOSE
         Constructor.
      '''
      self.label    = ['frameID', 'totalQ', 'bonded', 'two-body', 'coulomb', 'pV', 'psi']
      self.ncolumn  = len(self.label)

      self.fname    = ''
      self.nskip    = 0
      self.count    = 0
      self.energy   = [0.0] * (self.ncolumn - 1)
      self.variance = [0.0] * (self.ncolumn - 1)

   def parse(self, datafile, nskipIn = 0):
      '''
      PURPOSE
         Parse the energy data file and do the average.
      '''
      self.fname = datafile
      self.nskip = nskipIn
      self.count = 0
      for i in xrange(len(self.energy)):
         self.energy[i]   = 0.0
         self.variance[i] = 0.0

      nline, x = 0, 0.0
      for line in open(self.fname, 'r'):
         nline += 1
         if nline <= self.nskip: continue
         line = line.strip().split()
         for i in xrange(self.ncolumn-1):
            x = float(line[i+1])
            self.energy[i]   += x
            self.variance[i] += x*x
         self.count += 1

      for i in xrange(self.ncolumn-1):
         self.energy[i]   /= float(self.count)
         self.variance[i] /= float(self.count)
         self.variance[i] -= self.energy[i]**2

   def output(self, o = sys.stdout):
      '''
      PURPOSE
         Output average.
      '''
      energy   = self.energy
      variance = self.variance
      print 'data file: ', self.fname
      for i in xrange(self.ncolumn-1):
         o.write('{0:12}: {1:8.4f} {2:8.4f}\n'.format(self.label[i+1], energy[i], variance[i]))
      o.write('{0:12}: {1:8d}\n'.format('data counts', self.count))


   def outlabel(self, sel, o = sys.stdout):
      '''
      PURPOSE
         Output labels.
      '''
      o.write(reduce(lambda x, y: x + ' ' + y, map(lambda i: '{0:>10s}'.format(self.label[i+1]), sel)) + '\n')


   def seloutput(self, sel, o = sys.stdout):
      '''
      PURPOSE
         Output average for a selected set of energy terms.
      '''
      o.write(reduce(lambda x, y: x + ' ' + y, map(lambda i: '{0:>10.4f}'.format(self.energy[i]), sel)) + '\n')


if __name__ == '__main__':
   if len(sys.argv) < 2:
      print 'No data file provided. Do nothing.'
   else:
      engfile = EnergyData()
      engfile.parse(sys.argv[1])
      engfile.seloutput([0, 1, 2, 4] if len(sys.argv) < 3 else map(int, sys.argv[2:]))

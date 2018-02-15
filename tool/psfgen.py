#!/usr/bin/env python
import sys

f = open('param' if len(sys.argv) < 2 else sys.argv[1], 'r').readlines()
if len(f) < 8: sys.exit('Param file corrupted')

nBlock   = map(int, f[4].strip().split()[:2])
nPolymer = int(f[5].strip().split()[0])
nNeutral = int(f[6].strip().split()[0])
chargeCode = int(f[7].strip().split()[0])
nSalt = 0
if chargeCode > 0 and (chargeCode != 4 or chargeCode != 7):
   if chargeCode < 10:
      print 'The charging code is not implemented'
   else:
      nSalt = chargeCode - 10

# Calculate the total number of beads
nBeadsPerChain = sum(nBlock)
nBeads  = nBeadsPerChain * nPolymer
nBeads += nNeutral
nBeads += nSalt
if chargeCode < 0: nBeads += -chargeCode

# Output file
outf = open('Polymer.psf' if len(sys.argv) < 3 else sys.argv[2], 'w')

# Output the total number of beads
outf.write('{0:8d} !NATOMS\n'.format(nBeads))

# Output bead information
iAtom = 0
if chargeCode < 0:
   for i in xrange(nPolymer):
      for j in xrange(nBlock[0]):
         iAtom += 1
         outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'POL', 'O', 0, 0.0, 1.0, 0))
      for j in xrange(nBlock[1]):
         iAtom += 1
         outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'POL', 'N', 1, 0.0, 1.0, 0))
elif chargeCode == 4 or chargeCode >= 10:
   for i in xrange(nPolymer):
      q = 1.0 if i < nPolymer/2 else -1.0
      atype = 0 if i < nPolymer/2 else 1
      for j in xrange(sum(nBlock)):
         outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'POL', 'O', atype, q, 1.0, 0))
elif chargeCode == 7:
   for i in xrange(nPolymer):
      q = 1.0 if i < nPolymer/2 else -1.0
      atype = 0 if i < nPolymer/2 else 1
      for j in xrange(sum(nBlock)):
         if j%2 == 0:
            outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'POL', 'O', atype, 0.0, 1.0, 0))
         else:
            outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'POL', 'O', atype, q, 1.0, 0))

if chargeCode < 0:
   for i in xrange(-chargeCode/2):
      iAtom += 1
      outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'SOL', 'S', 2, 1.0, 1.0, 0))
   for i in xrange(-chargeCode/2):
      iAtom += 1
      outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'SOL', 'S', 3, -1.0, 1.0, 0))
elif chargeCode == 4 or chargeCode == 7:
   for i in xrange(nNeutral):
      iAtom += 1
      outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'SOL', 'S', 2, 0.0, 1.0, 0))
elif chargeCode >= 10:
   for i in xrange(nSalt):
      iAtom += 1
      q = 1.0 if i < nSalt/2 else -1.0
      atype = 2 if i < nSalt/2 else 3
      outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'SAL', 'S', atype, q, 1.0, 0))
   for i in xrange(nNeutral):
      iAtom += 1
      outf.write('{0:8d}{1:5d}{2:2d}{3:>7}{4:>3}{5:8d}{6:15.6e}{7:10.3f}{8:12d}\n'.format(iAtom, 0, 1, 'SOL', 'S', 4, 0.0, 1.0, 0))

# Output bond information
nBonds = (nBeadsPerChain - 1) * nPolymer
outf.write('\n{0:8d} !NBONDS: bonds\n'.format(nBonds))
bstr = ''
iBead = 0
for i in xrange(nPolymer):
   for j in xrange(nBeadsPerChain):
      iBead += 1
      bstr += '{0:8d}'.format(iBead)
      if len(bstr) == 64:
         outf.write(bstr + '\n')
         bstr = ''
      if j > 0 and j < nBeadsPerChain - 1:
         bstr += '{0:8d}'.format(iBead)
         if len(bstr) == 64:
            outf.write(bstr + '\n')
            bstr = ''
if len(bstr) > 0:
   outf.write(bstr + '\n')

# Output angle and dihedral, and other information
outf.write('\n{0:8d} !NTHETA: angles\n'.format(0))
outf.write('\n{0:8d} !NPHI: dihedrals\n'.format(0))
outf.write('\n{0:8d} !NIMPR: impropers\n'.format(0))
outf.write('\n{0:8d} !HDON\n'.format(0))
outf.write('\n{0:8d} !HACC\n'.format(0))
outf.write('\n{0:8d} !HNB\n'.format(0))

outf.close()

#!/usr/bin/env python
'''
PURPOSE
   A GridMC parameter file.
USAGE
   (1) ./ParamFile filename
   (2) Write another script to create a ParamFile object, modify the attributes, and save as follows.
           prmobj = ParamFile()
           prmobj.nGrid[2] = 9
           prmobj.resetSystem()
           prmfile = open('filename', 'w')
           prmobj.write(prmfile)
           prmfile.close() 
'''
import sys

class ParamFile:

   def __init__(self):
      '''
      PURPOSE
         Constructor.
      '''
      self.ncolumn = 36

      self.phiNeutral = 0
      self.phiSalt    = 0
      self.nChainUnit = 128.0
      self.nGridRe    = 6.0

      self.systemHeader  = '----- System parameters -----'
      self.nGrid               = [6, 6, 6]
      self.nType               = 3
      self.nBlock              = [16, 16]
      self.useEwald            = 1
      self.ewald               = [4.0, 5, 5]
      self.configFile          = '0'
      self.seed                = -19791980

      self.thermodynamicHeader = '----- Thermodynamics parameters -----'
      self.kappaN              = 50.0
      self.sqrtNbar            = 128.0
      self.esStrength          = 0.0
      self.isIsobaric          = 0
      self.barostatPressure    = 4103.49825
      self.oprmFlag            = 2
      self.qstar               = 3.8
      self.waveIndex           = [0, 0, 1]
      self.dosFlag             = 1
      self.psi                 = [0.001, 0.50, 0.0001]
      self.flatnessTol         = 0.2
      self.nIterationTol       = 25
      self.sTolerance          = 0.8
      self.nMaxIteration       = 40

      self.simulationHeader    = '----- Simulation parameters -----'
      self.nMCmoves            = 1000
      self.pBeadMove           = 0.2
      self.beadMoveSize        = 0.8
      self.pChainFlipMove      = 0.4
      self.pReptationMove      = 0.4
      self.nReptationMax       = 5
      self.volumeMoveInterval  = 10
      self.volumeMoveSize      = 0.08
      self.flatnessInterval    = 10
      self.energyInterval      = 1000
      self.configInterval      = 5000
      self.nDiagnosis          = 0
      self.baseInterval        = 100
      self.sfactorName         = 'StructureFactor'
      self.sfactorInterval     = 100
      self.correlationType     = [0, 1]
      self.maxWaveIndex        = 3

      # Set derived parameter
      self.resetSystem()

   def resetSystem(self):
      '''
      PURPOSE
         Calculate derived parameters: number of chains, solvents, and salts, etc.
      '''
      self.nGridRe  = 6.0 * pow(self.nChainUnit / 128.0, 1.0/3.0)
      self.boxL     = [ float(i) / self.nGridRe for i in self.nGrid ]
      self.nPolymer = int(round(reduce(lambda x, y: x*y, self.boxL) * self.nChainUnit))

      if not self.phiNeutral == 0:
         self.nNeutralSolvent = int(round(float(self.nPolymer * sum(self.nBlock)) * self.phiNeutral))
      else:
         self.nNeutralSolvent = 0

      if not self.phiSalt == 0:
         self.chargeCode = -int(round(float(self.nPolymer * sum(self.nBlock)) * self.phiSalt))
         if not self.chargeCode % 2 == 0:
            self.chargeCode -= 1
         self.nType = 4
         if self.nNeutralSolvent > 0:
            self.nType += 1
      else:
         self.chargeCode = 0

      self.nPolymer *= sum(self.nBlock)
      self.nPolymer -= self.nNeutralSolvent
      if self.chargeCode < 0:
         self.nPolymer += self.chargeCode
      self.nPolymer = int(round(float(self.nPolymer) / float(sum(self.nBlock))))

      self.chiN = [ [0.0]*i for i in range(1, self.nType) ]

 
   def writeSystem(self, o = sys.stdout):
      '''
      PURPOSE
         Echo system parameters.
      '''
      o.write(self.systemHeader + '\n' )

      put_var = self.put_var
      put_var(o, '9.6f', self.boxL, 'boxL')
      put_var(o, '3d',   self.nGrid, 'nGrid')
      put_var(o, '3d',   self.nType, 'nType')
      put_var(o, '3d',   self.nBlock, 'nBlock')
      put_var(o, '10d',  self.nPolymer, 'nPolymer')
      put_var(o, '10d',  self.nNeutralSolvent, 'nNeutralSolvent')
      put_var(o, '10d',  self.chargeCode, 'chargeCode')
      put_var(o, '10d',  self.useEwald, 'useEwald')
      if self.useEwald > 0:
         o.write('{0:<7.3f} {1:<3d} {2:<3d}'.format(self.ewald[0], self.ewald[1], self.ewald[2]).ljust(self.ncolumn) + 'ewald\n')
      put_var(o, '',     self.configFile, 'configFile')
      put_var(o, '12d',  self.seed, 'seed')


   def writeThermodynamic(self, o = sys.stdout):
      '''
      PURPOSE
         Echo thermodyanmic parameters.
      '''
      o.write(self.thermodynamicHeader + '\n' )

      put_var = self.put_var

      put_var(o, '6.2f', self.chiN, 'chiN')
      put_var(o, '8.2f', self.kappaN, 'kappaN')
      put_var(o, '12.6f', self.sqrtNbar, 'sqrtNbar')
      put_var(o, '8.2f', self.esStrength, 'esStrength')
      put_var(o, '1d', self.isIsobaric, 'isIsobaric')

      if not self.isIsobaric == 0:
         put_var(o, '13.6f', self.barostatPressure, 'barostatPressure')

      put_var(o, '3d', self.oprmFlag, 'oprmFlag (order parameter flag)')
      if self.oprmFlag == 1 or self.oprmFlag == 2:
         put_var(o, '3d', self.waveIndex, 'waveIndex')
      elif self.oprmFlag < 0:
         put_var(o, '8.5f', self.qstar, 'qstar')

      put_var(o, '1d', self.dosFlag, 'dosFlag')
      if self.dosFlag > 0:
         put_var(o, '8.5f', self.psi, 'psi (order parameter range)')
         put_var(o, '6.3g', self.flatnessTol, 'flatnessTol')
         hist = [self.nIterationTol, self.sTolerance, self.nMaxIteration]
         put_var(o, '6.3g', hist, 'nIterationTol sTolerance nMaxIteration')


   def writeSimulation(self, o = sys.stdout):
      '''
      PURPOSE
         Echo simulation parameters.
      '''
      o.write(self.simulationHeader + '\n' )

      put_var = self.put_var
      put_var(o, '15d', self.nMCmoves, 'nMCmoves')
      put_var(o, '6.3f', self.pBeadMove, 'pBeadMove')
      put_var(o, '6.3f', self.beadMoveSize, 'beadMoveSize (in unit of bond length)')
      put_var(o, '6.3f', self.pChainFlipMove, 'pChainFlipMove')
      put_var(o, '6.3f', self.pReptationMove, 'pReptationMove')
      put_var(o, '10d', self.nReptationMax, 'nReptationMax')

      if not self.isIsobaric == 0:
         put_var(o, '10d', self.volumeMoveInterval, 'volumeMoveInterval')
         put_var(o, '6.3f', self.volumeMoveSize, 'volumeMoveSize (in Re^3)')

      if self.dosFlag > 0:
         put_var(o, '10d', self.flatnessInterval, 'flatnessInterval')

      put_var(o, '10d', self.energyInterval, 'energyInterval')
      put_var(o, '10d', self.configInterval, 'configInterval')
      put_var(o, '10d', self.nDiagnosis, 'nDiagnosis')

      if self.nDiagnosis > 0:
         put_var(o, '10d', self.baseInterval, 'baseInterval')
         put_var(o, '', self.sfactorName, 'sfactorName')
         put_var(o, '10d', self.sfactorInterval, 'sfactorInterval')
         put_var(o, '4d', self.correlationType, 'correlationType')
         put_var(o, '4d', self.maxWaveIndex, 'maxWaveIndex')


   def write(self, o = sys.stdout):
      '''
      PURPOSE
         Formatted outputs.
      '''
      self.writeSystem(o)
      o.write('\n')
      self.writeThermodynamic(o)
      o.write('\n')
      self.writeSimulation(o)


   def put_var(self, o, code, var, name):
      '''
      PURPOSE
         Formated outputs: scalar, vector & array.
      '''
      if isinstance(var, list):
         if isinstance(var[0], list):
            map(lambda i: self.put_var(o, code, var[i], name if i == 0 else ''), range(len(var)))
         else:
            fmt = reduce(lambda x, y: x + ' ' + y,
                     map(lambda x: ('{0:<' + code + '}').format(x), var))
            o.write(fmt.ljust(self.ncolumn) + name + '\n')
      else:
         self.put_var(o, code, [var], name)


if __name__ == '__main__':
   '''
   PURPOSE
      Echo parameters.
   '''
   if len(sys.argv) > 1:
      ParamFile().write(open(sys.argv[1], 'w'))
      print 'Parameters written to ' + sys.argv[1]
   else:
      ParamFile().write()


#ifndef CHAINFLIPMOVE_CPP
#define CHAINFLIPMOVE_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <math.h>
#include "../util/global.h"
#include "../util/Vector.h"

#include "ChainFlipMove.h"
#include "../simulation/McSystem.h"

namespace GridMC
{
   using namespace Util;

   /*
   * Constructor.
   */
   ChainFlipMove::ChainFlipMove(McSystem& sysIn) : Move(sysIn)
   { name_.assign("ChainFlipMove"); }

   /*
   * Default destructor.
   */
   ChainFlipMove::~ChainFlipMove()
   {}

   /*
   * Read step size.
   */
   void ChainFlipMove::readParam(std::istream& in)
   {
      #if 0
      char   comment[200];
      std::string line;

      getline(in, line);
      if (line.size() > 0) {
         sscanf(line.c_str(), "%lf %s", &stepSize_, comment);
      } else {
         UTIL_THROW("reading error: bead move step size (in units of bond length)");
      }
      #endif
   }

   /*
   * Write parameters, pairing with readParam.
   */
   void ChainFlipMove::writeParam(std::ostream& out) const
   {
      #if 0
      out << "step size in bond length    " << stepSize_ << std::endl;
      #endif
   }

   /*
   * Select one molecule, and flip chain head and tail.
   */
   void ChainFlipMove::move()
   {
      double         deltaE;
      bool           result(false);
      vector<Vector> rTrial;
      int            id, iBead;
      int            nAB = system_.nAB_;
      int            iChain = system_.random_.IRandom(0, system_.nPolymers_-1);
      Particle       *ptrChain = system_.polymer_[iChain];

      // Backup system status variables.
      //system_.backupStatusVariable();

      // Backup bead positions.
      rTrial.resize(nAB, Vector(0.0));
      for (id = 0; id < nAB; ++id)
         rTrial[id] = ptrChain[id].r;

      // Calculate chain flip move weight.
      deltaE = 0.0;
      for (id = 0; id < nAB; ++id) {
         iBead = nAB * iChain + id;
         deltaE += system_.calculateMoveWeight(iBead, rTrial[nAB - 1 - id]);
         system_.updateStatusVariable(iBead, rTrial[nAB - 1 - id]);
      }

      // Metropolis test.
      if (deltaE <= 0.0) {
         result = true;
      } else if (system_.random_.Random() <= exp(-deltaE)) { 
         result = true;
      } else {
         result = false;
      }

      ++nAttempt_;
      if (result) { 
         ++nAccept_;
      } else {
         for (id = 0; id < nAB; ++id) {
            //ptrChain[id].r = rTrial[id];
            iBead = nAB * iChain + id;
            system_.calculateMoveWeight(iBead, rTrial[id]);
            system_.updateStatusVariable(iBead, rTrial[id]);
         }
         //system_.restoreStatusVariable();
      }

      if (system_.doDOS())
         system_.updateDOS();
   }
}

#endif

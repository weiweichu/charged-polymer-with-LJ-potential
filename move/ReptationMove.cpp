#ifndef REPTATIONMOVE_CPP
#define REPTATIONMOVE_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <math.h>

#include "ReptationMove.h"
#include "../simulation/McSystem.h"

namespace GridMC
{
   using namespace Util;

   /*
   * Constructor.
   */
   ReptationMove::ReptationMove(McSystem& sysIn) :
      Move(sysIn),
      nStepMax_(0)
   { name_.assign("ReptationMove"); }

   /*
   * Default destructor.
   */
   ReptationMove::~ReptationMove()
   {}

   /*
   * Read step size.
   */
   void ReptationMove::readParam(std::istream& in)
   {
      char   comment[200];
      std::string line;

      getline(in, line);
      if (line.size() > 0) {
         sscanf(line.c_str(), "%d %s", &nStepMax_, comment);
      } else {
         UTIL_THROW("reading error: max number of moving steps in reptation move");
      }
   }

   /*
   * Write parameters, pairing with readParam.
   */
   void ReptationMove::writeParam(std::ostream& out) const
   {
      out << "max no. of reptation steps  " << nStepMax_ << std::endl;
   }

   /*
   * Move each bead in system randomly.
   */
   void ReptationMove::move()
   {
      double         deltaE;
      bool           result(false);
      int            id, iBead;
      int            nAB = system_.nAB_;
      int            iChain = system_.random_.IRandom(0, system_.nPolymers_-1);
      int            nSteps = system_.random_.IRandom(1, nStepMax_);
      int            direction = system_.random_.IRandom(0, 1);
      Particle       *ptrChain = system_.polymer_[iChain];
      vector<Vector> rTrial(nAB, Vector(0.0));
      vector<Vector> rBonds(nSteps, Vector(0.0));
      Vector         rNew;
      double         length;

      // Save bead positions.
      for (id = 0; id < nAB; ++id)
         rTrial[id] = ptrChain[id].r;

      // Generate nSteps bonds using the old bond length and new orientations.
      if (direction == 1) {
         for (id = nSteps - 1; id >= 0; --id) {
            rNew.subtract(ptrChain[id].r, ptrChain[id+1].r);
            system_.pbcShift(rNew);
            length = rNew.abs();
            unitVector(rNew);
            rNew *= length;
            rBonds[(nSteps - 1) - id] = rNew;
         }
      } else {
         for (id = nAB - nSteps; id < nAB; ++id) {
            rNew.subtract(ptrChain[id].r, ptrChain[id-1].r);
            system_.pbcShift(rNew);
            length = rNew.abs();
            unitVector(rNew);
            rNew *= length;
            rBonds[id - (nAB - nSteps)] = rNew;
         }
      }

      // Calculate reptation move weight.
      deltaE = 0.0;

      if (direction == 1) { // Forward reptatin.

         for (id = 0; id < nAB - nSteps; ++id) {
            iBead = nAB * iChain + id;
            deltaE += system_.calculateMoveWeight(iBead, rTrial[id + nSteps]);
            system_.updateStatusVariable(iBead, rTrial[id + nSteps]);
         }
         rNew = rTrial[nAB - 1];
         for (id = nAB - nSteps; id < nAB; ++id) {
            iBead = nAB * iChain + id;
            rNew += rBonds[id - (nAB - nSteps)];
            system_.toPrimaryCell(rNew);
            deltaE += system_.calculateMoveWeight(iBead, rNew);
            system_.updateStatusVariable(iBead, rNew);
         }

      } else { // Backward reptation.

         for (id = nSteps; id < nAB; ++id) {
            iBead = nAB * iChain + id;
            deltaE += system_.calculateMoveWeight(iBead, rTrial[id - nSteps]);
            system_.updateStatusVariable(iBead, rTrial[id - nSteps]);
         }
         rNew = rTrial[0];
         for (id = nSteps - 1; id >= 0; --id) {
            iBead = nAB * iChain + id;
            rNew += rBonds[nSteps - 1 - id];
            system_.toPrimaryCell(rNew);
            deltaE += system_.calculateMoveWeight(iBead, rNew);
            system_.updateStatusVariable(iBead, rNew);
         }

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
            iBead = nAB * iChain + id;
            system_.calculateMoveWeight(iBead, rTrial[id]);
            system_.updateStatusVariable(iBead, rTrial[id]);
         }
      }

      if (system_.doDOS())
         system_.updateDOS();
   }


   /*
   * Generate a uniformly distributed unit vector.
   */
   void ReptationMove::unitVector(Vector& v) const
   {
      double ran1=0, ran2=0, ranh;
      double ransq;

      ransq=2.0;
      while (ransq >= 1.0) {
         ran1  = 1.0 - 2.0*system_.random_.Random();
         ran2  = 1.0 - 2.0*system_.random_.Random();
         ransq = ran1*ran1 + ran2*ran2;
      }
      ranh= 2.0*sqrt(1.0-ransq);
      v[0] = ran1*ranh;
      v[1] = ran2*ranh;
      v[2] = 1.0 - 2.0*ransq;
   }

}

#endif

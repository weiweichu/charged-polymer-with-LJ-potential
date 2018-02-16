#ifndef BEADMOVE_CPP
#define BEADMOVE_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <math.h>
#include "../util/global.h"
#include "../util/Vector.h"

#include "BeadMove.h"
#include "../simulation/McSystem.h"

namespace GridMC
{
   using namespace Util;

   /*
   * Constructor.
   */
   BeadMove::BeadMove(McSystem& sysIn) :
      Move(sysIn),
      stepSize_(1.0) 
   { name_.assign("BeadMove"); }

   /*
   * Default destructor.
   */
   BeadMove::~BeadMove()
   {}

   /*
   * Read step size.
   */
   void BeadMove::readParam(std::istream& in)
   {
      char   comment[200];
      std::string line;

      getline(in, line);
      if (line.size() > 0) {
         sscanf(line.c_str(), "%lf %s", &stepSize_, comment);
      } else {
         UTIL_THROW("reading error: bead move step size (in units of bond length)");
      }
   }

   /*
   * Write parameters, pairing with readParam.
   */
   void BeadMove::writeParam(std::ostream& out) const
   {
      out << "step size in bond length    " << stepSize_ << std::endl;
   }

   /*
   * Move each bead in system randomly.
   */
   void BeadMove::move()
   {
      double    deltaE, randomNum, maxStep(stepSize_*system_.bondL_);
      int       id, j;
      bool      result(false);
      Vector    rTrial;

      for (unsigned i = 0; i < system_.beads_.size(); ++i) {

         // Generate random bead id.
         do {
            randomNum = system_.random_.Random();
         } while (randomNum == 0);
         id = ceil(randomNum * double(system_.beads_.size())) - 1;    

         // Generate trial bead position.
         rTrial = system_.beads_[id].r;
         for (j = 0; j < Dimension; ++j)
            rTrial[j] += (2.0*system_.random_.Random()-1.0)*maxStep;
         system_.toPrimaryCell(rTrial);

         // Get the energy change.
         deltaE = system_.calculateMoveWeight(id, rTrial);

         // Metropolis test.
         if (deltaE <= 0.0) {
            result = true;
         } else if (system_.random_.Random() <= exp(-deltaE)) { 
            result = true;
         } else {
            //Log::file() << "deltaE " << deltaE << "   " <<exp(-deltaE) << endl;
            result = false;
         }

         ++nAttempt_;
         if (result) { 
           //Log::file() << "success bead move" << endl;
            system_.updateStatusVariable(id, rTrial);
            ++nAccept_;
         }
         if (system_.doDOS())
            system_.updateDOS();
      }
   }

}

#endif

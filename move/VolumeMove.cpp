#ifndef VOLUMEMOVE_CPP
#define VOLUMEMOVE_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <math.h>
#include "../util/global.h"
#include "../util/Vector.h"
#include "../simulation/McSystem.h"
#include "VolumeMove.h"

namespace GridMC
{
   using namespace Util;

   /*
   * Constructor.
   */
   VolumeMove::VolumeMove(McSystem& sysIn) :
      Move(sysIn),
      maxDeltaV_(0.0)
   { name_.assign("VolumeMove"); }

   /*
   * Default destructor.
   */
   VolumeMove::~VolumeMove()
   {}

   /*
   * Read step size.
   */
   void VolumeMove::readParam(std::istream& in)
   {
      char   comment[200];
      std::string line;

      // Read magnitude of volume change.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: maximum volume move step");
      sscanf(line.c_str(), "%lf %s", &maxDeltaV_, comment);
   }

   /*
   * Write parameters, pairing with readParam.
   */
   void VolumeMove::writeParam(std::ostream& out) const
   {
      out << "max volume move step size   " << maxDeltaV_ << endl;
   }

   /*
   * Change system volume within a given range.
   */
   void VolumeMove::move()
   {
      double  deltaV, deltaE;
      bool    result(false);

      deltaV = (2.0*system_.random_.Random() - 1.0) * maxDeltaV_;
      deltaE = system_.volumeChangeWeight(deltaV);

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
         system_.updateVolume(deltaV);
         ++nAccept_;
      }
      if (system_.doDOS())
         system_.updateDOS();
   }

}

#endif

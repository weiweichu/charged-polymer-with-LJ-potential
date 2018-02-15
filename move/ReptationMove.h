#ifndef REPTATIONMOVE_H
#define REPTATIONMOVE_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <string>
#include "../util/global.h"
#include "../util/Vector.h"
#include "Move.h"

namespace GridMC
{
   using namespace Util;

   /**
   * Reptation move.
   *
   * \ingroup Move_Module
   */
   class ReptationMove : public Move
   {
   
   public:
  
      /**
      * Default constructor.
      */
      ReptationMove(McSystem& sysIn);

      /**
      * Default destructor.
      */
      ~ReptationMove();

      /**
      * Read parameters.
      */
      virtual void readParam(std::istream& in);

      /**
      * Write parameters, pairing with readParam.
      */
      virtual void writeParam(std::ostream& out) const;

      /**
      * Apply the move.
      */
      virtual void move();

   protected:

      /// Number of move steps.
      int  nStepMax_;

      /**
      * Generate a uniformly distributed unit vector.
      */
      void unitVector(Vector& v) const;

   };

}
#endif

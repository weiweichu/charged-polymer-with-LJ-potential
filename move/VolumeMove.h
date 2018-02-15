#ifndef VOLUMEMOVE_H
#define VOLUMEMOVE_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <string>

#include "Move.h"

namespace GridMC
{

   /**
   * Volumeric change under a given barostat pressure.
   *
   * \ingroup Move_Module
   */
   class VolumeMove : public Move
   {
   
   public:
  
      /**
      * Default constructor.
      */
      VolumeMove(McSystem& sysIn);

      /**
      * Default destructor.
      */
      ~VolumeMove();

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

      /// Maximum volume move step.
      double maxDeltaV_;

   };

}
#endif

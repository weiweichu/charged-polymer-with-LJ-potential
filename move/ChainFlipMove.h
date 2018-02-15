#ifndef CHAINFLIPMOVE_H
#define CHAINFLIPMOVE_H

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
   * Head-tail flipping move.
   *
   * \ingroup Move_Module
   */
   class ChainFlipMove : public Move
   {
   
   public:
  
      /**
      * Default constructor.
      */
      ChainFlipMove(McSystem& sysIn);

      /**
      * Default destructor.
      */
      ~ChainFlipMove();

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

   };

}
#endif

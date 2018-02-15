#ifndef SEMIGRANDSALTMOVE_H
#define SEMIGRANDSALTMOVE_H

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
   * Semi-grand move: salts and solvent exchange identities.
   *
   * \ingroup Move_Module
   */
   class SemigrandSaltMove : public Move
   {
   
   public:
  
      /**
      * Default constructor.
      */
      SemigrandSaltMove(McSystem& sysIn);

      /**
      * Default destructor.
      */
      ~SemigrandSaltMove();

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

      /// Magnitude of move, in unit of bond length.
      double  dMu_;

   };

}
#endif

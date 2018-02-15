#ifndef PARTICLE_H
#define PARTICLE_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "../util/Vector.h"

namespace GridMC
{
   using namespace Util;

   /**
   * A point particle.
   *
   * \ingroup Simulation_Module
   */
   class Particle 
   {
   
   public:

      /**
      * Constructor.
      */
      Particle();

      /**
      * Default destructor.
      */
      ~Particle();

      /// Coordinates.
      Vector  r;
      int     t;
      double  q;

   };

}
#endif

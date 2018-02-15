#ifndef PARTICLE_CPP
#define PARTICLE_CPP

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "Particle.h"

namespace GridMC
{
   /*
   * Default constructor.
   */
   Particle::Particle() :
      r(0.0),
      t(0),
      q(0.0)
   {}

   /*
   * Default destructor.
   */
   Particle::~Particle()
   {}

}
#endif

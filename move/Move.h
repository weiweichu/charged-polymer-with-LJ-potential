#ifndef MOVE_H
#define MOVE_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <string>

namespace GridMC
{

   class  McSystem;

   /**
   * Base class for Monte Carlo type moves applied to system.
   *
   * \ingroup Move_Module
   */
   class Move
   {
   
   public:
  
      /**
      * Default constructor.
      */
      Move(McSystem& sysIn);

      /**
      * Default destructor.
      */
      virtual ~Move();

      /**
      * Read parameters.
      */
      virtual void readParam(std::istream& in) = 0;

      /**
      * Write parameters, pairing with readParam.
      */
      virtual void writeParam(std::ostream& out) const = 0;

      /**
      * Apply the move.
      */
      virtual void move() = 0;

      /**
      * Output results and statistics.
      */
      virtual void outputStatistics(std::ostream& out) const;

   protected:

      /// Interval.
      int nInterval_;

      /// Number attempted.
      int nAttempt_;

      /// Number accepted.
      int nAccept_;

      /// Reference to system.
      McSystem& system_;

      /// Name.
      std::string name_;

   };

}
#endif

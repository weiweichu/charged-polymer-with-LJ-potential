#ifndef DIAGNOSIS_H
#define DIAGNOSIS_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <vector>
#include <util/global.h>
#include <util/FileMaster.h>
#include <grid/GridMassCharge.h>

#include <simulation/System.h>

namespace GridMC
{

   /**
   * This class provide a generic interface for system diagnosis.
   *
   * \ingroup Diag_Module
   */
   class Diagnosis
   {
   
   public:
  
      /**
      * Default constructor.
      */
      Diagnosis(System& sysIn);

      /**
      * Default destructor.
      */
      virtual ~Diagnosis();

      /**
      * Read parameter.
      */
      virtual void readParam(std::istream& in);

      /**
      * Write parameter.
      */
      virtual void writeParam(std::ostream& out) const;

      /**
      * Return true if is at interval.
      */
      virtual bool isAtInterval(const long iStep);

      /**
      * Calculate order parameter.
      */
      virtual void sample(const long iStep) = 0;

      /**
      * Output result parameter.
      */
      virtual void output() const = 0;

   protected:

      /// Reference to system.
      System& system_;

      /// Number of bead types in the system.
      const Vector& boxL_;

      /// Reference to parent system.
      const IntVector& nGrid_;

      /// Reference to parent system.
      GridMassCharge& grid_;

      /// Reference to the beads array.
      const vector<Particle>& beads_;

      /// Reference to the file master.
      const FileMaster& fileMaster_;

      /// Sampling frequency.
      int   nInterval_;

      /// Number of samples.
      long  nSamples_;

      /// Name.
      std::string name_;

   };

   inline bool Diagnosis::isAtInterval(const long iStep)
   { return (iStep%nInterval_ == 0 ? true : false); }

}
#endif

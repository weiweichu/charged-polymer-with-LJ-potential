#ifndef IONNUMBERFRACTION_H
#define IONNUMBERFRACTION_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "Diagnosis.h"

namespace GridMC
{

   /**
   * Sample the number fraction of ions.
   *
   * \ingroup Diag_Module
   */
   class IonNumberFraction: public Diagnosis
   {
   
   public:
  
      /**
      * Default constructor.
      */
      IonNumberFraction(System& sysIn);

      /**
      * Default destructor.
      */
      virtual ~IonNumberFraction();

      /**
      * Read parameter.
      */
      virtual void readParam(std::istream& in);

      /**
      * Write parameter.
      */
      virtual void writeParam(std::ostream& out) const;

      /**
      * Calculate order parameter.
      */
      virtual void sample(const long iStep);

      /**
      * Write parameter.
      */
      virtual void output() const;

   private:

      /// Accumulated first and second moments of ion number fraction.
      double psi1_, psi2_;

      /// Data file.
      std::ofstream dataFile_;

   };

}
#endif

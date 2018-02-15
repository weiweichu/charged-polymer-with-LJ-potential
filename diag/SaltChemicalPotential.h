#ifndef SALTCHEMICALPOTENTIAL_H
#define SALTCHEMICALPOTENTIAL_H

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
   * Evaluate the chemical potential of a neutral solvent.
   *
   * \ingroup Diag_Module
   */
   class SaltChemicalPotential: public Diagnosis
   {
   
   public:
  
      /**
      * Default constructor.
      */
      SaltChemicalPotential(System& sysIn);

      /**
      * Default destructor.
      */
      virtual ~SaltChemicalPotential();

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

      /// Type of species.
      int  type_;

      /// First and second moments of fugacity.
      double f1_, f2_;

      /// Point to data file.
      std::ofstream dataFile_;

   };

}
#endif

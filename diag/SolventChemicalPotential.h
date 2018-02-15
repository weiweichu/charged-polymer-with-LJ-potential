#ifndef SOLVENTCHEMICALPOTENTIAL_H
#define SOLVENTCHEMICALPOTENTIAL_H

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
   class SolventChemicalPotential: public Diagnosis
   {
   
   public:
  
      /**
      * Default constructor.
      */
      SolventChemicalPotential(System& sysIn);

      /**
      * Default destructor.
      */
      virtual ~SolventChemicalPotential();

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

      /// Chemical potential.
      double mu1_, mu2_;

      /// Chemical potential.
      std::ofstream dataFile_;

   };

}
#endif

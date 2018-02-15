#ifndef CHAINCHEMICALPOTENTIAL_H
#define CHAINCHEMICALPOTENTIAL_H

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
   * Evaluate the chain chemical potential by deleting monomers at the end.
   *
   * \ingroup Diag_Module
   */
   class ChainChemicalPotential: public Diagnosis
   {
   
   public:
  
      /**
      * Default constructor.
      */
      ChainChemicalPotential(System& sysIn);

      /**
      * Default destructor.
      */
      virtual ~ChainChemicalPotential();

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

      /// First and second moments of fugacity.
      double f1_, f2_;

      /// Point to data file.
      std::ofstream dataFile_;

   };

}
#endif

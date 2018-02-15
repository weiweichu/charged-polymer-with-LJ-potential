#ifndef STRUCTUREFACTOR_H
#define STRUCTUREFACTOR_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <set>
#include <vector>
#include <util/IntVector.h>
#include <util/Vector.h>

#include "Diagnosis.h"

namespace GridMC
{

   typedef vector<double>     ValueV;
   typedef vector<Vector>     WaveV;
   typedef vector<IntVector>  WaveIndexV;

   /**
   * This class calculate the structure factor.
   *
   * \ingroup Diag_Module
   */
   class StructureFactor: public Diagnosis
   {
   
   public:
  
      /**
      * Default constructor.
      */
      StructureFactor(System& sysIn);

      /**
      * Default destructor.
      */
      virtual ~StructureFactor();

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

      /// Equlibration steps.
      int                 nEquilibration_;

      /// Type index.
      int                 tA_, tB_;

      /// maxIndex.
      int                 maxIndex_;

      /// Number of distinct waves.
      int                 nWaves_;

      /// List of wave indices squared.
      set<int>            waveIndexSq_;

      /// List of wave number.
      ValueV              waveNumber_;

      /// List of wave indices.
      vector<WaveIndexV>  waveIndex_;

      /// Wave vectors; one-to-one correspondence to waveIndex_.
      vector<WaveV>       waves_;

      /// Cosine and sine components of fourier modes.
      vector<ValueV>      cA_, sA_;
      vector<ValueV>      cB_, sB_;

      /// Correlation function.
      ValueV              S_;

      /**
      * Generate quivalent waves using wave indices.
      */
      void generateWaves();

   };

}
#endif

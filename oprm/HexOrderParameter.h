#ifndef HEXORDERPARAMETER_H
#define HEXORDERPARAMETER_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <vector>
#include <util/IntVector.h>

#include "OrderParameter.h"

namespace GridMC
{

   /**
   * This class represent the order parameter for a modulation wave of hexagonal symmetry.
   *
   * \ingroup Oprm_Module
   */
   class HexOrderParameter : public OrderParameter
   {
   
   public:
  
      /**
      * Default constructor.
      */
      HexOrderParameter(const System& sysIn);

      /**
      * Default destructor.
      */
      virtual ~HexOrderParameter();

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
      virtual void calculate();

      /**
      * Calculate orde parameter change due to bead move.
      */
      virtual double getChange(const Particle& bead, const Vector& rTrial);

      /**
      * Update order parameter by dpsi, or dpcos and dpsin.
      */
      virtual void updateByChange();

      /**
      * Reset parameters after box dimension is rescaled.
      */
      virtual void updateAfterBoxMove();

   private:

      /// Type flag for order parameter.
      int       oprmFlag_;

      /// Number of waves.
      unsigned          nWaves_;

      /// Wave vector index.
      vector<IntVector> waveIndex_;

      /// Wave vector.
      vector<Vector>    waveVector_;

      /// Parameter needed by sine and cosine modes.
      vector<double>    pcos_, psin_, dpcos_, dpsin_;

      /**
      * Calculate orde parameter using both cosine and sine mode.
      */
      void oprmFromWave();

      /**
      * Calculate orde parameter change using both cosine and sine wave mode.
      */
      void oprmChangeFromWave(const Particle& bead, const Vector& rTrial);

   };

   /*
   * Accessor to the order parameter.
   */
   inline void HexOrderParameter::updateByChange()
   {
      psi_ += dpsi_;
      if (oprmFlag_ > 1) {
         for (unsigned j = 0; j < nWaves_; ++j) {
            pcos_[j] += dpcos_[j];
            psin_[j] += dpsin_[j];
         }
      }
   }

}
#endif

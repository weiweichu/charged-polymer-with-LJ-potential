#ifndef LAMORDERPARAMETER_H
#define LAMORDERPARAMETER_H

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
   * This class represent the order parameter for a 1D modulation wave.
   *
   * \ingroup Oprm_Module
   */
   class LamOrderParameter : public OrderParameter
   {
   
   public:
  
      /**
      * Default constructor.
      */
      LamOrderParameter(const System& sysIn);

      /**
      * Default destructor.
      */
      virtual ~LamOrderParameter();

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

      /**
      * Backup status variables.
      */
      virtual void backup();

      /**
      * Restore status variables.
      */
      virtual void restore();



   private:

      /// Type flag for order parameter.
      int       oprmFlag_;

      /// Wave vector index.
      IntVector waveIndex_;

      /// Wave vector.
      Vector    waveVector_;

      /// Parameter needed by sine and cosine modes.
      double pcos_, psin_, dpcos_, dpsin_;

      /// Wave vector index.
      vector<IntVector> tripleWaveIndex_;

      /// Wave vector.
      vector<Vector>    tripleWaveVector_;

      /// Parameter needed by sine and cosine modes.
      vector<double>    triplepcos_, triplepsin_, dtriplepcos_, dtriplepsin_;

      /// Wave vector index.
      vector<IntVector> waveArrayIndex_;

      /// Wave vector.
      vector<Vector>    waveArrayVector_;

      /// Cosine mode.
      vector<double>    waveArrayPcos_, waveArrayDPcos_;

      /// Sine mode.
      vector<double>    waveArrayPsin_, waveArrayDPsin_;
 
      /// Bias weight.
      vector<double>    waveArrayWeight_;

      /// Number of waves.
      unsigned nWaves_;

      /// Max wave index number.
      int maxIndexSq_;

      /// Peak wave vector.
      double qstar_, qstarSq_;

      /**
      * Calculate orde parameter using cosine wave mode.
      */
      void oprmFromCosineWave();

      /**
      * Calculate orde parameter change using cosine wave mode.
      */
      void oprmChangeFromCosineWave(const Particle& bead, const Vector& rTrial);

      /**
      * Calculate orde parameter using both cosine and sine mode.
      */
      void oprmFromWave();

      /**
      * Calculate orde parameter change using both cosine and sine wave mode.
      */
      void oprmChangeFromWave(const Particle& bead, const Vector& rTrial);

      /**
      * Calculate orde parameter using <001> wave family.
      */
      void oprmFromTripleWave();

      /**
      * Calculate orde parameter change using <001> wave family.
      */
      void oprmChangeFromTripleWave(const Particle& bead, const Vector& rTrial);

      /**
      * Calculate orde parameter using biased wave array.
      */
      void oprmFromWaveArray();

      /**
      * Calculate orde parameter change using biased wave array.
      */
      void oprmChangeFromWaveArray(const Particle& bead, const Vector& rTrial);

      /**
      * Generate waves according to the maximum wave index.
      */
      void generateWaves();

   };

   /*
   * Accessor to the order parameter.
   */
   inline void LamOrderParameter::updateByChange()
   {
      psi_ += dpsi_;
      if (oprmFlag_ == 2) {
         pcos_ += dpcos_;
         psin_ += dpsin_;
      } else if (oprmFlag_ == 3) {
         for (unsigned j = 0; j < 3; ++j) {
            triplepcos_[j] += dtriplepcos_[j];
            triplepsin_[j] += dtriplepsin_[j];
         }
      } else if (oprmFlag_ < 0) {
         for (unsigned j = 0; j < nWaves_; ++j) {
            waveArrayPcos_[j] += waveArrayDPcos_[j];
            waveArrayPsin_[j] += waveArrayDPsin_[j];
         }
      }
   }

   /*
   * Backup status variables.
   */
   inline void LamOrderParameter::backup()
   {
      OrderParameter::backup();
   }

   /*
   * Restore status variables.
   */
   inline void LamOrderParameter::restore()
   {
      OrderParameter::restore();
   }

}
#endif

#ifndef LAMORDERPARAMETER_CPP
#define LAMORDERPARAMETER_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <math.h>
#include <util/Constants.h>
#include "LamOrderParameter.h"

namespace GridMC
{

   /*
   * Constructor.
   */
   LamOrderParameter::LamOrderParameter(const System& system) :
      OrderParameter(system),
      oprmFlag_(0),
      waveIndex_(0),
      waveVector_(0.0),
      pcos_(0.0),
      psin_(0.0),
      dpcos_(0.0),
      dpsin_(0.0),
      // triple wave
      tripleWaveIndex_(3, IntVector(0)),
      tripleWaveVector_(3, Vector(0.0)),
      triplepcos_(3, 0.0),
      triplepsin_(3, 0.0),
      dtriplepcos_(3, 0.0),
      dtriplepsin_(3, 0.0),
      // wave array
      waveArrayIndex_(),
      waveArrayVector_(),
      waveArrayPcos_(),
      waveArrayDPcos_(),
      waveArrayPsin_(),
      waveArrayDPsin_(),
      waveArrayWeight_(),
      nWaves_(0),
      maxIndexSq_(0),
      qstar_(0.0)
   {}

   /*
   * Default destructor.
   */
   LamOrderParameter::~LamOrderParameter()
   {}

   /*
   * Read parameter.
   */
   void LamOrderParameter::readParam(std::istream& in)
   {
      char        comment[200];
      std::string line;

      // Order parameter type flag.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: order parameter flag");
      sscanf(line.c_str(), "%d %s", &oprmFlag_, comment);

      // Read wave vector if needed.
      if (oprmFlag_ == 1 || oprmFlag_ == 2) {

         getline(in, line);
         if (line.size() <= 0)
            UTIL_THROW("reading error: wave index for lamellar order parameter");
         sscanf(line.c_str(), "%d %d %d %s", &waveIndex_[0], &waveIndex_[1], &waveIndex_[2], comment);

         // Calculate the wave vector.
         for (int i = 0; i < Dimension; ++i) 
            waveVector_[i] = 2.0 * Constants::Pi / boxL_[i] * double(waveIndex_[i]);

      } else if (oprmFlag_ == 3) {

         // Specify wave indices.
         tripleWaveIndex_[0] = IntVector(1, 0, 0);
         tripleWaveIndex_[1] = IntVector(0, 1, 0);
         tripleWaveIndex_[2] = IntVector(0, 0, 1);

         // Calculate the wave vector.
         for (int j = 0; j < 3; ++j)
            for (int i = 0; i < Dimension; ++i) 
               tripleWaveVector_[j][i] = 2.0 * Constants::Pi / boxL_[i] * double(tripleWaveIndex_[j][i]);

      } else if (oprmFlag_ < 0) {

         // Read qstar.
         getline(in, line);
         if (line.size() <= 0)
            UTIL_THROW("reading error: q* for lamellar order parameter");
         sscanf(line.c_str(), "%lf %s", &qstar_, comment);
         qstarSq_ = qstar_ * qstar_;

         // Use -oprmFlag_ as maximum wave index.
         maxIndexSq_ = -oprmFlag_;
         generateWaves();

      } else if (oprmFlag_ != 0) {

         UTIL_THROW("Unrecognized oprmFlag");

      }

      // Initialize number of bead.
      nBead_ = beads_.size();
   }

   /*
   * Write parameter.
   */
   void LamOrderParameter::writeParam(std::ostream& out) const
   {
      out << "Order parameter type flag   " << oprmFlag_ << std::endl;
      if (oprmFlag_ == 1 || oprmFlag_ == 2) {
         out << "wave index for lam oprm     ";
         out << waveIndex_[0] << "  ";
         out << waveIndex_[1] << "  ";
         out << waveIndex_[2] << std::endl;
      } else if (oprmFlag_ < 0) {
         out << "number of waves             " << nWaves_ << std::endl;
         for (unsigned i = 0; i < waveArrayIndex_.size(); ++i) {
            if (i == 0)
               out << "wave index for lam oprm     ";
            else
               out << "                            ";
            out << waveArrayIndex_[i][0] << "  ";
            out << waveArrayIndex_[i][1] << "  ";
            out << waveArrayIndex_[i][2] << std::endl;
         }
      }
   }

   /*
   * Calculate order parameter.
   */
   void LamOrderParameter::calculate()
   {
      if (oprmFlag_ == 0)
         oprmFromGrid();
      else if (oprmFlag_ == 1)
         oprmFromCosineWave();
      else if (oprmFlag_ == 2)
         oprmFromWave();
      else if (oprmFlag_ == 3)
         oprmFromTripleWave();
      else
         oprmFromWaveArray();
   }

   /*
   * Calculate order parameter change dure to bead move.
   */
   double LamOrderParameter::getChange(const Particle& bead, const Vector& rTrial)
   {
      if (oprmFlag_ == 0)
         oprmChangeFromGrid(bead.t);
      else if (oprmFlag_ == 1)
         oprmChangeFromCosineWave(bead, rTrial);
      else if (oprmFlag_ == 2)
         oprmChangeFromWave(bead, rTrial);
      else if (oprmFlag_ == 3)
         oprmChangeFromTripleWave(bead, rTrial);
      else
         oprmChangeFromWaveArray(bead, rTrial);
      return dpsi_;
   }

   /*
   * Calculate order parameter.
   */
   void LamOrderParameter::oprmFromWave()
   {
      double rdotq;

      pcos_ = 0.0;
      psin_ = 0.0;

      for (unsigned i = 0; i < nBead_; ++i) {
         rdotq = beads_[i].r.dot(waveVector_);
         if (beads_[i].t == 0) {
            pcos_ += cos(rdotq);
            psin_ += sin(rdotq);
         } else if (beads_[i].t == 1) {
            pcos_ -= cos(rdotq);
            psin_ -= sin(rdotq);
         }
      }

      pcos_ /= double(nBead_);
      psin_ /= double(nBead_);

      psi_ = sqrt(pcos_*pcos_ + psin_*psin_);
   }

   /*
   * Calculate order parameter change dure to bead move.
   */
   void LamOrderParameter::oprmChangeFromWave(const Particle& bead, const Vector& rTrial)
   {
      if (bead.t > 1) {
         dpcos_ = 0.0;
         dpsin_ = 0.0;
         dpsi_  = 0.0;
      } else {
         double dot1(bead.r.dot(waveVector_));
         double dot2(rTrial.dot(waveVector_));

         if (bead.t == 0) {
            dpcos_ = cos(dot2) - cos(dot1);
            dpsin_ = sin(dot2) - sin(dot1);
         } else {
            dpcos_ = cos(dot1) - cos(dot2);
            dpsin_ = sin(dot1) - sin(dot2);
         }
         dpcos_ /= double(nBead_);
         dpsin_ /= double(nBead_);

         dpsi_  = sqrt((pcos_+dpcos_)*(pcos_+dpcos_) + (psin_+dpsin_)*(psin_+dpsin_));
         dpsi_ -= psi_;
      }
   }

   /*
   * Calculate order parameter using <001> wave family.
   */
   void LamOrderParameter::oprmFromTripleWave()
   {
      unsigned i, j;
      double rdotq;

      for (j = 0; j < 3; ++j) {
         triplepcos_[j] = 0.0;
         triplepsin_[j] = 0.0;
      }

      for (i = 0; i < nBead_; ++i) {
         for (j = 0; j < 3; ++j) {
            rdotq = beads_[i].r.dot(tripleWaveVector_[j]);
            if (beads_[i].t == 0) {
               triplepcos_[j] += cos(rdotq);
               triplepsin_[j] += sin(rdotq);
            } else if (beads_[i].t == 1) {
               triplepcos_[j] -= cos(rdotq);
               triplepsin_[j] -= sin(rdotq);
            }
         }
      }

      for (j = 0; j < 3; ++j) {
         triplepcos_[j] /= double(nBead_);
         triplepsin_[j] /= double(nBead_);
      }

      psi_ = 0.0;
      for (j = 0; j < 3; ++j)
         psi_ += sqrt(triplepcos_[j]*triplepcos_[j] + triplepsin_[j]*triplepsin_[j]);
      psi_ /= sqrt(3.0);

   }

   /*
   * Calculate order parameter change dure to bead move.
   */
   void LamOrderParameter::oprmChangeFromTripleWave(const Particle& bead, const Vector& rTrial)
   {
      unsigned j;
      if (bead.t > 1) {
         for (j = 0; j < 3; ++j) {
            dtriplepcos_[j] = 0.0;
            dtriplepsin_[j] = 0.0;
         }
         dpsi_  = 0.0;
      } else {
         double dot1, dot2;

         if (bead.t == 0) {

            for (j = 0; j < 3; ++j) {
               dot1 = bead.r.dot(tripleWaveVector_[j]);
               dot2 = rTrial.dot(tripleWaveVector_[j]);

               dtriplepcos_[j] = (cos(dot2) - cos(dot1)) / double(nBead_);
               dtriplepsin_[j] = (sin(dot2) - sin(dot1)) / double(nBead_);
            }

         } else {

            for (j = 0; j < 3; ++j) {
               dot1 = bead.r.dot(tripleWaveVector_[j]);
               dot2 = rTrial.dot(tripleWaveVector_[j]);

               dtriplepcos_[j] = (cos(dot1) - cos(dot2)) / double(nBead_);
               dtriplepsin_[j] = (sin(dot1) - sin(dot2)) / double(nBead_);
            }

         }

         dpsi_ = 0.0;
         for (j = 0; j < 3; ++j)
            dpsi_ += sqrt( (triplepcos_[j] + dtriplepcos_[j]) * (triplepcos_[j] + dtriplepcos_[j]) + 
                           (triplepsin_[j] + dtriplepsin_[j]) * (triplepsin_[j] + dtriplepsin_[j]) );
         dpsi_ /= sqrt(3.0);
         dpsi_ -= psi_;
      }
      //Log::file() << "dpsi = " << dpsi_ << std::endl;
   }


   /*
   * Calculate order parameter using prescribed wave array.
   */
   void LamOrderParameter::oprmFromWaveArray()
   {
      unsigned i, j;
      double rdotq;

      for (j = 0; j < nWaves_; ++j) {
         waveArrayPcos_[j] = 0.0;
         waveArrayPsin_[j] = 0.0;
      }

      for (i = 0; i < nBead_; ++i) {
         for (j = 0; j < nWaves_; ++j) {
            rdotq = beads_[i].r.dot(waveArrayVector_[j]);
            if (beads_[i].t == 0) {
               waveArrayPcos_[j] += cos(rdotq);
               waveArrayPsin_[j] += sin(rdotq);
            } else if (beads_[i].t == 1) {
               waveArrayPcos_[j] -= cos(rdotq);
               waveArrayPsin_[j] -= sin(rdotq);
            }
         }
      }

      for (j = 0; j < nWaves_; ++j) {
         waveArrayPcos_[j] /= double(nBead_);
         waveArrayPsin_[j] /= double(nBead_);
      }

      psi_ = 0.0;
      for (j = 0; j < nWaves_; ++j)
         psi_ += pow(pow(waveArrayPcos_[j],2) + pow(waveArrayPsin_[j],2), 2) * waveArrayWeight_[j];
      psi_ = pow(psi_, 1.0/4.0);

   }

   /*
   * Calculate order parameter change dure to bead move using prescribed wave array.
   */
   void LamOrderParameter::oprmChangeFromWaveArray(const Particle& bead, const Vector& rTrial)
   {
      unsigned j;
      if (bead.t > 1) {
         for (j = 0; j < nWaves_; ++j) {
            waveArrayDPcos_[j] = 0.0;
            waveArrayDPsin_[j] = 0.0;
         }
         dpsi_  = 0.0;
      } else {
         double dot1, dot2;

         if (bead.t == 0) {

            for (j = 0; j < nWaves_; ++j) {
               dot1 = bead.r.dot(waveArrayVector_[j]);
               dot2 = rTrial.dot(waveArrayVector_[j]);

               waveArrayDPcos_[j] = (cos(dot2) - cos(dot1)) / double(nBead_);
               waveArrayDPsin_[j] = (sin(dot2) - sin(dot1)) / double(nBead_);
            }

         } else if (bead.t == 1) {

            for (j = 0; j < nWaves_; ++j) {
               dot1 = bead.r.dot(waveArrayVector_[j]);
               dot2 = rTrial.dot(waveArrayVector_[j]);

               waveArrayDPcos_[j] = (cos(dot1) - cos(dot2)) / double(nBead_);
               waveArrayDPsin_[j] = (sin(dot1) - sin(dot2)) / double(nBead_);
            }

         }

         dpsi_ = 0.0;
         for (j = 0; j < nWaves_; ++j)
            dpsi_ += pow( pow(waveArrayPcos_[j] + waveArrayDPcos_[j], 2) + 
                         +pow(waveArrayPsin_[j] + waveArrayDPsin_[j], 2), 2)
                     * waveArrayWeight_[j];
         dpsi_  = pow(dpsi_, 1.0/4.0);
         dpsi_ -= psi_;
      }
   }


   /*
   * Calculate order parameter using cosine mode.
   */
   void LamOrderParameter::oprmFromCosineWave()
   {
      double rdotq;

      psi_ = 0.0;
      for (unsigned i = 0; i < nBead_; ++i) {
         rdotq = beads_[i].r.dot(waveVector_);
         if (beads_[i].t == 0) {
            psi_ += cos(rdotq);
         } else if (beads_[i].t == 1) {
            psi_ -= cos(rdotq);
         }
      }
      psi_ /= double(nBead_);
   }

   /*
   * Calculate order parameter change due to bead move using consine mode.
   */
   void LamOrderParameter::oprmChangeFromCosineWave(const Particle& bead, const Vector& rTrial)
   {
      if (bead.t > 1) {
         dpsi_ = 0.0;
      } else {
         double dot1(bead.r.dot(waveVector_));
         double dot2(rTrial.dot(waveVector_));

         if (bead.t == 0) {
            dpsi_ = cos(dot2) - cos(dot1);
         } else {
            dpsi_ = cos(dot1) - cos(dot2);
         }
         dpsi_ /= double(nBead_);
      }
   }

   /*
   * Reset the wave vector after the box dimension change.
   */
   void LamOrderParameter::updateAfterBoxMove()
   {
      if (oprmFlag_ != 0)
         for (int i = 0; i < Dimension; ++i) 
            waveVector_[i] = 2.0 * Constants::Pi / boxL_[i] * double(waveIndex_[i]);
   }

   /*
   * Generate all equivalent waves from the wave indices.
   *
   * Assuming a cubic simulation box.
   */
   void LamOrderParameter::generateWaves()
   {
      // Generate the set of wave index magnitudes.
      if (fabs(boxL_[0] - boxL_[1]) + fabs(boxL_[0] - boxL_[2]) > Constants::Epsilon)
         UTIL_THROW("simulation box is cubic");

      // Generate the set of wave index magnitudes.
      int maxIndex, ix, iy, iz, sqsum;
      double basis = 2.0 * Constants::Pi / boxL_[0], weight, wNormal(0.0);
      Vector wave;

      maxIndex = ceil(sqrt(double(maxIndexSq_)));
      for (ix = -maxIndex; ix <= maxIndex; ++ix) { 
         for (iy = -maxIndex; iy <= maxIndex; ++iy) { 
            for (iz = 0; iz <= maxIndex; ++iz) { 
               if (iz > 0 || (iz == 0 && iy > 0) || (iz == 0 && iy == 0 && ix > 0)) {
                  sqsum = ix * ix + iy * iy + iz * iz;
                  if (sqsum <= maxIndexSq_) {
                     waveArrayIndex_.push_back(IntVector(ix, iy,iz));
                     wave = Vector(double(ix)*basis, double(iy)*basis, double(iz)*basis);
                     waveArrayVector_.push_back(wave);
                     weight = exp(- 0.5 * double(sqsum)*basis*basis / qstarSq_);
                     waveArrayWeight_.push_back(weight);
                     wNormal += weight;
                  }
               }
            }
         }
      }
      nWaves_ = waveArrayIndex_.size();
      for (unsigned j = 0; j < nWaves_; ++j)
         waveArrayWeight_[j] /= wNormal;

      waveArrayPcos_.resize(nWaves_, 0.0);
      waveArrayPsin_.resize(nWaves_, 0.0);
      waveArrayDPcos_.resize(nWaves_, 0.0);
      waveArrayDPsin_.resize(nWaves_, 0.0);
   }

}

#endif

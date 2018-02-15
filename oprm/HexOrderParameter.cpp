#ifndef HEXORDERPARAMETER_CPP
#define HEXORDERPARAMETER_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "HexOrderParameter.h"
#include <util/Constants.h>

namespace GridMC
{

   /*
   * Constructor.
   */
   HexOrderParameter::HexOrderParameter(const System& system) :
      OrderParameter(system),
      oprmFlag_(0),
      nWaves_(0),
      waveIndex_(),
      waveVector_(),
      pcos_(),
      psin_(),
      dpcos_(),
      dpsin_()
   {}

   /*
   * Default destructor.
   */
   HexOrderParameter::~HexOrderParameter()
   {}

   /*
   * Read parameter.
   */
   void HexOrderParameter::readParam(std::istream& in)
   {
      char        comment[200];
      std::string line;

      // Order parameter type flag.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: order parameter flag");
      sscanf(line.c_str(), "%d %s", &oprmFlag_, comment);

      if (oprmFlag_ < 0)
         UTIL_THROW("order parameter flag  can't be negative");

      // Set wave vector if needed.
      if (oprmFlag_ != 0) {
         waveIndex_.push_back(IntVector(1, -1, 0));
         waveIndex_.push_back(IntVector(0, 1, -1));
         waveIndex_.push_back(IntVector(-1, 0, 1));

         nWaves_ = 3;

         unsigned j;
         int i;
         for (j = 0; j < nWaves_; ++j) {
            waveVector_.push_back(Vector(0.0));
            for (i = 0; i < Dimension; ++i) 
               waveVector_[j][i] = 2.0 * Constants::Pi / boxL_[i] * double(waveIndex_[j][i]);
         }

         pcos_.resize(nWaves_, 0.0);
         psin_.resize(nWaves_, 0.0);
         dpcos_.resize(nWaves_, 0.0);
         dpsin_.resize(nWaves_, 0.0);
      }

      // Initialize number of bead.
      nBead_ = beads_.size();
   }

   /*
   * Write parameter.
   */
   void HexOrderParameter::writeParam(std::ostream& out) const
   {
      out << "Order parameter type flag   " << oprmFlag_ << std::endl;
   }

   /*
   * Calculate order parameter.
   */
   void HexOrderParameter::calculate()
   {
      if (oprmFlag_ == 0)
         oprmFromGrid();
      else
         oprmFromWave();
   }

   /*
   * Calculate order parameter change dure to bead move.
   */
   double HexOrderParameter::getChange(const Particle& bead, const Vector& rTrial)
   {
      if (oprmFlag_ == 0)
         oprmChangeFromGrid(bead.t);
      else
         oprmChangeFromWave(bead, rTrial);
      return dpsi_;
   }

   /*
   * Calculate order parameter.
   */
   void HexOrderParameter::oprmFromWave()
   {
      unsigned i, j;
      double rdotq;

      for (j = 0; j < nWaves_; ++j) {
         pcos_[j] = 0.0;
         psin_[j] = 0.0;
      }

      for (i = 0; i < nBead_; ++i) {
         for (j = 0; j < nWaves_; ++j) {
            rdotq = beads_[i].r.dot(waveVector_[j]);
            if (beads_[i].t == 0) {
               pcos_[j] += cos(rdotq);
               psin_[j] += sin(rdotq);
            } else if (beads_[i].t == 1) {
               pcos_[j] -= cos(rdotq);
               psin_[j] -= sin(rdotq);
            }
         }
      }

      for (j = 0; j < nWaves_; ++j) {
         pcos_[j] /= double(nBead_);
         psin_[j] /= double(nBead_);
      }

      psi_ = 0.0;
      for (j = 0; j < nWaves_; ++j)
         psi_ += sqrt(pcos_[j]*pcos_[j] + psin_[j]*psin_[j]);
      psi_ /= sqrt( double(nWaves_) );
   }

   /*
   * Calculate order parameter change dure to bead move.
   */
   void HexOrderParameter::oprmChangeFromWave(const Particle& bead, const Vector& rTrial)
   {
      unsigned j;
      if (bead.t > 1) {
         for (j = 0; j < nWaves_; ++j) {
            dpcos_[j] = 0.0;
            dpsin_[j] = 0.0;
         }
         dpsi_  = 0.0;
      } else {
         double dot1, dot2;

         if (bead.t == 0) {

            for (j = 0; j < nWaves_; ++j) {
               dot1 = bead.r.dot(waveVector_[j]);
               dot2 = rTrial.dot(waveVector_[j]);

               dpcos_[j] = (cos(dot2) - cos(dot1)) / double(nBead_);
               dpsin_[j] = (sin(dot2) - sin(dot1)) / double(nBead_);
            }

         } else {

            for (j = 0; j < nWaves_; ++j) {
               dot1 = bead.r.dot(waveVector_[j]);
               dot2 = rTrial.dot(waveVector_[j]);

               dpcos_[j] = (cos(dot1) - cos(dot2)) / double(nBead_);
               dpsin_[j] = (sin(dot1) - sin(dot2)) / double(nBead_);
            }

         }

         dpsi_ = 0.0;
         for (j = 0; j < nWaves_; ++j)
            dpsi_ += sqrt((pcos_[j]+dpcos_[j])*(pcos_[j]+dpcos_[j]) + (psin_[j]+dpsin_[j])*(psin_[j]+dpsin_[j]));
         dpsi_ /= sqrt( double(nWaves_) );
         dpsi_ -= psi_;
      }
   }

   /*
   * Reset the wave vector after the box dimension change.
   */
   void HexOrderParameter::updateAfterBoxMove()
   {
      if (oprmFlag_ != 0) {
         unsigned j;
         int i;
         for (j = 0; j < nWaves_; ++j) {
            for (i = 0; i < Dimension; ++i) 
               waveVector_[j][i] = 2.0 * Constants::Pi / boxL_[i] * double(waveIndex_[j][i]);
         }
      }
   }

}

#endif

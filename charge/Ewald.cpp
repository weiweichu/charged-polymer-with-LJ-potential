#ifndef EWALD_CPP
#define EWALD_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "Ewald.h"
#include "../util/Constants.h"
#include "../util/global.h"

#include <string>

namespace GridMC
{
   using namespace Util;

   /*
   * Constructor. Assume cell to be orthogonal.
   */
   Ewald::Ewald() :
      boxL_(1.0),
      boxV_(1.0),
      esStrength_(0.0),
      alpha_(1.0),
      alphaSqrt_(1.0),
      rMax_(4),
      kMax_(4),
      rvec_(),
      kvec_(),
      ksq_(),
      cosQ_(),
      sinQ_(),
      dcosQ_(),
      dsinQ_()
   {}


   /*
   * Default destructor.
   */
   Ewald::~Ewald()
   {}


   /**
   * Read parameters.
   */
   void Ewald::readParam(std::istream& in)
   {
      char comment[200];
      std::string line;

      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: Ewald summation parameters");
      sscanf(line.c_str(), "%lf %d %d %s", &alpha_, &rMax_, &kMax_, comment);
      alphaSqrt_ = sqrt(alpha_);
   }


   /**
   * Write parameters, pairing with readParam.
   */
   void Ewald::writeParam(std::ostream& out) const
   {
      out << "Ewald summation parameter   " << alpha_
          << "  " << rMax_ << "  " << kMax_ << std::endl;
   }


   /*
   * Set box dimension and Ewald parameters.
   */
   void Ewald::setBox(const Vector& boxLIn)
   {
      boxL_ = boxLIn;
      boxV_ = 1.0;
      for (int i = 0; i < Dimension; ++i) boxV_ *= boxL_[i];

      // Allocate real space shift arrays.
      int k1, k2, k3;
      Vector rv;

      // Use primary box (k1 = k2 = k3 = 0) as the first element.
      rvec_.push_back(Vector(0.0));

      // k1 = 0; k2 = 0; k3 > 0.
      k1 = 0;
      rv[0] = 0.0;

      k2 = 0;
      rv[1] = 0.0;
      for (k3 = 1; k3 <= rMax_; ++k3) {
         rv[2] = boxL_[2] * double(k3);
         rvec_.push_back(rv);
      }

      // k1 = 0; k2 > 0; k3 arbitrary.
      for (k2 = 1; k2 <= rMax_; ++k2) {
         rv[1] = boxL_[1] * double(k2);
         for (k3 = -rMax_; k3 <= rMax_; ++k3) {
            rv[2] = boxL_[2] * double(k3);
            rvec_.push_back(rv);
         }
      }

      // k1 > 0; k2 & k3 arbitrary.
      for (k1 = 1; k1 <= rMax_; ++k1) {
         rv[0] = boxL_[0] * double(k1);
         for (k2 = -rMax_; k2 <= rMax_; ++k2) {
            rv[1] = boxL_[1] * double(k2);
            for (k3 = -rMax_; k3 <= rMax_; ++k3) {
               rv[2] = boxL_[2] * double(k3);
               rvec_.push_back(rv);
            }
         }
      }

      // Allocate wave vector arrays. Assuming box to be orthogonal.
      Vector kv;
      for (k1 = -kMax_; k1 <= kMax_; ++k1) {
         kv[0] = 2.0 * Constants::Pi / boxL_[0] * double(k1);
         for (k2 = -kMax_; k2 <= kMax_; ++k2) {
            kv[1] = 2.0 * Constants::Pi / boxL_[1] * double(k2);
            for (k3 = -kMax_; k3 <= kMax_; ++k3) {
               kv[2] = 2.0 * Constants::Pi / boxL_[2] * double(k3);
               if (abs(k1) + abs(k2) + abs(k3) > 0) {
                  kvec_.push_back(kv);
                  ksq_.push_back(kv.square());
               }
            }
         }
      }

      // Allocate Fourier mode coefficient arrays.
      cosQ_.resize(kvec_.size(), 0.0);
      sinQ_.resize(kvec_.size(), 0.0);

      dcosQ_.resize(kvec_.size(), 0.0);
      dsinQ_.resize(kvec_.size(), 0.0);
   }


   /*
   * Rescale box dimension, interaction strength, and Green's function.
   */
   void Ewald::rescaleBox(const double s1)
   {
      boxL_ *= s1;
      boxV_ *= s1*s1*s1;

      // Rescale real space shift wave vectors.
      for (unsigned i = 0; i < rvec_.size(); ++i)
         rvec_[i] *= s1;

      // Rescale wave vectors.
      for (unsigned i = 0; i < kvec_.size(); ++i) {
         kvec_[i] /= s1;
         ksq_[i] /= (s1*s1);
      }
   }

   /*
   * Read parameter.
   */
   void Ewald::setCoulombStrength(const double strengthIn)
   { esStrength_ = strengthIn; }


   /*
   * Get the totalt coulomb energy.
   */
   double Ewald::getCoulombEnergy(const vector<Particle>& beads)
   {
      double    Eshort(0.0), Elong(0.0), Eself(0.0);
      double    qi, qj, kdotr, d, eij;
      unsigned  i, j, k;
      Vector    ri, rji, dr;

      //Log::file() << rvec_.size() << std::endl;
      //for (k = 0; k < rvec_.size(); ++k)
      //   Log::file() << rvec_[k] << std::endl;

      // Short-range contribution.
      for (i = 0; i < beads.size(); ++i) {
         qi = beads[i].q;
         if (fabs(qi) > Constants::Epsilon) {

            // Interaction between i and its images (a constant depending on rMax).
            eij = 0.0;
            for (k = 1; k < rvec_.size(); ++k) {
               d = rvec_[k].abs();
               eij += erfc(alphaSqrt_*d)/d;
            }
            Eshort += qi*qi*eij;

            // Interaction between i and j (j /= i).
            ri = beads[i].r;
            for (j = i + 1; j < beads.size(); ++j) {
               qj = beads[j].q;
               if (fabs(qj) > Constants::Epsilon) {
                  rji.subtract(beads[j].r, ri);
                  eij = 0.0;

                  // Both i and j are in the primary cell.
                  d = rji.abs();
                  eij += erfc(alphaSqrt_ * d) / d;

                  // Only i or j is in the primary cell.
                  for (k = 1; k < rvec_.size(); ++k) {
                     dr.add(rji, rvec_[k]);
                     d = dr.abs();
                     eij += erfc(alphaSqrt_ * d) / d;

                     dr.subtract(rji, rvec_[k]);
                     d = dr.abs();
                     eij += erfc(alphaSqrt_ * d) / d;
                  }

                  Eshort += qi * qj * eij;
               }
            }
         }
      }


      // Long-range contribution; first evaluate Fourier components.
      for (k = 0; k < kvec_.size(); ++k) {
         cosQ_[k] = 0.0;
         sinQ_[k] = 0.0;
      }

      for (i = 0; i < beads.size(); ++i) {
         qi = beads[i].q;
         if (fabs(qi) > Constants::Epsilon) {
            ri = beads[i].r;
            for (k = 0; k < kvec_.size(); ++k) {
               kdotr = kvec_[k].dot(ri);
               cosQ_[k] += qi * cos(kdotr);
               sinQ_[k] += qi * sin(kdotr);
            }
         }
      }

      for (k = 0; k < kvec_.size(); ++k)
         Elong += (pow(cosQ_[k],2) + pow(sinQ_[k],2))
                  * exp(-0.25*ksq_[k]/alpha_) / ksq_[k];
      Elong *= 2.0 * Constants::Pi / boxV_;

      // Self energy.
      for (i = 0; i < beads.size(); ++i) {
         qi = beads[i].q;
         if (fabs(qi) > Constants::Epsilon)
            Eself += qi*qi;
      }
      Eself *= sqrt(alpha_ / Constants::Pi);

      /*
      Log::file() << Elong / 128.0 * esStrength_;
      Log::file() << "  " << Eself / 128.0 * esStrength_;
      Log::file() << "  " << Eshort / 128.0 * esStrength_;
      Log::file() << "  " << (Elong - Eself + Eshort) / 128.0 * esStrength_ << std::endl;
      */

      return (Elong - Eself + Eshort) * esStrength_;
   }


   /*
   * Get coulomb energy of a particle move.
   */
   double Ewald::getCoulombEnergyChange(const vector<Particle>& beads,
      const int ibead, const Vector& rTrial)
   {
      double    dEshort(0.0), dElong(0.0), dEself(0.0);
      double    qi, qj, kdotr, d, dNew, eij;
      unsigned  j, k;
      Vector    ri, rji, rjiNew, dr;

      // Particle position and charge.
      qi = beads[ibead].q;

      if (fabs(qi) > Constants::Epsilon) {
         ri = beads[ibead].r;

         // Loop over all beads (j /= i).
         for (j = 0; j < beads.size(); ++j) {
            qj = beads[j].q;

            if (fabs(qj) > Constants::Epsilon && j != unsigned(ibead)) {
               rji.subtract(beads[j].r, ri);
               rjiNew.subtract(beads[j].r, rTrial);
               eij = 0.0;

               // Both i and j are in primary cell.
               d = rji.abs();
               dNew = rjiNew.abs();
               eij += erfc(alphaSqrt_ * dNew) / dNew - erfc(alphaSqrt_ * d) / d;

               // i or j are in image cell.
               for (k = 1; k < rvec_.size(); ++k) {
                  dr.add(rji, rvec_[k]);
                  d = dr.abs();
                  dr.add(rjiNew, rvec_[k]);
                  dNew = dr.abs();
                  eij += erfc(alphaSqrt_ * dNew) / dNew - erfc(alphaSqrt_ * d) / d;

                  dr.subtract(rji, rvec_[k]);
                  d = dr.abs();
                  dr.subtract(rjiNew, rvec_[k]);
                  dNew = dr.abs();
                  eij += erfc(alphaSqrt_ * dNew) / dNew - erfc(alphaSqrt_ * d)/ d;
               }

               dEshort += qi * qj * eij;
            }
         }

         // Long-range contribution; first evaluate Fourier component changes.
         for (k = 0; k < kvec_.size(); ++k) {
            kdotr = kvec_[k].dot(ri);
            dcosQ_[k] = -qi * cos(kdotr);
            dsinQ_[k] = -qi * sin(kdotr);
   
            kdotr = kvec_[k].dot(rTrial);
            dcosQ_[k] += qi * cos(kdotr);
            dsinQ_[k] += qi * sin(kdotr);
         }
   
         for (k = 0; k < kvec_.size(); ++k)
            dElong += (pow(cosQ_[k] + dcosQ_[k], 2) + pow(sinQ_[k] + dsinQ_[k], 2)
                      -pow(cosQ_[k], 2)             - pow(sinQ_[k], 2))
                     * exp(-0.25*ksq_[k]/alpha_) / ksq_[k];

         dElong *= 2.0 * Constants::Pi / boxV_;
   
         // Self energy.
         dEself = 0.0;
      }
 
      return (dElong - dEself + dEshort) * esStrength_;
   }

   /*
   * Update Fourier modes.
   */
   void Ewald::updateFourierModes()
   {
      for (unsigned k = 0; k < kvec_.size(); ++k) {
         cosQ_[k] += dcosQ_[k];
         sinQ_[k] += dsinQ_[k];
      }
   }
 
}

#endif

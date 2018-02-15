#ifndef EWALD_H
#define EWALD_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "../simulation/Particle.h"
#include "../util/Constants.h"

#include <vector>

namespace GridMC
{

   using namespace Util;
   using namespace std;

   /**
   * Ewald summation.
   *
   * \ingroup Grid_Module
   */
   class Ewald
   {

   public:
  
      /**
      * Default constructor.
      */
      Ewald();

      /**
      * Default destructor.
      */
      ~Ewald();

      /**
      * Set box dimension.
      */
      void setBox(const Vector& boxLIn);

      /**
      * Read parameters.
      */
      void readParam(std::istream& in);

      /**
      * Write parameters, pairing with readParam.
      */
      void writeParam(std::ostream& out) const;

      /**
      * Rescale box dimension.
      *
      * \param s1  scaling factor for box edge length.
      */
      virtual void rescaleBox(const double s1);

      /**
      * Set Coulomb interaction strength.
      */
      void setCoulombStrength(const double strengthIn);

      /**
      * Get the total Coulomb energy.
      */
      double getCoulombEnergy(const vector<Particle>& beads);

      /**
      * Get the energy associated with a bead move.
      */
      double getCoulombEnergyChange(const vector<Particle>& beads,
                const int ibead, const Vector& rTrial);

      /**
      * Update fourier modes.
      */
      void updateFourierModes();

   protected:

      /// Box dimension.
      Vector boxL_;

      /// Box volume.
      double boxV_;

      /// Strength of coulomb interaction (inverse dielectric constant).
      double  esStrength_;

      /// Parameters for Ewald summation.
      double  alpha_, alphaSqrt_;
      int     rMax_;
      int     kMax_;

      /// Fourier components of charges.
      std::vector<Vector> rvec_;
      std::vector<Vector> kvec_;
      std::vector<double> ksq_;
      std::vector<double> cosQ_;
      std::vector<double> sinQ_;
      std::vector<double> dcosQ_;
      std::vector<double> dsinQ_;

   };

}
#endif

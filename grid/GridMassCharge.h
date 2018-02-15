#ifndef GRIDMASSCHARGE_H
#define GRIDMASSCHARGE_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "GridMass.h"
#include "../simulation/Particle.h"
#include "../util/Constants.h"

#include <vector>

namespace GridMC
{

   using namespace Util;
   using namespace std;

   /**
   * The grid charges defined from TICG simulation.
   *
   * \ingroup Grid_Module
   */
   class GridMassCharge : public GridMass
   {

   public:
  
      /**
      * Default constructor.
      */
      GridMassCharge();

      /**
      * Default destructor.
      */
      ~GridMassCharge();

      /**
      * Rescale box dimension.
      *
      * \param s1  scaling factor for box edge length.
      */
      virtual void rescaleBox(const double s1);

      /**
      * Allocate mass and charge array.
      */
      void allocate(const int nTypeIn);

      /**
      * Set Coulomb interaction strength.
      */
      void setCoulombStrength(const double strengthIn, const double ratioIn);

      /**
      * Project the bead mass and charge onto grid.
      */
      virtual void insertBead(const Particle& bead);

      /**
      * Update density grid from sorted merged list.
      */
      virtual void updateFromList(const Particle& bead);

      /**
      * Update charge grid from sorted merged list.
      */
      void updateChargeGridFromList(const double q);

      /**
      * Get the total Coulomb energy.
      */
      double getCoulombEnergy() const;

      /**
      * Get the self Coulomb energy.
      */
      double getSelfEnergy() const;

      /**
      * Get the energy associated with a bead move.
      */
      double getCoulombEnergyChange(const double q);

      /**
      * Get the total charge on the grid.
      */
      double getTotalCharge() const;

      /**
      * Get the dielectric constant at given position.
      */
      double getDielectricPermittivity(const Vector& r);

      /**
      * Get the total charge on the grid.
      */
      void solvePoisson();

   protected:

      /// Strength of coulomb interaction (inverse dielectric constant).
      double  esStrength_;
      double  esRatio_;

      /// Self energy constant for a unit cube (UC).
      static const double UCself_;

      /// Green's function array.
      vector<double> gridG_;

      /// The dielectric constant array.
      double* dielectric_real_;

      /// The calculated Green's function array.
      double** green_real_;

      /// Charge array.
      std::vector< std::vector< std::vector<double> > > Q_;

      /**
      * Accumulate bead charge to grid charge density.
      */
      void insertCharge(const double q);

      /**
      * Build periodic Green's function.
      */
      void buildGridGreen();

      /*
      * Get the Green's function using real space position.
      */
      double calGreen(const Vector& r) const;

      /*
      * Normalization coefficients; related to symmetry of indices. Eq(6) in paper1.
      *
      * Designed for Dimension = 3
      */
      double gamma(const int k, const int el, const int m) const;

      /*
      * Compact index accessor to the Green's function; homogeneous medium.
      */
      double getGreen(const IntVector& ind) const;

      /*
      * Compact index accessor to the Green's function; inhomogeneous medium.
      */
      double getGreen(const IntVector& ind1, const IntVector& ind2) const;

   };

   /*
   * Get the total charge on the grid.
   */
   inline double GridMassCharge::getTotalCharge() const
   {
      int     i, j, k;
      double  q(0.0);
      for (i = 0; i < nGrid_[0]; ++i)
         for (j = 0; j < nGrid_[1]; ++j)
            for (k = 0; k < nGrid_[2]; ++k)
               q += Q_[i][j][k];
      return q;
   }

   /*
   * Normalization coefficients; related to symmetry of indices. Eq(6) in paper1.
   *
   * Designed for Dimension = 3
   */
   inline double GridMassCharge::gamma(const int k, const int el, const int m) const
   {
      if ( (k*el*m) != 0 )
         return 1.0;
      else if ( (k==0) && (el==0) && (m==0) )
         return 0.0;
      else if ( ( k==0 && el==0 ) || ( k==0 && m==0 ) || ( el==0 && m==0 ) )
         return 0.25;
      else
         return 0.5;
   }

   /*
   * Green's function 3-index accessor.
   *
   * Assume Dimension = 3.
   */
   inline double GridMassCharge::getGreen(const IntVector& ind) const
   {
      int index(0);
      index += abs(ind[0]) * nGridHalf_[1] * nGridHalf_[2];
      index += abs(ind[1]) * nGridHalf_[2];
      index += abs(ind[2]);
      return gridG_[index];
   }

   /*
   * Green's function 3-index pair accessor.
   *
   * Assume Dimension = 3.
   */
   inline double GridMassCharge::getGreen(const IntVector& ind1, const IntVector& ind2) const
   {
      int i = wrap(ind1);
      int j = wrap(ind2);
      return green_real_[i][j];
   }


   /*
   * Update mass density from merged sort list.
   */
   inline void GridMassCharge::updateFromList(const Particle& bead)
   {
      // Update grid mass.
      GridMass::updateFromList(bead);

      // Update grid charge.
      if (fabs(bead.q) > Constants::Epsilon) {
         IntVector ind;
         double q(bead.q);
         for (int i = 0; i < nChanged_; ++i) {
            unwrap(changeId_[i].first, ind);
            Q_[ind[0]][ind[1]][ind[2]] += changeId_[i].second * q;
         }
      }
   }

   /*
   * Update charge grid from merged sort list.
   */
   inline void GridMassCharge::updateChargeGridFromList(const double q)
   {
      // Update grid charge.
      if (fabs(q) > Constants::Epsilon) {
         IntVector ind;
         for (int i = 0; i < nChanged_; ++i) {
            unwrap(changeId_[i].first, ind);
            Q_[ind[0]][ind[1]][ind[2]] += changeId_[i].second * q;
         }
      }
   }

}
#endif

#ifndef GRIDMASS_H
#define GRIDMASS_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <vector>
#include "Grid.h"
#include "../simulation/Particle.h"

namespace GridMC
{

   using namespace std;
   using namespace Util;

   typedef vector<double>            MassV;
   typedef vector< vector<double> >  MassM;

   /**
   * The grid charges defined from TICG simulation.
   *
   * \ingroup Grid_Module
   */
   class GridMass : public Grid
   {
   public:
  
      /**
      * Default constructor.
      */
      GridMass();

      /**
      * Default destructor.
      */
      virtual ~GridMass();

      /**
      * Rescale box dimension.
      *
      * \param s1  scaling factor for box edge length.
      */
      virtual void rescaleBox(const double s1);

      /**
      * Set the interaction strength; for convenient evaluation of 2-body overlap enrgy.
       */
      void allocateMassArray(const int nTypeIn);

      /**
      * Map bead mass and charge onto grid.
      */
      virtual void insertBead(const Particle& bead);

      /**
      * Update grid mass array from sorted merged list.
      */
      virtual void updateFromList(const Particle& bead);

      /**
       * Update neighbor array.
       */
      void updateNeighbor(const int id, const Vector& r, const Vector& rTrial);

      /**
      * Set two-body interaction strength.
      */
      void setTwobodyInteraction(const MassM& chiN, const double kappaN, const double sqrtNbar, const int polymerN);


      /**
      * Get the total two body interaction energy.
      */
      double getTwobodyEnergy() const;

      /**
      * Calculate the two body energy change associated with a bead move.
      */
      double getTwobodyEnergyChange(const int type);

      /**
       * Calculate lj energy change associated with a bead move.
       */
      //double getLjEnergyChange();
      /**
      * Get the total two body interaction energy.
      */
      double getSiteEnergy(const MassV& mass) const;

      /**
      * Get the total masses on the grid.
      */
      MassV getTotalMass() const;

      /**
      * Return mass vector on a grid site.
      */
      const MassV& operator[](int i) const;

      /**
      * Return mass vector on grid site.
      */
      MassV& operator[](int i);

      /**
      * Initialize VTK array.
      */
      void initializeVTK();

      /**
      * Update VTK array.
      */
      void updateVTK();

      /**
      * Reset VTK array.
      */
      void resetVTK();


      /**
      * Output VTK array.
      */
      void writeVTK(ostream& out) const;

   protected:
 
      /// Number of mass types.
      int    nType_;

      /// Mass interaction array.
      MassM  uPair_;

      /// Bead density arrays.
      MassM  phi_;

      /// VTK arrays.
      MassM  vtk_;

      /// Number of frames in VTK.
      int    nFramesVTK_;
      
      /**
      * Add mass to grid.
      */
      void insertMass(const int type);

   };

   /*
   * Get the mass on a grid site.
   */
   inline const MassV& GridMass::operator[](int i) const
   { return phi_[i]; }

   /*
   * Get the mass on a grid site.
   */
   inline MassV& GridMass::operator[](int i)
   { return phi_[i]; }

   /*
   * Get the total masses on the grid.
   */
   inline MassV GridMass::getTotalMass() const
   {
      MassV mass(nType_, 0.0);
      for (unsigned i = 0; i < phi_.size(); ++i)
         for (int j = 0; j < nType_; ++j)
            mass[j] += phi_[i][j];
      return mass;
   }

   /*
   * Get the total two body interaction energy.
   */
   inline double GridMass::getSiteEnergy(const MassV& mass) const
   {
      double energy(0.0);
      int    i, j;
      for (i = 0; i < nType_; ++i)
         for (j = i; j < nType_; ++j)
            energy += mass[i] * uPair_[i][j] * mass[j];
      return energy;
   }

   /*
   * Get two body interaction energy.
   */
   inline double GridMass::getTwobodyEnergy() const
   {
      double energy(0.0);
      for (unsigned i = 0; i < phi_.size(); ++i)
         energy += getSiteEnergy(phi_[i]);
      return energy;
   }

   /*
   * Accumulate bead mass to density field.
   */
   inline void GridMass::insertMass(const int type)
   {
      for (unsigned i = 0; i < insertId_.size(); ++i)
         phi_[insertId_[i].first][type] += insertId_[i].second;
   }

   /*
   * Update mass density from merged sort list.
   */
   inline void GridMass::updateFromList(const Particle& bead)
   {
      for (int i = 0; i < nChanged_; ++i)
         phi_[changeId_[i].first][bead.t] += changeId_[i].second;
   }

}
#endif

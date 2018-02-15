#ifndef MCSYSTEM_H
#define MCSYSTEM_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <vector>
#include <oprm/LamOrderParameter.h>
#include "HistogramWeight.h"
#include "System.h"

namespace GridMC
{

   using namespace std;
   using namespace Util;

   class Simulation;

   /**
   * A collection of energy evaluators, needed by Monte Carlo moves.
   *
   * \ingroup Simulation_Module
   */
   class McSystem : public System
   {
   
   public:

      /**
      * Constructor.
      */
      McSystem(const Simulation& simulationIn);

      /**
      * Default destructor.
      */
      ~McSystem();

      /**
      * Read parameters.
      */
      virtual void readParam(istream& in);

      /**
      * Write parameters.
      */
      virtual void writeParam(ostream& out);

      /**
      * Evaluate system energy and order parameter.
      */
      void calculateStatusVariable();

      /**
      * Update Coulomb energy from input.
      */
      void updateCoulombEnergy(const double deltaCoulombE);

      /**
      * Update bead position, grid density, and system energies.
      */
      void updateStatusVariable(const int id, const Vector& rTrial);

      /**
      * Backup status variables.
      */
      void backupStatusVariable();

      /**
      * Restore status variables.
      */
      void restoreStatusVariable();

      /**
      * Update DOS based weight.
      */
      void updateDOS();

      /**
      * Calculate the statistical weight associated with a bead move.
      */
      double calculateMoveWeight(const int id, const Vector& rTrial);

      /**
      * Calculate the statistical weight associated with volume change.
      */
      double volumeChangeWeight(const double deltaV);

      /**
      * Affinely update the system volume.
      */
      void updateVolume(const double deltaV);

      /**
      * Get total bond energy.
      */
      double getBondEnergy() const;

      /**
      * Get total nonbond energy.
      */
      double getTwobodyEnergy() const;

      /**
      * Get Coulomb energy.
      */
      double getCoulombEnergy() const;

      /**
      * Get solvation energy.
      */
      double getSolvationEnergy() const;

      /**
       * Get lj energy.
       */
      double getLjEnergy() const;

      /**
      * Get pressure.
      */
      double getPV() const;

      /**
      * Bead bond energy change.
      */
      double getBondEnergyChange(const int id, const Vector& rTrial) const;

      /**
       * LJ energy change.
       */
      double getLjEnergyChange(const int id, const Vector& rTrial);
      /**
      * Get order parameter.
      */
      double getOrderParameter() const;

      /**
      * Return doDOS_.
      */
      bool doDOS() const;

      /**
      * Return doDOS_.
      */
      bool doneDOS() const;

      /**
      * Check the histogram flatness.
      */
      void checkDOSflatness();

      /**
      * Get histogram roughness.
      */
      double getHistogramRoughness() const;

      /**
      * Get the random mixing two body energy.
      */
      double getRandomMixingEnergy() const;

   private:

      friend class VolumeMove;
      friend class BeadMove;
      friend class SemigrandSaltMove;
      friend class Simulation;

      /// Bare flory-huggins parameter.
      vector< vector<double> >   chiN_;

      /// Incompressibility parameter.
      double kappaN_;

      /// Overlap constant; effectively, bead volume.
      double sqrtNbar_;

      /// Coulomb interaction strenght; inversely proportional to dielectric constant.
      double esStrength_;
      double esRatio_;
      /// Ion radius
      double aBorn_;
      /// LJ parameters
      double cutoff_;
      double sigma_;
      double epsi_;

      /// If is isobaric system.
      bool   isIsobaric_;
      double barostatPressure_;

      /// Density of states (DOS) sampling parameters.
      bool             doDOS_;
      HistogramWeight  dosWeight_;

      /// Status variables: energies and change.
      double bondEnergy_;
      double twobodyEnergy_;
      double coulombEnergy_;
      double solvationEnergy_;
      double ljEnergy_;
      double pV_;

      double dBondEnergy_;
      double dTwobodyEnergy_;
      double dCoulombEnergy_;
      double dSolvationEnergy_;
      double dLjEnergy_;

      double bondEnergyBak_;
      double twobodyEnergyBak_;
      double coulombEnergyBak_;
      double solvationEnergyBak_;
      double ljEnergyBak_;
      double pVBak_;

      /// Status variables: order parameter.
      LamOrderParameter psi_;

      /// Auxiliary integer used for generating trial configurations of specified psi value.
      int    i030, i001;
   };

   /**
   * Get total bond energy.
   */
   inline double McSystem::getBondEnergy() const
   { return bondEnergy_; }

   /*
   * Get total twobody energy.
   */
   inline double McSystem::getTwobodyEnergy() const
   { return twobodyEnergy_; }

   /**
   * Get coulomb energy.
   */
   inline double McSystem::getCoulombEnergy() const
   { return coulombEnergy_; }

   /**
   * Get solvation energy.
   */
   inline double McSystem::getSolvationEnergy() const
   { return solvationEnergy_; }

   /**
   * Get lj energy.
   */
   inline double McSystem::getLjEnergy() const
   { return ljEnergy_;}
   /**
   * Get coulomb energy.
   */
   inline double McSystem::getPV() const
   { return pV_; }

   /**
   * Get order parameter.
   */
   inline double McSystem::getOrderParameter() const
   { return psi_.get(); }

   /*
   * Return doDOS_.
   */
   inline bool McSystem::doDOS() const
   { return doDOS_; }

   /*
   * Return true if histogram is sufficiently flat.
   */
   inline bool McSystem::doneDOS() const
   { return doDOS_ && dosWeight_.doneIteration(); }

   /*
   * Return histogram roughness.
   */
   inline double McSystem::getHistogramRoughness() const
   { return dosWeight_.getRoughness(); }

   /*
   * Get the random mixing two body potential.
   */
   inline double McSystem::getRandomMixingEnergy() const
   {
      double nSites = double(nGrid_[0]*nGrid_[1]*nGrid_[2]);
      vector<double> mass(nType_, 0.0);
      double energy;
      mass[0] = double(nA_ * nPolymers_) / nSites;
      mass[1] = double(nB_ * nPolymers_) / nSites;
      mass[2] = double(nSolvents_) / nSites;
      energy = grid_.getSiteEnergy(mass);
      return (energy * nSites);
   }

   /*
   * Backup status variable.
   */
   inline void McSystem::backupStatusVariable()
   {
      bondEnergyBak_      = bondEnergy_;
      twobodyEnergyBak_   = twobodyEnergy_;
      coulombEnergyBak_   = coulombEnergy_;
      solvationEnergyBak_ = solvationEnergy_;
      ljEnergyBak_        = ljEnergy_;
      pVBak_              = pV_;
      //if (useEwald_ > 0)
      //   ewald_.backupFourierModes();
      //psi_.backup();
   }

   /*
   * Restore status variable.
   */
   inline void McSystem::restoreStatusVariable()
   {
      bondEnergy_      = bondEnergyBak_;
      twobodyEnergy_   = twobodyEnergyBak_;
      coulombEnergy_   = coulombEnergyBak_;
      solvationEnergy_ = solvationEnergyBak_;
      ljEnergy_        = ljEnergyBak_;
      pV_            = pVBak_;
      //if (useEwald_ > 0)
      //   ewald_.restoreFourierModes();
      //psi_.restore();
   }

}

#endif

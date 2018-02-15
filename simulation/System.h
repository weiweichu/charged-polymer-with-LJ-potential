#ifndef SYSTEM_H
#define SYSTEM_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <vector>
#include <math.h>

#include <util/FileMaster.h>
#include <util/Mersenne.h>
#include <util/Vector.h>
#include <util/IntVector.h>
#include "../grid/GridMassCharge.h"
#include "../charge/Ewald.h"
#include "Particle.h"

namespace GridMC
{

   using namespace std;
   using namespace Util;

   class Simulation;

   /**
   * A collection of particles.
   *
   * \ingroup Simulation_Module
   */
   class System 
   {
   
   public:

      /**
      * Constructor.
      */
      System(const Simulation& simulationIn);

      /**
      * Default destructor.
      */
      virtual ~System();

      /**
      * Read parameters.
      */
      virtual void readParam(std::istream& in);

      /**
      * Write parameters.
      */
      virtual void writeParam(std::ostream& out);

      /**
      * Set pointer to file master.
      */
      void setFileMaster(FileMaster& fileMaster);

      /**
      * Generate system configuration using random walks.
      */
      void generateRandomConfig();

      /**
      * Read configuration.
      */
      void readConfig(istream& in);

      /**
      * Output configuration.
      */
      void writeConfig(ostream& out) const;

      /**
      * Initialize VTK array.
      */
      void initializeVTK();

      /**
      * Update VTK array.
      */
      void updateVTK(const int iStep, const int flag = -1);

      /**
      * Output VTK array.
      */
      void writeVTK(ostream& out) const;

      /**
      * Solve inhomogeneous Poisson's equation.
      */
      void solvePoisson();

      /**
      * Constant reference to the box dimension.
      */
      const Vector& getBoxL() const;

      /**
      * Constant reference to nGrid.
      */
      const IntVector& getNGrid() const;

      /**
      * Constant reference to the grid object.
      */
      const GridMassCharge& getGrid() const;

      /**
      * Accesor to the number of bead types.
      */
      int getNType() const;

      /**
      * Constant reference to beads array.
      */
      const vector<Particle>& getBeads() const;

      /**
      * Constant reference to random number generator.
      */
      CRandomMersenne& getRandom();

      /**
      * Constant reference to file master.
      */
      const FileMaster& getFileMaster() const;

      /**
      * Generate random bond length.
      */
      double randomBondLength();

      /**
      * Return reference to the simulation object.
      */
      const Simulation& simulation() const;

      /**
      * Shift particle position into the primary cell.
      */
      void toPrimaryCell(Vector& r) const;

   protected:

      friend class BeadMove;
      friend class ChainFlipMove;
      friend class ReptationMove;
      friend class SemigrandSaltMove;
      friend class Diagnosis;
      friend class SolventChemicalPotential;
      friend class SaltChemicalPotential;
      friend class ChainChemicalPotential;
      friend class ChainChemicalPotentialInsertion;
      friend class IonNumberFraction;
      friend class Simulation;

      /// Reference to simulation (to access simulation parameters).
      const Simulation* simulationPtr_;

      /// Box geometry.
      Vector    boxL_;
      double    boxV_;
      IntVector nGrid_;

      /// Number of bead types.
      int       nType_;

      /// Molecule parameter.
      int       nA_, nB_, nAB_;
      int       nPolymers_, nPolymerBeads_;
      int       nIons_, nNeutralSolvents_, nSolvents_;
      int       micellePolymerNQ_, micelleNIons_, micelleIonQ_;
      int       qCode_, useEwald_;
      int       nFreeIons_;
      int       chargeCount_ = 0;
      vector<double> chargeDistribution_;
      double    chargeDensity_;
      /// Initial configuration parameter.
      char      configFileName_[200];

      /// Initial configuration parameter.
      double    bondL_;

      /// Beads array.
      vector<Particle> beads_;

      /// Pointer to polymers's first beads.
      vector<Particle*> polymer_;

      /// Pointer to first solvent bead.
      Particle* solvent_;

      /// Grid.
      GridMassCharge grid_;

      /// Ewald.
      Ewald ewald_;

      /// Random number generator.
      int             seed_;
      CRandomMersenne random_;

      /// File master.
      FileMaster*     fileMasterPtr_;

      /**
      * Project bead position and charge onto grid.
      */
      void makeGridDensity();

      /**
      * Create polymers and solvents.
      */
      void allocate();

      /**
      * Shift particle bond vector using the minimum image convention.
      */
      void pbcShift(Vector& dr) const;

   };

   /*
   * Return the simulation pointer.
   */
   inline const Simulation& System::simulation() const
   { return *simulationPtr_; }

   /*
   * Set pointer to file master.
   */
   inline void System::setFileMaster(FileMaster& fileMaster)
   { fileMasterPtr_ = &fileMaster; }

   /*
   * Constant reference to the box dimension.
   */
   inline const Vector& System::getBoxL() const
   { return boxL_; }

   /*
   * Output configuration.
   */
   inline const IntVector& System::getNGrid() const
   { return nGrid_; }

   /*
   * Output configuration.
   */
   inline const GridMassCharge& System::getGrid() const
   { return grid_; }

   /*
   * Get number of bead types.
   */
   inline int System::getNType() const
   { return nType_; }

   /*
   * Return reference to beads array.
   */
   inline const vector<Particle>& System::getBeads() const
   { return beads_; }

   /**
   * Constant reference to random number generator.
   */
   inline CRandomMersenne& System::getRandom()
   { return random_; }

   /**
   * Constant reference to file master.
   */
   inline const FileMaster& System::getFileMaster() const
   { return *fileMasterPtr_; }

   /*
   * Shift particle position into the primary cell.
   */
   inline void System::toPrimaryCell(Vector& r) const
   {
      for (int i = 0; i < Dimension; ++i)
         r[i] -= floor(r[i] / boxL_[i]) * boxL_[i];
   }


   /*
   * Shift particle bond vector using the minimum image convention.
   */
   inline void System::pbcShift(Vector& dr) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if (dr[i] > 0.5 * boxL_[i])
            dr[i] -= boxL_[i];
         else if (dr[i] < -0.5 * boxL_[i])
            dr[i] += boxL_[i];
      }
   }

   /*
   * Generate a random bond length. The desired distribution is proportional to
   *    x * x * Gaussian(x)
   */
   inline double System::randomBondLength()
   {
      double sd  = bondL_/sqrt(3.0);
      double xm2 = 7.0*sd;
             xm2 = xm2 * xm2;
      double x;

      //Generate trials until one is accepted
      for (;;) {
         x = sd*random_.GaussianRandom();
         if (x < 0) continue;
         if ( x*x/xm2 >= random_.Random() ) {
            return x;
         }
      }
   }

}

#endif

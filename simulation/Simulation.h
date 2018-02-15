#ifndef SIMULATION_H
#define SIMULATION_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <vector>
#include <util/FileMaster.h>
#include <move/BeadMove.h>
#include <move/VolumeMove.h>
#include <move/ChainFlipMove.h>
#include <move/ReptationMove.h>
#include <move/SemigrandSaltMove.h>
#include <diag/Diagnosis.h>

#include "McSystem.h"

namespace GridMC
{
   using namespace Util;

   /**
   * The simulation protocol.
   *
   * \ingroup Simulation_Module
   */
   class Simulation
   {
   
   public:

      /**
      * Constructor.
      */
      Simulation(const char* fileName);

      /**
      * Default destructor.
      */
      ~Simulation();

      /**
      * Run the simulation.
      */
      void run();

      /*
      * Input parameters.
      */
      void readParam(std::ifstream& prmfile);

      /*
      * Output parameters.
      */
      void writeParam(std::ostream& out) const;

   private:

      friend class McSystem;

      /// A system object contains box dimension, particle positions, etc.
      McSystem       system_;

      /// Number of Monte Carlo steps.
      long int       nMCsteps_;

      /// Monte Carlo moves.
      double         pBeadMove_;
      BeadMove       beadMove_;

      double         pChainFlipMove_;
      ChainFlipMove  chainFlipMove_;

      double         pReptationMove_;
      ReptationMove  reptationMove_;

      double             pSemigrandSaltMove_;
      SemigrandSaltMove  semigrandSaltMove_;

      int            nVolumeMoveInterval_;
      VolumeMove     volumeMove_;

      /// Diagnostic intervals.
      long int       iStep_;
      int            nFlatnessInterval_;
      int            nDielectricInterval_;
      int            nEnergyInterval_;
      int            nConfigInterval_;

      /// Diagnosis objects.
      int                baseInterval_;
      int                nDiagnosis_;
      vector<Diagnosis*> diagnosis_;

      /// Data I/O.
      FileMaster     fileMaster_;
      std::ofstream  logFile_;

      #ifdef UTIL_MPI

      /// Pointer to the simulation communicator.
      MPI::Intracomm* communicatorPtr_;

      #endif

      /**
      * Output system status variables.
      */
      void echoSystemStatus(std::ostream& out) const;

      /**
      * Read parameters for diagnosis objects.
      */
      void readDiagnosis(std::istream& in);

      /**
      * Write diagnosis object parameter.
      */
      void writeDiagnosis(std::ostream& out) const;

   };

}
#endif

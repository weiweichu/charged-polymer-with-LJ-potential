#ifndef SIMULATION_CPP
#define SIMULATION_CPP

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

#include "../util/global.h"
#include <diag/StructureFactor.h>
#include <diag/SolventChemicalPotential.h>
#include <diag/SaltChemicalPotential.h>
#include <diag/ChainChemicalPotential.h>
#include <diag/ChainChemicalPotentialInsertion.h>
#include <diag/IonNumberFraction.h>
#include "Simulation.h"

namespace GridMC
{
   using namespace Util;
   using namespace std;

   /*
   * Constructor.
   */
   Simulation::Simulation(const char* prmfileName) :
      system_(*this),
      nMCsteps_(0),
      //
      pBeadMove_(0.0),
      beadMove_(system_),
      //
      pChainFlipMove_(0.0),
      chainFlipMove_(system_),
      //
      pReptationMove_(0.0),
      reptationMove_(system_),
      //
      pSemigrandSaltMove_(0.0),
      semigrandSaltMove_(system_),
      //
      nVolumeMoveInterval_(0),
      volumeMove_(system_),
      //
      nFlatnessInterval_(0),
      nDielectricInterval_(0),
      nEnergyInterval_(0),
      nConfigInterval_(0),
      baseInterval_(0),
      nDiagnosis_(0),
      diagnosis_(),
      fileMaster_()
   {

      #ifdef UTIL_MPI

      if (!MPI::Is_initialized())
         MPI::Init();

      // Initialize MPI, set communicator to COMM_WORLD
      communicatorPtr_ = &(MPI::COMM_WORLD);

      // Set I/O prefix for file using MPI processor rank.
      std::string filePrefix;
      std::stringstream sMyId;
      sMyId << communicatorPtr_->Get_rank();
      filePrefix = sMyId.str();
      filePrefix += "/";
      fileMaster_.setFilePrefix(filePrefix);

      // Set log file for processor n to a new file named "n/log"
      fileMaster_.openOutputFile("log", logFile_);
      Log::setFile(logFile_);

      #endif

      // Set file master for system.
      system_.setFileMaster(fileMaster_);

      // Read paramters.
      std::ifstream prmfile;
      fileMaster_.openInputFile(std::string(prmfileName), prmfile);
      readParam(prmfile);
      prmfile.close();

      // Echo paramters.
      system_.writeParam(Log::file());
      Log::file() << std::endl;
      writeParam(Log::file());
      Log::file() << std::endl;
   }

   /*
   * Default destructor.
   */
   Simulation::~Simulation()
   {
      #ifdef UTIL_MPI

      if (logFile_.is_open())
        logFile_.close();

      MPI::Finalize();

      #endif
   }

   /*
   * Execute the simulation.
   */
   void Simulation::run()
   {
      iStep_ = 0;

      std::ofstream statusfile;
      fileMaster_.openOutputFile("status.dat", statusfile);

      #ifndef UTIL_MPI
      echoSystemStatus(Log::file());
      #endif
      echoSystemStatus(statusfile);

      // Prepare the move selection probability: beadMove, chainFlipMove, reptationMove, semigrandSaltMove
      double pMove[4];
      pMove[0] = pBeadMove_;
      pMove[1] = pMove[0] + pChainFlipMove_;
      pMove[2] = pMove[1] + pReptationMove_;
      pMove[3] = pMove[2] + pSemigrandSaltMove_;
      for (int i = 0; i < 4; ++i)
         pMove[i] /= pMove[3];
 
      // Run the simulation.
      double randomNum;
      for (iStep_ = 1; iStep_ <= nMCsteps_ && !system_.doneDOS(); ++iStep_) {

         randomNum = system_.random_.Random();
         if (randomNum < pMove[0])
            beadMove_.move();
         else if (randomNum < pMove[1])
            chainFlipMove_.move();
         else if (randomNum < pMove[2])
            reptationMove_.move();
         else
            semigrandSaltMove_.move();

         // Volume move.
         if (system_.isIsobaric_ && nVolumeMoveInterval_ > 0)
            if (iStep_ % nVolumeMoveInterval_ == 0)
               volumeMove_.move();

         // Energy output.
         if (nEnergyInterval_ > 0) {
            if (iStep_ % nEnergyInterval_ == 0) {
               #ifndef UTIL_MPI
               echoSystemStatus(Log::file());
               #endif
               echoSystemStatus(statusfile);
            }
         }

         // Configuration output.
         if (nConfigInterval_ > 0) {
            if (iStep_ % nConfigInterval_ == 0)
               system_.updateVTK(iStep_, 1);
         } else if (nConfigInterval_ < 0) {
            if (iStep_ % -nConfigInterval_ == 0)
               system_.updateVTK(iStep_, -1);
         }

         // Density of states sampling.
         if (system_.doDOS() && nFlatnessInterval_ > 0) {
            if (iStep_ % nFlatnessInterval_ == 0)
               system_.checkDOSflatness();
         }

         // Diagnosis.
         if (nDiagnosis_ > 0 && iStep_ % baseInterval_ == 0) {
            for (int i = 0; i < nDiagnosis_; ++i)
               diagnosis_[i]->sample(iStep_);
         }

         // Dielectric constant.
         if (nDielectricInterval_ > 0) {
            if (iStep_ % nDielectricInterval_ == 0) {
               //system_.calculateStatusVariable();
               system_.solvePoisson();
               system_.calculateStatusVariable();
            }
         }

      }
      statusfile.close();

      // Output system configuration.
      std::ofstream outconfig;
      fileMaster_.openOutputFile("config.out", outconfig);
      system_.writeConfig(outconfig);
      outconfig.close();

      if (nConfigInterval_ != 0) {
         fileMaster_.openOutputFile("config.vtk", outconfig);
         system_.writeVTK(outconfig);
         outconfig.close();
      }

      // Output move statistics.
      Log::file() << std::endl;
      beadMove_.outputStatistics(Log::file());
      chainFlipMove_.outputStatistics(Log::file());
      reptationMove_.outputStatistics(Log::file());
      if (system_.isIsobaric_)
         volumeMove_.outputStatistics(Log::file());

      // Output diagnostic average.
      for (int i = 0; i < nDiagnosis_; ++i)
         diagnosis_[i]->output();

      // Output histogram and weight if needed.
      if (system_.doDOS()) {
         std::ofstream dosFile;
         fileMaster_.openOutputFile("histogram.out", dosFile);
         system_.dosWeight_.outputHistogram(dosFile);
         dosFile.close();

         if (nFlatnessInterval_ > 0) {
            fileMaster_.openOutputFile("weight.out", dosFile);
            system_.dosWeight_.outputWeight(dosFile);
            dosFile.close();
         }
      }
   }

   /*
   * Output parameters.
   */
   void Simulation::readParam(std::ifstream& prmfile)
   {
      char comment[200];
      string   line;

      // Read system parameters and allocate arrays.
      system_.readParam(prmfile);

      // Parse simulation header.
      getline(prmfile, line);
      getline(prmfile, line);
      if (line.find("Simulation") == string::npos)
         UTIL_THROW("reading error: simulation header.");

      // Read number of MC moves.
      getline(prmfile, line);
      if (line.size() > 0) {
         sscanf(line.c_str(), "%ld  %s", &nMCsteps_, comment);
      } else {
         UTIL_THROW("reading error: nMCmoves");
      }

      // Read bead move parameter.
      getline(prmfile, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: pBeadMove");
      sscanf(line.c_str(), "%lf  %s", &pBeadMove_, comment);
      if (pBeadMove_ < 0.)
         UTIL_THROW("reading error: invalid beadmove probability");
      beadMove_.readParam(prmfile);

      // Read chain flip move parameter.
      getline(prmfile, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: pChainFlipMove");
      sscanf(line.c_str(), "%lf  %s", &pChainFlipMove_, comment);
      if (pChainFlipMove_ < 0.)
         UTIL_THROW("reading error: invalid chain flip move probability");
      chainFlipMove_.readParam(prmfile);

      // Read reptation move parameter.
      getline(prmfile, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: pReptationMoveInterval");
      sscanf(line.c_str(), "%lf  %s", &pReptationMove_, comment);
      if (pReptationMove_ < 0.)
         UTIL_THROW("reading error: invalid reptation move probability");
      reptationMove_.readParam(prmfile);

      // Read semigrand salt move parameter.
      getline(prmfile, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: pSemigrandSaltMoveInterval");
      sscanf(line.c_str(), "%lf  %s", &pSemigrandSaltMove_, comment);
      if (pSemigrandSaltMove_ < 0.)
         UTIL_THROW("reading error: invalid semigrand salt move probability");
      semigrandSaltMove_.readParam(prmfile);

      // Read volume move parameter.
      if (system_.isIsobaric_) {
         getline(prmfile, line);
         if (line.size() <= 0)
            UTIL_THROW("reading error: nVolumeMoveInterval");
         sscanf(line.c_str(), "%d  %s", &nVolumeMoveInterval_, comment);

         volumeMove_.readParam(prmfile);
      }

      // Read interval for checking histogram flatness.
      if (system_.doDOS()) {
         getline(prmfile, line);
         if (line.size() <= 0)
            UTIL_THROW("reading error: nFlatnessInterval");
         sscanf(line.c_str(), "%d  %s", &nFlatnessInterval_, comment);
      }

      // Read dielectric array update interval.
      getline(prmfile, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: nDielectricInterval");
      sscanf(line.c_str(), "%d  %s", &nDielectricInterval_, comment);

      // Read energy output interval.
      getline(prmfile, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: nEnergyInterval");
      sscanf(line.c_str(), "%d  %s", &nEnergyInterval_, comment);

      // Read configuration output interval.
      getline(prmfile, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: nConfigInterval");
      sscanf(line.c_str(), "%d  %s", &nConfigInterval_, comment);
      if (nConfigInterval_ > nDielectricInterval_)
         UTIL_THROW("Improper configuration interval: it should be smaller than dielectric interval");

      // Read base diagnosis interval.
      getline(prmfile, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: number of diagnosis objects");
      sscanf(line.c_str(), "%d  %s", &nDiagnosis_, comment);

      if (nDiagnosis_ > 0) {
         // Read base diagnosis interval.
         getline(prmfile, line);
         if (line.size() <= 0)
            UTIL_THROW("reading error: base diagnosis interval_");
         sscanf(line.c_str(), "%d  %s", &baseInterval_, comment);

         // Read diagnosis object.
         readDiagnosis(prmfile);
      }

   }


   /*
   * Output parameters.
   */
   void Simulation::writeParam(std::ostream& out) const
   {
      out << "-----  Simulation parameters -----" << endl;
      out << "number of MC moves          " << nMCsteps_ << endl;
      out << "bead move probability       " << pBeadMove_ << endl;
      beadMove_.writeParam(out);
      out << "chain flip move probability " << pChainFlipMove_ << endl;
      chainFlipMove_.writeParam(out);
      out << "reptation move probability  " << pReptationMove_ << endl;
      reptationMove_.writeParam(out);
      if (system_.isIsobaric_) {
         out << "volume move interval        " << nVolumeMoveInterval_ << endl;
         volumeMove_.writeParam(out);
      }

      if (system_.doDOS())
         out << "histogram flatness interval " << nFlatnessInterval_ << endl;
      out << "dielectric update interval  " << nDielectricInterval_ << endl;
      out << "energy output interval      " << nEnergyInterval_ << endl;
      out << "config output interval      " << nConfigInterval_ << endl;
      out << "number of diagnosis objects " << nDiagnosis_ << endl;

      // Write diagnosis prameters.
      if (nDiagnosis_ > 0) {
         out << "diagnosis base interval     " << baseInterval_ << endl;
         writeDiagnosis(out);
      }
   }

   /*
   * Read parameters for diagnosis objects.
   */
   void Simulation::readDiagnosis(std::istream& in)
   {
      char comment[200];
      string   line;

      for (int i = 0; i < nDiagnosis_; ++i) {
         // Read name of diagnosis.
         char name[200];
         getline(in, line);
         if (line.size() <= 0)
            UTIL_THROW("reading error: name of diagnosis object");
         sscanf(line.c_str(), "%s %s", name, comment);

         // Create diagnosis objects.
         if (std::string(name) == "StructureFactor") {
            diagnosis_.push_back(new StructureFactor(system_));
            diagnosis_[i]->readParam(in);
         } else if (std::string(name) == "SolventChemicalPotential") {
            diagnosis_.push_back(new SolventChemicalPotential(system_));
            diagnosis_[i]->readParam(in);
         } else if (std::string(name) == "SaltChemicalPotential") {
            diagnosis_.push_back(new SaltChemicalPotential(system_));
            diagnosis_[i]->readParam(in);
         } else if (std::string(name) == "ChainChemicalPotential") {
            diagnosis_.push_back(new ChainChemicalPotential(system_));
            diagnosis_[i]->readParam(in);
         } else if (std::string(name) == "ChainChemicalPotentialInsertion") {
            diagnosis_.push_back(new ChainChemicalPotentialInsertion(system_));
            diagnosis_[i]->readParam(in);
         } else if (std::string(name) == "IonNumberFraction") {
            diagnosis_.push_back(new IonNumberFraction(system_));
            diagnosis_[i]->readParam(in);
         } else {
            cout << name << endl;
            UTIL_THROW("Invalid diagnosis name");
         }
      }

   }

   /*
   * Write diagnosis object parameter.
   */
   void Simulation::writeDiagnosis(std::ostream& out) const
   {
      for (int i = 0; i < nDiagnosis_; ++i)
         diagnosis_[i]->writeParam(out);
   }


   /*
   * Output system status variables.
   */
   void Simulation::echoSystemStatus(std::ostream& out) const
   {
      double          Eb, Enb, Ecoulomb, Esolvation, Elj, pVorV, psi;
      double          normal_constant;
      vector<double>  mass(system_.grid_.getTotalMass());

      if (system_.nPolymers_ > 0)
         normal_constant = double(system_.nPolymers_);
      else
         normal_constant = double(system_.beads_.size());

      Eb         = system_.getBondEnergy() / normal_constant;
      Enb        = system_.getTwobodyEnergy() / normal_constant;
      Ecoulomb   = system_.getCoulombEnergy() / normal_constant;
      Esolvation = system_.getSolvationEnergy() / normal_constant;
      Elj        = system_.getLjEnergy() / normal_constant;
      if (system_.isIsobaric_)
         pVorV = system_.boxV_;
      else
         pVorV = system_.getPV() / normal_constant;
      if (system_.doDOS())
         psi   = system_.getHistogramRoughness();
      else
         psi   = system_.getOrderParameter();

      out.setf(ios::left, ios::adjustfield);
      out << setw(8) << iStep_;
      out.unsetf(ios::adjustfield);

      // uncomment the lines below if grid mass and charge needs to be print
      //for (int i = 0; i < system_.nType_; ++i)
      //   out << "  " << mass[i];
      out << setw(12) << system_.grid_.getTotalCharge();

      out << fixed;
      out << setw(12) << Eb;
      out << setw(12) << Enb;
      out << setw(12) << Ecoulomb;
      out << setw(12) << Esolvation;
      out << setw(12) << Elj;
      out << setw(12) << pVorV;
      out << setw(12) << psi;
      out << endl;
      out.unsetf(ios::floatfield);
   }
}

#endif

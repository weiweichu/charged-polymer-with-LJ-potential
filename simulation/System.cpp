#ifndef SYSTEM_CPP
#define SYSTEM_CPP

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include "../util/global.h"
#include <time.h>
#include <cstring>
#include <fstream>
#include <iomanip>

namespace GridMC
{
   using namespace Util;
   using namespace std;

   /*
   * Constructor. Assume cell to be orthogonal.
   */
   System::System(const Simulation& simulationIn) :
      simulationPtr_(&simulationIn),
      boxL_(1.0),
      boxV_(1.0),
      nGrid_(1),
      nType_(0),
      nA_(0),
      nB_(0),
      nAB_(0),
      nPolymers_(0),
      nPolymerBeads_(0),
      // Parameters only used by charged polymer simulation.
      nNeutralSolvents_(0),
      nSolvents_(0),
      // Parameters only used by charged micelle simulation.
      micellePolymerNQ_(0),  // Number of charged per charged block on polymer
      micelleNIons_(0),      // Number of macro-ions
      micelleIonQ_(0.0),     // Charge on the macro-ions
      //
      qCode_(0),
      useEwald_(0),
      seed_(1),
      random_(seed_)
   {}

   /*
   * Default destructor.
   */
   System::~System()
   {}

   /*
   * Read parameters.
   */
   void System::readParam(std::istream& in)
   {
      char   comment[200];
      string line;

      // Parse header.
      getline(in, line);
      if (line.find("System") == string::npos)
         UTIL_THROW("reading error: System header.");

      // Read box dimension.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: box dimension");
      sscanf(line.c_str(), "%lf %lf %lf %s", &boxL_[0], &boxL_[1], &boxL_[2], comment);
      boxV_ = boxL_[0] * boxL_[1] * boxL_[2];

      // Read grid dimension.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: grid dimension");
      sscanf(line.c_str(), "%d %d %d %s", &nGrid_[0], &nGrid_[1], &nGrid_[2], comment);

      // Read nType.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: number of bead types");
      sscanf(line.c_str(), "%d %s", &nType_, comment);
      if (nType_ < 2)
         UTIL_THROW("reading error: number of bead type must be greater than 2");

      // Read bond length.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: bond length dimension");
      sscanf(line.c_str(), "%lf %s", &bondL_, comment);
 
      // Read block lengths.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: block length");
      sscanf(line.c_str(), "%d %d %s", &nA_, &nB_, comment);
      nAB_ = nA_ + nB_;

      if (bondL_ < Constants::Epsilon) {
         if (nAB_ > 1)
            bondL_ = 1.0 / sqrt(double(nAB_ - 1));
         else
            bondL_ = 1.0;
      }

      // Read number of chains.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: number of chains");
      sscanf(line.c_str(), "%d %s", &nPolymers_, comment);
      nPolymerBeads_ = nAB_ * nPolymers_;

      // Read number of neutral solvents.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: number of neutral solvents");
      sscanf(line.c_str(), "%d %s", &nNeutralSolvents_, comment);
      nSolvents_ = nNeutralSolvents_;

      // Read charge code.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: charge code");
      sscanf(line.c_str(), "%d %s", &qCode_, comment);

      if (qCode_ == 6) { // Read micelle parameters if charged micelle simulatin.
         getline(in, line);
         if (line.size() <= 0)
            UTIL_THROW("reading error: charge code");
         sscanf(line.c_str(), "%d %d %d %s", &micellePolymerNQ_, &micelleNIons_, &micelleIonQ_, comment);
         nSolvents_ += micelleNIons_;
      }

      if (qCode_ == 8) {
        getline(in, line);
        if (line.size() <= 0)
          UTIL_THROW("reading error: charge code 8");
        sscanf(line.c_str(), "%lf %d %s", &chargeDensity_, &nFreeIons_, comment);
        nSolvents_ += nFreeIons_;
      }
      // Read Ewald code.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: Ewald code");
      sscanf(line.c_str(), "%d %s", &useEwald_, comment);
      if (useEwald_ > 0) {
         ewald_.readParam(in);
         ewald_.setBox(boxL_);
      }

      // Read configuration code.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: input configuration file name");
      sscanf(line.c_str(), "%s %s", configFileName_, comment);

      // Read random number generator seed.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: random number seed");
      sscanf(line.c_str(), "%d %s", &seed_, comment);
      if (seed_ < 0)
         seed_ = -seed_;
      else
         seed_ = (int)time(NULL);
      random_.RandomInit(seed_);

      // Allocate polymer and solvent arrays.
      allocate();
      Log::file() <<"after allocate" << endl;
      // Set grid, read thermodynamic parameters, and allocate grid's memory.
      grid_.setBox(boxL_, nGrid_);
      grid_.allocate(nType_);

      // Generate orr read configuration and initialize grid arrays.
      if (strlen(configFileName_) == 1 && configFileName_[0] == '0') {
         Log::file() << "Generate initial configuration from random walks." << endl;
         generateRandomConfig();
      } else {
         Log::file() << "Read initial configuration from file: " << configFileName_ << endl;
         ifstream  cfgfile;
         fileMasterPtr_->openInputFile(configFileName_, cfgfile);
         readConfig(cfgfile);
         cfgfile.close();
      }

      // Project beads to grids.
      makeGridDensity();
   }

   /*
   * Write parameters.
   */
   void System::writeParam(std::ostream& out)
   {
      out << "-----  System parameters -----" << endl;
      out << "box dimension               " << boxL_[0] << "  " << boxL_[1] << "  " << boxL_[2] << endl;
      out << "grid dimension              " << nGrid_[0] << "  " << nGrid_[1] << "  " << nGrid_[2] << endl;
      out << "type of masses              " << nType_ << endl;
      out << "bond length                 " << bondL_ << endl;
      out << "block length                " << nA_ << "  " << nB_ << endl;
      out << "number of chains            " << nPolymers_ << endl;
      out << "number of neutral solvents  " << nNeutralSolvents_ << endl;
      out << "charge code                 " << qCode_ << endl;
      if (qCode_ == 6)
         out << "micelle charge parameter    " << micellePolymerNQ_ << "  " << micelleNIons_ << "  " << micelleIonQ_ << endl;
      if (qCode_ == 8)
         out << "charged polymer parameter    " << chargeDensity_ << "  "  << nFreeIons_ << endl;
      out << "Ewald code                  " << useEwald_ << endl;
      if (useEwald_ > 0) ewald_.writeParam(out);
      out << "input config file name      " << configFileName_ << endl;
      out << "random generator seed       " << seed_ << endl;

      out << "(derived)" << endl;
      out << "polymer length              " << nAB_ << endl;
      out << "number of polymer beads     " << nPolymerBeads_ << endl;
      out << "number of solvents          " << nSolvents_ << endl;

      out << endl;
   }

   /*
   * Create polymer and solvents.
   */
   void System::allocate()
   {
      int        i, j, k;
      double     qA, qB;
      Particle*  ptr;

      // Calculate number of solvent beads by charge code.
      switch (qCode_) {
         case 0:
            break;
         case 1:
            nSolvents_ += nA_ * nPolymers_;
            break;
         case 2:
            nSolvents_ += nAB_ * nPolymers_;
            break;
         case 3:
            if (nA_ % 2 == 0)
               nSolvents_ += (nA_ / 2) * nPolymers_;
            else
               nSolvents_ += ((nA_ + 1) / 2) * nPolymers_;
            break;
         case 4: // Polymer chains are polycation/polyanion pairs
            if (nPolymers_ % 2 != 0)
               UTIL_THROW("Number of polymers uneven for polycation/polyanion pairs");
            break;
         case 5: // Divide the neutral solvents into two groups
            if (nNeutralSolvents_ % 2 != 0)
               UTIL_THROW("Number of neutral solvents is not an even number");
            if (nType_ < 4)
               UTIL_THROW("Number of particle types is less than 4");
            break;
         case 6: // Multiblock charged micelle.
            if (nNeutralSolvents_ < 0)
               UTIL_THROW("Number of counterions less than 0");
            if (micelleNIons_ * micelleIonQ_ % (micellePolymerNQ_ * nPolymers_) != 0)
               UTIL_THROW("Counterion charge doesn't match polymer charge");
            if (nType_ < 3)
               UTIL_THROW("Number of particle types is less than 3");
            break;
         case 7: // Alternating charged polymers.
            if (nPolymers_ % 2 != 0)
               UTIL_THROW("Number of polymers uneven");
            break;
         case 8:
            if (nType_ < 4 && nFreeIons_ > 0)
             	 UTIL_THROW("Number of particle types is less than 4");
        	  if (nFreeIons_ % 2 != 0)
        	   	 UTIL_THROW("Number of free ions is not even");
        	  if (chargeDensity_ > 1)
        	     UTIL_THROW("Charge density should no greater than one");
            //Log::file() << "charge distribution for polymers:  ";
            for (int i = 0; i < nB_ * nPolymers_; i++) {
	          	 if(random_.Random()<chargeDensity_){
                  chargeDistribution_.push_back(-1.0);
                  chargeCount_++;
              //    Log::file() << -1.0;
               }
               else{
                  chargeDistribution_.push_back(0.0);
                //   Log::file() << 0.0;
               }
            }
           // Log::file() << endl;
            Log::file() << "System.cpp: number of charged -1.0 type B beads:  " << chargeCount_  << endl;
            Log::file() << "System.cpp: charge density: real  " << chargeCount_ * 1.0 / (nB_ * nPolymers_) <<"  set  " << chargeDensity_ << endl;
        	  nSolvents_ += chargeCount_;
            break;
         default:  // Explicit ions.
            if (qCode_ >= 10) {
               nIons_ = qCode_ - 10;
               if (nIons_ % 2 != 0)
                  UTIL_THROW("Number of ions uneven");
               if (nPolymers_ % 2 != 0)
                  UTIL_THROW("Number of polymers uneven");
               if (nType_ < 5)
                  UTIL_THROW("Number of bead type (<5) contradict with charge code");
               nSolvents_ += nIons_;
            } else if (qCode_ < 0) {
               if (-qCode_ % 2 != 0)
                  UTIL_THROW("Number of salt beads is not even");
               if (nType_ < 4)
                  UTIL_THROW("Number of bead type contradict with charge code");
               if (nNeutralSolvents_ > 0 && nType_ < 5)
                  UTIL_THROW("Number of bead type contradict with charge code");
               nIons_ = -qCode_;
               nSolvents_ += nIons_;
            } else {
               if (qCode_ > 0)
                  UTIL_THROW("invalid charge code");
            }
      }
     
     // cout << chargeCount_ << "  " << nPolymerBeads_ << "   " <<nSolvents_ << endl; 
      // Allocate bead array and set molecule pointers.
      beads_.resize(nPolymerBeads_ + nSolvents_);
      for (i = 0; i < nPolymers_; ++i) polymer_.push_back(&(beads_[i*nAB_]));
      solvent_ = &(beads_[nPolymerBeads_]);

      // Assign default types.
      for (i = 0; i < nPolymers_; ++i) {
         ptr = polymer_[i];
         for (j = 0; j < nA_; ++j)
            ptr[j].t = 0;
         for (j = nA_; j < nAB_; ++j)
            ptr[j].t = 1;
      }
      for (i = 0; i < nSolvents_; ++i)
         solvent_[i].t = 2;

      // Assign charges.
      switch (qCode_) {

         case 0: // Charge neutral.

            Log::file() << "Bead charges are set to zero." << endl;
            for (i = 0; i < nPolymerBeads_; ++i) beads_[i].q = 0.0;
            for (i = 0; i < nSolvents_; ++i) solvent_[i].q = 0.0;
            break;

         case 1: // Use all solvents as counter ions to A block.

            Log::file() << "Create counter ion charges for A block." << endl;
            qA = -1.0;
            qB = 0.0;

            for (i = 0; i < nPolymers_; ++i) {
               ptr = polymer_[i];
               for (j = 0; j < nA_; ++j)
                  ptr[j].q = qA;
               for (j = nA_; j < nAB_; ++j)
                  ptr[j].q = qB;
            }

            for(i = 0; i < nA_*nPolymers_; ++i)
               solvent_[i].q = -qA;
            for(i = nA_*nPolymers_; i < nSolvents_; ++i)
               solvent_[i].q = 0.0;

            break;

         // The first half of chain or solvents qA; the scond half of chain or solvent qB.
         case 2:

            Log::file() << "Create counter ion charges for polymercations and polyanions." << endl;

            qA = -1.0;
            qB =  1.0;

            for (i = 0; i < nPolymerBeads_/2; ++i) {
               beads_[i].q = qA;
               solvent_[i].q = -qA;
            }

            for (i = nPolymerBeads_/2; i < nPolymerBeads_; ++i) {
               beads_[i].q = qB;
               solvent_[i].q = -qB;
            }

            for (i = nPolymerBeads_; i < nSolvents_; ++i)
               solvent_[i].q = 0.0;

            break;

         // Half of A block is charged, with counter ions.
         case 3:

            int nAcharged;
            nAcharged = ( nA_ % 2 == 0 ? nA_/2 : (nA_+1)/2 );

            Log::file() << "Create counter ion charges for polymercations and polyanions." << endl;
            qA = -1.0;
            qB = 0.0;

            for (i = 0; i < nPolymers_; ++i) {
               ptr = polymer_[i];
               for (j = 0; j < nA_; j += 2)
                  ptr[j].q = qA;
               for (j = 1; j < nA_; j += 2 )
                  ptr[j].q = 0.0;
               for (j = nA_; j < nAB_; ++j)
                  ptr[j].q = qB;
            }

            for(i = 0; i < nPolymers_*nAcharged; ++i)
               solvent_[i].q = -qA;
            for(i = nPolymers_*nAcharged; i < nSolvents_; ++i)
               solvent_[i].q = 0.0;

            break;

         // Polymers are either polycations or polyanions. Solvents are neutral.
         case 4:

            Log::file() << "Assign charges to polycations and polyanions." << endl;

            for (i = 0; i < nPolymers_; ++i) {
               qA = (i < nPolymers_/2 ? 1.0 : -1.0);
               ptr = polymer_[i];
               for (j = 0; j < nAB_; ++j) ptr[j].q = qA;
            }

            break;

         // Split the neutral solvents into two groups.
         case 5:

            Log::file() << "Split the neutral solvents into two groups." << endl;

            for (i = 0; i < nSolvents_ / 2; ++i) {
               solvent_[i].t = 2;
               solvent_[i].q = 0.0;
            }
            for (i = nSolvents_ / 2; i < nSolvents_; ++i) {
               solvent_[i].t = 3;
               solvent_[i].q = 0.0;
            }

            break;

         // Assign particle types and charges for "counter ion + alternating charged polymer".
         case 6:

            Log::file() << "Assign charges and bead types to counterions and polymer beads." << endl;

            int nChargeBlock, neutralBlockLength, nHead, nTail;

            // No. of charged blocks per chain.
            nChargeBlock = (micelleNIons_ * micelleIonQ_) / (micellePolymerNQ_ * nPolymers_);

            // No. of beads per neutral blocks.
            neutralBlockLength = nAB_ - micellePolymerNQ_ * nChargeBlock;
            if (neutralBlockLength % (nChargeBlock + 1) != 0)
               UTIL_THROW("Neutral block beads cannot be evently distributed");
            neutralBlockLength /= (nChargeBlock + 1);

            // Number of beads in the head and tail (neutral) blocks.
            nHead = neutralBlockLength / 2;
            nTail = neutralBlockLength - nHead;

            for (i = 0; i < nPolymers_; ++i) {
               ptr = polymer_[i];
               for (j = 0; j < nAB_; ++j) {
                  ptr[j].t = 0;

                  if (j < nHead || j >= nAB_ - nTail)
                     ptr[j].q = 0.0;
                  else if ( (j - nHead) % (micellePolymerNQ_ + neutralBlockLength) < micellePolymerNQ_)
                     ptr[j].q = (micelleIonQ_ > 0 ? -1.0 : 1.0);
                  else
                     ptr[j].q = 0.0;
               }
            }

            for (i = 0; i < nSolvents_; ++i) {
               solvent_[i].t = 2;
               if (i < nNeutralSolvents_)
                  solvent_[i].q = 0.0;
               else
                  solvent_[i].q = double(micelleIonQ_);
            }

            break;

         // Alternatingly charged polymers.
         case 7:

            Log::file() << "Assign charges to alternatingly charged polymers." << endl;

            for (i = 0; i < nPolymers_; ++i) {
               qA = (i < nPolymers_/2 ? 1.0 : -1.0);
               ptr = polymer_[i];
               for (j = 1; j < nAB_; j += 2) ptr[j].q = qA;
            }
            break;
         // Charged polymer
         case 8:

            Log::file() << "Assign charges to charged polymers." << endl;
            
            k = 0;
            for (i = 0; i < nPolymers_; ++i) {
              ptr = polymer_[i];
              for (j = 0; j < nA_; ++j ){
                ptr[j].q = 0.0;
              }
              for (j = nA_; j < nA_ + nB_; ++j) {
                ptr[j].q = chargeDistribution_[k];
            //    Log::file() << ptr[j].q;
                k++;
              }

            }
            Log::file() << endl;
            for (i = 0; i < chargeCount_; i++) {
              solvent_[i].q = 1.0;
            }
            for (i = chargeCount_; i < nFreeIons_ / 2; i++) {
              solvent_[i].t = 3;
              solvent_[i].q = -1.0;
            }
            for (i = chargeCount_ + nFreeIons_ / 2; i < chargeCount_ + nFreeIons_; i++) {
              solvent_[i].t = 2;
              solvent_[i].q = 1.0;
            }

            for ( i = chargeCount_ + nFreeIons_; i < nSolvents_; i++) {
              solvent_[i].t = 4;
              solvent_[i].q = 0.0;
            }
            break;

         default: // Explicitly specified ions qCode_: > 10 or < 0.

            Log::file() << "Assign charges to ions and polyelectrolytes." << endl;

            for (i = 0; i < nIons_ / 2; ++i) {
               solvent_[i].t = 2;
               solvent_[i].q = 1.0;
            }
            for (i = nIons_ / 2; i < nIons_; ++i) {
               solvent_[i].t = 3;
               solvent_[i].q = -1.0;
            }
            for (i = nIons_; i < nSolvents_; ++i) {
               solvent_[i].t = 4;
               solvent_[i].q = 0.0;
            }

            if (qCode_ >= 10) {

               for (i = 0; i < nPolymers_; ++i) {
                  qA = (i < nPolymers_/2 ? 1.0 : -1.0);
                  ptr = polymer_[i];
                  for (j = 0; j < nAB_; ++j) ptr[j].q = qA;
               }

            }
      }
 
   }

   /*
   * Project bead position and charge onto grid.
   */
   void System::makeGridDensity()
   {
      Particle *ptr = &(beads_[0]);
      for (int i = 0; i < nPolymerBeads_ + nSolvents_; ++i) {
         grid_.insertBead(*ptr);
         if (fabs(ptr->q) > Constants::Epsilon) {
           grid_.fillNeighborList(i, ptr->r, 1);
         }
         ++ptr;
      }
   }

   /*
   * Generate system configuration using random walks.
   */
   void System::generateRandomConfig()
   {
      int       iChain, iBead, idim;
      double    R1, R2, R3;
      Vector    r0, bond;
      Particle  *ptr;

      // Generate polymer bead positions.
      for (iChain = 0; iChain < nPolymers_; ++iChain) {
         ptr = &beads_[iChain*nAB_];

         // The first bead.
         iBead = 0;
         for (idim = 0; idim < Dimension; ++idim)
            r0[idim] = boxL_[idim] * random_.Random();
         toPrimaryCell(r0);
         ptr[iBead].r = r0;

         // The other n - 1 beads.
         for (iBead = 1; iBead < nAB_; ++iBead) {
            do { // Generate a randomly oriented vector.
               R1 = 2.0 * random_.Random() - 1.0;
               R2 = 2.0 * random_.Random() - 1.0;
               R3 = R1 * R1 + R2 * R2;
            } while (R3 >= 1.0);

            bond[0] = 2.0 * sqrt(1.0 - R3) * R1 * bondL_;
            bond[1] = 2.0 * sqrt(1.0 - R3) * R2 * bondL_;
            bond[2] = (1.0 - 2.0 * R3) * bondL_;

            r0 += bond;
            toPrimaryCell(r0);
            ptr[iBead].r = r0;
         }
      }

      // Generate solvent bead positions.
      ptr = &beads_[nPolymerBeads_];
      for (iBead = 0; iBead < nSolvents_; ++iBead) {
         for (idim = 0; idim < Dimension; ++idim)
            r0[idim] = boxL_[idim] * random_.Random();
         toPrimaryCell(r0);
         ptr[iBead].r = r0;
      }
   }

   /*
   * Read system configuration (bead position, type, et al.).
   */
   void System::readConfig(istream& in)
   {
      char      comment[200];
      string    line;
      int       nTotal, nPolymers, nAB, nSolvents, i, type;
      Particle  *ptr;

      // Parse config header.
      getline(in, line);
      sscanf(line.c_str(), "%d %s", &nTotal, comment);
      if (nTotal != nPolymerBeads_ + nSolvents_)
         UTIL_THROW("config file error: bead total number mismatch");

      getline(in, line);
      sscanf(line.c_str(), "%d  %d  %d %s", &nPolymers, &nAB, &nSolvents, comment);
      if (nPolymers != nPolymers_)
         UTIL_THROW("config file error: nPolymers mismatch");
      if (nAB != nAB_)
         UTIL_THROW("config file error: nAB mismatch");
      if (nSolvents != nSolvents_)
         UTIL_THROW("config file error: nSolvents mismatch");
      if (nTotal != nPolymers * nAB + nSolvents)
         UTIL_THROW("config file error: corrupt config file");

      // Read polymer and solvent bead position.
      ptr = &(beads_[0]);
      for (i = 0; i < nPolymerBeads_ + nSolvents_; ++i) {
         getline(in, line);
         sscanf(line.c_str(), "%d  %lf  %lf  %lf  %lf  %s", &type, &(ptr->r[0]), &(ptr->r[1]), &(ptr->r[2]), &(ptr->q), comment);
         if (ptr->t != type) {
            Log::file() << "bead id:            " << i << endl;
            Log::file() << "bead type:          " << ptr->t << endl;
            Log::file() << "bead type in file:  " << type << endl;
            UTIL_THROW("config file error: bead type mismatch");
         }
         toPrimaryCell(ptr->r);
         ++ptr;
      }
   }

   /*
   * Write system configuration (bead position, type, et al.).
   */
   void System::writeConfig(ostream& out) const
   {
      // Output config header.
      out << nPolymerBeads_ + nSolvents_ << endl;
      out << nPolymers_ << "  " << nAB_ << "  " << nSolvents_ << "   nPolymer nPolymerLength nSolvent" << endl;

      // Save and reset precision.
      streamsize oldprec = out.precision(10);
      out.setf(ios::fixed, ios::floatfield);

      // Write polymer bead positions.
      vector< Vector >  rChain(nAB_, Vector(0.0));

      int      i, j, k;
      const Particle *ptr;
      Vector   dr, rcm, rshift;

      for (i = 0; i < nPolymers_; ++i) {
         ptr = &(beads_[i*nAB_]);

         rChain[0] = ptr[0].r;
         rcm = rChain[0];
         for (j = 1; j < nAB_; ++j) {
            dr.subtract(ptr[j].r, ptr[j-1].r);
            pbcShift(dr);
            rChain[j].add(rChain[j-1], dr);
            rcm += rChain[j];
         }
         rcm /= double(nAB_);
         for (k = 0; k < Dimension; ++k)
            rshift[k] = floor(rcm[k] / boxL_[k]) * boxL_[k];

         for (j = 0; j < nAB_; ++j) {
            out.setf(ios::left, ios::adjustfield);
            out << setw(2) << ptr[j].t;

            rChain[j] -= rshift;
            out.unsetf(ios::adjustfield);
            out << setw(15) << rChain[j][0];
            out << setw(15) << rChain[j][1];
            out << setw(15) << rChain[j][2];
            out << setw(15) << ptr[j].q << endl;
         }
      }

      // Write solvent position.
      ptr = &(beads_[nPolymerBeads_]);
      for (i = 0; i < nSolvents_; ++i) {
         out.setf(ios::left, ios::adjustfield);
         out << setw(2) << ptr->t;
         out.unsetf(ios::adjustfield);
         out << setw(15) << ptr->r[0];
         out << setw(15) << ptr->r[1];
         out << setw(15) << ptr->r[2];
         out << setw(15) << ptr->q << endl;
         ++ptr;
      }

      // Restore I/O formats.
      out.unsetf(ios::floatfield);
      out.precision(oldprec);

   }

   /*
   * Initialize VTK array.
   */
   void System::initializeVTK()
   { grid_.initializeVTK(); }

   /*
   * Update VTK array.
   */
   void System::updateVTK(const int iStep, const int flag)
   {
      grid_.updateVTK();

      if (flag > 0) {
         std::ofstream outconfig;
         std::string name("config.");
         std::ostringstream s;
         s << iStep;
         name.append(s.str()); 
         fileMasterPtr_->openOutputFile(name.c_str(), outconfig);
         writeConfig(outconfig);
         outconfig.close(); 
      }
   }

   /*
   * Output VTK array.
   */
   void System::writeVTK(ostream& out) const
   { grid_.writeVTK(out); }

   /*
   * Output VTK array.
   */
   void System::solvePoisson()
   {
      grid_.solvePoisson();
      grid_.resetVTK();
   }

}

#endif

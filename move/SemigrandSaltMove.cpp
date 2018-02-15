#ifndef SEMIGRANDSALTMOVE_CPP
#define SEMIGRANDSALTMOVE_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <math.h>
#include "../util/global.h"
#include "../util/Vector.h"

#include "SemigrandSaltMove.h"
#include "../simulation/McSystem.h"

namespace GridMC
{
   using namespace Util;

   /*
   * Constructor.
   */
   SemigrandSaltMove::SemigrandSaltMove(McSystem& sysIn) :
      Move(sysIn),
      dMu_(0.0) 
   { name_.assign("SemigrandSaltMove"); }

   /*
   * Default destructor.
   */
   SemigrandSaltMove::~SemigrandSaltMove()
   {}

   /*
   * Read chemical potential difference: mu(salt) - mu(solvent).
   */
   void SemigrandSaltMove::readParam(std::istream& in)
   {
      char   comment[200];
      std::string line;

      getline(in, line);
      if (line.size() > 0) {
         sscanf(line.c_str(), "%lf %s", &dMu_, comment);
      } else {
         UTIL_THROW("reading error: chemical potential difference (salt - solvent)");
      }
   }

   /*
   * Write parameters, pairing with readParam.
   */
   void SemigrandSaltMove::writeParam(std::ostream& out) const
   {
      out << "mu(salt) - mu(solvent)      " << dMu_ << std::endl;
   }

   /*
   * Semi-grand solvent-salt switching move.
   *
   * Assumption: solvent beads and salt beads are identical as far as two-body interactions are considered.
   */
   void SemigrandSaltMove::move()
   {
      int       nNeutralSolvents(system_.nNeutralSolvents_);
      int       nIons(system_.nIons_), nCations(nIons / 2);
      int       nSolvents(nNeutralSolvents + nIons);
      int       id, idCation, idAnion, neutral1, neutral2;
      double    deltaE, semigrandWeight;
      Particle  *solvents(system_.solvent_);
      Vector    rCation, rAnion;

      // Choose a bead (neutral solvent or ion) at random.
      id = system_.random_.IRandom(0, nSolvents - 1);

      if (id >= nIons) { // If it's a solvent; select another solvent and propose a solvent -> ion-pair move.
         neutral1 = id;   
         neutral2 = system_.random_.IRandom(nIons, nSolvents - 2);
         if (neutral2 >= id) neutral2 += 1;
         semigrandWeight = exp(dMu_) * double(nNeutralSolvents) / double(nIons + 2);

         idCation = nCations;  // Index of the new cation.
         idAnion = nIons + 1; // Index of the new anion.

         rCation = solvents[neutral1].r;
         rAnion  = solvents[neutral2].r;
         system_.grid_.findAffectedSites(rAnion, rCation);
      } else { // If it's an ion; select an ion of different type, then propose a ion-pair -> solvent move.
         if (id < nCations) {
            idCation = id;
            idAnion = system_.random_.IRandom(nCations, nIons - 1);
         } else {
            idAnion = id;
            idCation = system_.random_.IRandom(0, nCations - 1);
         }
         semigrandWeight = exp(-dMu_) * double(nIons) / double(nNeutralSolvents + 2);

         neutral1 = nIons - 1; // Index of the new neutral solvent.
         neutral2 = nIons; // Index of the new neutral solvent.

         rCation = solvents[idCation].r;
         rAnion  = solvents[idAnion].r;
         system_.grid_.findAffectedSites(rCation, rAnion);
      }
      deltaE = system_.grid_.getCoulombEnergyChange(1.0);

      // Metropolis test.
      ++nAttempt_;

      if (system_.random_.Random() <= exp(-deltaE) * semigrandWeight) { 

         // Update bead indices, nIons, and nNeutralSolvents.
         if (id >= nIons) {

            // Save positions of the first two neutral beads: algorithm still OK when neutral- 1 or 2 equals 0 or 1.
            solvents[neutral1].r = solvents[nIons].r;
            solvents[neutral2].r = solvents[nIons+1].r;

            // Mutate the first two vacant neutral beads into anions: created anion and the original first anion.
            solvents[nIons].r = rAnion;
            solvents[nIons].q = -1.0;
            solvents[nIons].t = 3;
            solvents[nIons+1].r = solvents[nCations].r;
            solvents[nIons+1].q = -1.0;
            solvents[nIons+1].t = 3;

            // Mutate the original first anion into cation.
            solvents[nCations].r = rCation;
            solvents[nCations].q = 1.0;
            solvents[nCations].t = 2;

            // Update bead count.
            system_.nIons_ += 2;
            system_.nNeutralSolvents_ -= 2;

         } else {

            // Save position of the last cation: algorithm still OK when idCation = nCations-1.
            solvents[idCation].r = solvents[nCations-1].r;

            // Mutate the last two anions into the last cation and the deleted anion.
            if (idAnion < nIons - 2) {
               solvents[nCations - 1].r = solvents[nIons - 1].r;
               solvents[nCations - 1].q = -1.0;
               solvents[nCations - 1].t = 3;
               solvents[idAnion].r = solvents[nIons - 2].r;
            } else {
               if (idAnion == nIons - 2) {
                  solvents[nCations - 1].r = solvents[nIons - 1].r;
                  solvents[nCations - 1].q = -1.0;
                  solvents[nCations - 1].t = 3;
               } else {
                  solvents[nCations - 1].r = solvents[nIons - 2].r;
                  solvents[nCations - 1].q = -1.0;
                  solvents[nCations - 1].t = 3;
               }
            }

            // Change the last two anions to neutral beads.
            solvents[nIons - 2].r = rCation;
            solvents[nIons - 2].q = 0.0;
            solvents[nIons - 2].t = 4;
            solvents[nIons - 1].r = rAnion;
            solvents[nIons - 1].q = 0.0;
            solvents[nIons - 1].t = 4;

            // Update beads count.
            system_.nIons_ -= 2;
            system_.nNeutralSolvents_ += 2;

         }

         system_.grid_.updateChargeGridFromList(1.0);
         system_.updateCoulombEnergy(deltaE);

         ++nAccept_;
      }

      //if (system_.doDOS())
      //   system_.updateDOS();

   }

}

#endif

#ifndef CHAINCHEMICALPOTENTIALINSERTION_CPP
#define CHAINCHEMICALPOTENTIALINSERTION_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iomanip>
#include <map>
#include "ChainChemicalPotentialInsertion.h"

namespace GridMC
{

   /*
   * Constructor.
   */
   ChainChemicalPotentialInsertion::ChainChemicalPotentialInsertion(System& system) :
      Diagnosis(system),
      f1_(0.0),
      f2_(0.0),
      dataFile_(NULL)
   { name_.assign("ChainChemicalPotentialInsertion"); }

   /*
   * Default destructor.
   */
   ChainChemicalPotentialInsertion::~ChainChemicalPotentialInsertion()
   { dataFile_.close(); }

   /*
   * Read sampling interval.
   */
   void ChainChemicalPotentialInsertion::readParam(std::istream& in)
   {
      Diagnosis::readParam(in);
      fileMaster_.openOutputFile("dEChainInsert.dat", dataFile_);
   }

   /*
   * Write parameter.
   */
   void ChainChemicalPotentialInsertion::writeParam(std::ostream& out) const
   {
      Diagnosis::writeParam(out);
   }

   /*
   * Sample.
   */
   void ChainChemicalPotentialInsertion::sample(const long iStep)
   {
      if (isAtInterval(iStep) == false) return;

      int       nPolymers(system_.nPolymers_);
      int       nPolycations(nPolymers/2);
      int       nAB(system_.nAB_);
      int       idCation, idAnion, idim;
      Particle  *ptrPC, *ptrPA;  // Pointer to polycation and polyanion
      Vector    r1, r2, bond;
      double    R1, R2;
      double    f(0.0), dEbond(0.0), dENeutral(0.0), dECharge(0.0);
      Particle  widomBead;    // Virtual bead to insert.

      // Pickup a polycation at random.
      idCation = system_.random_.IRandom(0, nPolycations - 1);
      ptrPC = system_.polymer_[idCation]; 

      // Generate a new bond.
      for (idim = 0; idim < 3; ++idim) {
         R1 = system_.random_.Random();
         R2 = system_.random_.Random();
         bond[idim] = sqrt(-2.0 * log(R1)) * cos(2.0 * Constants::Pi * R2) / 3.0;
      }
      dEbond = bond.square();

      if (system_.random_.Random() > 0.5) {
         r1 = ptrPC[nAB-1].r;
      } else {
         r1 = ptrPC[0].r;
      }
      r1 += bond;
      system_.pbcShift(r1);

      // Pickup a polyanion at random.
      idAnion = system_.random_.IRandom(nPolycations, nPolymers - 1);
      ptrPA = system_.polymer_[idAnion]; 

      // Generate the second new bond.
      for (idim = 0; idim < 3; ++idim) {
         R1 = system_.random_.Random();
         R2 = system_.random_.Random();
         bond[idim] = sqrt(-2.0 * log(R1)) * sin(2.0 * Constants::Pi * R2) / 3.0;
      }
      dEbond += bond.square();

      if (system_.random_.Random() > 0.5) {
         r2 = ptrPA[nAB-1].r;
      } else {
         r2 = ptrPA[0].r;
      }
      r2 += bond;
      system_.pbcShift(r2);

      // Energy change due to the addition of two new bonds.
      dEbond *= 1.5;

      // Get the neutral energy change.
      grid_.findAffectedSites(r1);
      dENeutral = grid_.getTwobodyEnergyChange(0);

      // Update the mass field to get the neutral energy chagne properly.
      widomBead.t = 0;
      widomBead.r = r1;
      widomBead.q = 0.0;
      grid_.updateFromList(widomBead);

      grid_.findAffectedSites(r2);
      dENeutral += grid_.getTwobodyEnergyChange(0);

      // Restore the mass field.
      grid_.findAffectedSites(r1, -2);
      grid_.updateFromList(widomBead);

      // Get the charge energy change.
      grid_.findAffectedSites(r2, r1);
      dECharge = grid_.getCoulombEnergyChange(1.0);

      // Output data and accumulate the statistics.
      //f = exp(- (dEbond + dENeutral + dECharge) );
      f = dEbond + dENeutral + dECharge;

      streamsize oldprec = dataFile_.precision(6);
      dataFile_.setf(ios::fixed, ios::floatfield);
      dataFile_ << std::setw(15) << f << std::endl;
      dataFile_.precision(oldprec);
      dataFile_.unsetf(ios::floatfield);

      f1_ += f;
      f2_ += f * f;
      nSamples_ += 1;
   }


   /*
   * Output result.
   */
   void ChainChemicalPotentialInsertion::output() const
   {
      double ave = f1_ / double(nSamples_);
      double var = f2_ / double(nSamples_) - ave*ave;

      std::ofstream out;
      fileMaster_.openOutputFile("dEChainInsert.out", out);

      streamsize oldprec = out.precision(6);
      out.setf(ios::fixed, ios::floatfield);
      out.setf(ios::left, ios::adjustfield);

      out << setw(25) << "Average: ";
      out << setw(15) << ave << std::endl;
      out << setw(25) << "Standard deviation: ";
      out << setw(15) << sqrt(var) << std::endl;

      out.precision(oldprec);
      out.setf(ios::right, ios::adjustfield);
      out.unsetf(ios::floatfield);

      out.close();
   }

}

#endif

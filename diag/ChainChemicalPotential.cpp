#ifndef CHAINCHEMICALPOTENTIAL_CPP
#define CHAINCHEMICALPOTENTIAL_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iomanip>
#include <map>
#include "ChainChemicalPotential.h"

namespace GridMC
{

   /*
   * Constructor.
   */
   ChainChemicalPotential::ChainChemicalPotential(System& system) :
      Diagnosis(system),
      f1_(0.0),
      f2_(0.0),
      dataFile_(NULL)
   { name_.assign("ChainChemicalPotential"); }

   /*
   * Default destructor.
   */
   ChainChemicalPotential::~ChainChemicalPotential()
   { dataFile_.close(); }

   /*
   * Read sampling interval.
   */
   void ChainChemicalPotential::readParam(std::istream& in)
   {
      Diagnosis::readParam(in);
      fileMaster_.openOutputFile("dEChain.dat", dataFile_);
   }

   /*
   * Write parameter.
   */
   void ChainChemicalPotential::writeParam(std::ostream& out) const
   {
      Diagnosis::writeParam(out);
   }

   /*
   * Sample.
   */
   void ChainChemicalPotential::sample(const long iStep)
   {
      if (isAtInterval(iStep) == false) return;

      int       nPolymers(system_.nPolymers_);
      int       nPolycations(nPolymers/2);
      int       nAB(system_.nAB_);
      int       idCation, idAnion;
      Particle  *ptrPC, *ptrPA;  // Pointer to polycation and polyanion
      Vector    r1, r2, bond;
      double    f(0.0), dEbond(0.0), dENeutral(0.0), dECharge(0.0);
      Particle  widomBead;    // Virtual bead to insert.

      // Pickup a polycation at random.
      idCation = system_.random_.IRandom(0, nPolycations - 1);
      ptrPC = system_.polymer_[idCation]; 
      if (system_.random_.Random() > 0.5) {
         r1 = ptrPC[nAB-1].r;
         bond.subtract(r1, ptrPC[nAB-2].r);
      } else {
         r1 = ptrPC[0].r;
         bond.subtract(r1, ptrPC[1].r);
      }
      system_.pbcShift(bond);
      dEbond = bond.square();

      // Pickup a polyanion at random.
      idAnion = system_.random_.IRandom(nPolycations, nPolymers - 1);
      ptrPA = system_.polymer_[idAnion]; 
      if (system_.random_.Random() > 0.5) {
         r2 = ptrPA[nAB-1].r;
         bond.subtract(r2, ptrPA[nAB-2].r);
      } else {
         r2 = ptrPA[0].r;
         bond.subtract(r2, ptrPA[1].r);
      }
      system_.pbcShift(bond);
      dEbond += bond.square();

      // Energy change due to the absence of two new bonds.
      dEbond *= -(1.5 / system_.bondL_ / system_.bondL_);

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

      // Convert to mass deletion convention.
      dENeutral *= -1.0;

      // Get the charge energy change.
      grid_.findAffectedSites(r1, r2);
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
   void ChainChemicalPotential::output() const
   {
      double ave = f1_ / double(nSamples_);
      double var = f2_ / double(nSamples_) - ave*ave;

      std::ofstream out;
      fileMaster_.openOutputFile("dEChain.out", out);

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

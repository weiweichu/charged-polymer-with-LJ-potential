#ifndef SALTCHEMICALPOTENTIAL_CPP
#define SALTCHEMICALPOTENTIAL_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iomanip>
#include <map>
#include "SaltChemicalPotential.h"

namespace GridMC
{

   /*
   * Constructor.
   */
   SaltChemicalPotential::SaltChemicalPotential(System& system) :
      Diagnosis(system),
      type_(-1),  // Invalid type ID
      f1_(0.0),
      f2_(0.0),
      dataFile_(NULL)
   { name_.assign("SaltChemicalPotential"); }

   /*
   * Default destructor.
   */
   SaltChemicalPotential::~SaltChemicalPotential()
   { dataFile_.close(); }

   /*
   * Read sampling interval.
   */
   void SaltChemicalPotential::readParam(std::istream& in)
   {
      Diagnosis::readParam(in);

      // Read number of wave and wave indices.
      char        comment[200];
      std::string line;

      // Maximum wave index.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: species type");
      sscanf(line.c_str(), "%d %s", &type_, comment);
      f1_ = 0.0;
      f2_ = 0.0;

      fileMaster_.openOutputFile("dESalt.dat", dataFile_);
   }

   /*
   * Write parameter.
   */
   void SaltChemicalPotential::writeParam(std::ostream& out) const
   {
      Diagnosis::writeParam(out);
      out << "Species type                " << type_ << std::endl;
   }

   /*
   * Sample.
   */
   void SaltChemicalPotential::sample(const long iStep)
   {
      if (isAtInterval(iStep) == false) return;

      double    f(0.0), dENeutral(0.0), dECharge(0.0);
      Vector    r1, r2;
      Particle  widomBead;    // Virtual bead to insert.

      // Generate trial bead position.
      for (int j = 0; j < Dimension; ++j) {
         r1[j] = boxL_[j] * system_.getRandom().Random();
         r2[j] = boxL_[j] * system_.getRandom().Random();
      }
      system_.toPrimaryCell(r1);
      system_.toPrimaryCell(r2);

      // Get the neutral energy change.
      grid_.findAffectedSites(r1);
      dENeutral = grid_.getTwobodyEnergyChange(type_);

      // Update the mass field to get the neutral energy chagne properly.
      widomBead.t = type_;
      widomBead.r = r1;
      widomBead.q = 0.0;
      grid_.updateFromList(widomBead);

      grid_.findAffectedSites(r2);
      dENeutral += grid_.getTwobodyEnergyChange(type_);

      // Restore the mass field.
      grid_.findAffectedSites(r1, -2);
      grid_.updateFromList(widomBead);

      // Get the charge energy change.
      grid_.findAffectedSites(r2, r1);
      dECharge = grid_.getCoulombEnergyChange(1.0);

      // Output data and accumulate the statistics.
      //f = exp(- (dENeutral + dECharge) );
      f = dENeutral + dECharge;

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
   void SaltChemicalPotential::output() const
   {
      double ave = f1_ / double(nSamples_);
      double var = f2_ / double(nSamples_) - ave*ave;

      std::ofstream out;
      fileMaster_.openOutputFile("dESalt.out", out);

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

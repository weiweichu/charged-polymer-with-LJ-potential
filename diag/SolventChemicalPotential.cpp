#ifndef SOLVENTCHEMICALPOTENTIAL_CPP
#define SOLVENTCHEMICALPOTENTIAL_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iomanip>
#include <map>
#include "SolventChemicalPotential.h"

namespace GridMC
{

   /*
   * Constructor.
   */
   SolventChemicalPotential::SolventChemicalPotential(System& system) :
      Diagnosis(system),
      type_(-1),  // Invalid type ID
      mu1_(0.0),
      mu2_(0.0),
      dataFile_(NULL)
   { name_.assign("SolventChemicalPotential"); }

   /*
   * Default destructor.
   */
   SolventChemicalPotential::~SolventChemicalPotential()
   { dataFile_.close(); }

   /*
   * Read sampling interval.
   */
   void SolventChemicalPotential::readParam(std::istream& in)
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
      mu1_ = 0.0;
      mu2_ = 0.0;

      fileMaster_.openOutputFile("dESolvent.dat", dataFile_);
   }

   /*
   * Write parameter.
   */
   void SolventChemicalPotential::writeParam(std::ostream& out) const
   {
      Diagnosis::writeParam(out);
      out << "Species type                " << type_ << std::endl;
   }

   /*
   * Sample.
   */
   void SolventChemicalPotential::sample(const long iStep)
   {
      if (isAtInterval(iStep) == false) return;

      double mu(0.0), deltaE;
      Vector rTrial;

      // Generate trial bead position.
      for (int j = 0; j < Dimension; ++j)
         rTrial[j] = boxL_[j] * system_.getRandom().Random();
      system_.toPrimaryCell(rTrial);

      // Get the system energy change.
      grid_.findAffectedSites(rTrial);
      deltaE = grid_.getTwobodyEnergyChange(type_);

      //mu = exp(-deltaE);
      mu = deltaE;

      // Output data and accumulate the statistics.
      streamsize oldprec = dataFile_.precision(6);
      dataFile_.setf(ios::fixed, ios::floatfield);
      dataFile_ << std::setw(15) << mu << std::endl;
      dataFile_.precision(oldprec);
      dataFile_.unsetf(ios::floatfield);

      mu1_ += mu;
      mu2_ += mu * mu;
      nSamples_ += 1;
   }


   /*
   * Output result.
   */
   void SolventChemicalPotential::output() const
   {
      double ave = mu1_ / double(nSamples_);
      double var = mu2_ / double(nSamples_) - ave*ave;

      std::ofstream out;
      fileMaster_.openOutputFile("dESolvent.out", out);

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

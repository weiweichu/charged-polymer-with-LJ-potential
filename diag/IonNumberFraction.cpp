#ifndef IONNUMBERFRACTION_CPP
#define IONNUMBERFRACTION_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iomanip>
#include <map>
#include "IonNumberFraction.h"

namespace GridMC
{

   /*
   * Constructor.
   */
   IonNumberFraction::IonNumberFraction(System& system) :
      Diagnosis(system),
      psi1_(0.0),
      psi2_(0.0),
      dataFile_(NULL)
   { name_.assign("IonNumberFraction"); }

   /*
   * Default destructor.
   */
   IonNumberFraction::~IonNumberFraction()
   { dataFile_.close(); }

   /*
   * Read sampling interval.
   */
   void IonNumberFraction::readParam(std::istream& in)
   {
      Diagnosis::readParam(in);
      fileMaster_.openOutputFile("psi.dat", dataFile_);
   }

   /*
   * Write parameter.
   */
   void IonNumberFraction::writeParam(std::ostream& out) const
   {
      Diagnosis::writeParam(out);
   }

   /*
   * Sample.
   */
   void IonNumberFraction::sample(const long iStep)
   {
      if (isAtInterval(iStep) == false) return;

      double psi = double(system_.nIons_) / double(system_.nIons_ + system_.nNeutralSolvents_);

      // Output data and accumulate the statistics.
      streamsize oldprec = dataFile_.precision(6);
      dataFile_.setf(ios::fixed, ios::floatfield);
      dataFile_ << std::setw(15) << psi << std::endl;
      dataFile_.precision(oldprec);
      dataFile_.unsetf(ios::floatfield);

      psi1_ += psi;
      psi2_ += psi * psi;
      nSamples_ += 1;
   }


   /*
   * Output result.
   */
   void IonNumberFraction::output() const
   {
      double ave = psi1_ / double(nSamples_);
      double var = psi2_ / double(nSamples_) - ave * ave;

      std::ofstream out;
      fileMaster_.openOutputFile("psi.out", out);

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

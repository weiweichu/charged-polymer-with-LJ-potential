#ifndef HISTOGRAMWEIGHT_CPP
#define HISTOGRAMWEIGHT_CPP

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <math.h>
#include <iomanip>
#include "HistogramWeight.h"

namespace GridMC
{

   using namespace std;
   using namespace Util;

   /*
   * Default constructor.
   */
   HistogramWeight::HistogramWeight() :
      nIteration_(0),
      deltaW_(1.0),
      weight_(),
      nIterationTol_(25),
      sTolerance_(0.5),
      nMaxIteration_(40)
   {}

   /*
   * Default destructor.
   */
   HistogramWeight::~HistogramWeight()
   {}

   /*
   * Read parameters.
   */
   void HistogramWeight::readParam(istream& in)
   {
      Histogram::readParam(in);

      nIteration_ = 1;
      deltaW_ = 1.0;
      weight_.resize(nBin_, 0.0);

      // Iteration parameters. Algorithm parameter.
      char   comment[200];
      string line;

      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: histogram iterating parameters");
      sscanf(line.c_str(), "%d %lf %d %s", &nIterationTol_, &sTolerance_, &nMaxIteration_, comment);
   }

   /*
   * Write parameters.
   */
   void HistogramWeight::writeParam(ostream& out) const
   {
      Histogram::writeParam(out);
      out << "iteration parameters        ";
      out << nIterationTol_ << "  "  << sTolerance_ << "  " << nMaxIteration_ << endl;
   }

   /*
   * Adopt a newly sample value.
   */
   void HistogramWeight::sample(const double value, const int intervalIn)
   {
      if (value > min_ && value < max_) {
         int i = binIndex(value);
         histogram_[i] += 1;
         nSample_ += 1;
 
         if (intervalIn > 0)
            weight_[i] += deltaW_;
      } else {
         Log::file() << "value:  " << value << endl;
         Log::file() << "weight: " << getWeight(value) << endl;
         UTIL_THROW("Invalid sample value");
      }
   }

   /*
   * Read the prescribed bias weight.
   */
   void HistogramWeight::readWeight(istream& in, const int itr)
   {
      string  line;
      double  x;
      char    comment[200];
      for (int i = 0; i < nBin_; ++i) {
         getline(in, line);
         if (line.size() <= 0)
            UTIL_THROW("reading error: histogram weight");
         sscanf(line.c_str(), "%lf %lf %s", &x, &weight_[i], comment);
         weight_[i] *= -1.0;
      }
      deltaW_ = pow(0.5, itr - 1);
      nIteration_ = itr;
   }

   /*
   * Cut the weight step size by half.
   */
   void HistogramWeight::resetHistogram()
   {
      deltaW_ *= 0.5;
      nIteration_ += 1;
      clearHistogram();
      if (nIteration_ > nIterationTol_)
         rescaleTolerance(sTolerance_);
   }

   /*
   * Clear the histogram, decrease deltaW, and shift the weight such that the
   * maximum (free energy minimum) sits at zero.
   */
   void HistogramWeight::shiftWeight()
   {
      double maxW(-1.0E14);
      int    i;
      for (i = 0; i < nBin_; ++i)
         if (weight_[i] > maxW) maxW = weight_[i];

      for (i = 0; i < nBin_; ++i)
         weight_[i] -= maxW;
   }

   /*
   * Output the bias weight.
   */
   void HistogramWeight::outputWeight(ostream& out) const
   {
      // Output.
      out.setf(ios::fixed, ios::floatfield);
      streamsize oldprec = out.precision();
      double x = min_ - 0.5 * binWidth_;
      for (int i = 0; i < nBin_; ++i) {
         x += binWidth_;
         out.setf(ios::left, ios::adjustfield);
         out.precision(5);
         out << setw(10) << x;
         out.unsetf(ios::adjustfield);

         out.precision(8);
         out << -weight_[i] << endl;
      }
      out.precision(oldprec);
      out.unsetf(ios::floatfield);
   }

}
#endif

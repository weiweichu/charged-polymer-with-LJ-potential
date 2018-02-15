#ifndef HISTOGRAM_CPP
#define HISTOGRAM_CPP

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <math.h>
#include <string>
#include <iomanip>
#include "Histogram.h"

namespace GridMC
{

   using namespace std;
   using namespace Util;

   /*
   * Default constructor.
   */
   Histogram::Histogram() :
      min_(0.0),
      max_(0.0),
      binWidth_(0.0),
      nBin_(0),
      nSample_(0),
      histogram_(),
      flatnessTol_(0.0)
   {}

   /*
   * Default destructor.
   */
   Histogram::~Histogram()
   {}

   /*
   * Read parameters.
   */
   void Histogram::readParam(istream& in)
   {
      char   comment[200];
      string line;

      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: range of histogram");
      sscanf(line.c_str(), "%lf %lf %lf %s", &min_, &max_, &binWidth_, comment);

      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: flatness tolerance");
      sscanf(line.c_str(), "%lf %s", &flatnessTol_, comment);

      nBin_ = int((max_ - min_) / binWidth_);
      binWidth_ = (max_ - min_) / double(nBin_);
      histogram_.resize(nBin_, 0);
   }

   /*
   * Write parameters.
   */
   void Histogram::writeParam(ostream& out) const
   {
      out << "range of histogram          ";
      out << min_ << "  "  << max_ << "  " << binWidth_ << endl;

      out << "tolerance of flatness       ";
      out << flatnessTol_ << endl;
   }

   /*
   * Output the histogram normalized by mean value.
   */
   void Histogram::outputHistogram(ostream& out) const
   {
      // Calculate the mean value.
      double mean(0.0);
      for (int i = 0; i < nBin_; ++i)
         mean += double(histogram_[i]);
      mean /= double(nBin_);

      // Output.
      out.setf(ios::fixed, ios::floatfield);
      streamsize oldprec = out.precision(3);
      double x = min_ - 0.5 * binWidth_;
      for (int i = 0; i < nBin_; ++i) {
         x += binWidth_;
         out.setf(ios::left, ios::adjustfield);
         out << setw(8) << x;
         out.unsetf(ios::adjustfield);
         out << setw(6) << double(histogram_[i]) / mean;
         out << endl;
      }
      out.precision(oldprec);
      out.unsetf(ios::floatfield);
   }

   /*
   * Resclae the tolerance factor.
   */
   void Histogram::rescaleTolerance(const double s)
   { flatnessTol_ *= s; }

   /* 
   * Zero all accumulators.
   */
   void Histogram::clearHistogram() 
   {
      for (int i=0; i < nBin_; ++i)
         histogram_[i] = 0;
      nSample_ = 0;
   }

   /* 
   * Add a value to the histogram
   */
   void Histogram::sample(double value)
   {
      if (value > min_ && value < max_) {
         int i = binIndex(value);
         histogram_[i] += 1;
         nSample_ += 1;
      } else {
         Log::file() << "value: " << value << endl;
         UTIL_THROW("Invalid sample value");
      }
   }

   /*
   * Roughness measures the deviation from a flat histogram. What's implemented
   * below is the maximum fractional deviation from the mean.
   */
   double Histogram::getRoughness() const
   {
      double dx, dev(0.0);
      double mean = double(nSample_) / double(nBin_);
      for (int i = 0; i < nBin_; ++i) {
         dx = fabs(double(histogram_[i]) - mean);
         if (dx > dev) dev = dx;
      }
      return (dev / mean);
   }

}
#endif

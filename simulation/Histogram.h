#ifndef HISTOGRAM_H
#define HISTOGRAM_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <vector>
#include <iostream>
#include <util/global.h>

namespace GridMC
{

   /**
   * Maintain the histogram of a double variable and evaluate the weight.
   *
   * \ingroup Simulation_Module
   */
   class Histogram 
   {
   
   public:

      /**
      * Constructor.
      */
      Histogram();

      /**
      * Default destructor.
      */
      virtual ~Histogram();

      /**
      * Read parameters.
      */
      virtual void readParam(std::istream& in);

      /**
      * Write parameters.
      */
      virtual void writeParam(std::ostream& out) const;

      /**
      * Output histogram.
      */
      virtual void outputHistogram(std::ostream& out) const;

      /**
      * Accumulate a new sample value to histogram.
      */
      virtual void sample(const double value);

      /**
      * If the histogram is sufficiently flat.
      */
      bool isFlat() const;

      /**
      * Return the roughness of histogram as measured by the deviation from the mean.
      */
      double getRoughness() const;

      /*
      * Rescale tolerance.
      */
      void rescaleTolerance(const double s = 0.5);

      /**
      * Add an element to histogram.
      */
      void clearHistogram();

      /**
      * Return the index of the bin for a value.
      *
      * \param value sampled value
      */
      int binIndex(double value) const;

   protected:

      /// Histogram parameters and variables.
      double         min_;
      double         max_;
      double         binWidth_;
      int            nBin_;
      unsigned long  nSample_;
      std::vector<unsigned long> histogram_;

      /// Flatness tolerance.
      double         flatnessTol_;

   };

   /*
   * Return the index of the bin for a value.
   */
   inline bool Histogram::isFlat() const
   { return (getRoughness() < flatnessTol_); }

   /*
   * Return the index of the bin for a value.
   */
   inline int Histogram::binIndex(double value) const
   { return int((value - min_)/binWidth_); }

}
#endif

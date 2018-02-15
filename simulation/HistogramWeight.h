#ifndef HISTOGRAMWEIGHT_H
#define HISTOGRAMWEIGHT_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "Histogram.h"

namespace GridMC
{

   /**
   * Maintain the histogram of a double variable and evaluate the weight.
   *
   * \ingroup Simulation_Module
   */
   class HistogramWeight : public Histogram
   {
   
   public:

      /**
      * Constructor.
      */
      HistogramWeight();

      /**
      * Default destructor.
      */
      virtual ~HistogramWeight();

      /**
      * Read parameters.
      */
      virtual void readParam(std::istream& in);

      /**
      * Write parameters.
      */
      virtual void writeParam(std::ostream& out) const;

      /**
      * Accumulate a new value to histogram and bias weight.
      */
      virtual void sample(const double value, const int intervalIn);

      /**
      * Read prescribed weight function.
      */
      void readWeight(std::istream& in, const int itr);

      /**
      * Get the weight by bin index.
      */
      double getWeight(const int id) const;

      /**
      * Get the weight by value.
      */
      double getWeight(const double value) const;

      /**
      * Get the count of iteration rounds.
      */
      int getNIteration() const;

      /**
      * Clear histogram and decrease the weight step size.
      */
      void resetHistogram();

      /**
      * Shift the bias weight such that the maximum sits at zero.
      */
      void shiftWeight();

      /**
      * Output the bias weight.
      */
      virtual void outputWeight(std::ostream& out) const;

      /**
      * Return true the maximum iteration number if reached.
      */
      bool doneIteration() const;

   private:

      /// Weighting variables.
      int                 nIteration_;
      double              deltaW_;
      std::vector<double> weight_;

      int                 nIterationTol_;
      double              sTolerance_;
      int                 nMaxIteration_;
   };


   /*
   * Get the count of iteration rounds.
   */
   inline int HistogramWeight::getNIteration() const
   { return nIteration_; }


   /*
   * Return true if the maximum iteration number is reached.
   */
   inline bool HistogramWeight::doneIteration() const
   { return nIteration_ > nMaxIteration_; }


   /*
   * Get bias weight.
   */
   inline double HistogramWeight::getWeight(const int id) const
   {
      if (id >= 0 && id < nBin_)
         return weight_[id];
      else
         return 1.0E14;
   }

   /*
   * Get bias weight.
   */
   inline double HistogramWeight::getWeight(const double value) const
   {
      if (value > min_ && value < max_) {
         int i = binIndex(value);
         return weight_[i];
      } else {
         return 1.0E14;
      }
   }

}
#endif

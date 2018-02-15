#ifndef POISSONGREENFUNC_H
#define POISSONGREENFUNC_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <vector>
#include <utility>
#include <fftw3.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_complex_math.h>

#include "../util/global.h"
#include "../util/IntVector.h"
#include "../util/Vector.h"

namespace GridMC
{

   using namespace Util;

   /**
   * An PoissonGreenFunc object solves the Poisson's equation on lattice.
   *
   * \ingroup Grid_Module
   */
   class PoissonGreenFunc
   {
   
   public:
  
      /**
      * Default constructor.
      */
      PoissonGreenFunc();

      /**
      * Default destructor.
      */
      virtual ~PoissonGreenFunc();

      /**
      * Set Grid length and number of grid.
      */
      void setBox(const Vector& boxIn, const IntVector& nGridIn);

      /**
      * Rescale box dimension.
      *
      * \param s1  scaling factor for box edge length.
      */
      virtual void rescaleBox(const double s1);

      /**
      * Get the pointer to dielectric array.
      */
      double* getDielectricArray();

      /**
      * Get the pointer to green's function array.
      */
      double* getGreenFuncArray();

      /**
      * The main method: calculate the Green's function.
      */
      void calculateGreenFunc();

   protected:

      /// Box dimension.
      Vector boxL_;

      /// Box volume.
      double boxV_;

      /// Number of grid points.
      IntVector nGrid_;

      /// Half number of grid points along each direction.
      IntVector nGridHalf_;

      /// Grid size along each direction.
      Vector drGrid_;

      /// Number of waves.
      int nReal_, nFourier_;

      /// Dielectric constant arrays and fftw plan.
      double        *dielectric_real;
      fftw_complex  *dielectric_fourier;
      fftw_plan     p;

      /// Poisson's operator, Green's function and fftw plan.
      gsl_matrix_complex  *poisson_matrix;
      fftw_complex        **green_fourier;
      double              **green_real;
      fftw_plan           pinv;

      /**
      * Wrap 3d indices to 1d index.
      */
      int wrap(const IntVector& ind) const;

      /**
      * Unwrap 1d indices to 3d index.
      */
      void unwrap(const int id, IntVector& ind) const;

      /**
      * Wrap 3d indices to 1d  FFTW-cmplex index.
      */
      int wrap_c(const IntVector& ind) const;

      /**
      * Unwrap FFTW-c 1d indices to 3d index.
      */
      void unwrap_c(const int id, IntVector& ind) const;

      /**
      * Shift index to the primary cell.
      */
      void to_primary(IntVector& ind) const;

      /**
      * Shift index to the first Brillouin zone.
      */
      void to_FBZ(IntVector& ind) const;

      /**
      * Check if the Poisson matrix is Hermitian: G(i,j) = G(j,i).
      */
      void check_poisson_hermitian() const;

      /**
      * Check if the Poisson matrix is "central symmetric": G(-i,-j)=G(i,j)*.
      */
      void check_poisson_duality() const;

      /**
      * Check if the Green's matrix is symmetric: G_real(i,j) = G_real(j,i).
      */
      void check_green_symmetry() const;

      /**
      * Output Green's matrix.
      */
      void output_green(std::ostream& out) const;

   };

   /*
   * Get the pointer to green's function array.
   */
   inline double* PoissonGreenFunc::getDielectricArray()
   { return dielectric_real; }

   /*
   * Get the pointer to green's function array.
   */
   inline double* PoissonGreenFunc::getGreenFuncArray()
   { return green_real[0]; }

   /*
   * Wrap 3d indices to 1d index.
   */
   inline int PoissonGreenFunc::wrap(const IntVector& ind) const
   { return ((ind[0]*nGrid_[1] + ind[1])*nGrid_[2] + ind[2]); }

   /*
   * Unwrap 1d indices to 3d index.
   */
   inline void PoissonGreenFunc::unwrap(const int id, IntVector& ind) const
   {
      ind[2] = id % nGrid_[2];
      ind[1] = ((id - ind[2]) / nGrid_[2]) % nGrid_[1];
      ind[0] = (id - ind[2] - ind[1] * nGrid_[2]) / nGrid_[1] / nGrid_[2];
   }

   /*
   * Wrap 3d indices to 1d FFTW-complex index.
   *
   * ASSUMPTION
   *    ind[2] < nGridHalf_[2]
   */
   inline int PoissonGreenFunc::wrap_c(const IntVector& ind) const
   { return ((ind[0]*nGrid_[1] + ind[1])*nGridHalf_[2] + ind[2]); }

   /*
   * Unwrap FFTW-c 1d indices to 3d index.
   */
   inline void PoissonGreenFunc::unwrap_c(const int id, IntVector& ind) const
   {
      ind[2] = id % nGridHalf_[2];
      ind[1] = ((id - ind[2]) / nGridHalf_[2]) % nGrid_[1];
      ind[0] = (id - ind[2] - ind[1] * nGridHalf_[2]) / nGrid_[1] / nGridHalf_[2];
   }

   /*
   * Shift index to the primary cell.
   */
   inline void PoissonGreenFunc::to_primary(IntVector& ind) const
   {
      for (int i = 0; i < Dimension; ++i)
         if (ind[i] < 0)
            ind[i] += nGrid_[i];
   }

   /*
   * Shift index to the first Brillouin zone.
   */
   inline void PoissonGreenFunc::to_FBZ(IntVector& ind) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if (ind[i] >= nGridHalf_[i])
            ind[i] -= nGrid_[i];
         else {
            if (nGrid_[i]%2 == 0) {
               if (ind[i] <= -nGridHalf_[i] + 1)
                  ind[i] += nGrid_[i];
            } else {
               if (ind[i] <= -nGridHalf_[i])
                  ind[i] += nGrid_[i];
            }
         }
      }
   }

}
#endif

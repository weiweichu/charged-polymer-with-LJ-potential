#ifndef POISSONGREENFUNC_CPP
#define POISSONGREENFUNC_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "PoissonGreenFunc.h"
#include "../util/Constants.h"

#include <math.h>
#include <iomanip>

namespace GridMC
{
   using namespace std;

   /*
   * Constructor.
   *
   * Assumptions:
   *    Unit cell shape is orthogonal.
   */
   PoissonGreenFunc::PoissonGreenFunc():
      boxL_(1.0),
      boxV_(1.0),
      nGrid_(2),
      nGridHalf_(1),
      drGrid_(0.5),
      nReal_(0),
      nFourier_(0)
   {}

   /*
   * Default destructor.
   */
   PoissonGreenFunc::~PoissonGreenFunc()
   {
      if (dielectric_real) fftw_free(dielectric_real);
      if (dielectric_fourier) fftw_free(dielectric_fourier);
      if (p) fftw_destroy_plan(p);

      if (poisson_matrix) gsl_matrix_complex_free(poisson_matrix);
      if (green_fourier[0]) fftw_free(green_fourier[0]);
      if (green_fourier) fftw_free(green_fourier);
      if (green_real[0]) fftw_free(green_real[0]);
      if (green_real) fftw_free(green_real);
      if (pinv) fftw_destroy_plan(pinv);
   }

   /*
   * Set unit cell length and grid dimension and allocate arrays.
   */
   void PoissonGreenFunc::setBox(const Vector& boxLIn, const IntVector& nGridIn)
   {
      int i;

      boxL_ = boxLIn;
      boxV_ = 1.0;
      for (i = 0; i < Dimension; ++i) boxV_ *= boxL_[i];

      nGrid_ = nGridIn;
      for (i = 0; i < Dimension; ++i) {
         if (nGrid_[i]%2 == 0)
            nGridHalf_[i] = nGrid_[i] / 2 + 1;
         else
            nGridHalf_[i] = (nGrid_[i] + 1) / 2;
      }

      for (i = 0; i < Dimension; ++i)
         drGrid_[i] = boxL_[i] / double(nGrid_[i]);

      // Calculate array size.
      nReal_ = 1;
      for (i = 0; i < Dimension; ++i)
         nReal_ *= nGrid_[i];

      nFourier_ = 1;
      for (i = 0; i < Dimension - 1; ++i)
         nFourier_ *= nGrid_[i];
      nFourier_ *= nGridHalf_[Dimension - 1];

      // Allocate dielectric constant array; create plan.
      dielectric_real = (double*) fftw_malloc(sizeof(double) * nReal_);
      dielectric_fourier = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nFourier_);

      int dim[Dimension];
      for (i = 0; i < Dimension; ++i)
         dim[i] = nGrid_[i];
      p = fftw_plan_dft_r2c(Dimension, dim, dielectric_real, dielectric_fourier, FFTW_MEASURE);

      // Allocate Poisson's operator, and green's functions; create plan.
      poisson_matrix = gsl_matrix_complex_calloc(nReal_ - 1, nReal_ - 1);

      double* foor = (double*) fftw_malloc(sizeof(double) * nReal_ * nReal_);
      green_real = (double**) fftw_malloc(sizeof(double*) * nReal_);
      for (i = 0; i < nReal_; ++i)
         green_real[i] = foor + i * nReal_;

      fftw_complex* fooc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nReal_ * nFourier_);
      green_fourier = (fftw_complex**) fftw_malloc(sizeof(fftw_complex*) * nReal_);
      for (i = 0; i < nReal_; ++i)
         green_fourier[i] = fooc + i * nFourier_;

      int dim2[2*Dimension];
      for (i = 0; i < Dimension; ++i) {
         dim2[i] = nGrid_[i];
         dim2[i+Dimension] = nGrid_[i];
      }
      pinv = fftw_plan_dft_c2r(2*Dimension, dim2, green_fourier[0], green_real[0], FFTW_MEASURE);
   }

   /*
   * Reset unit cell length and grid dimension.
   */
   void PoissonGreenFunc::rescaleBox(const double s1)
   {
      boxL_   *= s1;
      boxV_   *= s1 * s1 * s1;
      drGrid_ *= s1;
   }

   /*
   * Calculate the Green's function by matrix inversion.
   *
   * ALGORITHM
   *    (1) Fourier-transform the dielectric constant.
   *    (2) Construct the Poisson's operator.
   *    (3) Remove the zero mode and invert the operator --> green's function, G(q, q').
   *    (4) Inverse-Fourier-transform the green function.
   *    (5) Change notation from G(r, -r') to G(r, r').
   */
   void PoissonGreenFunc::calculateGreenFunc()
   {
      int           i, j, j2, k;               /// Wave index.
      int           qdim;                          /// Rank index.
      IntVector     vi, vj, vk;                    /// Wave index on grid.
      IntVector     viFBZ, vjFBZ;                  /// Wave index on grid.
      Vector        wave_unit;                     /// Unit length squared of waves.
      fftw_complex  *elek;                         /// Pointer to k vector.
      int           sign;                          /// If less than Nyquist frequency.
      double        dot_ij;                        /// Dot product of waves.
      double        re, im;                        /// Real and Imaginary components.
      double        scale(1.0);                    /// FFT normalization.
      double        swap;                          /// Temporary swap variable.

      // Wave vector units.
      for (qdim = 0; qdim < Dimension; ++qdim)
         wave_unit[qdim] = pow(2.0 * Constants::Pi / boxL_[qdim], 2);

      // Calculate the Fourier coefficients for the dielectric constants.
      fftw_execute(p);

      // Build the Poisson's operator; entries with zero frequency are eliminated.
      for (i = 0; i < nReal_ - 1; ++i) {
         unwrap(i+1, vi);
         viFBZ = vi;
         to_FBZ(viFBZ);

         for (j = 0; j < nReal_ - 1; ++j) {
            unwrap(j+1, vj);
            vjFBZ = vj;
            to_FBZ(vjFBZ);

            // Dot product of i and j wave vectors. "Same sign" convention is used for boundary waves.
            dot_ij = 0.0;
            for (qdim = 0; qdim < Dimension; ++qdim) {
               if (nGrid_[qdim] % 2 == 0 && (viFBZ[qdim] == nGrid_[qdim] / 2 || vjFBZ[qdim] == nGrid_[qdim] / 2) ) {
                  dot_ij += wave_unit[qdim] * abs(viFBZ[qdim] * vjFBZ[qdim]);
               } else {
                  dot_ij += wave_unit[qdim] * viFBZ[qdim] * vjFBZ[qdim];
               }
            }

            // Wave vector for the dieletric constant.
            vk.subtract(vi, vj);
            to_primary(vk);

            sign = 1;
            if (vk[Dimension - 1] >= nGridHalf_[Dimension-1]) {
               sign = -1;
               vk *= -1;
               to_primary(vk);
            }

            k = wrap_c(vk);
            elek = dielectric_fourier + k;

            re = dot_ij * (*elek)[0];
            im = dot_ij * (*elek)[1] * sign;

            gsl_matrix_complex_set(poisson_matrix, i, j, gsl_complex_rect(re, im));
         }
      }

      #ifdef GREEN_SYMMETRY
      cout << "----------------------------------------------" << endl;

      // Check that the symmetry of Poisson's matrix is preserved.
      check_poisson_hermitian();
      check_poisson_duality();

      cout << "In-place inversion" << endl;
      #endif
 
      // Invert the poisson matrix.
      gsl_linalg_complex_cholesky_decomp(poisson_matrix);
      gsl_linalg_complex_cholesky_invert(poisson_matrix);

      // Check that the symmetry of Green's matrix is preserved.
      #ifdef GREEN_SYMMETRY
      check_poisson_hermitian();
      check_poisson_duality();
      #endif
 
      // Store the inverted Poisson matrix onto the green's function array.
      for (i = 0; i < nReal_; ++i) {
         for (j = 0; j < nFourier_; ++j) {
            if (i == 0 || j == 0) {
               green_fourier[i][j][0] = 0.0;
               green_fourier[i][j][1] = 0.0;
            } else {
               unwrap_c(j, vj);
               j2 = wrap(vj);
               gsl_complex tmp = gsl_matrix_complex_get(poisson_matrix, i-1, j2-1);
               green_fourier[i][j][0] = GSL_REAL(tmp);
               green_fourier[i][j][1] = GSL_IMAG(tmp);
            }
         }
      }

      // Inverse Fourier transform of the Green's operator.
      fftw_execute(pinv);

      // Flip the sign of arguments for the real space Green's function: G(r, -r') --> G(r, r').
      scale = 4.0 * Constants::Pi; // The factor left for the Bjerrum length.
      for (i = 0; i < Dimension; ++i)
         scale *= double(nGrid_[i]) / boxL_[i];

      for (i = 0; i < nReal_; ++i) {
         for (j = 0; j < nReal_; ++j) {
            unwrap(j, vj);

            vj *= -1;
            to_primary(vj);
            j2 = wrap(vj);

            if (j < j2) {
               swap = green_real[i][j];
               green_real[i][j] = green_real[i][j2] * scale;
               green_real[i][j2] = swap * scale;
            } else if (j == j2) {
               green_real[i][j] *= scale;
            }
         }
      }

      #ifdef GREEN_SYMMETRY
      check_green_symmetry();
      cout << "----------------------------------------------" << endl;
      #endif

   }

   /*
   * Check if the poisson's matrix is Hermitian.
   */
   void PoissonGreenFunc::check_poisson_hermitian() const
   {
      cout << "Hermitian check: G(i,j)   = G(j,i)*" << endl;

      int i, j;
      for (i = 0; i < nReal_-1; ++i) {
         for (j = 0; j <= i; ++j) {
            gsl_complex tmp1 = gsl_matrix_complex_get(poisson_matrix, i, j);
            gsl_complex tmp2 = gsl_matrix_complex_get(poisson_matrix, j, i);

            if (fabs(GSL_IMAG(tmp1) + GSL_IMAG(tmp2)) > Constants::Epsilon ||
                fabs(GSL_REAL(tmp1) - GSL_REAL(tmp2)) > Constants::Epsilon) {
               cout << "(" << i << ", " << j << "): " << GSL_REAL(tmp1) << "  " << GSL_IMAG(tmp1) << endl;
               cout << "(" << j << ", " << i << "): " << GSL_REAL(tmp2) << "  " << GSL_IMAG(tmp2) << endl;
            }
         }
      }
   }


   /*
   * Check if the poisson's matrix is "central symmetric": G(-i,-j)=G(i,j)*.
   */
   void PoissonGreenFunc::check_poisson_duality() const
   {
      cout << "Duality   check: G(-i,-j) = G(i,j)*" << endl;

      // Calculate auxiliary array sizes.
      int        i, j, i2, j2;
      IntVector  vi, vj, vi2, vj2;

      for (i = 0; i < nReal_-1; ++i) {
         unwrap(i+1, vi);

         vi2 = vi;
         vi2 *= -1;
         to_primary(vi2);
         i2 = wrap(vi2) - 1;

         for (j = 0; j < nReal_-1; ++j) {
            gsl_complex tmp1 = gsl_matrix_complex_get(poisson_matrix, i, j);

            unwrap(j+1, vj);

            vj2 = vj;
            vj2 *= -1;
            to_primary(vj2);
            j2 = wrap(vj2) - 1;

            gsl_complex tmp2 = gsl_matrix_complex_get(poisson_matrix, i2, j2);

            if (fabs(GSL_IMAG(tmp1) + GSL_IMAG(tmp2)) > Constants::Epsilon ||
                fabs(GSL_REAL(tmp1) - GSL_REAL(tmp2)) > Constants::Epsilon) {
               cout << "(" << vi[0] << " " << vi[1] << " " << vi[2] << ",   "
                           << vj[0] << " " << vj[1] << " " << vj[2] << "):  ";
               cout << GSL_REAL(tmp1) << "  " << GSL_IMAG(tmp1) << endl;

               cout << "(" << vi2[0] << " " << vi2[1] << " " << vi2[2] << ",  "
                           << vj2[0] << " " << vj2[1] << " " << vj2[2] << "):  ";
               cout << GSL_REAL(tmp2) << "  " << GSL_IMAG(tmp2) << endl << endl;
            }
         }

      }
   }


   /*
   * Check if the green's matrix is symmetric.
   */
   void PoissonGreenFunc::check_green_symmetry() const
   {
      cout << "Symmetry  check: G_real(i,j) = G_real(j,i)" << endl;

      for (int i = 0; i < nReal_; ++i) {
         for (int j = 0; j < i; ++j) {
            if (fabs(green_real[i][j] - green_real[j][i]) > Constants::Epsilon) {
               cout << "(" << i << ", " << j << "): " << setw(12) << green_real[i][j] << "  ";
               cout << "(" << j << ", " << i << "): " << setw(12) << green_real[j][i] << endl;
            }
         }
      }
   }


   /*
   * Output the Green's matrix.
   */
   void PoissonGreenFunc::output_green(std::ostream& out) const
   {
      out.setf(ios::fixed, ios::floatfield);
      for (int i = 0; i < nReal_; ++i) {
         for (int j = 0; j < nReal_; ++j) {
            out << setw(12) << green_real[i][j];
         }
         out << endl;
      }
      out.unsetf(ios::floatfield);
   }

}

#endif

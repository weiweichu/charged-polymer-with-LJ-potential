#ifndef GRIDMASSCHARGE_CPP
#define GRIDMASSCHARGE_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "GridMassCharge.h"
#include "../util/Constants.h"
#include "../util/global.h"

#include <string>
#include <iomanip>

namespace GridMC
{
   using namespace Util;

   const double GridMassCharge::UCself_ = 1.882312645;

   /*
   * Constructor. Assume cell to be orthogonal.
   */
   GridMassCharge::GridMassCharge() :
      GridMass(),
      esStrength_(0.0),
      esRatio_(1.0)
   {}

   /*
   * Default destructor.
   */
   GridMassCharge::~GridMassCharge()
   {
      if (green_real_) fftw_free(green_real_);
   }

   /*
   * Rescale box dimension, interaction strength, and Green's function.
   */
   void GridMassCharge::rescaleBox(const double s1)
   {
      GridMass::rescaleBox(s1);
      double invs = 1.0 / s1;
      for (unsigned i = 0; i < gridG_.size(); ++i)
         gridG_[i] *= invs;
   }

   /*
   * Set Coulomb interaction strength.
   */
   void GridMassCharge::setCoulombStrength(const double strengthIn, const double ratioIn)
   {
      esStrength_ = strengthIn;
      esRatio_ = ratioIn;
   }

   /*
   * Allocate memory.
   */
   void GridMassCharge::allocate(const int nTypeIn)
   {
      // Allocate mass array.
      allocateMassArray(nTypeIn);

      // Allocate charge array.
      vector<double> zero1(nGrid_[2], 0.0);
      vector< vector<double> > zero2(nGrid_[1], zero1);
      Q_.resize(nGrid_[0], zero2);

      // Allocate and calcualte green's array.
      buildGridGreen();

      // Map pointers for dielectric constant and green's function.
      dielectric_real_ = poissonGreenFunc_.getDielectricArray();

      int nReal(1);
      for (int i = 0; i < Dimension; ++i)
         nReal *= nGrid_[i];

      double* foo = poissonGreenFunc_.getGreenFuncArray();
      green_real_ = (double**) fftw_malloc(sizeof(double*) * nReal);
      for (int i = 0; i < nReal; ++i)
         green_real_[i] = foo + i * nReal;
   }

   /*
   * Add bead's mass and charge contributions.
   */
   void GridMassCharge::insertBead(const Particle& bead)
   {
      GridMass::insertBead(bead);
      insertCharge(bead.q);
   }

   /*
   * Accumulate bead charge to density field.
   */
   void GridMassCharge::insertCharge(const double q)
   {
      IntVector ind;
      for (unsigned i = 0; i < insertId_.size(); ++i) {
         unwrap(insertId_[i].first, ind);
         Q_[ind[0]][ind[1]][ind[2]] += q * insertId_[i].second;
      }
   }

   /*
   * Get coulomb energy of a particle move.
   */
   double GridMassCharge::getCoulombEnergyChange(const double q)
   {

      double deltaE(0.0);

      if (fabs(q) > Constants::Epsilon) {
         int       i, j;
         double    dQ, dQ2;
         IntVector ia, ib;
         #ifndef GREEN_NEW
         IntVector dv;
         #endif

         // Linear terms.
         for (i = 0; i < nChanged_; ++i) {
            unwrap(changeId_[i].first, ia); 
            dQ = q * changeId_[i].second;

            for (ib[0] = 0; ib[0] < nGrid_[0]; ++ib[0]) {
            for (ib[1] = 0; ib[1] < nGrid_[1]; ++ib[1]) {
            for (ib[2] = 0; ib[2] < nGrid_[2]; ++ib[2]) {
               #ifdef GREEN_NEW
               deltaE += dQ * getGreen(ia, ib) * Q_[ib[0]][ib[1]][ib[2]];
               #else
               dv.subtract(ia, ib);
               shiftDIndex(dv);
               deltaE += dQ * getGreen(dv) * Q_[ib[0]][ib[1]][ib[2]];
               #endif
            } } }
         }

         // Quadratic terms.
         for (i = 0; i < nChanged_; ++i) {
            unwrap(changeId_[i].first, ia); 
            dQ = q * changeId_[i].second;

            // Self term.
            #ifdef GREEN_NEW
            deltaE += dQ * getGreen(ia, ia) * dQ * 0.5;
            #else
            deltaE += dQ * getGreen(IntVector(0)) * dQ * 0.5;
            #endif

            // Cross terms.
            for (j = i + 1; j < nChanged_; ++j) {
               unwrap(changeId_[j].first, ib); 
               dQ2 = q * changeId_[j].second;

               #ifdef GREEN_NEW
               deltaE += dQ * getGreen(ia, ib) * dQ2;
               #else
               dv.subtract(ia, ib);
               shiftDIndex(dv);
               deltaE += dQ * getGreen(dv) * dQ2;
               #endif
            }
         }

         deltaE *= esStrength_;
      }
 
      return deltaE;
   }


   /*
   * Get the totalt coulomb energy.
   */
   double GridMassCharge::getCoulombEnergy() const
   {
      double ese(0.0);
      IntVector ia, ib;
      #ifndef GREEN_NEW
      IntVector dv;
      #endif

      for (ia[0] = 0; ia[0] < nGrid_[0]; ++ia[0]) {
      for (ia[1] = 0; ia[1] < nGrid_[1]; ++ia[1]) {
      for (ia[2] = 0; ia[2] < nGrid_[2]; ++ia[2]) {

         for (ib[0] = 0; ib[0] < nGrid_[0]; ++ib[0]) {
         for (ib[1] = 0; ib[1] < nGrid_[1]; ++ib[1]) {
         for (ib[2] = 0; ib[2] < nGrid_[2]; ++ib[2]) {

            #ifdef GREEN_NEW
            ese += Q_[ia[0]][ia[1]][ia[2]] * getGreen(ia, ib) * Q_[ib[0]][ib[1]][ib[2]];
            #else
            dv.subtract(ia, ib);
            shiftDIndex(dv);
            ese += Q_[ia[0]][ia[1]][ia[2]] * getGreen(dv) * Q_[ib[0]][ib[1]][ib[2]];
            #endif

         } } }

      } } }

      return (ese * 0.5 * esStrength_);
   }

   /*
   * Return the self interaction energy.
   */
   double GridMassCharge::getSelfEnergy() const
   {
      double self(0.0), q, selfCoeff;
      IntVector ia;

      selfCoeff = getGreen(IntVector(0));

      for (ia[0] = 0; ia[0] < nGrid_[0]; ++ia[0]) {
      for (ia[1] = 0; ia[1] < nGrid_[1]; ++ia[1]) {
      for (ia[2] = 0; ia[2] < nGrid_[2]; ++ia[2]) {
         q = Q_[ia[0]][ia[1]][ia[2]];
         #ifdef GREEN_NEW
         self += q * getGreen(ia, ia) * q;
         #else
         self += q * selfCoeff * q;
         #endif
      } } } 

      return (0.5 * self * esStrength_);
   }


   /*
   * Periodic Green's function. Eq(5) in paper1.
   *
   * Assumed Dimension = 3.
   */
   void GridMassCharge::buildGridGreen()
   {
      Vector r;
      IntVector ind;
      int counter(0);
      unsigned nfG(1);

      // allocate Green's function array
      for (int i = 0; i < Dimension; ++i) nfG *= nGridHalf_[i];
      gridG_.assign(nfG, 0.0);

      for (ind[0] = 0; ind[0] < nGridHalf_[0]; ++ind[0]) {
      for (ind[1] = 0; ind[1] < nGridHalf_[1]; ++ind[1]) {
      for (ind[2] = 0; ind[2] < nGridHalf_[2]; ++ind[2]) {

         // Caution: think about the case to be excluded.
         if (ind[0] + ind[1] + ind[2] > 0) {
            r[0] = double(ind[0]) / double(nGrid_[0]) * boxL_[0];
            r[1] = double(ind[1]) / double(nGrid_[1]) * boxL_[1];
            r[2] = double(ind[2]) / double(nGrid_[2]) * boxL_[2];
            gridG_[counter++] = calGreen(r);
         } else {
            //gridG_[counter++] = 0.0;
            gridG_[counter++] = UCself_ / drGrid_[0];
         }

      } } } // loop over position indices
   }


   /*
   * Get the Green's function by real-time calculation.
   *
   * Assumpsion: 0 <= r[i] <= boxL_[i] for i = 0, 1, 2
   */
   double GridMassCharge::calGreen(const Vector& r) const
   {
      IntVector qind;
      double tmp1, tmp2, Pi2(2.0*Constants::Pi), g3D(0.0); 

      // Think carefully about the case to be excluded.
      if (fabs(r[0]) + fabs(r[1]) + fabs(r[2]) > 0) {

         for (qind[0] = 0; qind[0] < nGridHalf_[0]; ++qind[0]) {
         for (qind[1] = 0; qind[1] < nGridHalf_[1]; ++qind[1]) {
         for (qind[2] = 0; qind[2] < nGridHalf_[2]; ++qind[2]) {

            tmp1 = gamma(qind[0], qind[1], qind[2]);
            if (tmp1 > 0.0) {
               tmp1 *= cos(Pi2 * double(qind[0]) * r[0] / boxL_[0]);
               tmp1 *= cos(Pi2 * double(qind[1]) * r[1] / boxL_[1]);
               tmp1 *= cos(Pi2 * double(qind[2]) * r[2] / boxL_[2]);

               tmp2  = double(qind[0]*qind[0]) / (boxL_[0]*boxL_[0]);
               tmp2 += double(qind[1]*qind[1]) / (boxL_[1]*boxL_[1]);
               tmp2 += double(qind[2]*qind[2]) / (boxL_[2]*boxL_[2]);

               g3D += tmp1 / tmp2;
            }

         } } } // loop over wave vector indices
      }
      return (g3D * 8.0 / Constants::Pi / boxV_);
   }

   /*
   * Calculate the array for dielectric constant and solve the Poisson's equation.
   */
   void GridMassCharge::solvePoisson()
   {
      IntVector   ind;
      int         i;

      /// Set the dielectric coefficients.
      for (ind[0] = 0; ind[0] < nGrid_[0]; ++ind[0]) {
         for (ind[1] = 0; ind[1] < nGrid_[1]; ++ind[1]) {
            for (ind[2] = 0; ind[2] < nGrid_[2]; ++ind[2]) {
               i = wrap(ind);
               dielectric_real_[i] = (vtk_[i][0] + vtk_[i][1]*esRatio_) / (vtk_[i][0] + vtk_[i][1]);
            }
         }
      }
      cout<<"dielectic permittivity: "<<dielectric_real_[1]<<endl;
      /// Calculate the Green's function.
      poissonGreenFunc_.calculateGreenFunc();

      #ifdef GREEN_CONSISTENCE
      /// Check the consistency of the two definition.
      IntVector ia, ib, dv;
      cout.setf(ios::fixed, ios::floatfield);
      for (ia[0] = 0; ia[0] < nGrid_[0]; ++ia[0]) {
      for (ia[1] = 0; ia[1] < nGrid_[1]; ++ia[1]) {
      for (ia[2] = 0; ia[2] < nGrid_[2]; ++ia[2]) {

         for (ib[0] = 0; ib[0] < nGrid_[0]; ++ib[0]) {
         for (ib[1] = 0; ib[1] < nGrid_[1]; ++ib[1]) {
         for (ib[2] = 0; ib[2] < nGrid_[2]; ++ib[2]) {

            dv.subtract(ia, ib);
            shiftDIndex(dv);
            if (fabs(getGreen(dv) - getGreen(ia, ib) > Constants::Epsilon)) {
               cout << setw(3) << ia[0];
               cout << setw(3) << ia[1];
               cout << setw(3) << ia[2];
               cout << setw(7) << ib[0];
               cout << setw(3) << ib[1];
               cout << setw(3) << ib[2];
               cout << setw(12) << getGreen(dv) / getGreen(ia, ib) << endl;
            }

         } } }
      } } }
      cout.unsetf(ios::floatfield);
      #endif

   }


   /*
   * Get the dielectric permittivity at given position.
   */
   double GridMassCharge::getDielectricPermittivity(const Vector& r)
   {
      double epsilon(0.0);
      fillSiteList(r, +1);
      for (int ic = 0; ic < 8; ++ic) {
         epsilon += dielectric_real_[ insertId_[ic].first ] * insertId_[ic].second;
      }
      return epsilon;
   }

}

#endif

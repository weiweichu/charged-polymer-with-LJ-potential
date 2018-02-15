#ifndef GRIDMASS_CPP
#define GRIDMASS_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "GridMass.h"
#include "../util/global.h"

#include <string>
#include <iomanip>

namespace GridMC
{
   using namespace Util;

   /*
   * Constructor. Assume cell to be orthogonal.
   */
   GridMass::GridMass() :
      Grid(),
      nType_(0),
      uPair_(),
      phi_(),
      vtk_(),
      nFramesVTK_(0)
   {}

   /*
   * Default destructor.
   */
   GridMass::~GridMass()
   {}

   /*
   * Reset box dimension and update interaction strength.
   */
   void GridMass::rescaleBox(const double s1)
   {
      Grid::rescaleBox(s1);

      double invs = 1.0 / (s1 * s1 * s1);
      for (int i = 0; i < nType_; ++i)
         for (int j = 0; j < nType_; ++j)
            uPair_[i][j] *= invs;
   }

   /*
   * Set the interaction strength; for convenient evaluation of 2-body overlap enrgy.
   */
   void GridMass::allocateMassArray(const int nTypeIn)
   {
      nType_ = nTypeIn;

      MassV zero(nType_, 0.0);
      uPair_.resize(nType_, zero);
      phi_.resize(nGrid_[0]*nGrid_[1]*nGrid_[2], zero);
   }

   /*
   * Set the interaction strength; for convenient evaluation of 2-body overlap enrgy.
   */
   void GridMass::setTwobodyInteraction(const MassM& chiN,
          const double kappaN, const double sqrtNbar, const int polymerN)
   {
      int    i, j;
      double normalization(1.0);

      normalization /= (sqrtNbar * double(polymerN * polymerN));
      normalization /= (drGrid_[0] * drGrid_[1] * drGrid_[2]);

      for (i = 0; i < nType_; ++i)
         uPair_[i][i] = 0.5 * kappaN * normalization;

      for (i = 0; i < nType_; ++i) {
         for (j = i + 1; j < nType_; ++j)
            uPair_[i][j] = (chiN[i][j] + kappaN) * normalization;
         for (j = 0; j < i; ++j)
            uPair_[i][j] = uPair_[j][i];
      }
   }


   void GridMass::updateNeighbor(const int id, const Vector& r, const Vector& rTrial) {
     fillNeighborList(id, r, -1);
     fillNeighborList(id, rTrial, 1);
//     int total = 0;
//     for (int i = 0; i < int(neighbor_.size()); i++) {
//       total += int(neighbor_[i].size());
//     }
//     Log::file() << "number of charged beads after neighbor update   " <<total << std::endl;
   }
   /*
   * Insert bead to the grid.
   */
   void GridMass::insertBead(const Particle& bead)
   {
      fillSiteList(bead.r, +1);
      insertMass(bead.t);
   }

   /*
   * Get the two body energy change of a proposed bead move using the affected site list.
   */
   double GridMass::getTwobodyEnergyChange(const int type)
   {
      MassV  mass(nType_, 0);
      double deltaPhi, deltaE(0.0);
      int    i, j;

      for (i = 0; i < nChanged_; ++i) {
         mass = phi_[changeId_[i].first];
         deltaPhi = changeId_[i].second;
         for (j = 0; j < nType_; ++j) {
            if (j == type)
               deltaE += deltaPhi * uPair_[j][j] * (2.0 * mass[j] + deltaPhi);
            else
               deltaE += deltaPhi * uPair_[type][j] * mass[j];
         }
      }
      return deltaE;
   }

   /*
   * Initialize VTK array.
   */
   void GridMass::initializeVTK()
   {
      MassV zero(nType_, 0.0);
      vtk_.resize(nGrid_[0]*nGrid_[1]*nGrid_[2], zero);
      nFramesVTK_ = 0;
   }

   /*
   * Update VTK array.
   */
   void GridMass::updateVTK()
   {
      unsigned i, j;
      for (i = 0; i < vtk_.size(); ++i)
         for (j = 0; j < unsigned(nType_); ++j)
            vtk_[i][j] += phi_[i][j];
      nFramesVTK_ += 1;
   }

   /*
   * Reset VTK array.
   */
   void GridMass::resetVTK()
   {
      unsigned i, j;
      for (i = 0; i < vtk_.size(); ++i)
         for (j = 0; j < unsigned(nType_); ++j)
            vtk_[i][j] = 0.0;
      nFramesVTK_ = 0;
   }

   /*
   * Output normalized VTK array.
   */
   void GridMass::writeVTK(ostream& out) const
   {
      // Output VTK header.
      out << "# vtk DataFile Version 2.0" << endl;
      out << "GridMC configuration" << endl;
      out << "ASCII" << endl;
      out << "DATASET STRUCTURED_POINTS" << endl;
      out << "DIMENSIONS " << nGrid_[0] << " " << nGrid_[1] << " " << nGrid_[2] << endl;
      out << "ORIGIN 0 0 0" << endl;
      out << "SPACING " << drGrid_[0] << " " << drGrid_[1] << " " << drGrid_[2] << endl << endl;

      // Output VTK data header.
      int nj;
      out << "POINT_DATA " << nGrid_[0]*nGrid_[1]*nGrid_[2] << endl;
      out << "SCALARS volume_fraction float ";
      if (nType_ <= 5) {
         out << nType_ << endl;
         nj = nType_;
      } else {
         out << 4 << endl;
         Log::file() << "Number components greater than 4. Use 4 as required by VTK." << endl;
         nj = nType_;
      }
      out << "LOOKUP_TABLE default" << endl;

      // Save and reset precision.
      streamsize oldprec = out.precision(3);
      out.setf(ios::fixed, ios::floatfield);

      IntVector ivec;
      int i, j;
      for (ivec[0] = 0; ivec[0] < nGrid_[0]; ++ivec[0]) {
      for (ivec[1] = 0; ivec[1] < nGrid_[1]; ++ivec[1]) {
      for (ivec[2] = 0; ivec[2] < nGrid_[2]; ++ivec[2]) {
         i = wrap(ivec);
         out << setw(8) << vtk_[i][0] / double(nFramesVTK_);
         for (j = 1; j < nj; ++j) {
            out << setw(9) << vtk_[i][j] / double(nFramesVTK_);
         }
         out << endl;
      } } }

      // Restore I/O formats.
      out.unsetf(ios::floatfield);
      out.precision(oldprec);

   }

}

#endif

#ifndef ORDERPARAMETER_CPP
#define ORDERPARAMETER_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <math.h>
#include <vector>
#include "OrderParameter.h"

namespace GridMC
{
   /*
   * Constructor.
   */
   OrderParameter::OrderParameter(const System& system) :
      boxL_(system.getBoxL()),
      nGrid_(system.getNGrid()),
      grid_(system.getGrid()),
      beads_(system.getBeads()),
      nBead_(0),
      psi_(0.0),
      dpsi_(0.0)
   {}

   /*
   * Default destructor.
   */
   OrderParameter::~OrderParameter()
   {}

   /*
   * Order parameter based on grid density.
   *
   * Note that nBead_ is initialized in the readParam method of the derived class.
   */
   void OrderParameter::oprmFromGrid()
   {
      psi_ = 0.0;
      for (int i = 0; i < nGrid_[0]*nGrid_[1]*nGrid_[2]; ++i)
         psi_ += pow((grid_[i][0] - grid_[i][1]), 2);
      psi_ /= double(nBead_);
   }

   /*
   * Order parameter change based on grid mass density.
   *
   * (1) Must be called after non-bonded energy change be calculated.
   * (2) nBead_ is initialized in the readParam method of the derived class.
   */
   void OrderParameter::oprmChangeFromGrid(const int type)
   {
      dpsi_ = 0.0;

      if (type == 0 || type == 1) {
         double  phi1, phi2;
         int     nChanged = grid_.getNChanged();
         const   vector< pair<int,double> >& siteList = grid_.getChangedSiteList();

         for (int i = 0; i < nChanged; ++i) {
            phi1 = grid_[siteList[i].first][0];
            phi2 = grid_[siteList[i].first][1];
            dpsi_ -= (phi2 - phi1)*(phi2 - phi1);

            if (type == 0)
               phi1 += siteList[i].second;
            else
               phi2 += siteList[i].second;
            dpsi_ += (phi2 - phi1)*(phi2 - phi1);
         }
         dpsi_ /= double(nBead_);
      }
   }

}

#endif

#ifndef MOVE_CPP
#define MOVE_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iomanip>
#include "Move.h"

namespace GridMC
{
   /*
   * Constructor.
   */
   Move::Move(McSystem& sysIn) :
      nInterval_(0),
      nAttempt_(0),
      nAccept_(0),
      system_(sysIn),
      name_("")
   {}

   /*
   * Default destructor.
   */
   Move::~Move()
   {}

   /*
   * Output statistics.
   */
   void Move::outputStatistics(std::ostream& out) const
   {
      using namespace std;
      out.setf(ios::left, ios::adjustfield);
      out << setw(15) << name_;
      out.unsetf(ios::adjustfield);

      out << setw(15) << nAttempt_;
      out << setw(15) << nAccept_;
      out << setw(15) << (nAttempt_ == 0 ? 0.0 : double(nAccept_)/double(nAttempt_));
      out << endl;
   }

}

#endif

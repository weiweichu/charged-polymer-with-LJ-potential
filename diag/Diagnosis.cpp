#ifndef DIAGNOSIS_CPP
#define DIAGNOSIS_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "Diagnosis.h"

namespace GridMC
{

   /*
   * Constructor.
   */
   Diagnosis::Diagnosis(System& system) :
      system_(system),
      boxL_(system.boxL_),
      nGrid_(system.nGrid_),
      grid_(system.grid_),
      beads_(system.beads_),
      fileMaster_(*(system.fileMasterPtr_)),
      nInterval_(0),
      nSamples_(0),
      name_()
   {}

   /*
   * Default destructor.
   */
   Diagnosis::~Diagnosis()
   {}

   /*
   * Read sampling interval.
   */
   void Diagnosis::readParam(std::istream& in)
   {
      char        comment[200];
      std::string line;

      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: diagnosing interval");
      sscanf(line.c_str(), "%d %s", &nInterval_, comment);

      if (nInterval_ <= 0)
         UTIL_THROW("Invalid diagnosing inteval");
   }

   /*
   * Write parameter.
   */
   void Diagnosis::writeParam(std::ostream& out) const
   {
      out << "diagnosis name              " << name_ << std::endl;
      out << "nInterval                   " << nInterval_ << std::endl;
   }

}

#endif

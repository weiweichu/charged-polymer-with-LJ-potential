#ifndef FILEMASTER_CPP
#define FILEMASTER_CPP

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "FileMaster.h"

namespace Util
{

   /*
   * default constructor.
   */
   FileMaster::FileMaster() : filePrefix_("")
   {}

   /*
   * Default destructor.
   */
   FileMaster::~FileMaster()
   {}

   /*
   * Set prefix for file IO.
   */
   void FileMaster::setFilePrefix(const std::string& prefixIn)
   { filePrefix_ = prefixIn; }

   /*
   * Set prefix for file IO.
   */
   void FileMaster::prependPath(std::string& name) const
   {
      // Construct filename = inputPrefix_ + name
      std::string filename;
      filename += filePrefix_;
      filename += name;
      name = filename;
   }

   /*
   * Open and return an input file named filePrefix_ + name.
   */
   void FileMaster::openInputFile(const std::string& name, std::ifstream& in) const
   {
      // Construct filename = inputPrefix_ + name
      std::string filename;
      filename += filePrefix_;
      filename += name;

      in.open(filename.c_str());

      // Check for error opening file
      if (in.fail()) {
         std::string message = "Error opening input file. Filename: ";
         message += filename;
         UTIL_THROW(message.c_str());
      }


   }

   /*
   * Open and return an output file named filePrefix_ + name.
   */
   void FileMaster::openOutputFile(const std::string& name, std::ofstream& out) const
   {
      // Construct filename = outputPrefix_ + name
      std::string filename;
      filename += filePrefix_;
      filename += name;

      out.open(filename.c_str());

      // Check for error opening file
      if (out.fail()) {
         std::string message = "Error opening output file. Filename: ";
         message += filename;
         UTIL_THROW(message.c_str());
      }
   }

}
#endif

#ifndef FILEMASTER_H
#define FILEMASTER_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <fstream>
#include "util/global.h"

namespace Util
{

   /**
   * A class for storing I/O paths and generating files.
   *
   * \ingroup Util_Module
   */
   class FileMaster 
   {
   
   public:

      /**
      * default constructor.
      */
      FileMaster();

      /**
      * Default destructor.
      */
      ~FileMaster();

      /**
      * Set prefix for file IO.
      */
      void setFilePrefix(const std::string& prefixIn);

      /**
      * Get full file name.
      */
      void prependPath(std::string& name) const;

      /**
      * Open and return an input file named filePrefix_ + name.
      */
      void openInputFile(const std::string& name, std::ifstream& in) const;

      /**
      * Open and return an output file named filePrefix_ + name.
      */
      void openOutputFile(const std::string& name, std::ofstream& out) const;

   private:

      /// Common prefix to 
      std::string    filePrefix_;
 
   };

}
#endif

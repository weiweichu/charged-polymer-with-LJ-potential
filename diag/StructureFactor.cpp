#ifndef STRUCTUREFACTOR_CPP
#define STRUCTUREFACTOR_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iomanip>
#include <map>
#include "StructureFactor.h"

namespace GridMC
{

   /*
   * Constructor.
   */
   StructureFactor::StructureFactor(System& system) :
      Diagnosis(system),
      nEquilibration_(0),
      tA_(0),
      tB_(0),
      maxIndex_(0),
      nWaves_(0),
      waveIndexSq_(),
      waveNumber_(),
      waveIndex_(),
      waves_(),
      cA_(),
      sA_(),
      cB_(),
      sB_(),
      S_()
   { name_.assign("StructureFactor"); }

   /*
   * Default destructor.
   */
   StructureFactor::~StructureFactor()
   {}

   /*
   * Read sampling interval.
   */
   void StructureFactor::readParam(std::istream& in)
   {
      Diagnosis::readParam(in);

      // Read number of wave and wave indices.
      char        comment[200];
      std::string line;

      // Number of equilibration steps.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: equilibration steps");
      sscanf(line.c_str(), "%d %s", &nEquilibration_, comment);

      // Maximum wave index.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: maximum wave index");
      sscanf(line.c_str(), "%d %d %s", &tA_, &tB_, comment);

      // Maximum wave index.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: maximum wave index");
      sscanf(line.c_str(), "%d %s", &maxIndex_, comment);

      if (maxIndex_ <= 0)
         UTIL_THROW("maximum wave index must be greater than 0");

      // Generate all equivalent waves.
      generateWaves();

      // Resize data arrays.
      cA_.resize(nWaves_);
      sA_.resize(nWaves_);
      cB_.resize(nWaves_);
      sB_.resize(nWaves_);
      int i, j;
      for (i = 0; i < nWaves_; ++i) {
         j = waveIndex_[i].size();
         cA_[i].resize(j, 0.0);
         sA_[i].resize(j, 0.0);
         cB_[i].resize(j, 0.0);
         sB_[i].resize(j, 0.0);
      }
      S_.resize(nWaves_, 0.0);
   }

   /*
   * Write parameter.
   */
   void StructureFactor::writeParam(std::ostream& out) const
   {
      Diagnosis::writeParam(out);

      // Number of waves.
      out << "Equilibration steps         " << nEquilibration_ << std::endl;
      out << "Mode types                  " << tA_ << "  " << tB_ << std::endl;
      out << "Max wave index              " << maxIndex_ << std::endl;
      out << "Number of distinct waves    " << nWaves_ << std::endl;

      // Wave indices.
      std::set<int>::iterator itr;
      int i(-1);
      for (itr = waveIndexSq_.begin(); itr != waveIndexSq_.end(); ++itr) {
         out << "                            ";
         out << "wave index squared: " << *itr;
         out << " (" << waveIndex_[++i].size() << " waves)" << std::endl;
         for (unsigned j = 0; j < waveIndex_[i].size(); ++j) {
            out << "                            ";
            out << setw(4) << waveIndex_[i][j][0] << "  ";
            out << setw(4) << waveIndex_[i][j][1] << "  ";
            out << setw(4) << waveIndex_[i][j][2] << std::endl;
         }
      }

   }

   /*
   * Sample.
   */
   void StructureFactor::sample(const long iStep)
   {
      if (isAtInterval(iStep) == false || iStep <= nEquilibration_) return;

      // Reset fourier mode values.
      int      i;
      unsigned j;
      for (i = 0; i < nWaves_; ++i) {
         for (j = 0; j < waves_[i].size(); ++j) {
            cA_[i][j] = 0.0;
            sA_[i][j] = 0.0;
            cB_[i][j] = 0.0;
            sB_[i][j] = 0.0;
         }
      }

      // Calculate fourier modes.
      unsigned k;
      int      t;
      Vector   R;
      double   rdotq;
      int      tA1(tA_), tA2(-1), tB1(tB_), tB2(-1);

      if (tA_ < 0) {
         if (tA_ == -1) {
            tA1 = 0;
            tA2 = 1;
         } else if (tA_ == -2) {
            tA1 = 2;
            tA2 = 3;
         } else {
            UTIL_THROW("Invalid type index in structure factor calculation");
         }
      }

      if (tB_ < 0) {
         if (tB_ == -1) {
            tB1 = 0;
            tB2 = 1;
         } else if (tB_ == -2) {
            tB1 = 2;
            tB2 = 3;
         } else {
            UTIL_THROW("Invalid type index in structure factor calculation");
         }
      }


      for (k = 0; k < beads_.size(); ++k) {
         t = beads_[k].t; 
         R = beads_[k].r; 

         // First component for structure factor.
         if (t == tA1 || t == tA2) {

            for (i = 0; i < nWaves_; ++i) {
               for (j = 0; j < waves_[i].size(); ++j) {
                  rdotq = R.dot(waves_[i][j]);
                  cA_[i][j] += cos(rdotq);
                  sA_[i][j] += sin(rdotq);
               }
            }

         // Second component for structure factor.
         } else if (t == tB1 || t == tB2) {

            for (i = 0; i < nWaves_; ++i) {
               for (j = 0; j < waves_[i].size(); ++j) {
                  rdotq = R.dot(waves_[i][j]);
                  cB_[i][j] += cos(rdotq);
                  sB_[i][j] += sin(rdotq);
               }
            }

         }
      }

      // Calculate the structure factor.
      for (i = 0; i < nWaves_; ++i) {
         for (j = 0; j < waves_[i].size(); ++j) {
            if (tB_ != tA_) {
               S_[i] += (cA_[i][j]*cB_[i][j] + sA_[i][j]*sB_[i][j]) / double(waves_[i].size()*beads_.size());
            } else {
               S_[i] += (cA_[i][j]*cA_[i][j] + sA_[i][j]*sA_[i][j]) / double(waves_[i].size()*beads_.size());
            }
         }
      }

      nSamples_ += 1;
   }


   /*
   * Output result.
   */
   void StructureFactor::output() const
   {
      std::ofstream out;
      fileMaster_.openOutputFile("Sq.dat", out);

      streamsize oldprec = out.precision(6);
      out.setf(ios::fixed, ios::floatfield);
      for (int i = 0; i < nWaves_; ++i) {
         out.setf(ios::left, ios::adjustfield);
         out << std::setw(4) << waveIndex_[i][0].square();
         out.unsetf(ios::adjustfield);
         out << std::setw(12) << waveNumber_[i];
         out << std::setw(15) << S_[i] / double(nSamples_) << std::endl;
      }
      out.precision(oldprec);
      out.unsetf(ios::floatfield);

      out.close();
   }

   /*
   * Generate all equivalent waves from the wave indices.
   *
   * Assuming a cubic simulation box.
   */
   void StructureFactor::generateWaves()
   {
      // Generate the set of wave index magnitudes.
      if (fabs(boxL_[0] - boxL_[1]) + fabs(boxL_[0] - boxL_[2]) > Constants::Epsilon)
         UTIL_THROW("simulation box is not of cubic shape");

      // Generate the set of wave index magnitudes.
      int ix, iy, iz, sqsum;

      for (ix = -maxIndex_; ix <= maxIndex_; ++ix) { 
         for (iy = -maxIndex_; iy <= maxIndex_; ++iy) { 
            for (iz = 0; iz <= maxIndex_; ++iz) { 
               sqsum = ix * ix + iy * iy + iz * iz;
               if (sqsum > 0)
                  waveIndexSq_.insert(sqsum);
            }
         }
      }

      // Map set to ordinal number.
      map<int,int> ids;
      std::set<int>::iterator itr;
      int id(-1);
      for (itr = waveIndexSq_.begin(); itr != waveIndexSq_.end(); ++itr)
         ids.insert(pair<int,int>(*itr, ++id));

      // Allocate containers.
      nWaves_ = waveIndexSq_.size();
      waveNumber_.resize(nWaves_, 0.0);
      waveIndex_.resize(nWaves_);
      waves_.resize(nWaves_);

      // Fill the wave number array.
      double basis = 2.0 * Constants::Pi / boxL_[0];
      for (itr = waveIndexSq_.begin(); itr != waveIndexSq_.end(); ++itr)
         waveNumber_[ids[*itr]] = basis * sqrt(double(*itr));

      // Fill the wave index array.
      for (ix = -maxIndex_; ix <= maxIndex_; ++ix) { 
         for (iy = -maxIndex_; iy <= maxIndex_; ++iy) { 

            // To remove degeneracy; consider iz == 0 separately.
            sqsum = ix * ix + iy * iy;
            if (sqsum > 0) {
               if (iy == 0) {
                  if (ix > 0)
                     waveIndex_[ids[sqsum]].push_back(IntVector(ix,0,0));
               } else if (iy > 0) {
                  waveIndex_[ids[sqsum]].push_back(IntVector(ix,iy,0));
               }
            }

            // iz == 0;
            for (iz = 1; iz <= maxIndex_; ++iz) { 
               sqsum = ix * ix + iy * iy + iz * iz;
               if (sqsum > 0)
                  waveIndex_[ids[sqsum]].push_back(IntVector(ix,iy,iz));
            }
         }
      }

      // Fill the wave array.
      Vector wave;
      for (int i = 0; i < nWaves_; ++i) {
         if (waveIndex_[i].size() <= 0)
            UTIL_THROW("wave index not constructed properly");
         for (unsigned j = 0; j < waveIndex_[i].size(); ++j) {
            wave[0] = basis * double(waveIndex_[i][j][0]);
            wave[1] = basis * double(waveIndex_[i][j][1]);
            wave[2] = basis * double(waveIndex_[i][j][2]);
            waves_[i].push_back(wave);
         }
      }
   }

}

#endif

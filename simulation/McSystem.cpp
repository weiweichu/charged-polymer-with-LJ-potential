#ifndef MCSYSTEM_CPP
#define MCSYSTEM_CPP

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <utility>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

#include "../util/global.h"
#include "McSystem.h"
#include "Simulation.h"
#include <iostream>

namespace GridMC
{
   /*
   * Constructor.
   */
   McSystem::McSystem(const Simulation& simulationIn) :
      System(simulationIn),
      // interaction parameter
      chiN_(),
      kappaN_(1.0),
      sqrtNbar_(1.0),
      esStrength_(0.0),
      esRatio_(1.0),
      aBorn_(0.01),
      cutoff_(1.0),
      sigma_(1.0),
      epsi_(1.0),
      // barostat parameter
      isIsobaric_(false),
      barostatPressure_(0.0),
      // density of states parameter
      doDOS_(false),
      dosWeight_(),
      // system status parameter
      bondEnergy_(0.0),
      twobodyEnergy_(0.0),
      coulombEnergy_(0.0),
      solvationEnergy_(0.0),
      ljEnergy_(0.0),
      pV_(0.0),
      // system status parameter change
      dBondEnergy_(0.0),
      dTwobodyEnergy_(0.0),
      dCoulombEnergy_(0.0),
      dSolvationEnergy_(0.0),
      dLjEnergy_(0.0),
      // order parameter
      psi_(*this),
      i030(0),
      i001(0)
   {}

   /*
   * Default destructor.
   */
   McSystem::~McSystem()
   {}

   /*
   * Read parameters.
   */
   void McSystem::readParam(istream& in)
   {
      // Molecular topology parameters.
      System::readParam(in);

      // Auxiliary varialbe for reading file.
      char   comment[200];
      string line;

      // Parse header.
      getline(in, line);
      getline(in, line);
      if (line.find("Thermodynamic") == string::npos)
         UTIL_THROW("reading error: Thermodynamic parameter header.");

      // Read chi-parameter.
      vector<double> zero(nType_, 0.0);
      chiN_.resize(nType_, zero);

      int i, j;
      for (i = 1; i < nType_; ++i) {
         for (j = 0; j < i; ++j) {
            in >> chiN_[i][j];
            chiN_[j][i] = chiN_[i][j];
         }
         getline(in, line);
      }

      // Read kappaN.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: kappaN");
      sscanf(line.c_str(), "%lf %s", &kappaN_, comment);

      // Read sqrt Nbar.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: sqrtNbar");
      sscanf(line.c_str(), "%lf %s", &sqrtNbar_, comment);
      //Log::file() << "sqrtNbar  " << sqrtNbar_ << endl;
      // Read esStrength (electrostatic interaction strength).
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: esStrength");
      sscanf(line.c_str(), "%lf %lf %lf %s", &esStrength_, &esRatio_, &aBorn_, comment);
      // Read lj parameter(cutoff radius, sigma, epsi)
      getline(in, line);
      if (line.size() <= 0)
        UTIL_THROW("reading error: lj parameter");
      sscanf(line.c_str(), "%lf %lf %lf %s", &cutoff_, &sigma_, &epsi_, comment);
      // Read flag for barostat.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: barostat flag");
      int baroflag;
      sscanf(line.c_str(), "%d %s", &baroflag, comment);
      isIsobaric_ = (baroflag == 0 ? false : true); 
      // Read barostat pressure.
      if (isIsobaric_) {
         getline(in, line);
         if (line.size() <= 0)
            UTIL_THROW("reading error: barostat pressure");
         sscanf(line.c_str(), "%lf %s", &barostatPressure_, comment);
      }

      // Read parameters for order parameter.
      psi_.readParam(in);

      // Read flag for DOS sampling.
      getline(in, line);
      if (line.size() <= 0)
         UTIL_THROW("reading error: density of states flag");

      int dosflag;
      sscanf(line.c_str(), "%d %s", &dosflag, comment);
      doDOS_ = (dosflag <= 0 ? false : true); 

      // Read range of order parameter for DOS sampling and/or prescribed bias weight.
      if (doDOS_) {
         dosWeight_.readParam(in);

         if (dosflag > 1) {
            string   weightFileName;
            ifstream weightFile;
            stringstream itrstr;
            itrstr << dosflag - 1;

            weightFileName = string("weight.");
            weightFileName += itrstr.str();
            fileMasterPtr_->openInputFile(weightFileName, weightFile);
            dosWeight_.readWeight(weightFile, dosflag - 1);
            weightFile.close();

            dosWeight_.shiftWeight();
            dosWeight_.resetHistogram();
         }
      }

      // Set grid interaction and allocate arrays.
      grid_.setTwobodyInteraction(chiN_, kappaN_, sqrtNbar_, nAB_);
      grid_.setCoulombStrength(esStrength_, esRatio_);
      grid_.setLjInteraction(cutoff_);

      // Set interaction strength for Ewald.
      if (useEwald_ > 0)
         ewald_.setCoulombStrength(esStrength_);

      initializeVTK();
      updateVTK(0);
      solvePoisson();

      // Calculate status variable and update DOS.
      calculateStatusVariable();
      if (doDOS_)
         updateDOS();
   }

   /*
   * Write parameters.
   */
   void McSystem::writeParam(ostream& out)
   {
      System::writeParam(out);
      out << "----- Thermodynamic parameters ----- " << endl;

      out.setf(ios::left, ios::adjustfield);
      out << "chiN parameters             " << setw(10) << chiN_[0][1] << endl;
      for (int i = 2; i < nType_; ++i) {
         out << "                            ";
         for (int j = 0; j < i; ++j)
            out << setw(10) << chiN_[i][j];
         out << endl;
      }
      out.unsetf(ios::adjustfield);

      out << "incompressibility           " << kappaN_ << endl;
      out << "sqrtNbar                    " << sqrtNbar_ << endl;
      out << "inverse dielectric const    " << esStrength_ << "  " << esRatio_ << endl;
      out << "ion radius                  " << aBorn_ << endl;
      out << "cutoff radius               " << cutoff_ << endl;
      out << "lj sigma                    " << sigma_ << endl;
      out << "lj epsilon                  " << epsi_ << endl;
      out << "isIsobaric                  " << isIsobaric_ << endl;
      if (isIsobaric_)
         out << "barostat pressure           " << barostatPressure_ << endl;
      psi_.writeParam(out);
      out << "doDOS                       " << doDOS_ << endl;
      if (doDOS_)
         dosWeight_.writeParam(out);
   }


   /*
   * Check the flatness of the order parameter histogram.
   */
   void McSystem::checkDOSflatness()
   {
      if (dosWeight_.isFlat()) {
         int itr = dosWeight_.getNIteration();

         // Check if the folder exists at the first iteration.
         if (itr > 0) {
            struct stat buf;
            string name("doshistory");
            fileMasterPtr_->prependPath(name);
            if (stat(name.c_str(), &buf) != 0 || !S_ISDIR(buf.st_mode)) {
               if (mkdir(name.c_str(), S_IRWXU | S_IRGRP | S_IROTH) != 0)
                  UTIL_THROW("Folder doshistory doesn't exist or can't be created.");
            }
         }

         // Output histogram.
         ofstream dosFile;
         stringstream itrstr;
         itrstr << itr;

         string fileName("doshistory/histogram.");
         fileName += itrstr.str();
         fileMasterPtr_->openOutputFile(fileName, dosFile);
         dosWeight_.outputHistogram(dosFile);
         dosFile.close();

         // Output weight. 
         dosWeight_.shiftWeight();
         fileName = string("doshistory/weight.");
         fileName += itrstr.str();
         fileMasterPtr_->openOutputFile(fileName, dosFile);
         dosWeight_.outputWeight(dosFile);
         dosFile.close();

         // Reset histogram and weight if needed.
         dosWeight_.resetHistogram();
      }
   }

   /*
   * Evaluate system energies.
   */
   void McSystem::updateDOS()
   { dosWeight_.sample(psi_.get(), simulation().nFlatnessInterval_); }

   /*
   * Evaluate system energies.
   */
   void McSystem::calculateStatusVariable()
   {
      // Bond energy.
      int       i, j;
      Vector    r1, r2, dr;
      Particle  *ptr;

      bondEnergy_ = 0.0;
      for (i = 0; i < nPolymers_; ++i) {
         ptr = polymer_[i];
         r1 = ptr[0].r;
         for (j = 1; j < nAB_; ++j) {
            r2 = ptr[j].r;
            dr.subtract(r2, r1);
            pbcShift(dr);
            bondEnergy_ += dr.square();
            r1 = r2;
         }
      }
      bondEnergy_ *= (1.5 / bondL_ / bondL_);
      bondEnergy_ -= 1.5*double(nPolymerBeads_ - nPolymers_);

      // Two-body energy.
      twobodyEnergy_  = grid_.getTwobodyEnergy();
      twobodyEnergy_ -= getRandomMixingEnergy();

      // Coulomb energy.
      if (useEwald_ > 0)
         coulombEnergy_ = ewald_.getCoulombEnergy(beads_);
      else
         coulombEnergy_ = grid_.getCoulombEnergy();

      // Born solvation energy.
      double Qi, permittivity;


      solvationEnergy_ = 0.0;
      for (i = 0; i < int(beads_.size()); ++i) {
         Qi = beads_[i].q;
         if (fabs(Qi) > Constants::Epsilon) {
            permittivity = grid_.getDielectricPermittivity(beads_[i].r);
            solvationEnergy_ += Qi * Qi / permittivity;
         }
      }
      solvationEnergy_ *= esStrength_ / (2.0 * aBorn_);

      // LJ energy
      int total = 0;
      double distotal = 0.0;
      ljEnergy_ = 0.0;
      for (i = 0; i < int(beads_.size()); ++i) {
        for (j = 0; j < i; ++j) {
          if (fabs(beads_[i].q) > Constants::Epsilon && fabs(beads_[j].q) > Constants::Epsilon) {
            total ++;
            r1 = beads_[i].r;
            r2 = beads_[j].r;
            toPrimaryCell(r1);
            toPrimaryCell(r2);
            for (int m = 0; m < 3; m++) {
              dr[m] = fabs(r1[m] - r2[m]);
              if (dr[m] > boxL_[m] / 2.0) {
                dr[m] = boxL_[m] - dr[m];
                if (dr[m] > boxL_[m] / 2.0) {
                  Log::file() << "distance over half box" << endl;
                }
              }
            }
            double r = dr.abs();
            distotal += r;
            if (r < cutoff_) {
              ljEnergy_ += 4 * epsi_ * (pow(sigma_ / r, 12) - pow(sigma_ / r, 6));
            }
          }
        }
      }
      // Pressure.
      pV_ = -2.0*bondEnergy_/3.0 + twobodyEnergy_ + coulombEnergy_/3.0 + ljEnergy_;

      // Order parameter.
      psi_.calculate();
   }

   /*
   * Update coulomb energy (does not apply to Ewald summation).
   */
   void McSystem::updateCoulombEnergy(const double deltaCoulombE)
   {
      coulombEnergy_ += deltaCoulombE;
   }

   /*
   * Update bead position, grid density, and system energy.
   */
   void McSystem::updateStatusVariable(const int id, const Vector& rTrial)
   {
      if (fabs(beads_[id].q) > Constants::Epsilon) {
        grid_.updateNeighbor(id, beads_[id].r, rTrial);
      }
      beads_[id].r = rTrial;
      grid_.updateFromList(beads_[id]);
      bondEnergy_      += dBondEnergy_;
      twobodyEnergy_   += dTwobodyEnergy_;
      coulombEnergy_   += dCoulombEnergy_;
      solvationEnergy_ += dSolvationEnergy_;
      ljEnergy_        += dLjEnergy_;
      if (useEwald_ > 0 && fabs(beads_[id].q) > Constants::Epsilon)
         ewald_.updateFourierModes();
      pV_   = -2.0*bondEnergy_/3.0 + twobodyEnergy_ + coulombEnergy_/3.0 + ljEnergy_;
      psi_.updateByChange();

      //Log::file() << "dE: b nb c s lj    " << dBondEnergy_<< "  "<< dTwobodyEnergy_ <<"  "<<dCoulombEnergy_ <<"  "<<dSolvationEnergy_ <<"  "<<dLjEnergy_ << endl;
      /*
      if (psi_.get() > 0.25 && psi_.get() < 0.35 && i030 == 0) {
         std::ofstream outconfig;
         fileMasterPtr_->openOutputFile("config.030", outconfig);
         writeConfig(outconfig);
         outconfig.close();
         i030 = 1;
      } else if (psi_.get() > 0.0 && psi_.get() < 0.01 && i001 == 0) {
         std::ofstream outconfig;
         fileMasterPtr_->openOutputFile("config.001", outconfig);
         writeConfig(outconfig);
         outconfig.close();
         i001 = 1;
      }
      */
   }

   /*
   * Calcualte the statistical weight associated with a bead move.
   */
   double McSystem::calculateMoveWeight(const int id, const Vector& rTrial)
   {
      // Bond energy change.
      if (id < nPolymerBeads_)
         dBondEnergy_ = getBondEnergyChange(id, rTrial);
      else
         dBondEnergy_ = 0.0;

      // Nonbonded energy change: two body and Coulomb.
      grid_.findAffectedSites(beads_[id].r, rTrial);
      dTwobodyEnergy_ = grid_.getTwobodyEnergyChange(beads_[id].t);
      if (useEwald_ > 0)
         dCoulombEnergy_ = ewald_.getCoulombEnergyChange(beads_, id, rTrial);
      else
         dCoulombEnergy_ = grid_.getCoulombEnergyChange(beads_[id].q);
 
      dLjEnergy_ = getLjEnergyChange(id, rTrial);
      double Qi = beads_[id].q;
      // Solvation energy change.
      if (fabs(Qi) > Constants::Epsilon) {
         double permittivity;
         permittivity = grid_.getDielectricPermittivity(beads_[id].r);
         dSolvationEnergy_ = -esStrength_ * Qi * Qi / (2.0 * permittivity * aBorn_);

         permittivity = grid_.getDielectricPermittivity(rTrial);
         dSolvationEnergy_ += esStrength_ * Qi * Qi / (2.0 * permittivity * aBorn_);
      }
      else{
           dSolvationEnergy_ = 0.0;
      } 

      // Order parameter change, which used the list of affected sites.
      double dpsi = psi_.getChange(beads_[id], rTrial);

      // Order parameter change, which used the list of affected sites.
      double dWeight(0.0);
      if (doDOS_) {
         double psi = psi_.get();
         dWeight  = dosWeight_.getWeight(psi + dpsi);
         dWeight -= dosWeight_.getWeight(psi);
      }

      return (dBondEnergy_ + dTwobodyEnergy_ + dCoulombEnergy_ + dSolvationEnergy_ + dLjEnergy_ + dWeight);
   }
 
   /*
   * Bond energy change.
   */
   double McSystem::getBondEnergyChange(const int id, const Vector& rTrial) const
   {
      if (id >= nPolymerBeads_)
         UTIL_THROW("A solvent bead can't have bond energy");

      double  deltaE(0.0);
      Vector  r1, r2, dr;

      // Identify bead's relative position on polymer.
      int beadId = id % nAB_;
      r1 = beads_[id].r;

      // Calculate old energy.
      if (beadId > 0) {
         r2 = beads_[id-1].r;
         dr.subtract(r2, r1);
         pbcShift(dr);
         deltaE -= dr.square();

         dr.subtract(r2, rTrial);
         pbcShift(dr);
         deltaE += dr.square();
      }

      if (beadId < nAB_ - 1) {
         r2 = beads_[id+1].r;
         dr.subtract(r2, r1);
         pbcShift(dr);
         deltaE -= dr.square();

         dr.subtract(r2, rTrial);
         pbcShift(dr);
         deltaE += dr.square();
      }

      return (deltaE * 1.5 / bondL_ / bondL_);
   }

   double McSystem::getLjEnergyChange(const int id, const Vector& rTrial)
   {
      double Qi = beads_[id].q;
      Vector r2, dr, r1;
      double dLjEnergy(0.0), dis(0.0);
      
//      double lj1(0.0), lj2(0.0);
//      int i, j;
//
//     for (i = 0; i < int(beads_.size()); i++) {
//        for (j = 0; j < i; j++) {
//          if (fabs(beads_[i].q) > Constants::Epsilon && fabs(beads_[j].q) > Constants::Epsilon) {
//            r1 = beads_[i].r;
//            r2 = beads_[j].r;
//            toPrimaryCell(r1);
//            toPrimaryCell(r2);
//            for (int m = 0; m < 3; m++) {
//              dr[m] = fabs(r1[m] - r2[m]);
//              if (dr[m] > boxL_[m] / 2.0) {
//                dr[m] = boxL_[m] - dr[m];
//              }
//            }
//            dis = dr.abs();
//            if (dis < cutoff_) {
//              lj1 += 4 * epsi_ * (pow(sigma_ / dis, 12) - pow(sigma_ / dis, 6));
//            }
//          }
//        }
//      }
//
//   //  Log::file() << "lj1   " <<lj1 << endl;
//
//      for (i = 0; i < int(beads_.size()); i++) {
//        for (j = 0; j < i; j++) {
//          if (fabs(beads_[i].q) > Constants::Epsilon && fabs(beads_[j].q) > Constants::Epsilon) {
//            if (i == id) {
//              r1 = rTrial;
//            }
//            else {
//              r1 = beads_[i].r;
//            }
//            if (j == id) {
//              r2 = rTrial;
//            }
//            else {
//              r2 = beads_[j].r;
//            }
//            toPrimaryCell(r1);
//            toPrimaryCell(r2);
//            for (int m = 0; m < 3; m++) {
//              dr[m] = fabs(r1[m] - r2[m]);
//              if (dr[m] > boxL_[m] / 2.0) {
//                dr[m] = boxL_[m] - dr[m];
//              }
//            }
//            dis = dr.abs();
//            if (dis < cutoff_) {
//              lj2 += 4 * epsi_ * (pow(sigma_ / dis, 12) - pow(sigma_ / dis, 6));
//            }
//          }
//        }
//      }
//     //Log::file() << "lj2   " <<lj2 << endl;
//     double temp = lj2 - lj1;
//      return temp;
//
      dLjEnergy = 0.0;
      r1 = beads_[id].r;
      if (fabs(Qi) > Constants::Epsilon) {
        grid_.findNeighborList(r1);
        vector<int> neighborList = grid_.getNeighborList();
        for (int k = 0; k < int(neighborList.size()); k++) {

          if (neighborList[k] != id) {
            r2 = beads_[neighborList[k]].r;
            toPrimaryCell(r1);
            toPrimaryCell(r2);
            for (int m = 0; m < 3; m++) {
              dr[m] = fabs(r1[m] - r2[m]);
              if (dr[m] > boxL_[m] / 2.0) {
                dr[m] = boxL_[m] - dr[m];
              }
            }
            dis = dr.abs();
            if (dis < cutoff_){
              dLjEnergy -= 4 * epsi_ * (pow(sigma_/dis, 12) - pow(sigma_/dis, 6));
            } 
          }
        }
      }

      r1 = rTrial;
      if (fabs(Qi) > Constants::Epsilon) {
        grid_.findNeighborList(r1);
        vector<int> neighborList = grid_.getNeighborList();
        for (int k = 0; k < int(neighborList.size()); k++) {
          if (neighborList[k] != id) {
            r2 = beads_[neighborList[k]].r;
            toPrimaryCell(r1);
            toPrimaryCell(r2);
            for (int m = 0; m < 3; m++) {
              dr[m] = fabs(r1[m] - r2[m]);
              if (dr[m] > boxL_[m] / 2.0) {
                dr[m] = boxL_[m] - dr[m];
              }
            }
            dis = sqrt(dr.square());

            if (dis < cutoff_){
              dLjEnergy += 4 * epsi_ * (pow(sigma_/dis, 12) - pow(sigma_/dis, 6));
            } 
          }
        }
      }
      return dLjEnergy;
   }
   /*
   * Calculate the statistical weight associated with volume change.
   *
   * Note: No contribution from the order parameter.
   */
   double McSystem::volumeChangeWeight(const double deltaV)
   {
      double weight(0.0);
      double s = (boxV_ + deltaV) / boxV_;
      weight  = (pow(s, 2.0/3.0) - 1.0) * bondEnergy_;
      weight += (1.0/s - 1.0) * (twobodyEnergy_ + ljEnergy_);
      weight += (pow(s,-1.0/3.0) - 1.0) * coulombEnergy_;
      weight += barostatPressure_ * deltaV;
      // ideal gas contribution
      weight -= double(nPolymerBeads_ + nSolvents_)*log(s);
      return weight;
   }

   /*
   * Affinely update the system volume, and update the status variables accordingly.
   */
   void McSystem::updateVolume(const double deltaV)
   {
      double s = (boxV_ + deltaV) / boxV_;
      double s1 = pow(s, 1.0/3.0);

      // Box dimension.
      boxV_ += deltaV;
      if (boxV_ <= 0.0)
         UTIL_THROW("Negative volume for volume change.");
      boxL_ *= s1;

      // Bead position.
      for (int j = 0; j < nPolymerBeads_ + nSolvents_; ++j)
         beads_[j].r *= s1;

      // System energy. Order parameter doesn't change.
      bondEnergy_    *= pow(s, 2.0/3.0);
      twobodyEnergy_ /= s;
      coulombEnergy_ *= pow(s, -1.0/3.0);
      solvationEnergy_ /= s1;
      ljEnergy_ /= s;
      pV_ = -2.0*bondEnergy_/3.0 + twobodyEnergy_ + coulombEnergy_/3.0 + ljEnergy_;

      // Grid.
      grid_.rescaleBox(s1);

      // Ewald.
      if (useEwald_ > 0)
         ewald_.rescaleBox(s1);

      // Update wave vector in the order parameter.
      psi_.updateAfterBoxMove();
   }

}

#endif

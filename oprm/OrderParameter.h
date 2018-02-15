#ifndef ORDERPARAMETER_H
#define ORDERPARAMETER_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <util/global.h>

#include <util/Vector.h>
#include <simulation/System.h>
#include <grid/GridMassCharge.h>

namespace GridMC
{

   /**
   * This class provide a generic interface for constructing system speciifc
   * order parameters (double number).
   *
   * \ingroup Oprm_Module
   */
   class OrderParameter
   {
   
   public:
  
      /**
      * Default constructor.
      */
      OrderParameter(const System& sysIn);

      /**
      * Default destructor.
      */
      virtual ~OrderParameter();

      /**
      * Read parameter.
      */
      virtual void readParam(std::istream& in) = 0;

      /**
      * Write parameter.
      */
      virtual void writeParam(std::ostream& out) const = 0;

      /**
      * Calculate order parameter.
      */
      virtual void calculate() = 0;

      /**
      * Calculate orde parameter change due to bead move.
      */
      virtual double getChange(const Particle& bead, const Vector& rTrial) = 0;

      /**
      * Reset parameters after box dimension is rescaled.
      */
      virtual void updateAfterBoxMove() = 0;

      /**
      * Update order parameter by dpsi.
      */
      virtual void updateByChange();

      /**
      * Set order parameter value.
      */
      void set(const double psiIn);

      /**
      * Accessor to order parameter.
      */
      double get() const;

      /**
      * Backup status variables.
      */
      virtual void backup();

      /**
      * Restore status variables.
      */
      virtual void restore();

   protected:

      /// Number of bead types in the system.
      const Vector& boxL_;

      /// Reference to parent system.
      const IntVector& nGrid_;

      /// Reference to parent system.
      const GridMassCharge& grid_;

      /// Pointer to the beads array.
      const vector<Particle>& beads_;

      /// Number of beads.
      unsigned nBead_;

      /// Order parameter.
      double psi_;
      double psiBak_;

      /// Change in order parameter (most useful in MC simulation).
      double dpsi_;

      /**
      * Calculate orde parameter using grid density.
      */
      void oprmFromGrid();

      /**
      * Calculate orde parameter change using grid density.
      */
      void oprmChangeFromGrid(const int type);

   };

   /*
   * Accessor to the order parameter.
   */
   inline void OrderParameter::updateByChange()
   { psi_ += dpsi_; }

   /*
   * Accessor to the order parameter.
   */
   inline void OrderParameter::set(const double psiIn)
   { psi_ = psiIn; }

   /*
   * Accessor to the order parameter.
   */
   inline double OrderParameter::get() const
   { return psi_; }

   /*
   * Backup status variables.
   */
   inline void OrderParameter::backup()
   {
      psiBak_  = psi_;
   }

   /*
   * Restore status variables.
   */
   inline void OrderParameter::restore()
   {
      psi_  = psiBak_;
   }

}
#endif

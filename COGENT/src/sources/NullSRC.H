#ifndef  _VIRTUALSRC_H_
#define  _VIRTUALSRC_H_

#include "KineticSpecies.H"
#include "REAL.H"
#include "SRCInterface.H"
#include "ParmParse.H"
#include <sstream>

#include "NamespaceHeader.H"

/**
 * NullSRC neutrals operator class.
*/
class NullSRC
   : public SRCInterface
{
   public:
      
      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      NullSRC() {;}
      
      /// Destructor.
      /**
       */
      virtual ~NullSRC() {;}
  
      /// Accumulates the RHS of the neutral operator.
      /**
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] species index of species to which neutral is applied.
       * @param[in] time current time.
       */
      virtual void evalSrcRHS( KineticSpecies&                 rhs_species,
                               const KineticSpeciesPtrVect&    kinetic_species_phys,
                               const CFG::FluidSpeciesPtrVect& fluid_species_phys,
                               const int                       species,
                               const Real                      time ) {;}
  
      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDtExplicitTI( const KineticSpeciesPtrVect& soln ) { return DBL_MAX; }
      Real computeDtImExTI( const KineticSpeciesPtrVect& soln ) { return DBL_MAX; }
};


#include "NamespaceFooter.H"

#endif

#ifndef _GKSOURCES_H_
#define _GKSOURCES_H_

#include "REAL.H"
#include "SRCInterface.H"
#include "KineticSpecies.H"
#include "ParmParse.H"

#include <map>

#include "NamespaceHeader.H"

/**
 * Gyrokinetic sources operator class.
*/
class GKSources
{
   public:

      /// Constructor with initialization.
      /**
       *  @param[in] pp the input database.
       */
      GKSources( const int verbose );

      /// Destructor.
      /**
       */
      virtual ~GKSources();

      /// Accumulates the RHS of the sources operator.
      /**
       * @param[out] rhs data holder for rhs (without ghost cells).
       * @param[in] soln current solution with which to evaluate RHS.
       * @param[in] time current time.
       */
      virtual void accumulateRHS( KineticSpeciesPtrVect&          rhs,
                                  const KineticSpeciesPtrVect&    kinetic_species_phys,
                                  const CFG::FluidSpeciesPtrVect& fluid_species_phys,
                                  const Real                      time );

      /// Compute a stable time step.
      /**
       * Computes and returns an estimate of the maximum stable time step.
       */
      Real computeDtExplicitTI( const KineticSpeciesPtrVect& soln );
      Real computeDtImExTI( const KineticSpeciesPtrVect& soln );

      /// Compute the time scale.
      /**
       * Computes and returns an estimate of the time scale.
       */
      Real computeTimeScale( const KineticSpeciesPtrVect& soln );

      /// returns the source model associated with the input name
      /**
       * @param[in] name String name of the species.
       */
      SRCInterface& sourceModel( const std::string& name );
  
      void preTimeStep(const KineticSpeciesPtrVect&, const Real, const KineticSpeciesPtrVect&);

   private:

      // prevent copying
      GKSources( const GKSources& );
      const GKSources& operator=( const GKSources& );

      bool m_verbose;
      std::map<std::string,int> m_species_map;
      std::vector<SRCInterface*> m_source_model;
      std::map<std::string,int> m_source_model_name;
};

#include "NamespaceFooter.H"

#endif

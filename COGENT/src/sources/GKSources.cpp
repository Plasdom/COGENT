#include "GKSources.H"

#include "SRCInterface.H"
// #include "FixedBckgr.H"
#include "PrescribedSources.H"
// #include "FluidNeutrals.H"
#include "NullSRC.H"

#include <float.h>
#include <sstream>

#include "NamespaceHeader.H"

GKSources::GKSources( const int a_verbose )
   : m_verbose(a_verbose)
{
   bool more_kinetic_species(true);
   int count(0);
   while (more_kinetic_species) {

      // look for data specifying another kinetic species
      std::stringstream s;
      s << "kinetic_species." << count+1;
      ParmParse ppspecies( s.str().c_str() );

      std::string species_name("Invalid");
      if (ppspecies.contains("name")) {
         ppspecies.get("name", species_name);
         
         std::string src_type("None");
         SRCInterface* src(NULL);
         
         if (ppspecies.contains( "src" )) {
            ppspecies.get( "src", src_type );
            const std::string prefix( "SRC." + species_name );
            ParmParse ppsrc( prefix.c_str() );
         
            if (src_type == "PrescribedSources") {
               src = new PrescribedSources( ppsrc, m_verbose );
            }
           
         }
         else {
            if (procID()==0) {
              cout << "Unrecognized source model; setting model to NULL" << endl;
            }
            src = new NullSRC();
         }
         
         m_source_model.push_back( src );
         typedef std::map<std::string,int>::value_type valType;
         m_species_map.insert( valType( species_name, count ) );
	      m_source_model_name.insert( valType( src_type, count ) );
         count++;
      }
      else {
         more_kinetic_species = false;
      }
   }
}


GKSources::~GKSources()
{
   for (int i(0); i<m_source_model.size(); i++ ) {
      delete m_source_model[i];
   }
}


SRCInterface& GKSources::sourceModel( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index((*it).second);
   return *(m_source_model[index]);
}


void GKSources::accumulateRHS( KineticSpeciesPtrVect&            a_rhs,
                                const KineticSpeciesPtrVect&      a_kinetic_species_phys,
                                const CFG::FluidSpeciesPtrVect&   a_fluid_species_phys,
                                const Real                        a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      SRCInterface& SRC( sourceModel( species_name ) );
      SRC.evalSrcRHS( rhs_species,
                      a_kinetic_species_phys,
                      a_fluid_species_phys,
                      species,
                      a_time );
   }
}


Real GKSources::computeDtExplicitTI( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
  int count(0);
  std::map<std::string,int>::iterator it;
  for (it=m_source_model_name.begin(); it!=m_source_model_name.end(); ++it) {
    Real tmp = m_source_model[it->second]->computeDtExplicitTI(soln,it->second);
    dt = (tmp < dt ? tmp : dt);
    count++;
  }
  return (count ? dt : -1);

}

Real GKSources::computeDtImExTI( const KineticSpeciesPtrVect& soln )
{
  Real dt(DBL_MAX);
  int count(0);
  std::map<std::string,int>::iterator it;
  for (it=m_source_model_name.begin(); it!=m_source_model_name.end(); ++it) {
    Real tmp = m_source_model[it->second]->computeDtImExTI(soln,it->second);
    dt = (tmp < dt ? tmp : dt);
    count++;
  }
  return (count ? dt : -1);

}

Real GKSources::computeTimeScale( const KineticSpeciesPtrVect& soln )
{
  std::map<std::string,int>::iterator it;
  Real scale = DBL_MAX;
  int count = 0;
  for (it=m_source_model_name.begin(); it!=m_source_model_name.end(); ++it) {
    Real tmp = m_source_model[it->second]->TimeScale(soln,it->second);
    scale = (tmp < scale ? tmp : scale);
    count++;
  }
  return (count ? scale : -1);

}

void GKSources::preTimeStep( const KineticSpeciesPtrVect& a_soln,
                              const Real a_time,
                              const KineticSpeciesPtrVect& a_soln_physical )

{
  for (int species(0); species<a_soln.size(); species++) {
    KineticSpecies&           soln_species(*(a_soln[species]));
    const std::string         species_name(soln_species.name());
    SRCInterface& SRC( sourceModel( species_name ) );
    SRC.preTimeStep(a_soln, species, a_time, a_soln_physical);
  }
}


#include "NamespaceFooter.H"

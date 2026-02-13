#include <math.h>

#include "FluidNeutrals.H"

#include "Directions.H"
#include "PhaseGeom.H"
#include "PhaseBlockCoordSys.H"
#include "NeutralsF_F.H"
#include "KineticFunctionLibrary.H"
#include "MaxwellianKineticFunctionF_F.H"
#include "KineticFunctionUtils.H"
#include "CollisionsF_F.H"
#include "MomentOp.H"
#include "ConstFact.H"
#include "Kernels.H"
#include "inspect.H"
#include "GKDiagnostics.H"

#include "NamespaceHeader.H" 


FluidNeutrals::FluidNeutrals( ParmParse& a_ppntr, const int a_verbosity )
   : m_verbosity(a_verbosity),
     m_neutr_vel(NULL),
     m_neutr_temp(NULL),
     m_first_call_tscale(true)
{
   parseParameters( a_ppntr );

   if (m_verbosity>0) {
      printParameters();
   }
   m_first_call = true;
}


FluidNeutrals::~FluidNeutrals()
{
}


void FluidNeutrals::evalNtrRHS(KineticSpecies&                   a_rhs_species,
                               const KineticSpeciesPtrVect&      a_soln,
                               const CFG::FluidSpeciesPtrVect&   a_fluid_species_phys,
                               const int                         a_species,
                               const Real                        a_time )
{
   /*
    Evaluates the ionization and charge-exchange contribution to the RHS as
    if m_iz_fixed_source = false:
       d(J*fB)/dt = J * (<sigmaV_ionization> * neutr_dens * fB
    if m_iz_fixed_source = true:
      d(J*fB)/dt = J * (<sigmaV_ionization> * neutr_dens * fB_neutral
 
    if m_chx_model_friction = false:
                  + <sigmaV_chx> (n_ion * fB_neutral - n_neutral * fB)
    if m_chx_model_friction = true:
                  + <sigmaV_chx> * n_neutral * vpar * upar * mass * fB_unshifted_Maxw / Ti)
    */
   
   // Get physical solution distribution function (Bstar_par*dfn) for the current species
   const KineticSpecies& soln_species( *(a_soln[a_species]) );
   const LevelData<FArrayBox>& soln_dfn( soln_species.distributionFunction() );

   //Get geometry
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );

   //Get neutral species
   int count(0);
   int neutral_species_index;
   for (int species(0); species<a_fluid_species_phys.size(); species++) {
      //const CFG::FluidSpecies& fluid_species( *(a_fluid_species_phys[species]) );
      const CFG::FluidSpecies& fluid_species( static_cast<CFG::FluidSpecies&>(  *(a_fluid_species_phys[species])) );
      if (fluid_species.isSpecies("neutrals")) {
         neutral_species_index = species;
         count++;
      }
   }
   if ( count != 1) MayDay::Error("FluidNeutrals:: input file must contain a single neutrals species");
   
   //Get neutral species density
   //const CFG::FluidSpecies& neutral_species( static_cast<CFG::FluidSpecies&>(  *(a_fluid_species_phys[neutral_species_index])) );
   const CFG::CFGVars& neutral_species = *(a_fluid_species_phys[neutral_species_index]);
   const CFG::LevelData<CFG::FArrayBox>& neutral_density_cfg( neutral_species.cell_var("density") );
   LevelData<FArrayBox> neutral_density_inj;
   phase_geom.injectConfigurationToPhase( neutral_density_cfg, neutral_density_inj);
   
   if (m_first_call) {

      if (m_include_line_radiation)
      {
         computeLineradClsFreq(soln_species, phase_geom, soln_dfn.getBoxes(), m_linerad_cls_freq);
         // phase_geom.plot
         // GKDiagnostics::GKDiagnostics dgn = GKDiagnostics::
         // GKDiagnostics::plotPhaseVar(m_linerad_cls_freq, 
         //                             phase_geom, 
         //                             "linerad_cls_freq", 
         //                             a_time );
      }      
      double dens_norm, vel_norm, temp_norm;
      computeChxNormalization(m_ionization_norm, m_chx_norm, dens_norm, vel_norm, temp_norm, m_SI_input);
      
      // Get electron temperature
      m_Te_cfg.define(mag_geom.grids(), 1, CFG::IntVect::Zero);
      if (m_use_Te_in_rate_calculation || m_include_iz_energy_loss)
      {
         // Identify the electron KineticSpecies in the a_soln input argument
         {
            for (int species(0); species<a_soln.size(); species++) {
               int species_charge = (*(a_soln[species])).charge();
               if (species_charge == -1)
               {
                  electron_species_index = species;
               }
            }
         }

         // Calculate the electron temperature (in eV) from the electron distribution
         (*(a_soln[electron_species_index])).temperature(m_Te_cfg);
      }
      else 
      {
         // Use the prescribed value from the input file
         m_electron_temp->assign( m_Te_cfg, mag_geom, a_time);
      }

      // Create the source distribution for ionization which accounts for energy loss of colliding particles (Maxwellian at Tiz, ne)
      if (m_include_iz_energy_loss)
      {
         m_iz_source_maxw.define( phase_geom.gridsFull(), 1, IntVect::Zero );
         m_iz_source_dfn.define( phase_geom.gridsFull(), 1, IntVect::Zero );
         CFG::LevelData<CFG::FArrayBox> iz_source_velocity_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero );
         for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
            iz_source_velocity_cfg[dit].setVal(0.);
         }
         double T_norm;
         ParmParse ppunits( "units" );
         ppunits.get("temperature",T_norm);     //[eV]
         m_iz_source_temperature_cfg.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
         for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
            m_iz_source_temperature_cfg[dit].copy(m_Te_cfg[dit]);
            m_iz_source_temperature_cfg[dit].divide(2.0);
            m_iz_source_temperature_cfg[dit].plus(-(E_iz/3.0)/T_norm);
         }
         m_iz_source_density_cfg.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
         for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
            m_iz_source_density_cfg[dit].setVal(1.0);
         }
         // soln_species.numberDensity( m_electron_density_cfg );

         MaxwellianKernel<FArrayBox> maxwellian(m_iz_source_density_cfg, m_iz_source_temperature_cfg, iz_source_velocity_cfg);
         maxwellian.eval(m_iz_source_maxw,soln_species);
      }

      // Get ionization reaction coefficient (sigmaV)
      CFG::LevelData<CFG::FArrayBox> ionization_sigmaV( mag_geom.grids(), 1, CFG::IntVect::Zero );
      computeIonizationSigmaV(ionization_sigmaV, a_time);

      phase_geom.injectConfigurationToPhase( ionization_sigmaV, m_ionization_sigmaV);
      
      if (m_neutr_temp != NULL) {
           CFG::LevelData<CFG::FArrayBox> neutral_temperature_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero );
           m_neutr_temp->assign( neutral_temperature_cfg, mag_geom, a_time);

           CFG::LevelData<CFG::FArrayBox> neutral_velocity_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero );
           CFG::LevelData<CFG::FArrayBox> neutral_density_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero );
           for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
              neutral_velocity_cfg[dit].setVal(0.);
              neutral_density_cfg[dit].setVal(1.);
           }
           
           m_neutral_dfn.define( phase_geom.gridsFull(), 1, IntVect::Zero );
           m_neutral_maxw.define( phase_geom.gridsFull(), 1, IntVect::Zero );
           MaxwellianKernel<FArrayBox> maxwellian(neutral_density_cfg, neutral_temperature_cfg, neutral_velocity_cfg);
           maxwellian.eval(m_neutral_maxw,soln_species);
      }

      if (m_include_chx) {
         m_neutral_velocity_cfg.define( mag_geom.grids(), 1, CFG::IntVect::Zero );
         m_neutr_vel->assign( m_neutral_velocity_cfg, mag_geom, a_time);
         for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
            m_neutral_velocity_cfg[dit].divide(vel_norm);
         }
      
      
        if (m_neutr_temp != NULL) {
           CFG::LevelData<CFG::FArrayBox> neutral_temperature_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero );
           m_neutr_temp->assign( neutral_temperature_cfg, mag_geom, a_time);

           for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
              neutral_temperature_cfg[dit].divide(temp_norm);
           }
         
           MaxwellianKernel<FArrayBox> maxwellian(neutral_density_cfg, neutral_temperature_cfg, m_neutral_velocity_cfg);
           maxwellian.eval(m_neutral_dfn,soln_species);

        }
      }      
   }
  
 
   //Calculate ionization RHS
   if (m_use_Te_in_rate_calculation)
   {
      (*(a_soln[electron_species_index])).temperature(m_Te_cfg);
      CFG::LevelData<CFG::FArrayBox> ionization_sigmaV( mag_geom.grids(), 1, CFG::IntVect::Zero );
      computeIonizationSigmaV(ionization_sigmaV, a_time);
      phase_geom.injectConfigurationToPhase( ionization_sigmaV, m_ionization_sigmaV);
   }
   LevelData<FArrayBox> ionization_rhs( soln_dfn.getBoxes(), soln_dfn.nComp(), IntVect::Zero );
   if (m_include_iz_energy_loss)
   {
      updateIzSourceDfn(soln_species);
      computeIonizationWithEnLossRHS( ionization_rhs, m_iz_source_dfn, soln_dfn, neutral_density_inj );
   }
   else 
   {
      if (m_neutr_temp== NULL) {
         computeIonizationRHS( ionization_rhs,  soln_dfn, neutral_density_inj );
      }
      else {
         updateNeutralDfn(soln_species);
         computeIonizationRHS( ionization_rhs,  m_neutral_dfn, neutral_density_inj );
      }
   }

   //Calculate line radiation RHS
   LevelData<FArrayBox> linerad_cls_rhs( soln_dfn.getBoxes(), soln_dfn.nComp(), IntVect::Zero );
   if (m_include_line_radiation)
   {
      evalLineradClsRHS( linerad_cls_rhs,  *(a_soln[a_species]), neutral_density_inj, a_time);  
      // evalLineradClsRHS( linerad_cls_rhs,  soln_species, a_time);  
   }
   
   //Calculate charge-exchange RHS
   LevelData<FArrayBox> chx_rhs( soln_dfn.getBoxes(), soln_dfn.nComp(), IntVect::Zero );
   if (m_include_chx) {
      
      //Get the size of ghost vector
      IntVect ghost = soln_dfn.ghostVect();
      CFG::IntVect ghost_cfg;
      ghost_cfg[0] = ghost[0];
      ghost_cfg[1] = ghost[1];
     
      //Get ion species density
      CFG::LevelData<CFG::FArrayBox> ion_density_cfg( mag_geom.grids(), 1, ghost_cfg );
      soln_species.numberDensity( ion_density_cfg );
      LevelData<FArrayBox> ion_density;
      phase_geom.injectConfigurationToPhase(ion_density_cfg, ion_density);

      //Get ion species temperature
      CFG::LevelData<CFG::FArrayBox> ion_parallelVel_cfg( mag_geom.grids(), 1, ghost_cfg );
      soln_species.parallelParticleFlux( ion_parallelVel_cfg );
      for (CFG::DataIterator dit(ion_density_cfg.dataIterator()); dit.ok(); ++dit) {
         ion_parallelVel_cfg[dit].divide(ion_density_cfg[dit]);
      }
      CFG::LevelData<CFG::FArrayBox> ion_temperature_cfg( mag_geom.grids(), 1, ghost_cfg );
      soln_species.pressure(ion_temperature_cfg, ion_parallelVel_cfg);
      for (CFG::DataIterator dit(ion_density_cfg.dataIterator()); dit.ok(); ++dit) {
        ion_temperature_cfg[dit].divide(ion_density_cfg[dit]);
      }
      LevelData<FArrayBox> ion_temperature;
      phase_geom.injectConfigurationToPhase(ion_temperature_cfg, ion_temperature);

      if (!m_chx_model_friction) {
         //Compute neutral dfn using the ion temperature
         if (m_neutr_temp == NULL) {
            MaxwellianKernel<FArrayBox> maxwellian(neutral_density_cfg, ion_temperature_cfg, m_neutral_velocity_cfg);
            maxwellian.eval(m_neutral_dfn, soln_species);
         }

         computeChargeExchangeRHS(chx_rhs,
                                  soln_dfn,
                                  ion_density,
                                  ion_temperature,
                                  neutral_density_inj,
                                  m_neutral_dfn );
      }

      if (m_chx_model_friction) {
         //Create unshifted ion Maxwellian and compute the parallel velocity difference
         CFG::LevelData<CFG::FArrayBox> zero_vel_cfg( mag_geom.grids(), 1, ghost_cfg );
         CFG::LevelData<CFG::FArrayBox> u_par_cfg( mag_geom.grids(), 1, ghost_cfg );
         for (CFG::DataIterator dit(ion_density_cfg.dataIterator()); dit.ok(); ++dit) {
            zero_vel_cfg[dit].setVal(0.0);
            u_par_cfg[dit].copy( m_neutral_velocity_cfg[dit]);
            u_par_cfg[dit].minus( ion_parallelVel_cfg[dit]);
         }
      
         LevelData<FArrayBox> ion_unshifted_maxw( phase_geom.gridsFull(), 1, IntVect::Zero );
         MaxwellianKernel<FArrayBox> maxwellian(ion_density_cfg, ion_temperature_cfg, zero_vel_cfg);
         maxwellian.eval(ion_unshifted_maxw, soln_species);
         
         LevelData<FArrayBox> u_par;
         phase_geom.injectConfigurationToPhase(u_par_cfg, u_par);
         
         double mass = soln_species.mass();
         computeModelChargeExchangeRHS(chx_rhs,
                                       neutral_density_inj,
                                       ion_temperature,
                                       u_par,
                                       ion_unshifted_maxw,
                                       mass,
                                       phase_geom );

      }
   }
   
   //Multiply by J to convert to the computational space
   phase_geom.multJonValid(ionization_rhs);
   phase_geom.multJonValid(chx_rhs);
   if (m_include_line_radiation){
      //TODO: Find out why this is necessary, and why it isn't used in ConsDragDiff (or maybe it is?)
      phase_geom.multJonValid(linerad_cls_rhs);
   }
   //TODO: Do we need to do this for linerad_cls_rhs too?

   //Add neutral-model RHS to the total RHS
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   for (DataIterator rdit(soln_dfn.dataIterator()); rdit.ok(); ++rdit) {
       rhs_dfn[rdit].plus( ionization_rhs[rdit] );
       if (m_include_chx) {
          rhs_dfn[rdit].plus( chx_rhs[rdit] );
       }
       if (m_include_line_radiation)
       {
         //  rhs_dfn[rdit].plus( linerad_cls_rhs[rdit] ); 
       }
   }  

   m_first_call = false;
}

void FluidNeutrals::computeLineradClsFreq(const KineticSpecies&     a_soln_species,
                                          const PhaseGeom&          a_phase_geom,
                                          const DisjointBoxLayout&  a_grids,
                                          LevelData<FArrayBox>&     a_cls_freq)
{
  
  // define collision frequency
  if (!m_linerad_cls_freq.isDefined()) m_linerad_cls_freq.define(a_grids, 1, IntVect::Zero);

  // Define other constants used in calculating collision frequency
  double ech = Constants::ELEMENTARY_CHARGE; 
  double pi = Constants::PI;  
  double mass = a_soln_species.mass();
  double mass_real = mass * Constants::MASS_OF_PROTON;
  double ni_norm = 1.0e19;
  double C = (8 * ech) / ((pow(pi,0.5))*ni_norm) * pow(2*ech/mass_real,gamma/2.0);
  double AbarC = Abar / C;

  
  // Compute normalisation constants
  double m0, m_norm, v_norm, cls_norm, T_norm, B_norm, L;
  ParmParse ppunits( "units" );
  ppunits.get("mass",m0);
  ppunits.get("temperature",T_norm);
  ppunits.get("magnetic_field",B_norm);
  ppunits.get("length",L);
  m_norm = m0 * Constants::MASS_OF_PROTON;
  v_norm = sqrt(Constants::ELEMENTARY_CHARGE * T_norm / m_norm);
  cls_norm = L/v_norm;

  // get injected magnetic field
  const LevelData<FArrayBox>& inj_B = a_phase_geom.getBFieldMagnitude();
  
  // computeSelfConsistFreq(m_cls_freq, m_dens_inj, m_temp_inj, mass, charge);
  const DisjointBoxLayout& grids = m_linerad_cls_freq.disjointBoxLayout();
  DataIterator dit(grids.dataIterator());
  const LevelData<FArrayBox>& B_injected( a_phase_geom.getBFieldMagnitude() );
  for (dit.begin(); dit.ok(); ++dit) {

    //TODO: Implement non-uniform velocity normalisation in same way as done in ConsDragDiff operator
    FArrayBox& this_cls_freq = m_linerad_cls_freq[dit];
    const FArrayBox& this_B = B_injected[dit];
    const FArrayBox& this_velocity = a_phase_geom.getVelocityRealCoords()[dit];
    FORT_COMPUTE_LINERAD_CLS_FREQ(CHF_BOX(this_cls_freq.box()),
                                   CHF_FRA1(this_cls_freq,0),
                                   CHF_CONST_FRA(this_velocity),
                                   CHF_CONST_FRA1(this_B,0),
                                   CHF_CONST_REAL(v_norm),
                                   CHF_CONST_REAL(mass),
                                   CHF_CONST_REAL(mass_real),
                                   CHF_CONST_REAL(ech),
                                   CHF_CONST_REAL(AbarC),
                                   CHF_CONST_REAL(alpha),
                                   CHF_CONST_REAL(beta),
                                   CHF_CONST_REAL(gamma),
                                   CHF_CONST_REAL(V0));

    this_cls_freq.mult(cls_norm);

  }
  
}

void FluidNeutrals::evalLineradClsRHS( LevelData<FArrayBox>&         a_rhs,
                                       KineticSpecies&         a_soln_species,
                                       const LevelData<FArrayBox>& a_neutral_density,
                                       const Real                    a_time )
{
  /*
    Evaluates a model collision operator for the cumulative effect of inelastic  
    excitation/de-excitation collisions between electrons and an atomic ion species 
    The model operator is (J. Roeltgen et al., Nuclear Fusion 65 (2025)):
    df/dt_coll=div(flux_coll), where flux_coll = [coll_freq*v*f].
    The implementation is based on the conservative drag-diffusion model operator, 
    implemented in COGENT as "ConsDragDiff" and outlined in Justin's thesis (Chapter 4).
  */

  CH_TIME("FluidNeutrals::evalLineradClsRHS");
  
  // get vlasov RHS for the current species
  LevelData<FArrayBox>& rhs_dfn = a_soln_species.distributionFunction();

  // get solution distribution function (f*J*Bstarpar) for the current species
  LevelData<FArrayBox>& soln_fBJ = a_soln_species.distributionFunction();
  double mass = a_soln_species.mass();
  double charge = a_soln_species.charge();

  // get coordinate system parameters
  const PhaseGeom& phase_geom = a_soln_species.phaseSpaceGeometry();
  const CFG::MagGeom & mag_geom = phase_geom.magGeom();
  const DisjointBoxLayout& dbl = soln_fBJ.getBoxes();
  
  // copy soln_fBJ so can perform exchange to ghost cells
  // we need to do so because we pass here a computational dfn
  // that does not have ghost cells filled
  if (!m_fBJ_vel_ghost.isDefined()) {
      IntVect velGhost = 2*IntVect::Unit;
      for (int dir=0; dir<CFG_DIM; dir++) {
         velGhost[dir] = 0;
      }
      m_fBJ_vel_ghost.define(dbl, 1, velGhost);
  }

  for (DataIterator dit(soln_fBJ.dataIterator()); dit.ok(); ++dit) {
      m_fBJ_vel_ghost[dit].copy(soln_fBJ[dit], dbl[dit]);
  }
  //since we only need ghost information in velocity space only
  // use simple exchange, instead of fillInternalGhosts
  m_fBJ_vel_ghost.exchange();
  
  // compute the divergence of individual velocity space fluxes
  // i.e., Psi quantities from Justin's thesis Chapter 4 (multiplied by J * nu)
  if (!m_Jpsi_linerad.isDefined()) m_Jpsi_linerad.define(dbl, 1, IntVect::Zero);
//   computeLineradVelFluxesDiv(a_rhs, m_fBJ_vel_ghost, phase_geom, mass);
  computeLineradVelFluxesDiv(m_Jpsi_linerad, m_fBJ_vel_ghost, phase_geom, mass);

  DataIterator cdit = a_rhs.dataIterator();
  for (cdit.begin(); cdit.ok(); ++cdit)
  {
    const FArrayBox& this_Jpsi = m_Jpsi_linerad[cdit];
    FORT_EVAL_LINERAD_RHS( CHF_BOX(a_rhs[cdit].box()),
                           CHF_FRA1(a_rhs[cdit],0),
                           CHF_CONST_FRA(this_Jpsi),
                           CHF_CONST_FRA1(a_neutral_density[cdit],0));
  }


  // check conservation properties
  if (m_first_call && m_linerad_diagnostics) {
    linerad_diagnostics(a_rhs, a_soln_species, a_time);
    exit(1);
  }
  
  m_first_call = false;
}

void FluidNeutrals::computeLineradVelFluxesDiv(LevelData<FArrayBox>&        a_Jpsi,
                                               LevelData<FArrayBox>&  a_fBJ,
                                               const PhaseGeom&             a_phase_geom,
                                               const double                 a_mass)
{
  /*This computes Jpsi divergence quantities.
   See Justin's thesis Chapter 4*/
  
  // get injected magnetic field
  const LevelData<FArrayBox>& inj_B = a_phase_geom.getBFieldMagnitude();
  
  // get injected cell-centered velocity spatial normalization
  const int use_spatial_vel_norm = a_phase_geom.spatialVelNorm();
  const LevelData<FArrayBox>& inj_vel_norm = a_phase_geom.getVelNorm();
      
  // get problem domain and number of vpar and mu cells
  const ProblemDomain& phase_domain = a_phase_geom.domain();
  const Box& domain_box = phase_domain.domainBox();
  int num_vp_cells = domain_box.size(VPARALLEL_DIR);
  int num_mu_cells = domain_box.size(MU_DIR);

  // create the temporary fluxes needed to compute the conservative
  // mean velocity and conservative temperature
  const DisjointBoxLayout& dbl = a_fBJ.getBoxes();
  if (!m_linerad_fluxes.isDefined())  m_linerad_fluxes.define(dbl, 1, IntVect::Zero);
  DataIterator dit = m_linerad_fluxes.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    // get phase space dx
    const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(dbl[dit]);
    const RealVect& phase_dx =  block_coord_sys.dx();
    
    const FArrayBox& nu_on_patch = m_linerad_cls_freq[dit];
    const FArrayBox& fBJnu_on_patch = a_fBJ[dit].mult(nu_on_patch);
    const FArrayBox& B_on_patch   = inj_B[dit];
        
    const FArrayBox* this_velnormptr;

   //  // Create a dummy to be used in Fortran if no velocity normalization
   //  FArrayBox dummy(B_on_patch.box(),2);
        
    for (int dir=0; dir<SpaceDim; dir++)
    {
      // we can use second_order = true to have consistency between
      // the PC and the collisional operator
       
      int is_second_order = (m_second_order_linerad_cls) ? 1 : 0;

      FORT_EVAL_LINERAD_FLUX(CHF_BOX(m_linerad_fluxes[dit][dir].box()),
                             CHF_FRA(m_linerad_fluxes[dit][dir]),
                             CHF_CONST_FRA1(fBJnu_on_patch,0),
                             CHF_CONST_FRA1(B_on_patch,0),
                             CHF_CONST_REAL(a_mass),
                             CHF_CONST_REALVECT(phase_dx),
                             CHF_CONST_INT(dir),
                             CHF_CONST_INT(num_vp_cells),
                             CHF_CONST_INT(num_mu_cells),
                             CHF_CONST_INT(is_second_order),
                             CHF_CONST_INT(use_spatial_vel_norm));

    }
  }

  // compute the divergence of the velocity space fluxes
  a_phase_geom.mappedGridDivergenceFromFluxNormals(a_Jpsi, m_linerad_fluxes);
  
}

void FluidNeutrals::computeIonizationRHS(LevelData<FArrayBox>& a_rhs,
                                         const LevelData<FArrayBox>& a_soln_dfn,
                                         const LevelData<FArrayBox>& a_neutral_density) const
{
   for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {

      FORT_COMPUTE_IONIZATION(CHF_BOX(a_rhs[dit].box()),
                              CHF_FRA1(a_rhs[dit],0),
                              CHF_CONST_FRA1(a_soln_dfn[dit],0),
                              CHF_CONST_FRA1(a_neutral_density[dit],0),
                              CHF_CONST_FRA1(m_ionization_sigmaV[dit],0));
   }
}

void FluidNeutrals::computeIonizationWithEnLossRHS(LevelData<FArrayBox>& a_rhs,
                                         const LevelData<FArrayBox>& a_source_dfn,
                                         const LevelData<FArrayBox>& a_soln_dfn,
                                         const LevelData<FArrayBox>& a_neutral_density) const
{
   for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {

      FORT_COMPUTE_IONIZATION_ENERGY_LOSS(CHF_BOX(a_rhs[dit].box()),
                              CHF_FRA1(a_rhs[dit],0),
                              CHF_CONST_FRA1(a_source_dfn[dit],0),
                              CHF_CONST_FRA1(a_soln_dfn[dit],0),
                              CHF_CONST_FRA1(a_neutral_density[dit],0),
                              CHF_CONST_FRA1(m_ionization_sigmaV[dit],0));
   }
}


void FluidNeutrals::computeChargeExchangeRHS(LevelData<FArrayBox>& a_rhs,
                                             const LevelData<FArrayBox>& a_soln_dfn,
                                             const LevelData<FArrayBox>& a_ion_density,
                                             const LevelData<FArrayBox>& a_ion_temperature,
                                             const LevelData<FArrayBox>& a_neutral_density,
                                             const LevelData<FArrayBox>& a_neutral_dfn) const
{

   for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {
      
      FORT_COMPUTE_CHARGE_EXCHANGE(CHF_BOX(a_rhs[dit].box()),
                                   CHF_FRA1(a_rhs[dit],0),
                                   CHF_CONST_FRA1(a_soln_dfn[dit],0),
                                   CHF_CONST_FRA1(a_ion_density[dit],0),
                                   CHF_CONST_FRA1(a_ion_temperature[dit],0),
                                   CHF_CONST_FRA1(a_neutral_density[dit],0),
                                   CHF_CONST_FRA1(a_neutral_dfn[dit],0));
      
      
      a_rhs[dit].mult(m_chx_norm);
   }
}


void FluidNeutrals::computeModelChargeExchangeRHS(LevelData<FArrayBox>&       a_rhs,
                                                  const LevelData<FArrayBox>& a_neutral_density,
                                                  const LevelData<FArrayBox>& a_ion_temperature,
                                                  const LevelData<FArrayBox>& a_par_vel_shift,
                                                  const LevelData<FArrayBox>& a_ion_unshifted_maxw,
                                                  const double                a_mass,
                                                  const PhaseGeom&            a_phase_geom) const
{
   
   const DisjointBoxLayout& grids = a_ion_unshifted_maxw.disjointBoxLayout();
   
   for (DataIterator dit( a_rhs.dataIterator() ); dit.ok(); ++dit) {
      
      // Get the physical velocity coordinates for this part of phase space
      const PhaseBlockCoordSys& block_coord_sys = a_phase_geom.getBlockCoordSys(grids[dit]);
      FArrayBox velocityRealCoords(a_rhs[dit].box(), VEL_DIM);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);
      
      FORT_COMPUTE_MODEL_CHARGE_EXCHANGE(CHF_BOX(a_rhs[dit].box()),
                                         CHF_FRA1(a_rhs[dit],0),
                                         CHF_CONST_FRA1(a_neutral_density[dit],0),
                                         CHF_CONST_FRA1(a_ion_temperature[dit],0),
                                         CHF_CONST_FRA1(a_par_vel_shift[dit],0),
                                         CHF_CONST_FRA1(a_ion_unshifted_maxw[dit],0),
                                         CHF_CONST_FRA1(velocityRealCoords,0),
                                         CHF_CONST_REAL(a_mass));
      
      
      a_rhs[dit].mult(m_chx_norm);
   }
}

void FluidNeutrals::computeIonizationSigmaV(CFG::LevelData<CFG::FArrayBox>&   a_ionization_sigmaV,
                                            const Real                        a_time) const
{
   
   //Universal constants (in CGS)
   double mp = 1.6726e-24;
   
   //Get normalization parameters (units)
   double N, T, L;
   ParmParse ppunits( "units" );
   ppunits.get("number_density",N);  //[m^{-3}]
   ppunits.get("temperature",T);     //[eV]
   ppunits.get("length",L);          //[m]
   
   double Tcgs = 1.602e-12 * T; //[erg]
   double Lcgs  = 1.0e2 * L;   //[cm]
   
   double time_norm = Lcgs / sqrt(Tcgs/mp); //[s]

   const CFG::DisjointBoxLayout& grids( a_ionization_sigmaV.getBoxes() );
   
   CFG::LevelData<CFG::FArrayBox> a2(grids, 1, CFG::IntVect::Zero);
   double fac = pow(T/10.0,2);
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      a2[dit].copy(m_Te_cfg[dit]);
      a2[dit].mult(m_Te_cfg[dit]);
      a2[dit].mult(fac);
   }
   
   //Get <sigmaV> in m^3/s
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      a_ionization_sigmaV[dit].copy(a2[dit]);
      a_ionization_sigmaV[dit].plus(3.0);
      a_ionization_sigmaV[dit].divide(a2[dit]);
      a_ionization_sigmaV[dit].invert(1.0);
      a_ionization_sigmaV[dit].mult(3.0e-14);
   }
   
   //Perform <sigmaV> normalization
   double norm = time_norm * N;
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      a_ionization_sigmaV[dit].mult(norm);
   }
}

void FluidNeutrals::updateIzSourceDfn(const KineticSpecies&  a_species)
{
   double T_norm;
   ParmParse ppunits( "units" );
   ppunits.get("temperature",T_norm);     //[eV]

   const PhaseGeom& phase_geom = a_species.phaseSpaceGeometry();
   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   const DisjointBoxLayout& grids = m_iz_source_maxw.disjointBoxLayout();
   
   CFG::LevelData<CFG::FArrayBox> iz_source_velocity_cfg( mag_geom.grids(), 1, CFG::IntVect::Zero );
   for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {
      iz_source_velocity_cfg[dit].setVal(0.);
   }
   (a_species).temperature(m_Te_cfg);
   // (a_species).numberDensity(m_electron_density_cfg);
   for (CFG::DataIterator dit(mag_geom.grids()); dit.ok(); ++dit) {

      double T_iz;
      CFG::BoxIterator bit(m_Te_cfg[dit].box());
      for (bit.begin(); bit.ok(); ++bit) {
         CFG::IntVect iv = bit();
         T_iz = 0.5 * m_Te_cfg[dit](iv,0) * T_norm - E_iz / 3.0;
         if (T_iz < T_iz_floor){
            T_iz = T_iz_floor;
         }
         // if (m_Te_cfg[dit](iv,0) * T_norm < (2/3) * E_iz){
         //    T_iz = 0.5 * m_Te_cfg[dit](iv,0) * T_norm;
         // }
         m_iz_source_temperature_cfg[dit](iv,0) = T_iz / T_norm;
      }

      // m_iz_source_temperature_cfg[dit].copy(T_iz);
      // m_iz_source_temperature_cfg[dit].divide(2.0);
      // m_iz_source_temperature_cfg[dit].plus(-(E_iz/3.0)/T_norm);

   }

   MaxwellianKernel<FArrayBox> maxwellian(m_iz_source_density_cfg, m_iz_source_temperature_cfg, iz_source_velocity_cfg);
   maxwellian.eval(m_iz_source_maxw,a_species);

   CFG::LevelData<CFG::FArrayBox> ne(mag_geom.grids(), 1, CFG::IntVect::Zero);
   a_species.numberDensity( ne );
   phase_geom.injectConfigurationToPhase(ne, m_ne_inj);
   
   for (DataIterator dit( m_iz_source_dfn.dataIterator() ); dit.ok(); ++dit) {
      
      FORT_MULT_NI(CHF_BOX(m_iz_source_dfn[dit].box()),
                   CHF_FRA1(m_iz_source_dfn[dit],0),
                   CHF_CONST_FRA1(m_iz_source_maxw[dit],0),
                   CHF_CONST_FRA1(m_ne_inj[dit],0));
   }

   // Compute density moment of this iz source dfn
   MomentOp& moment_op = MomentOp::instance();
   CFG::LevelData<CFG::FArrayBox> source_dfn_density_moment( mag_geom.grids(), 1, CFG::IntVect::Zero );
   moment_op.compute(source_dfn_density_moment, a_species, m_iz_source_dfn, DensityKernel<FArrayBox>());
   
   // Normalise this distribution to the density moment of the electron distribution to ensure particle conservation
   LevelData<FArrayBox> source_dfn_density_moment_inj;
   phase_geom.injectConfigurationToPhase(source_dfn_density_moment, source_dfn_density_moment_inj);
   for (DataIterator dit(grids.dataIterator() ); dit.ok(); ++dit) 
   {
      FORT_ENFORCE_INPUT_DENS_PROF(CHF_FRA(m_iz_source_dfn[dit]),
                                   CHF_BOX(grids[dit]),
                                   CHF_CONST_FRA1(source_dfn_density_moment_inj[dit],0),
                                   CHF_CONST_FRA1(m_ne_inj[dit],0));
   }

}

void FluidNeutrals::updateNeutralDfn(const KineticSpecies&  a_ion_species)
{
   const PhaseGeom& phase_geom = a_ion_species.phaseSpaceGeometry();
   const CFG::MagGeom& mag_geom( phase_geom.magGeom() );
   const CFG::DisjointBoxLayout& mag_grids = mag_geom.grids();
   
   CFG::LevelData<CFG::FArrayBox> ni(mag_grids, 1, CFG::IntVect::Zero);
   a_ion_species.numberDensity( ni );
   
   phase_geom.injectConfigurationToPhase(ni, m_ni_inj);
   
   for (DataIterator dit( m_neutral_dfn.dataIterator() ); dit.ok(); ++dit) {
      
      FORT_MULT_NI(CHF_BOX(m_neutral_dfn[dit].box()),
                   CHF_FRA1(m_neutral_dfn[dit],0),
                   CHF_CONST_FRA1(m_neutral_maxw[dit],0),
                   CHF_CONST_FRA1(m_ni_inj[dit],0));
   }

   // Compute density moment of this neutral dfn
   MomentOp& moment_op = MomentOp::instance();
   CFG::LevelData<CFG::FArrayBox> neutral_dfn_density_moment( mag_geom.grids(), 1, CFG::IntVect::Zero );
   moment_op.compute(neutral_dfn_density_moment, a_ion_species, m_neutral_dfn, DensityKernel<FArrayBox>());
   
   // Normalise this distribution to the density moment of the ion distribution to ensure particle conservation
   const DisjointBoxLayout& grids = m_neutral_dfn.disjointBoxLayout();
   LevelData<FArrayBox> neutral_dfn_density_moment_inj;
   phase_geom.injectConfigurationToPhase(neutral_dfn_density_moment, neutral_dfn_density_moment_inj);
   for (DataIterator dit(grids.dataIterator() ); dit.ok(); ++dit) 
   {
      FORT_ENFORCE_INPUT_DENS_PROF(CHF_FRA(m_neutral_dfn[dit]),
                                   CHF_BOX(grids[dit]),
                                   CHF_CONST_FRA1(neutral_dfn_density_moment_inj[dit],0),
                                   CHF_CONST_FRA1(m_ni_inj[dit],0));
   }
   
   
}

void FluidNeutrals::computeChxNormalization(double&       a_ioniz_norm,
                                            double&       a_chx_norm,
                                            double&       a_dens_norm,
                                            double&       a_vel_norm,
                                            double&       a_temp_norm,
                                            const bool    a_SI_input) const
{

   //If a_SI_input = true: expect that the physical parameters are provided as
   //neutral_density[1/m^3], <sigmaV>[m^3/s], neutral_Vpar[m/s], neutral_Temp[eV]
   //If a_SI_input = false: expect all input quantities are provided in COGENT units
   
   //Universal constants (in CGS)
   double mp = 1.6726e-24;

   //Get normalization parameters (units)
   double N, T, L;
   ParmParse ppunits( "units" );
   ppunits.get("number_density",N);  //[m^{-3}]
   ppunits.get("temperature",T);     //[eV]
   ppunits.get("length",L);          //[m]

   double Tcgs = 1.602e-12 * T; //[erg]
   double Lcgs  = 1.0e2 * L;   //[cm]
   
   double time_norm = Lcgs / sqrt(Tcgs/mp); //[s]

   if (a_SI_input) {
      a_ioniz_norm = time_norm * N;
      a_chx_norm = time_norm * 1.98e-14 * sqrt(T) * N;
      a_dens_norm = N;
      a_vel_norm = sqrt(Tcgs/mp);
      a_temp_norm = T;
   }
   
   else {
      a_ioniz_norm = 1.0;
      a_chx_norm = time_norm * 1.98e-14 * sqrt(T) * N;
      a_dens_norm = 1.0;
      a_vel_norm = 1.0;
      a_temp_norm = 1.0;

   }
   
}

inline
void FluidNeutrals::parseParameters( ParmParse& a_ppntr )
{

   if (a_ppntr.contains("include_charge_exchange")) {
      a_ppntr.get( "include_charge_exchange", m_include_chx );
   }
   else{
      m_include_chx = false;
   }

   if (a_ppntr.contains("charge_exchange_model_friction")) {
      a_ppntr.get( "charge_exchange_model_friction", m_chx_model_friction );
   }
   else{
      m_chx_model_friction = false;
   }

   if (a_ppntr.contains("SI_input")) {
      a_ppntr.get( "SI_input", m_SI_input );
   }
   else{
      m_SI_input = false;
   }

   if (a_ppntr.contains("use_Te_in_rate_calculation")) {
      a_ppntr.get( "use_Te_in_rate_calculation", m_use_Te_in_rate_calculation );
   }
   else 
   {
      m_use_Te_in_rate_calculation = false;
   }

   if (a_ppntr.contains("include_iz_energy_loss")) {
      a_ppntr.get( "include_iz_energy_loss", m_include_iz_energy_loss );
   }
   else 
   {
      m_include_iz_energy_loss = false;
   }

   if (a_ppntr.contains("include_line_radiation")) {
      a_ppntr.get( "include_line_radiation", m_include_line_radiation );
   }
   else 
   {
      m_include_line_radiation = false;
   }

   if (a_ppntr.contains("second_order_linerad_cls")) {
      a_ppntr.get( "second_order_linerad_cls", m_second_order_linerad_cls );
   }
   else 
   {
      m_second_order_linerad_cls = false;
   }

   if (a_ppntr.contains("line_radiation_diagnostics")) {
      a_ppntr.get( "line_radiation_diagnostics", m_linerad_diagnostics );
   }
   else 
   {
      m_linerad_diagnostics = false;
   }

   // KineticFunctionLibrary* library = KineticFunctionLibrary::getInstance();
   // std::string function_name;
  
   CFG::GridFunctionLibrary* grid_library = CFG::GridFunctionLibrary::getInstance();
   std::string grid_function_name;

   if (m_use_Te_in_rate_calculation == false)
   {
      if (a_ppntr.contains("electron_temperature")) {
         a_ppntr.get( "electron_temperature", grid_function_name );
         m_electron_temp = grid_library->find( grid_function_name );
      }
      else{
            MayDay::Warning("FluidNeutrals:: to compute ionization rate, either specify a prescribed electron temperature or use m_use_Te_in_rate_calculation = true");
      }
   }
   
   if (a_ppntr.contains("neutral_temperature")) {
      a_ppntr.get( "neutral_temperature", grid_function_name );
      m_neutr_temp = grid_library->find( grid_function_name );
   }
 
   if (m_include_chx) {
   
      if (a_ppntr.contains("parallel_velocity")) {
         a_ppntr.get( "parallel_velocity", grid_function_name );
         m_neutr_vel = grid_library->find( grid_function_name );
      }
      else{
         MayDay::Error("FluidNeutrals:: neutral parallel velocity must be specified ");
      }


      if (!a_ppntr.contains("neutral_temperature")) {
         MayDay::Warning("FluidNeutrals:: ion temperature is used for neutral temperature");
      }
   }
  
}


inline
void FluidNeutrals::printParameters()
{
   if (procID()==0) {
      std::cout << "FluidNeutrals neutral parameters:" << std::endl;
      std::cout << "  SI_input  =  " << m_SI_input;
      std::cout << "  include_charge_exchange  =  " << m_include_chx;

      std::cout << "  Neutral Density:" << std::endl;

      if (m_include_chx) {

         std::cout << "  Neutral parallel velocity:" << std::endl;
         m_neutr_vel->printParameters();
         
         if (m_neutr_temp != NULL) {
            std::cout << "  Neutral Temperature:" << std::endl;
            m_neutr_temp->printParameters();
         }
      }
   }
}

Real FluidNeutrals::computeDtExplicitTI(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  return TimeScale(a_soln, a_species);
}

Real FluidNeutrals::computeDtImExTI(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  return TimeScale(a_soln, a_species);
}

Real FluidNeutrals::TimeScale(const KineticSpeciesPtrVect& a_soln, const int a_species)
{
  //Simple calculation that assumes normalized neutral density the order of unity,
  //and that charge-exchange process is the stiffest process  


  if (m_first_call_tscale) {
    double dens_norm, vel_norm, temp_norm;
    computeChxNormalization(m_ionization_norm, m_chx_norm, dens_norm, vel_norm, temp_norm, m_SI_input);
  }

  Real time_scale;
  if (m_include_chx) {
    time_scale = 1.0/m_chx_norm;
  }
  else {
    time_scale = 1.0/m_ionization_norm;
  }

  m_first_call_tscale = false;

  return time_scale;
}

void FluidNeutrals::diagnostics(const LevelData<FArrayBox>& a_rhs,
                                const KineticSpecies&       a_rhs_species,
                                const double                a_time) const
{

  //Get geometry                                                                                                                                                            
  const PhaseGeom& phase_geom = a_rhs_species.phaseSpaceGeometry();
  const CFG::MultiBlockLevelGeom& mag_geom( phase_geom.magGeom() );

  //Plot particle source
  CFG::LevelData<CFG::FArrayBox> particle_src( mag_geom.grids(), 1, CFG::IntVect::Zero );                                                                                 
  a_rhs_species.numberDensity( particle_src );                                                                                                                              
  phase_geom.plotConfigurationData( "particle_src", particle_src, a_time );                                                                                                 
  //Plot parallel particle flux (nVpar) source
  CFG::LevelData<CFG::FArrayBox> parMom_src( mag_geom.grids(), 1, CFG::IntVect::Zero );
  a_rhs_species.parallelParticleFlux( parMom_src );
  phase_geom.plotConfigurationData( "parMom_src", parMom_src, a_time );

}

void FluidNeutrals::linerad_diagnostics(const LevelData<FArrayBox>& a_rhs,
                               const KineticSpecies&       a_rhs_species,
                               const double                a_time) const
{
  //Get geometry
  const PhaseGeom& phase_geom = a_rhs_species.phaseSpaceGeometry();
  const CFG::MagGeom& mag_geom( phase_geom.magGeom() );

  //Get moment operator
  MomentOp& moment_op = MomentOp::instance();
  
  //Plot particle source
  CFG::LevelData<CFG::FArrayBox> particle_src( mag_geom.grids(), 1, CFG::IntVect::Zero );
  moment_op.compute(particle_src, a_rhs_species, a_rhs, DensityKernel<FArrayBox>());
  phase_geom.plotConfigurationData( "plt_linerad_particle_src_plots/linerad_particle_src0000.", particle_src, a_time );
  
  //Plot parallel momentum source
  CFG::LevelData<CFG::FArrayBox> parMom_src( mag_geom.grids(), 1, CFG::IntVect::Zero );
  moment_op.compute(parMom_src, a_rhs_species, a_rhs, ParallelVelKernel<FArrayBox>());
  phase_geom.plotConfigurationData( "plt_linerad_parmom_src_plots/linerad_parmom_src0000.", parMom_src, a_time );
  
  //Plot energy source
  CFG::LevelData<CFG::FArrayBox> energy_src( mag_geom.grids(), 1, CFG::IntVect::Zero );
  moment_op.compute(energy_src, a_rhs_species, a_rhs, KineticEnergyKernel<FArrayBox>());
  phase_geom.plotConfigurationData( "plt_linerad_energy_src_plots/linerad_energy_src0000.", energy_src, a_time );

  //Plot collision frequency
  GKDiagnostics m_diagnostics_fg;
  m_diagnostics_fg.plotPhaseVar(m_linerad_cls_freq, phase_geom, "plt_linerad_cls_freq_plots/linerad_cls_freq0000.", a_time);
}

#include "NamespaceFooter.H"

#ifndef  _FLUIDNEUTRALS_H_
#define  _FLUIDNEUTRALS_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "KineticFunction.H"
#include "REAL.H"
#include "NTRInterface.H"
#include "ParmParse.H"
#include <sstream>

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "GridFunction.H"
#include "GridFunctionLibrary.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM



#include "NamespaceHeader.H"

/**
 * FluidNeutrals neutrals operator class.
*/
class FluidNeutrals
   : public NTRInterface
{
public:

  /// Constructor with initialization.
  /**
   *  @param[in] pp the input database.
   */
  FluidNeutrals( ParmParse&  a_ppntr,
        const int a_verbosity );

  /// Destructor.
  /**
   */
  virtual ~FluidNeutrals();

  /// Accumulates the RHS of the collision operator.
  /**
   * @param[out] rhs data holder for rhs (without ghost cells).
   * @param[in] soln current solution with which to evaluate RHS.
   * @param[in] species index of species to which collision is applied.
   * @param[in] time current time.
   */
  virtual void evalNtrRHS( KineticSpecies&                  rhs_species,
                           const KineticSpeciesPtrVect&     soln,
                           const CFG::FluidSpeciesPtrVect&  fluid_species_phys,
                           const int                        species,
                           const Real                       time );

  /// Compute a stable time step.
  /**
   * Computes and returns an estimate of the maximum stable time step.
   */
  Real computeDtExplicitTI(const KineticSpeciesPtrVect& a_soln, const int a_species);
  Real computeDtImExTI(const KineticSpeciesPtrVect& a_soln, const int a_species);
  Real TimeScale(const KineticSpeciesPtrVect& a_soln, const int a_species);
   
private:

   int m_verbosity;

   bool m_include_chx, m_chx_model_friction, m_first_call, m_use_Te_in_rate_calculation, m_include_iz_energy_loss, m_include_line_radiation;
   bool m_SI_input;
   int electron_species_index;

   double m_ionization_norm;
   double m_chx_norm;

   RefCountedPtr<KineticFunction>   m_neutr_func;   
   RefCountedPtr<CFG::GridFunction> m_neutr_vel;
   RefCountedPtr<CFG::GridFunction> m_neutr_temp;
   RefCountedPtr<CFG::GridFunction> m_electron_temp;

   CFG::LevelData<CFG::FArrayBox> m_neutral_density_cfg;
   CFG::LevelData<CFG::FArrayBox> m_neutral_velocity_cfg;
   CFG::LevelData<CFG::FArrayBox> m_iz_source_density_cfg;
   CFG::LevelData<CFG::FArrayBox> m_iz_source_temperature_cfg;
   CFG::LevelData<CFG::FArrayBox> m_Te_cfg;
   
   
   LevelData<FArrayBox> m_ionization_sigmaV;
   LevelData<FArrayBox> m_neutral_dfn;
   LevelData<FArrayBox> m_iz_source_maxw;
   LevelData<FArrayBox> m_iz_source_dfn;
   LevelData<FArrayBox> m_neutral_maxw;
   LevelData<FArrayBox> m_ni_inj;
   LevelData<FArrayBox> m_ne_inj;
   LevelData<FArrayBox> m_linerad_cls_freq;
   LevelData<FArrayBox> m_linerad_cls_freq_vel_ghost;
   LevelData<FArrayBox> m_fBJ_vel_ghost;
   LevelData<FluxBox> m_linerad_fluxes;
   LevelData<FArrayBox> m_Jpsi_linerad;
   bool m_second_order_linerad_cls;
   bool m_linerad_diagnostics;
   double E_iz = 13.6;
   double T_iz_floor = 3.0;

   bool m_first_call_tscale;

   // Define the fit parameters for line radiation (TODO: Get these from the input file instead)
  // H0
   double Abar = 5.5949129e-32;
   double alpha = 8.0000102e3;
   double beta = 7.9587517e-1;
   double gamma = -1.3919964e0;
   double V0 = 3.5201735e0;
  // Ar8+
  // double Abar = 5.8864840e-32;
  // double alpha = 1.0739759e9;
  // double beta = 1.4272764e0;
  // double gamma = -1.5696360e0;
  // double V0 = 1.7127784e1;
  // // Li1+
  // double Abar = 1.8391597e-32;
  // double alpha = 1.5550315e4;
  // double beta = 4.8963899e-1;
  // double gamma = -1.4801717e0;
  // double V0 = 7.7295814e0;
   
   /// Computes the ionization contribution to RHS.
   /*
      rhs = <sigmaV_ionization> * neutr_dens * fB
    */
   void computeIonizationRHS( LevelData<FArrayBox>& rhs,
                              const LevelData<FArrayBox>& soln_dfn,
                              const LevelData<FArrayBox>& neutral_density) const;

    void computeIonizationWithEnLossRHS( LevelData<FArrayBox>& rhs,
                              const LevelData<FArrayBox>& source_dfn,
                              const LevelData<FArrayBox>& soln_dfn,
                              const LevelData<FArrayBox>& neutral_density) const;
   
    /// Computes the line radiation as per J. Roeltgen et al., Nucl. Fus. 65 (2025)                    
    void computeLineradClsFreq(const KineticSpecies&     soln_species,
                        const PhaseGeom&          phase_geom,
                        const DisjointBoxLayout&  grids,
                        LevelData<FArrayBox>&     linerad_cls_freq);

    void computeLineradVelFluxesDiv(LevelData<FArrayBox>&        Jpsi,
                           LevelData<FArrayBox>&  fBJ,
                           const PhaseGeom&             phase_geom,
                           const double                 mass);   

    virtual void evalLineradClsRHS( LevelData<FArrayBox>&  rhs,
                                    KineticSpecies&        soln_species,
                                    const LevelData<FArrayBox>& neutral_density,
                                    const Real             time);

   /// Computes the charge-exchange contribution to RHS.
   /*
    rhs = <sigmaV_chx> (n_ion * fB_neutral - n_neutral * fB)
    */

   void computeChargeExchangeRHS(LevelData<FArrayBox>& rhs,
                                 const LevelData<FArrayBox>& soln_dfn,
                                 const LevelData<FArrayBox>& ion_density,
                                 const LevelData<FArrayBox>& ion_temperature,
                                 const LevelData<FArrayBox>& neutral_density,
                                 const LevelData<FArrayBox>& neutral_dfn) const;

   /// Computes the model friction charge-exchange contribution to RHS.
   /*
    rhs = <sigmaV_chx> * n_neutral * vpar * upar * mass_ion * fB_unshifted_Maxw / Ti)
    where 
      upar = neutral_par_vel - ion_par_vel
      fB_unshifted_Maxw is an unshifted Maxwellian with ion density and temperature
    */
   void computeModelChargeExchangeRHS(LevelData<FArrayBox>&       rhs,
                                      const LevelData<FArrayBox>& neutral_density,
                                      const LevelData<FArrayBox>& ion_temperature,
                                      const LevelData<FArrayBox>& par_vel_shift,
                                      const LevelData<FArrayBox>& ion_unshifted_maxw,
                                      const double                mass,
                                      const PhaseGeom&            phase_geom) const;

   
   
   /// Computes normalizations to handle the case of SI_input data
   void computeChxNormalization(double&       ioniz_norm,
                                double&       chx_norm,
                                double&       dens_norm,
                                double&       vel_norm,
                                double&       temp_norm,
                                const bool   SI_input) const;

  
   /// Computes ionization reaaction coefficient sigmaV
   void computeIonizationSigmaV(CFG::LevelData<CFG::FArrayBox>&   ionization_sigmaV,
                                const Real                        time) const;
   
   
   void updateNeutralDfn(const KineticSpecies&  a_ion_species);
   
   void updateIzSourceDfn(const KineticSpecies&  a_ion_species);
   
   
   /// Parse parameters.
   /**
    * Private method to obtain control parameters from "NTR.species" section
    * of ParmParse database, where "species" is replaced by the species name
    *
    * \param[in] ppntr NTR.species input database.
    */
   void parseParameters( ParmParse& pp );
   
   /// Print parameters.
   /**
    * Private method to print parameters.
    */
   void printParameters();
   
   //Plots particle, parallel momentum and energy sources
   void diagnostics(const LevelData<FArrayBox>& a_rhs,
		    const KineticSpecies&       a_rhs_species,
		    const double                a_time) const;
    void linerad_diagnostics(const LevelData<FArrayBox>& a_rhs,
		    const KineticSpecies&       a_rhs_species,
		    const double                a_time) const;


};


#include "NamespaceFooter.H"

#endif

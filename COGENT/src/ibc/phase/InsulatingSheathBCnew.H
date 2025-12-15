#ifndef _InsulatingSheathBCnew_H_
#define _InsulatingSheathBCnew_H_

#include <map>
#include "BoundaryBoxLayout.H.multidim"
#include "KineticSpeciesBC.H"

#ifdef CH_MULTIDIM

#include "NamespaceHeader.H"


/// Logical BC class
/**
 * Under "insulating_sheath_new", there are three options: "lsbc", "absorption", and "reflection".
 * "absorption" implies that there is zero ingoing flux on the boundary for every incoming velocity.
 * In this case, the ghost cells are filled with some values (not necessary zeros) to ensure that PDF on the boudary
 * faces is zero.
 * "reflection" is the opposite of "absorption", so that for every face value of the outgoing PDF there is a reflected
 * face value of the ingoing PDF: f_face[-v] = f_face[+v]. The ghost cells are filled to ensure this condition.
 * "lsbc" stands for the Logical Sheath Boudary Condition, where absorption BC is applied for the ions, and the electron
 * total boundary flux is equal to ion total boundary flux.
 *
 * Since BC types for every species are set individually in the input file, there might be some poorly bahaving scenarios.
 * Here, we lisred all possible combinations and the resulting behavior, with L, A, and R standing for lsbc, absorption, and reflection,
 * and the first letter standing for ions, the second for electrons.
 *
 * >> LL - Ion BC does nothing and relies on electron lsbc to be called. Ion absorption BC is implemented there.
 *         Electron lsbc is then executed.
 * >> LA - Ion BC does nothing and relies on electron lsbc to be called. Ion BC is NEVER called, causing unexpected behavior.
 *         Electron absorption BC is then executed.
 * >> LR - Ion BC does nothing and relies on electron lsbc to be called. Ion BC is NEVER called, causing unexpected behavior.
 *         Electron reflection BC is then executed.
 * >> AL - Ion absorption BC is called. Electron lsbc is called and performs ion absorption BC once again wasting resoures.
 * >> AA - Ion absorption BC is called. Electron absorption BC is called.
 * >> AR - Ion absorption BC is called. Electron reflection BC is called.
 * >> RL - Ion reflection BC is called. Electron lsbc is called and overwrites ion reflection BC to absorption one.
 * >> RA - Ion reflection BC is called. Electron absorption BC is called.
 * >> RR - Ion reflection BC is called. Electron reflection BC is called.
 *
 * As follows from the description above, LA, LR, AL, and RL combinations should be avoided in the input file.
 */
class InsulatingSheathBC_new
{
  public:

   /// Constructor
   /** Constructor
    * @param[in] coords  the boundary box layout
    */
   InsulatingSheathBC_new( const BoundaryBoxLayoutPtr& bdry_layout,
                    const ParmParse&            pp,
                    const int id,
                    const Real& a_time);

   /// Destructor
   /** Destructor
    */
   ~InsulatingSheathBC_new() {delete m_dv_par_e; delete  m_v_par_max_e; delete[] m_phi_val;}


   
   Real computeSheathBC(LevelData<FArrayBox>& phi_bc,
                        const KineticSpeciesPtrVect& species,
                        const int& species_index) const;
   
   void fillBoundaryData(LevelData<FArrayBox>& bdry_data,
                         const LevelData<FArrayBox>& dfn) const;
   void fillBoundaryData_VG(LevelData<FArrayBox>& a_bdry_data,
                                  const LevelData<FArrayBox>& a_dfn) const;
   
   void computeBoundaryIonCurrent(LevelData<FArrayBox>& bdry_ion_current,
                                  const KineticSpeciesPtrVect& species) const;
   void computeBoundaryElectronCurrent(LevelData<FArrayBox>& bdry_ion_current,
                                  const KineticSpeciesPtrVect& species,
                                  Real a_phi) const;                                  
   
   void applyBC(KineticSpeciesPtrVect& species,
                const int& species_index,
                const LevelData<FluxBox>& velocity,
                const CFG::LevelData<CFG::FArrayBox>& phi );
                
  void applyBC_iter(Real* p_res,
                KineticSpeciesPtrVect& species,
                const int& species_index,
                const LevelData<FluxBox>& velocity,
                const CFG::LevelData<CFG::FArrayBox>& phi );
                
  Real computeAdvection(const KineticSpeciesPtrVect& a_species,
                                         const int& a_species_index,
                                         const LevelData<FluxBox>& a_velocity) const;             

  void computeHeatFluxes( Real * p_res,
                          const KineticSpeciesPtrVect& a_species,
                          const int& a_species_index,
                          const LevelData<FluxBox>& a_velocity) const;
                          
   /// Sum the src array over velocity space
   /** Sum the src array over velocity space
    *
    * @param[in]  src  Data to be summed over velocity space
    * @param[out] dst  velocity-space sums
    */
   void sum( const LevelData<FArrayBox>& src,
             LevelData<FArrayBox>&       dst ) const;
   Real simple_sum( const LevelData<FArrayBox>& a_src, bool a_flag=false) const;
   
   static const string pp_name;
   static int print_ghosts;
   
private:

   void parseParameters(const ParmParse& pp);
   std::string m_advect_scheme;
   DisjointBoxLayout m_grids_inj;
   DisjointBoxLayout m_grids_full;
   const BoundaryBoxLayoutPtr m_bdry_layout;
   bool m_compute_potential_BC;
   std::string m_sheath_bc_type;
   Real* m_dv_par_e;
   Real* m_v_par_max_e;
   Real* m_phi_val;
   int m_iter_number;
   Real m_time;
   int m_id;
   Real m_tolerance;
   static map<int, Real> m_id_map;
   bool m_debug;
   bool m_newton;
   LevelData<FArrayBox> m_phi_injected;

   
};

#include "NamespaceFooter.H"

#endif /* CH_MULTIDIM */

#endif
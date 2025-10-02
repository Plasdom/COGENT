#include "InsulatingSheathBCnew.H"
#include "Directions.H"
#include "FlipGrids.H"
#include "LogicalSheathBCF_F.H"
#include <fstream>
#include <cstring>
#include <map>
#include <utility>
//#include "SpaceUtils.H"
//#include "SpaceUtilsF_F.H"

//#include "altFaceAverages.H"
//#include "altFaceAveragesF_F.H"

#include "SpaceUtils.H.multidim"
#include "SpaceUtilsF_F.H"


#undef CH_SPACEDIM
#include "ReductionOps.H.multidim"
#include "ReductionCopier.H.multidim"
#include "SpreadingCopier.H.multidim"
#ifdef CH_SPACEDIM
#undef CH_SPACEDIM
#endif

#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

const string InsulatingSheathBC_new::pp_name = "insulating_sheath_bc";
map<int, Real> InsulatingSheathBC_new::m_id_map = {}; 
int InsulatingSheathBC_new::print_ghosts = 0; // this counter will trigger whn to prit ghost cells

//------------------------------------------------------------------------------

double computeFaceValueWENO(double* a_ptr, double a_epsilon=10e-8)
{
  /** This method returns the cell-centered value on a phase between cell 2 and 3 using WENO scheme.
    * In particular, a_ptr has a length of 5, where a_ptr[0] = u[-2], a_ptr[1] = u[-1],... a_ptr[4] = u[2],
    * where increasing indexing shows advection direction.
    */
  //double alpha1, alpha2, alpha3;
  double beta1, beta2, beta3;
  double d1 = 0.1, d2 = 0.6, d3 = 0.3;
  double face1, face2, face3;
  face1 = (2*a_ptr[0] - 7*a_ptr[1] + 11*a_ptr[2]) / 6.0;
  face2 = (-1*a_ptr[1] + 5*a_ptr[2] + 2*a_ptr[3]) / 6.0;
  face3 = (2*a_ptr[2] + 5*a_ptr[3] - 1*a_ptr[4]) / 6.0;


  beta1 = 13.0/12.0 * (a_ptr[0] - 2*a_ptr[1] + a_ptr[2])*(a_ptr[0] - 2*a_ptr[1] + a_ptr[2]) + (a_ptr[0] - 4*a_ptr[1] + 3*a_ptr[2])*(a_ptr[0] - 4*a_ptr[1] + 3*a_ptr[2]) / 4;
  beta2 = 13.0/12.0 * (a_ptr[1] - 2*a_ptr[2] + a_ptr[3])*(a_ptr[1] - 2*a_ptr[2] + a_ptr[3]) + (a_ptr[1] - a_ptr[3])*(a_ptr[1] - a_ptr[3]) / 4;
  beta3 = 13.0/12.0 * (a_ptr[2] - 2*a_ptr[3] + a_ptr[4])*(a_ptr[2] - 2*a_ptr[3] + a_ptr[4]) + (3*a_ptr[2] - 4*a_ptr[3] + a_ptr[4])*(3*a_ptr[2] - 4*a_ptr[3] + a_ptr[4]) / 4;
//std::cout<<"beta: "<<beta1<<"  "<<beta2<<"  "<<beta3<<"\n";
  // Use the same beta to store alpha

  beta1 = d1/(a_epsilon+beta1)/(a_epsilon+beta1);
  beta2 = d2/(a_epsilon+beta2)/(a_epsilon+beta2);
  beta3 = d3/(a_epsilon+beta3)/(a_epsilon+beta3);

  double a_s = beta1 + beta2 + beta3;

  beta1 = beta1 / a_s;
  beta2 = beta2 / a_s;
  beta3 = beta3 / a_s;
  double face = beta1*face1 + beta2*face2 + beta3*face3;

  return face;
}

//------------------------------------------------------------------------------

InsulatingSheathBC_new::InsulatingSheathBC_new( const BoundaryBoxLayoutPtr&  a_bdry_layout,
                                  const ParmParse&             a_pp,
                                  const int id,
                                  const Real& a_time)

: m_bdry_layout(a_bdry_layout),
  m_compute_potential_BC(false),
  m_sheath_bc_type("second_order"),
  m_debug(false),
  m_tolerance(1.0e-12),
  m_advect_scheme("uw1"),
  m_iter_number(10),
  m_id(id),
  m_time(a_time),
  m_newton(false)
  // m_sheath_bc_type("first_order")
{
   // Create boundary dbl, which is one-cell-wide in vpar and mu
   const DisjointBoxLayout& grids_full( a_bdry_layout->disjointBoxLayout() );
   m_grids_full.define( a_bdry_layout->disjointBoxLayout() );
   DisjointBoxLayout grids_tmp;
   adjCellLo(grids_tmp, grids_full, VPARALLEL_DIR, -1);
   adjCellLo(m_grids_inj, grids_tmp, MU_DIR, -1);
    //grids_full.physDomain(); returns const ref to problem domain
   // VG add
   
    m_dv_par_e = new Real[1];
    m_v_par_max_e = new Real[1];
    m_phi_val = new Real[1];
    m_phi_val[0] = 1.7;

  // Open file for the potential history and clean it
  if(procID()==0 && !m_id_map.count(m_id)){
    string filename = "phi_history.txt";
    ofstream fout(filename);
    fout.close();
  }
  const int& dir( m_bdry_layout->dir() );
  const Side::LoHiSide& side( m_bdry_layout->side() );
  DisjointBoxLayout tmp_layout;
  Real tmp=0;
  
  if (side == Side::Lo) {adjCellHi(tmp_layout, m_grids_full, dir, 1);tmp=1;}
  if (side == Side::Hi) {adjCellLo(tmp_layout, m_grids_full, dir, 1);tmp=2;}
  Real tmp1;
  tmp1=tmp;
   parseParameters( a_pp );
}

Real
InsulatingSheathBC_new::computeSheathBC(LevelData<FArrayBox>& a_phi_bc,
                                 const KineticSpeciesPtrVect& a_species,
                                 const int& a_species_index) const
{

////check if both sides work
////const Side::LoHiSide& side( m_bdry_layout->side() );
////if(sign(side)==-1) {return 0.0;}
   //const Side::LoHiSide& side( m_bdry_layout->side() );
   //if(sign(side)==1) {MayDay::Error( "DONE" );}
   
   // Get coordinate system parameters
   //const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );
   //const MultiBlockCoordSys& coord_sys( *(geometry.coordSysPtr()) );
   //const PhaseGrid& phase_grid = geometry.phaseGrid();
   
   // Get dfn and create bdry_data container
   const KineticSpecies& species_physical( *(a_species[a_species_index]) );
   const DisjointBoxLayout& bdry_dbl( m_bdry_layout->disjointBoxLayout() );
   ////LevelData<FArrayBox> bdry_data(bdry_dbl, 1, IntVect::Zero);
   LevelData<FArrayBox> bdry_data(m_grids_inj, 1, IntVect::Zero);
   const LevelData<FArrayBox>& dfn = species_physical.distributionFunction();


   
   //Get ion parallel current or other data on the boundary,
   // UNCOMMENT THIS CODE ONCE IT'S TESTED
   
    ////LevelData<FArrayBox> bdry_ion_current(bdry_dbl, 1, IntVect::Zero);
    LevelData<FArrayBox> bdry_ion_current(m_grids_inj, 1, IntVect::Zero);
    computeBoundaryIonCurrent(bdry_ion_current, a_species);
    ////LevelData<FArrayBox> residual_current(bdry_dbl, 1, IntVect::Zero);
    LevelData<FArrayBox> residual_current(m_grids_inj, 1, IntVect::Zero);
   ////LevelData<FArrayBox> bdry_ele_current(bdry_dbl, 1, IntVect::Zero);
   LevelData<FArrayBox> bdry_ele_current(m_grids_inj, 1, IntVect::Zero);
   
    //computeBoundaryElectronCurrent(bdry_ele_current, a_species, 0.0);
  
  /* double x = 0;
   DataIterator dit = m_grids_inj.dataIterator();
  //for (DataIterator dit( m_grids_inj ); dit.ok(); ++dit) {
  for (dit.begin(); dit.ok(); ++dit) {
   const Box& this_box = bdry_ele_current[dit].box();
       for (BoxIterator bit(this_box); bit.ok(); ++bit) {
           IntVect iv = bit();
           x = 5+6*iv[0];
       }
   }

   cout<< procID()<<"   "<<x<<endl;
   

    DataIterator dit1 = m_grids_full.dataIterator();
    for (dit1.begin(); dit1.ok(); ++dit1) {
   const Box& this_box = bdry_ele_current[dit1].box();
       for (BoxIterator bit(this_box); bit.ok(); ++bit) {
           IntVect iv = bit();
           x = 5+6*iv[0];
       }
   }
    
       MayDay::Error( "DONE" ); 
   //MayDay::Error( "DONE\n" );
    */
    
   Real phi = 0.0;
   Real phi_hi = 100000;
   
   for (int s_index(0); s_index<a_species.size(); s_index++) {
     // Works for 2 species only
     const KineticSpecies& species_physical( *(a_species[s_index]) );
     const double charge = species_physical.charge();
     if ( charge > 0.0 ) {continue;} 
     
     const PhaseGeom& geometry( species_physical.phaseSpaceGeometry() );
     // Get dv_par
     const VEL::VelCoordSys& vel_coords = geometry.velSpaceCoordSys();
     const VEL::RealVect& vel_dx = vel_coords.dx();
     const Real dv_par = vel_dx[0];
     m_dv_par_e[0] = dv_par;
    
     const Box& domain_box = (geometry.domain()).domainBox();
     IntVect lo_end(domain_box.smallEnd());
     IntVect hi_end(domain_box.bigEnd());
     int v_min = lo_end[VPARALLEL_DIR];
     int v_max = hi_end[VPARALLEL_DIR];
     int N_v = (v_max - v_min + 1)/2;
     Real* arr_v = new Real[N_v];
     double vel_max = dv_par*N_v;
     m_v_par_max_e[0] = vel_max;
     
     const double mass = species_physical.mass();
     phi_hi = fabs(mass*vel_max*vel_max/2/charge);
   }
   
   
   
      
   Real phi_lo = 0.0;
   Real residual = DBL_MAX;
   Real res_der;
   //LevelData<FArrayBox> bdry_ele_current(bdry_dbl, 1, IntVect::Zero);
   //while (residual > 1.0e-6) {

/*           for(int i=0;i<5;i++){
             phi = m_phi_val[0]*1.1;
             Real local_residual = 0.0;
             computeBoundaryElectronCurrent(bdry_ele_current, a_species, phi);
             for (DataIterator dit( m_grids_inj ); dit.ok(); ++dit) {
               Real local_res_tmp = 0.0;
               residual_current[dit].setVal(0.0);
               residual_current[dit].plus(bdry_ion_current[dit]);
               residual_current[dit].plus(bdry_ele_current[dit]);
               const Box& this_box = residual_current[dit].box();
               for (BoxIterator bit(this_box); bit.ok(); ++bit) {
                 IntVect iv = bit();
                 local_res_tmp += residual_current[dit](iv, 0);
               }
               local_residual += local_res_tmp;
             }
#ifdef CH_MPI
             MPI_Allreduce(&local_residual, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else         
             residual = local_residual;
#endif
             res_der = residual;
             
             phi = m_phi_val[0];
             local_residual = 0.0;
             computeBoundaryElectronCurrent(bdry_ele_current, a_species, phi);
             for (DataIterator dit( m_grids_inj ); dit.ok(); ++dit) {
               Real local_res_tmp = 0.0;
               residual_current[dit].setVal(0.0);
               residual_current[dit].plus(bdry_ion_current[dit]);
               residual_current[dit].plus(bdry_ele_current[dit]);
               const Box& this_box = residual_current[dit].box();
               for (BoxIterator bit(this_box); bit.ok(); ++bit) {
                 IntVect iv = bit();
                 local_res_tmp += residual_current[dit](iv, 0);
               }
               local_residual += local_res_tmp;
             }
#ifdef CH_MPI
             MPI_Allreduce(&local_residual, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else         
             residual = local_residual;
#endif
             res_der = (res_der - residual)/m_phi_val[0]/0.1;
             
             phi -= residual/res_der;
             
             if (procID()==0) {
               cout<<"phi =  "<<phi<<endl;
               cout<<"residual =  "<<residual<<endl;
             }
             m_phi_val[0] = phi;
}    
            if (procID()==0) {
               cout<<endl;
             }
*/

   
   // Try only 20 iterations  for now
   for(int i=0;i<160;i++){
     Real local_residual = 0.0;
     computeBoundaryElectronCurrent(bdry_ele_current, a_species, phi);
     //for (DataIterator dit( bdry_dbl ); dit.ok(); ++dit) {
     for (DataIterator dit( m_grids_inj ); dit.ok(); ++dit) {
       Real local_res_tmp = 0.0;
       residual_current[dit].setVal(0.0);
       residual_current[dit].plus(bdry_ion_current[dit]);
       residual_current[dit].plus(bdry_ele_current[dit]);
       const Box& this_box = residual_current[dit].box();
       for (BoxIterator bit(this_box); bit.ok(); ++bit) {
           IntVect iv = bit();
           local_res_tmp += residual_current[dit](iv, 0);
       }
       local_residual += local_res_tmp;
       
     }
////cout<<procID()<<"  local_residual  "<< local_residual <<endl;  
     
#ifdef CH_MPI
      MPI_Allreduce(&local_residual, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      residual = local_residual;
#endif
    
    if(residual<0){
      phi_lo = phi;
      phi = 0.5*(phi + phi_hi);
    }
    else{
     phi_hi = phi;
     phi = 0.5*(phi_lo + phi);
    }
            if (procID()==0) {
           cout<<"phi =  "<<phi<<endl;
           cout<<"residual =  "<<residual<<endl;
           }
   }

   
   
////MayDay::Error( "DONE\n" ); 
   /*
   for (DataIterator dit( bdry_dbl ); dit.ok(); ++dit) {
     a_phi_bc[dit].setVal(phi);
   }*/
   return phi;
/*
   // Solve for potential BC
   Real residual = DBL_MAX;
   while (residual > 1.0e-6) {
      
      // Fill boundary data object with dfn
      fillBoundaryData(bdry_data, dfn);

      for (DataIterator dit(bdry_dbl); dit.ok(); ++dit) {
         //const Box& bdry_box = bdry_dbl[dit];
         //const FArrayBox& this_data = bdry_data[dit];
         //FArrayBox& this_phi = a_phi_bc[dit];
         
         //   Here, do something to the phase-space bdry_data using
          //  the current iteration of a_phi_bc. Note that this_phi is defined
         //   with m_grids_inj, which are flattened in vpar and mu.
         
      }
      
      // Integrate over velocity space
      LevelData<FArrayBox> bdry_data_summed(m_grids_inj, 1, IntVect::Zero);
      sum(bdry_data, bdry_data_summed);

      Real local_residual = DBL_MAX;
      for (DataIterator dit(m_grids_inj); dit.ok(); ++dit) {
         //const FArrayBox& this_data = bdry_data_summed[dit];
         //FArrayBox& this_phi = a_phi_bc[dit];
         FArrayBox this_residual(m_grids_inj[dit], 1);
         
         //   Here, compute residual and update a_phi_bc
         
         local_residual = this_residual.norm();
      }
      
#ifdef CH_MPI
      MPI_Allreduce(&local_residual, &residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
      residual = local_residual;
#endif
   }
   */
}

//----------------------------------------------------------------------------------------------------------

void
InsulatingSheathBC_new::fillBoundaryData_VG(LevelData<FArrayBox>& a_bdry_data,
                                  const LevelData<FArrayBox>& a_dfn) const
{
   /*
    Copy distribution function from the last valid cell
    into the bdry_data object that occupies one-cell-wide boundary layer.
    */
   const DisjointBoxLayout& bdry_grids( a_bdry_data.disjointBoxLayout() );
   for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

      const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
      const FArrayBox& this_dfn( a_dfn[interior_dit] );
      //FArrayBox& this_data( a_bdry_data[dit] );
      
      //FArrayBox tmp(this_dfn.box(), 1);
      //tmp.copy(this_dfn);
      const Box& tmp_box = a_bdry_data[dit].box();
      const Side::LoHiSide& side( m_bdry_layout->side() );
      if (sign(side)==1) {MayDay::Error( "SIDE 1\n" );}
      for (BoxIterator bit(tmp_box); bit.ok(); ++bit) {
        IntVect iv = bit();
        IntVect iv_sh = iv;
        iv_sh[1] = iv_sh[1] - sign(side);
        a_bdry_data[dit](iv,0) = this_dfn(iv_sh,0);
      }

   }
}

//----------------------------------------------------------------------------------------------------------

void
InsulatingSheathBC_new::fillBoundaryData(LevelData<FArrayBox>& a_bdry_data,
                                  const LevelData<FArrayBox>& a_dfn) const
{
   /*
    Copy distribution function from the last valid cell
    into the bdry_data object that occupies one-cell-wide boundary layer.
    */
   const DisjointBoxLayout& bdry_grids( a_bdry_data.disjointBoxLayout() );
   for (DataIterator dit( bdry_grids ); dit.ok(); ++dit) {

      const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
      const FArrayBox& this_dfn( a_dfn[interior_dit] );
      FArrayBox& this_data( a_bdry_data[dit] );
      
      FArrayBox tmp(this_dfn.box(), 1);
      tmp.copy(this_dfn);
      
      Real garbage;
      /*FORT_COMPUTE_TOTAL_CURRENT( CHF_CONST_FRA1(this_dfn,0),
                                 CHF_BOX(tmp),
                                 CHF_REAL(garbage ));
      */
      const int& dir( m_bdry_layout->dir() );
      const Side::LoHiSide& side( m_bdry_layout->side() );
      // VG why shift?
      tmp.shift(dir, sign(side));// VG CHANGE
      this_data.copy(tmp);
      /*FORT_COMPUTE_TOTAL_CURRENT( CHF_CONST_FRA1(this_dfn,0),
                                 CHF_BOX(tmp),
                                 CHF_REAL(garbage ));      */
      /*const Box& tmp_box = flippedData[fdit].box();
	    for (BoxIterator bit(tmp_box); bit.ok(); ++bit) {
           IntVect iv = bit();
           if (::isnan(flippedData[fdit](iv, 0))) {
             cout<< iv << "  -- PRE nan-- " << flippedData[fdit](iv, 0) << endl;
             MayDay::Error( "LogicalSheathBC::applyBC: nan f_r found" );
           }
         }
        */ 
         
   }
}

//----------------------------------------------------------------------------------------------------------

void
InsulatingSheathBC_new::computeBoundaryElectronCurrent(LevelData<FArrayBox>& a_bdry_ion_current,
                                           const KineticSpeciesPtrVect& a_species,
                                           Real a_phi) const
{
   LevelData<FArrayBox> bdry_kernel(m_grids_full, 1, IntVect::Zero);
   for (DataIterator dit( m_grids_full ); dit.ok(); ++dit) {
      bdry_kernel[dit].setVal(0.0);
   }
   
   // for test
   LevelData<FArrayBox> bdry_kernel1(m_grids_full, 1, IntVect::Zero);
   for (DataIterator dit( m_grids_full ); dit.ok(); ++dit) {
      bdry_kernel1[dit].setVal(0.0);
   }

   LevelData<FArrayBox> this_bdry_data(m_grids_full, 1, IntVect::Zero);
   
   for (int s_index(0); s_index<a_species.size(); s_index++) {

      const KineticSpecies& species_physical( *(a_species[s_index]) );
      const double mass = species_physical.mass();
      const double charge = species_physical.charge();
   
      // VG if ( species_physical.charge() > 0.0 ) continue;

      if ( charge > 0.0 ) {continue;}    
      

      
      /* VG extra stuff
      */
      //Start
      const LevelData<FArrayBox>& dfn = species_physical.distributionFunction();
      fillBoundaryData(this_bdry_data, dfn);
      const DisjointBoxLayout& grids = dfn.getBoxes();
      
      Real vel_esc = sqrt(fabs(2.0*charge*a_phi/mass));
      const PhaseGeom& geometry( species_physical.phaseSpaceGeometry() );
      // Finish
      for (DataIterator dit( m_grids_full ); dit.ok(); ++dit) {
         /* VG extra stuff
         */
         //Start
        const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
        const FArrayBox& this_dfn( dfn[interior_dit] );
        const Box& vel_box = this_dfn.box();
      
         //const Box& vel_box = m_grids_full[dit];
         FArrayBox velocityRealCoords(vel_box, VEL_DIM);
         const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[interior_dit]);
         block_coord_sys.getVelocityRealCoords(velocityRealCoords);
         
         //Finish
         const Side::LoHiSide& side( m_bdry_layout->side() );
         const int& dir( m_bdry_layout->dir() );
         velocityRealCoords.shift(dir, sign(side));// VG CHANGE
         const Box& tmp_box = this_bdry_data[dit].box();
         const int sn = sign(side);
         
         Real c_mass = species_physical.mass();
         c_mass = sqrt(c_mass*c_mass*c_mass);
                  
         // Assume only one species of electrons!!!
         //Real c_mass = species_physical.mass();
         //c_mass = sqrt(c_mass*c_mass*c_mass);
         Real factor = sign(side)*species_physical.charge()/c_mass;
         FORT_COMPUTE_BC_KERNEL( CHF_FRA1(bdry_kernel1[dit],0),
                                 CHF_BOX(tmp_box),
                                 CHF_CONST_FRA1(this_bdry_data[dit],0),
                                 CHF_CONST_FRA1(velocityRealCoords,0),
                                 CHF_CONST_REAL(factor),
                                 CHF_CONST_REAL(vel_esc),
                                 CHF_CONST_INT(sn),
                                 CHF_CONST_REAL(m_dv_par_e[0]) ); 
                                 
         
         this_bdry_data[dit].mult(sign(side)*species_physical.charge());
         this_bdry_data[dit].divide(c_mass);  // because we do not know how current is computed
         this_bdry_data[dit].mult(velocityRealCoords);
         //bdry_kernel[dit].plus(this_bdry_data[dit]);
         //const Box& tmp_box = this_bdry_data[dit].box();
         for (BoxIterator bit(tmp_box); bit.ok(); ++bit) {
	         IntVect iv = bit();
            if (sign(side)*(iv[2]+0.5)<0) {this_bdry_data[dit](iv,0)=0;}
            else{
              Real vel_tmp = fabs(velocityRealCoords(iv,0));
              Real vel_hi = vel_tmp + 0.5*m_dv_par_e[0];
              Real vel_lo = vel_tmp - 0.5*m_dv_par_e[0];
              if(vel_hi<vel_esc) {this_bdry_data[dit](iv, 0) = 0.0; continue;}
              if(vel_lo>vel_esc) {continue;}
              this_bdry_data[dit](iv, 0) *= (vel_hi-vel_esc)/m_dv_par_e[0];
              ////if (fabs(vel_tmp)<vel_esc) {this_bdry_data[dit](iv, 0) = 0.0;}
            }
            //if (sign(side)==-1) {this_bdry_data[dit](iv,0)=0;}
	       }

         bdry_kernel[dit].plus(this_bdry_data[dit]);
         
         

      } 
   }
   
   // for test
   Real xxxx = 0;
   for (DataIterator dit( m_grids_full ); dit.ok(); ++dit) {
     const Box& tmp_box = bdry_kernel[dit].box();
     Real val;
     FORT_COMPARE_KERNELS( CHF_FRA1(bdry_kernel[dit],0), CHF_FRA1(bdry_kernel1[dit],0), CHF_BOX(tmp_box), CHF_REAL(val));
     xxxx += val; 
   }
      
   /*
    Multiply by the proper kernel here (e.g., v_parallel)
   */
   
   //**************computeIntegratedMomentFluxNormals compare
   const DisjointBoxLayout& dbl_curr = a_bdry_ion_current.getBoxes();
   Real tmp_current = simple_sum(bdry_kernel);
   for (DataIterator dit( dbl_curr ); dit.ok(); ++dit) {
      a_bdry_ion_current[dit].setVal(tmp_current);
   }
      if (procID()==0) {
   cout<<"ELE CURRENT  "<<tmp_current<<endl;
   }
   
   //sum(bdry_kernel, a_bdry_ion_current);
}

//----------------------------------------------------------------------------------------------------------

void
InsulatingSheathBC_new::computeBoundaryIonCurrent(LevelData<FArrayBox>& a_bdry_ion_current,
                                           const KineticSpeciesPtrVect& a_species) const
{
   LevelData<FArrayBox> bdry_kernel(m_grids_full, 1, IntVect::Zero);
   for (DataIterator dit( m_grids_full ); dit.ok(); ++dit) {
      bdry_kernel[dit].setVal(0.0);
   }

   LevelData<FArrayBox> this_bdry_data(m_grids_full, 1, IntVect::Zero);
   
   for (int s_index(0); s_index<a_species.size(); s_index++) {

      const KineticSpecies& species_physical( *(a_species[s_index]) );
      
      // VG if ( species_physical.charge() > 0.0 ) continue;

      if ( species_physical.charge() < 0.0 ) {continue;}      
      
      const LevelData<FArrayBox>& dfn = species_physical.distributionFunction();
      //fillBoundaryData(this_bdry_data, dfn);
      fillBoundaryData(this_bdry_data, dfn);
      const DisjointBoxLayout& grids = dfn.getBoxes();
       
      /* VG extra stuff
      */
      //Start
      const PhaseGeom& geometry( species_physical.phaseSpaceGeometry() );
      // Finish
      for (DataIterator dit( m_grids_full ); dit.ok(); ++dit) {
         /* VG extra stuff
         */
         //Start
        const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
        const FArrayBox& this_dfn( dfn[interior_dit] );
        const Box& vel_box = this_dfn.box();
      
         //const Box& vel_box = m_grids_full[dit];
         FArrayBox velocityRealCoords(vel_box, VEL_DIM);
         const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(grids[interior_dit]);
         block_coord_sys.getVelocityRealCoords(velocityRealCoords);
         //Finish
         const Side::LoHiSide& side( m_bdry_layout->side() );
         const int& dir( m_bdry_layout->dir() );
         velocityRealCoords.shift(dir, sign(side));// VG CHANGE
        //inspect(velocityRealCoords);
         const Box& tmp_box = this_bdry_data[dit].box();
         /*Real garbage = 0; // for local test only
         for (BoxIterator bit(tmp_box); bit.ok(); ++bit) {
	         IntVect iv = bit();
                garbage =    velocityRealCoords(iv, 0);
                garbage = 0;
	       }*/
         /*      Real garbage;
         FORT_COMPUTE_TOTAL_CURRENT( CHF_CONST_FRA1(velocityRealCoords,0),
                                 CHF_BOX(tmp_box),
                                 CHF_REAL(garbage ));*/
                                 
         //this_bdry_data[dit].mult(species_physical.charge());
         Real c_mass = species_physical.mass();
         c_mass = sqrt(c_mass*c_mass*c_mass);
         this_bdry_data[dit].mult(sign(side)*species_physical.charge());
         this_bdry_data[dit].divide(c_mass);  // because we do not know how current is computed
         this_bdry_data[dit].mult(velocityRealCoords);
         for (BoxIterator bit(tmp_box); bit.ok(); ++bit) {
	         IntVect iv = bit();
            if (sign(side)*(iv[2]+0.5)<0) {this_bdry_data[dit](iv,0)=0;}
            //if (sign(side)==-1) {this_bdry_data[dit](iv,0)=0;}
	       }
         bdry_kernel[dit].plus(this_bdry_data[dit]);
      }
   }
   
   
      
   /*
    Multiply by the proper kernel here (e.g., v_parallel)
   */
   
   //sum(bdry_kernel, a_bdry_ion_current);
   
   const DisjointBoxLayout& dbl_curr = a_bdry_ion_current.getBoxes();
   Real tmp_current = simple_sum(bdry_kernel);
   for (DataIterator dit( dbl_curr ); dit.ok(); ++dit) {
      a_bdry_ion_current[dit].setVal(tmp_current);
   }
   if (procID()==0) {
   cout<<"ION CURRENT  "<<tmp_current<<endl;
   }
/*
   ---------------
   #include "SpaceUtils.H.multidim"

std::string advScheme ="uw3";
const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );

// Get mapped components of parallel velocity
LevelData<FluxBox> normal_vel( grids, 1, IntVect::Unit );
geometry.computeMetricTermProductAverage(normal_vel, a_velocity, false);
   
bool fourth_order = (geometry.secondOrder()) ? false : true;
LevelData<FluxBox> faceDist( grids, 1, IntVect::Unit );
      // Overwrite advection quantity BCs with extrapolation
      const KineticSpecies& species_physical( *(a_species[s_index]) );
      if ( species_physical.charge() > 0.0 ) continue;
      const LevelData<FArrayBox>& dfn = species_physical.distributionFunction();
      
      const IntVect& ghosts = dfn.ghostVect();
      LevelData<FArrayBox> tmp(grids, 1, ghosts );
      for (DataIterator dit( a_dist.dataIterator() ); dit.ok(); ++dit) {
         tmp[dit].copy(a_dist[dit]);
      }
      CFG::MagGeom& mag_geom = geometry.magGeom();
      mag_geom.extrapolateToPhysicalGhosts(tmp, fourth_order);
SpaceUtils::upWindToFaces(faceDist, tmp, normal_vel, m_advScheme, fourth_order);
// Compute flux
LevelData<FluxBox> flux( grids, CFG_DIM, IntVect::Unit );
computeFlux(flux, faceDist, a_velocity, fourth_order);
   
   void TwoFieldNeutralsOp::computeAdvectionRHS(LevelData<FArrayBox>&         a_div,
                                             const LevelData<FArrayBox>&   a_dist,
                                             const LevelData<FluxBox>&     a_velocity,
                                             const bool                    a_recycling_bc,
                                             const bool                    a_homogeneous_flux_bc,
                                             const Real                    a_time) const
{

   // Get boxes
   const DisjointBoxLayout& grids( a_div.getBoxes() );

   bool fourth_order = (m_geometry.secondOrder()) ? false : true;

   // Get mapped components of parallel velocity
   LevelData<FluxBox> normal_vel( grids, 1, IntVect::Unit );
   m_geometry.computeMetricTermProductAverage(normal_vel, a_velocity, false);

   // Get distribution on faces
   LevelData<FluxBox> faceDist( grids, 1, IntVect::Unit );
   if (m_include_advection_bc) {
      SpaceUtils::upWindToFaces(faceDist, a_dist, normal_vel, m_advScheme, fourth_order);
   }
   else {
      // Overwrite advection quantity BCs with extrapolation
      const IntVect& ghosts = a_dist.ghostVect();
      LevelData<FArrayBox> tmp(grids, 1, ghosts );
      for (DataIterator dit( a_dist.dataIterator() ); dit.ok(); ++dit) {
         tmp[dit].copy(a_dist[dit]);
      }
      m_geometry.extrapolateToPhysicalGhosts(tmp, fourth_order);
      SpaceUtils::upWindToFaces(faceDist, tmp, normal_vel, m_advScheme, fourth_order);
   }
// Compute flux
   LevelData<FluxBox> flux( grids, CFG_DIM, IntVect::Unit );
   computeFlux(flux, faceDist, a_velocity, fourth_order);

   // Compute integrated normal flux
   LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
   m_geometry.computeMetricTermProductAverage(NTF_normal, flux, fourth_order);

   // Average the normal component of NTF consistent across block interfaces
   m_geometry.averageAtBlockBoundaries(NTF_normal);

   // Apply recycling BC
   if (a_recycling_bc) {

      // Apply zero flux BC
      if (a_homogeneous_flux_bc) {
         LevelData<FluxBox> zero_flux(grids, 1, IntVect::Zero);
         setZero(zero_flux);
         m_fluid_bc.at(0)->applyRecyclingBC(NTF_normal, zero_flux, a_time);
      }

      // Apply ion flux BC. We should only do it for the perpendicular ion momentum advection.
      else {
        LevelData<FluxBox> tmp(grids, 1, IntVect::Zero);
        for (DataIterator dit(grids); dit.ok(); ++dit) {
          for (int dir=0; dir<SpaceDim; ++dir) {
            tmp[dit][dir].copy(m_ion_normal_flux[dit][dir]);
            tmp[dit][dir].mult(faceDist[dit][dir]);
          }
        }
        m_fluid_bc.at(0)->applyRecyclingBC(NTF_normal, tmp, a_time);
      }
   }

   ---------------
*/

}

//----------------------------------------------------------------------------------------------------------

void
InsulatingSheathBC_new::sum( const LevelData<FArrayBox>& a_src,
                      LevelData<FArrayBox>&       a_dst ) const
{
   /*
    Sum over the velocity space
    */
   
   CH_assert(a_src.nComp() == a_dst.nComp());

   const DisjointBoxLayout& dst_grids = a_dst.getBoxes();
   const DisjointBoxLayout& src_grids = a_src.getBoxes();
   const ProblemDomain& problem_domain = src_grids.physDomain();

   // Initialize the destination, since SumOp does not do that for us.
   DataIterator dit = a_dst.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
      a_dst[dit].setVal(0.);
   }
   
   DisjointBoxLayout grids_tmp;
   adjCellLo(grids_tmp, src_grids, VPARALLEL_DIR, -1);

   // Initialize the destination, since SumOp does not do that for us.
   LevelData<FArrayBox> tmp(grids_tmp, a_src.nComp(), IntVect::Zero);
   for (DataIterator dit(grids_tmp); dit.ok(); ++dit) {
      tmp[dit].setVal(0.);
   }
   
   // Define ReductionCopier to compute intersections (sum in the poloidal direction)
   ReductionCopier reduceCopierVp(src_grids, grids_tmp, problem_domain, VPARALLEL_DIR);
   
   SumOp opVp(VPARALLEL_DIR);
   opVp.scale = 1.0;

   // Do the summing operation -- sums data in src along the vparallel direction
   a_src.copyTo(a_src.interval(), tmp, tmp.interval(), reduceCopierVp, opVp);

   // Define ReductionCopier to compute intersections (sum in the mu  direction)
   ReductionCopier reduceCopierMu(grids_tmp, dst_grids, problem_domain, MU_DIR);
   
   SumOp opMu(MU_DIR);
   opMu.scale = 1.0;

   // Do the summing operation -- sums data in src along the nu direction
   tmp.copyTo(tmp.interval(), a_dst, a_dst.interval(), reduceCopierMu, opMu);
}


Real
InsulatingSheathBC_new::simple_sum( const LevelData<FArrayBox>& a_src, bool a_flag) const
{
  /*
    A sum over all x, y, v_par, and mu
  */
  const DisjointBoxLayout& grids = a_src.getBoxes();
  Real current_all = 0.0;
  DataIterator dit = grids.dataIterator();
  Real current_loc = 0.0;
  for (dit.begin(); dit.ok(); ++dit) {
     Real curr_loc_tmp = 0.0;
     const FArrayBox& this_dfn = a_src[dit];
     Box this_box = this_dfn.box();
     FORT_COMPUTE_TOTAL_CURRENT( CHF_CONST_FRA1(this_dfn,0),
                                 CHF_BOX(this_box),
                                 CHF_REAL(curr_loc_tmp) );
    current_loc += curr_loc_tmp;
   }
   if (a_flag) {
     current_loc=1.0;
     std::cout<<to_string(procID())<<"  "<<current_loc<<endl;
     //if (procID()==0) {std::cout<<"\n\n";}
   }
   
   // Sum over all ranks.
   // We could also use mpi_comm local to this boundary date (if needed for speed-up)
   MPI_Allreduce(&current_loc, &current_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return current_all;
}

//-----------------------------------------------------------------------------------

void InsulatingSheathBC_new::applyBC( KineticSpeciesPtrVect& a_species,
                               const int& a_species_index,
                               const LevelData<FluxBox>& a_velocity,
                               const CFG::LevelData<CFG::FArrayBox>& a_phi)
{
   //////// TMP degug of potential
   if (0&& print_ghosts==164)
   {
     KineticSpecies& species_physical_i( *(a_species[0]) );
     KineticSpecies& species_physical_e( *(a_species[1]) );
     const PhaseGeom& geometry_i( species_physical_i.phaseSpaceGeometry() );
     const CFG::MagGeom& mag_geom = geometry_i.magGeom();
     const CFG::DisjointBoxLayout& mag_grids = mag_geom.gridsFull();
     CFG::LevelData<CFG::FArrayBox> n_i( mag_grids, 1, CFG::IntVect::Zero );
     CFG::LevelData<CFG::FArrayBox> n_e( mag_grids, 1, CFG::IntVect::Zero );
     species_physical_i.chargeDensity(n_i);
     species_physical_e.chargeDensity(n_e);

     mag_geom.plotCellData( "phi_at_164", a_phi, 80.0);
     mag_geom.plotCellData( "n_deu_at_164.", n_i, 80.0);
     mag_geom.plotCellData( "n_ele_at_164.", n_e, 80.0);
     
   }
   
   //////////////// Diagnostics to print distributions fucntion and 2 ghost cells for each mu and v par
    if (0 && procID()==0) {
      KineticSpecies& species_physical_i( *(a_species[0]) );
      KineticSpecies& species_physical_e( *(a_species[1]) );
      const LevelData<FArrayBox>& dfn_i = species_physical_i.distributionFunction();
      const LevelData<FArrayBox>& dfn_e = species_physical_e.distributionFunction();
      
      string f_name = "output_i00.txt";
      ofstream fout_i0(f_name);
      f_name = "output_i63.txt";
      ofstream fout_i63(f_name);
      f_name = "output_e00.txt";
      ofstream fout_e0(f_name);
      f_name = "output_e63.txt";
      ofstream fout_e63(f_name);
      DataIterator dit = m_grids_full.dataIterator();   
      for (dit.begin(); dit.ok(); ++dit){
        const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
        //const FArrayBox& this_dfn_i( dfn_i[interior_dit] );
        //const FArrayBox& this_dfn_e( dfn_e[interior_dit] );
        for(int im=0;im<64;im++) {
          for(int iv=-64;iv<64;iv++){
            IntVect iv_tmp( IntVect::Unit );
            iv_tmp[0]=0;
            iv_tmp[1]=-2;
            iv_tmp[2]=iv;
            iv_tmp[3]=im;
            double tmp_2 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=-1;
            double tmp_1 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=0;
            double tmp0 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=1;
            double tmp1 = dfn_i[interior_dit](iv_tmp,0);
            fout_i0 << im << "  " << iv << "  " << tmp_2 << "  " << tmp_1 << "  " << tmp0 << "  " << tmp1 <<endl;
            // Do in a reverse order: domain cells first, then ghost cells
            iv_tmp[1]=62;
            tmp_2 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=63;
            tmp_1 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=64;
            tmp0 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=65;
            tmp1 = dfn_i[interior_dit](iv_tmp,0);
            fout_i63 << im << "  " << iv << "  " << tmp_2 << "  " << tmp_1 << "  " << tmp0 << "  " << tmp1 <<endl;
            // Repeat for electrons
            iv_tmp[1]=-2;
            tmp_2 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=-1;
            tmp_1 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=0;
            tmp0 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=1;
            tmp1 = dfn_e[interior_dit](iv_tmp,0);
            fout_e0 << im << "  " << iv << "  " << tmp_2 << "  " << tmp_1 << "  " << tmp0 << "  " << tmp1 <<endl;
            // Do in a reverse order: domain cells first, then ghost cells
            iv_tmp[1]=62;
            tmp_2 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=63;
            tmp_1 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=64;
            tmp0 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=65;
            tmp1 = dfn_e[interior_dit](iv_tmp,0);
            fout_e63 << im << "  " << iv << "  " << tmp_2 << "  " << tmp_1 << "  " << tmp0 << "  " << tmp1 <<endl;
            }
        }
      }
      fout_i0.close();
      fout_i63.close();
      fout_e0.close();
      fout_e63.close();
      MayDay::Error( "GHOST TEST DONE" );
   //////////////
    }
   
   
   if (0&& procID()==0) {
     cout<<"current print_ghost: "<<print_ghosts<<endl;
   }
   
   ////// Tmp output for ghost
   if (0 && print_ghosts==400) // means 4 RK steps, 2 species, 50 steps
   {////// diable for now
     KineticSpecies& species_physical_i( *(a_species[0]) );
     KineticSpecies& species_physical_e( *(a_species[1]) );
     const LevelData<FArrayBox>& dfn_i = species_physical_i.distributionFunction();
     const LevelData<FArrayBox>& dfn_e = species_physical_e.distributionFunction();
     const PhaseGeom& geometry_i( species_physical_i.phaseSpaceGeometry() );
     const PhaseGeom& geometry_e( species_physical_e.phaseSpaceGeometry() );
     geometry_i.plotData( "dfn_i_before", dfn_i, 0, 0);
     geometry_e.plotData( "dfn_e_before", dfn_e, 0, 0);
     //if (m_time>0.07) {MayDay::Error( "GHOST TEST DONE" );}
     
   }

     
    print_ghosts++;
   
   const int& dir( m_bdry_layout->dir() );
   const Side::LoHiSide& side( m_bdry_layout->side() );
   
   KineticSpecies& species_physical( *(a_species[a_species_index]) );
   LevelData<FArrayBox>& soln( species_physical.distributionFunction() );
   
   const double mass = species_physical.mass();
   const double charge = species_physical.charge();
   ////// if (charge>0) {return;} TMP FIX to set fully reflective BC
   if (charge>0) { return;
     Real pp[2] = {0.0, 0.0}; // this does not matter, as the following method will ignore pp and a_phi anyway
     applyBC_iter(pp, a_species, 0, a_velocity, a_phi);
     return;
   }
   
   
   const PhaseGeom& geometry( species_physical.phaseSpaceGeometry() );
   geometry.injectConfigurationToPhase(a_phi, m_phi_injected);     
   
   const VEL::VelCoordSys& vel_coords = geometry.velSpaceCoordSys();
   const VEL::RealVect& vel_dx = vel_coords.dx();
   const Real dv_par = vel_dx[0];
   m_dv_par_e[0] = dv_par;
   const Box& domain_box = (geometry.domain()).domainBox();
   IntVect lo_end(domain_box.smallEnd());
   IntVect hi_end(domain_box.bigEnd());
   int v_min = lo_end[VPARALLEL_DIR];
   int v_max = hi_end[VPARALLEL_DIR];
   int N_v = (v_max - v_min + 1)/2;
   double vel_max = dv_par*N_v;
   m_v_par_max_e[0] = vel_max;
   
   Real phi_hi = fabs(mass*vel_max*vel_max/2/charge);
   Real phi_MAX = phi_hi;
   Real phi_lo = 0.0;
   Real d_phi = (phi_hi-phi_lo)/1000; // Define a small d_phi like 1/1000

   
   // Total boundary currents I_i - ion, I_e - electron, I_s - sum of the two
   // Currently implemented for two species plasmas only
   Real I_i, I_e, I_s;
   // p_res is the array to sent to applyBC_iter that will get the species current [0] and send the potential [1]
   Real p_res[2];
   p_res[0] = 0; p_res[1] = (phi_hi+phi_lo)/2;
   Real phi_prev;
   // Check number of iterations
   int iters = m_iter_number;
        
   // For t=0 use bisect method, as we do not know initail value for the potential
   //if(m_time<m_tolerance) {
   ////// the following lines turn on anf off newton solver
   if (!m_id_map.count(m_id) || m_newton==false) 
   //////if (true)
   {
     // Apply for ions and get the ion current
     applyBC_iter(p_res, a_species, 0, a_velocity, a_phi);
     I_i = p_res[0];
   
     for (int i=0; i<m_iter_number; i++){
       applyBC_iter(p_res, a_species, 1, a_velocity, a_phi);
       I_e = p_res[0];
       I_s = sign(side)*(I_i+I_e);
       // Save previous iteration of phi
       phi_prev = p_res[1];
       if (I_s>0) {
         phi_hi = p_res[1];
         p_res[1] = (p_res[1] + phi_lo)/2;
       }
       else {
         phi_lo = p_res[1];
         p_res[1] = (p_res[1] + phi_hi)/2;
       }
     }
     ////// TMP fix v_esc = phi_high
     ////// p_res[1] = phi_MAX;
     
     // Update the id map with the solution for phi
     //if (m_id_map.count(m_id)) {MayDay::Error( "Boundary id map is corrupt!");}
     m_id_map.insert ( std::pair<int,Real>(m_id, p_res[1]) );
   }
   
   // For all other steps when the previous solution is known and the map is created
   else {
     // Check if map is ready map
     if (!m_id_map.count(m_id)) {MayDay::Error( "Boundary id map is not initialized!");}
     // Get the potential guess
     Real phi_0 = m_id_map[m_id];
     p_res[1] = phi_0;
     phi_prev = phi_0;
            
     // Apply for ions
     applyBC_iter(p_res, a_species, 0, a_velocity, a_phi);
     I_i = p_res[0];
     // Apply for electrons
     applyBC_iter(p_res, a_species, 1, a_velocity, a_phi);
     I_e = p_res[0];
     
     // Prepare variables for Newton iterations
     Real y0 = sign(side)*(I_e+I_i);
     Real a, b, x0, x1, y1;
     x1 = phi_0 + d_phi;
     
     // Begin Newton iterations
     // y0 = a*x0+b, y1 = a*x1_b, where x0 is phi_0, x1 is phi+d_phi, y0 = sign(side)*(I_e+I_i), y1 = the same at phi+d_phi
     for (int i=0; i<m_iter_number; i++){
       p_res[1] = x1;
       // Compute electron current once again
       applyBC_iter(p_res, a_species, 1, a_velocity, a_phi);
       I_e = p_res[0];
       y1 = sign(side)*(I_e+I_i);
       // Save Newton solver from division to zero
       if(fabs(y1-y0)<m_tolerance) {iters = i; break;}
       a = (y1-y0)/d_phi;
       b = y1 - a*x1;
       // Update x1 and x0 
       x0 = x1;
       x1 = -b/a;
       // Update y0
       y0 = y1;
       d_phi = x1-x0;
       // Save previous iteration of phi
       phi_prev = x0;
     }
     // Update the map
     m_id_map[m_id] = p_res[1];
     ////p_res[1] = x1;
   } 

   Real residual = fabs(p_res[1]-phi_prev);
   // Print convergence results
   if (procID()==0 && m_debug) {
     cout<<"boundary ID = "<<m_id<<"   dir = "<<dir<<"   side = "<<sign(side)<<endl;
     cout<<"phi = "<<p_res[1]<<"   I_i = " <<I_i<<"   I_e = "<<I_e<<"\n";
     cout<<"iterations = "<<iters<<"   residual = "<<residual<<endl;
     cout<<"tolerance = "<<m_tolerance<<endl;
   }
   
   if (procID()==0) {
     if(::isnan(residual) || residual>0.01) {
       cout<<"WARNING: residual of Insulating sheath failed!\n";
       cout<<"boundary ID = "<<m_id<<"   dir = "<<dir<<"   side = "<<sign(side)<<endl;
       cout<<"phi = "<<p_res[1]<<"   I_i = " <<I_i<<"   I_e = "<<I_e<<"\n";
       cout<<"iterations = "<<iters<<"   residual = "<<residual<<endl;
     }
   }
   

   if(0 && sign(side)>0)
   {
    if (procID()==0) {
      KineticSpecies& species_physical_i( *(a_species[0]) );
      KineticSpecies& species_physical_e( *(a_species[1]) );
      const LevelData<FArrayBox>& dfn_i = species_physical_i.distributionFunction();
      const LevelData<FArrayBox>& dfn_e = species_physical_e.distributionFunction();
      
      string f_name = "output_i00.txt";
      ofstream fout_i0(f_name);
      f_name = "output_i63.txt";
      ofstream fout_i63(f_name);
      f_name = "output_e00.txt";
      ofstream fout_e0(f_name);
      f_name = "output_e63.txt";
      ofstream fout_e63(f_name);
      DataIterator dit = m_grids_full.dataIterator();   
      for (dit.begin(); dit.ok(); ++dit){
        const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
        //const FArrayBox& this_dfn_i( dfn_i[interior_dit] );
        //const FArrayBox& this_dfn_e( dfn_e[interior_dit] );
        for(int im=0;im<64;im++) {
          for(int iv=-64;iv<64;iv++){
            IntVect iv_tmp( IntVect::Unit );
            iv_tmp[0]=0;
            iv_tmp[1]=-2;
            iv_tmp[2]=iv;
            iv_tmp[3]=im;
            double tmp_2 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=-1;
            double tmp_1 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=0;
            double tmp0 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=1;
            double tmp1 = dfn_i[interior_dit](iv_tmp,0);
            fout_i0 << im << "  " << iv << "  " << tmp_2 << "  " << tmp_1 << "  " << tmp0 << "  " << tmp1 <<endl;
            // Do in a reverse order: domain cells first, then ghost cells
            iv_tmp[1]=62;
            tmp_2 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=63;
            tmp_1 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=64;
            tmp0 = dfn_i[interior_dit](iv_tmp,0);
            iv_tmp[1]=65;
            tmp1 = dfn_i[interior_dit](iv_tmp,0);
            fout_i63 << im << "  " << iv << "  " << tmp_2 << "  " << tmp_1 << "  " << tmp0 << "  " << tmp1 <<endl;
            // Repeat for electrons
            iv_tmp[1]=-2;
            tmp_2 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=-1;
            tmp_1 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=0;
            tmp0 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=1;
            tmp1 = dfn_e[interior_dit](iv_tmp,0);
            fout_e0 << im << "  " << iv << "  " << tmp_2 << "  " << tmp_1 << "  " << tmp0 << "  " << tmp1 <<endl;
            // Do in a reverse order: domain cells first, then ghost cells
            iv_tmp[1]=62;
            tmp_2 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=63;
            tmp_1 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=64;
            tmp0 = dfn_e[interior_dit](iv_tmp,0);
            iv_tmp[1]=65;
            tmp1 = dfn_e[interior_dit](iv_tmp,0);
            fout_e63 << im << "  " << iv << "  " << tmp_2 << "  " << tmp_1 << "  " << tmp0 << "  " << tmp1 <<endl;
            }
        }
      }
      fout_i0.close();
      fout_i63.close();
      fout_e0.close();
      fout_e63.close();
      MayDay::Error( "GHOST TEST DONE" );
   //////////////
    }
   }
   // This it test why I_e is of wrong sigh
   if(I_i*I_e>0){////// TMP 0 && to check fully reflective ion BC
   
     const IntVect& ghost_vect( soln.ghostVect() );
     //VG if (side == Side::Lo) hi_end[dir] = lo_end[dir] + ghost_vect[dir];
     //VG if (side == Side::Hi) lo_end[dir] = hi_end[dir] - ghost_vect[dir];
     if (side == Side::Lo) hi_end[dir] = lo_end[dir] + 4;
     if (side == Side::Hi) lo_end[dir] = hi_end[dir] - 4;
                   
     Box refl_bnd_box(lo_end, hi_end);
     const DisjointBoxLayout& dbl = soln.getBoxes();
     int reflectDir = VPARALLEL_DIR;
     int reflectCoord = 0;
     Vector<Tuple<DataIndex, 2> > boxCorrelation;
   
     DisjointBoxLayout flippedGrids;
     getFlippedGrids(flippedGrids, boxCorrelation, dbl, refl_bnd_box, reflectDir, reflectCoord);
     LevelData<FArrayBox> flippedData(flippedGrids, 1);
     soln.copyTo(flippedData);
           
     // Iterate over patches of flipped data, a small subset of the full data
     DataIterator fdit = flippedData.dataIterator();
     for (fdit.begin(); fdit.ok(); ++fdit) {
            
     // find the iterator value for the UNFLIPPED data corresponding to this flipped data
     DataIndex regDataIndex;
     for (int n=0; n<boxCorrelation.size(); n++)
     {
        if (boxCorrelation[n][1] == fdit() )
        {
           regDataIndex = boxCorrelation[n][0];
        }
     }
            
     FArrayBox& this_dfn = soln[regDataIndex];
     const Box& this_box = dbl[regDataIndex];
     Box boundaryBox = adjCellBox(this_box, dir, side, ghost_vect[dir]);
      
     std::string str_out;
     int ind = fdit.i().intCode();
     string f_name = "output" + std::to_string(ind) + ".txt";
     ofstream fout(f_name);
     ////IntVect iv_lo = this_box.loVect();
     IntVect iv_lo = boundaryBox.smallEnd();
     IntVect iv_hi = boundaryBox.bigEnd();
     ////IntVect iv_hi = this_box.hiVect();
     int pID = flippedGrids.procID(fdit.i());
     fout<<"box_m_grids   "<<pID<<std::endl;
     str_out = "lo: "+std::to_string(iv_lo[0]) + " " + std::to_string(iv_lo[1]) + " " + std::to_string(iv_lo[2]) + " " + std::to_string(iv_lo[3]);
     fout<<str_out<<std::endl;
     str_out = "hi: "+std::to_string(iv_hi[0]) + " " + std::to_string(iv_hi[1]) + " " + std::to_string(iv_hi[2]) + " " + std::to_string(iv_hi[3]);
     fout<<str_out<<std::endl;
     fout<<std::endl;
      
     for (BoxIterator bit(boundaryBox); bit.ok(); ++bit) {
       IntVect iv = bit();
       str_out = std::to_string(iv[0]) + " " + std::to_string(iv[1]) + " " + std::to_string(iv[2]) + " " + std::to_string(iv[3]) + " " + std::to_string(this_dfn(iv, 0));
       fout<<str_out<<std::endl;
     }
     fout.close(); 
      
     }
     MayDay::Error( "Test of dfn is done");
   }
   

   
   // Save potential history
   if(procID()==0 && m_id==0) {
     string filename = "phi_history.txt";
     ofstream fout(filename, std::ios_base::app);
     fout<<m_time<<"   "<<p_res[1]<<endl;
     fout.close();
   }
      // Save heat fluxes
   if(m_id==1) {
     Real p_hf_i[3];// Two components for deu and two for ele: deu flux, ele flux, deu mv^2, deu mu b, ele mv^2, ele mv
     Real p_hf_e[3];
     computeHeatFluxes( p_hf_i, a_species, 0, a_velocity);
     computeHeatFluxes( p_hf_e, a_species, 1, a_velocity);
     if (procID()==0){
       string filename = "parallel_heat_flux_history.txt";
       ofstream fout(filename, std::ios_base::app);
       fout<<m_time<<"   "<<p_hf_i[0]<<"   "<<p_hf_i[1]<<"   "<<p_hf_i[2]<<"   "<<p_hf_e[0]<<"   "<<p_hf_e[1]<<"   "<<p_hf_e[2]<<endl;
       fout.close();
     }
   }
   
   /*
   // Test Newton solver
   
   // Check if map is ready map
   if (!m_id_map.count(m_id)) {MayDay::Error( "Boundary id map is not initialized!");}
   // Get the potential guess
   phi_0 = m_id_map[m_id];
   
   
   // Getting the value of the previos potential from file (yes, it is sick) and the value of the phi shift for the first Newton iteration
     std::map<char,int> mymap;

  // first insert function version (single parameter):
  mymap.insert ( std::pair<char,int>('a',100) );
  mymap.insert ( std::pair<char,int>('z',200) );
   
   double phi_0, d_phi;
   if(procID()==0){
     string filename = "phi_now.txt";
     ifstream fin(filename, std::ios_base::in);
     fin>>phi_0;
     fin>>d_phi;
     fin.close();
     if (d_phi<-1.0e6) {d_phi=phi_hi/1000;} // Define a small d_phi like 1/1000}
   }
   MPI_Bcast( &phi_0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Bcast( &d_phi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   
   // Get currents
   Real I_i, I_e;
   // p_res is the array to sent to applyBC_iter that will get the species current [0] and send the potential [1]
   p_res[0] = 0; p_res[1] = phi_0;
   // Apply for ions
   applyBC_iter(p_res, a_species, 0, a_velocity, a_phi);
   I_i = p_res[0];
   // Apply for electrons
   applyBC_iter(p_res, a_species, 1, a_velocity, a_phi);
   I_e = p_res[0];
   Real y0 = sign(side)*(I_e+I_i);
   Real a, b, x0, x1, y1;
   x1 = phi_0 + d_phi;
   
   // Begin Newton iterations
   // y0 = a*x0+b, y1 = a*x1_b, where x0 is phi_0, x1 is phi+d_phi, y0 = sign(side)*(I_e+I_i), y1 = the same at phi+d_phi
   for (int i=0; i<m_iter_number; i++){
     p_res[1] = x1;
     // Compute electron current once again
     applyBC_iter(p_res, a_species, 1, a_velocity, a_phi);
     I_e = p_res[0];
     y1 = sign(side)*(I_e+I_i);
     a = (y1-y0)/d_phi;
     b = y1 - a*x1;
     // Update x1 and x0 
     x0 = x1;
     x1 = -b/a;
     // Update y0
     y0 = y1;
     d_phi = x1-x0;
   }
   
   //The final update
   p_res[1] = x1;
   applyBC_iter(p_res, a_species, 1, a_velocity, a_phi);
   I_e = p_res[0];
   
   // Save changes back to the file
   if(procID()==0){
     string filename = "phi_now.txt";
     ofstream fout(filename);
     fout<<x1<<endl;
     fout<<(x1-phi_0);
     fout.close();
     // Save history too
     filename = "phi_history.txt";
     ofstream fout_h(filename, std::ios_base::app);
     fout_h<<x1<<endl;
     fout_h.close();
   }
   if (procID()==0)  {cout<<"dir = "<<dir<<"   "<<"side = "<<sign(side)<<endl;}
   if(procID()==0) {cout<<"phi = "<<p_res[1]<<"  I_i = " <<I_i<<"  I_e = "<<I_e<<"\n"<<"residual = "<<fabs(x1-x0)<<endl;}
   */
   
   ////// When we are done with electrons, save all pdf of electrons and ions
   //////if(a_species_index==1) {
   if(0 && a_species_index==1 && print_ghosts==400) {// disable for now
     KineticSpecies& species_physical_i( *(a_species[0]) );
     KineticSpecies& species_physical_e( *(a_species[1]) );
     const LevelData<FArrayBox>& dfn_i = species_physical_i.distributionFunction();
     const LevelData<FArrayBox>& dfn_e = species_physical_e.distributionFunction();
     const PhaseGeom& geometry_i( species_physical_i.phaseSpaceGeometry() );
     const PhaseGeom& geometry_e( species_physical_e.phaseSpaceGeometry() );
     geometry_i.plotData( "dfn_i_after", dfn_i, 0, 0);
     geometry_e.plotData( "dfn_e_after", dfn_e, 0, 0);
     if (m_time>0.07) {MayDay::Error( "GHOST TEST DONE" );}
   } 
}

//-----------------------------------------------------------------------------------

void InsulatingSheathBC_new::applyBC_iter( Real* p_res,
                                KineticSpeciesPtrVect& a_species,
                               const int& a_species_index,
                               const LevelData<FluxBox>& a_velocity,
                               const CFG::LevelData<CFG::FArrayBox>& a_phi)
{    
         
             
   
   /*
    This BC reflects particles that cannot make the potential barrier
    */
   const int& dir( m_bdry_layout->dir() );
   const Side::LoHiSide& side( m_bdry_layout->side() );
   
   KineticSpecies& species_physical( *(a_species[a_species_index]) );
   LevelData<FArrayBox>& soln( species_physical.distributionFunction() );

   const double mass = species_physical.mass();
   double charge = species_physical.charge();
            
   const PhaseGeom& geometry( species_physical.phaseSpaceGeometry() );
   ///geometry.injectConfigurationToPhase(a_phi, m_phi_injected);
   ///LevelData<FArrayBox> phi_injected; move this one from here to apply, since we do nto want to call this function many times 033123
   ///geometry.injectConfigurationToPhase(a_phi, phi_injected);
   
   Real phi_val = p_res[1];
   if (charge>0) {phi_val=0.0;}//////  TMP change to set fully reflective ion BC
   //////if (charge>0) {phi_val=10000000000.0; charge*=-1;}
   
   DataIterator dit_inj = m_phi_injected.dataIterator();
     for (dit_inj.begin(); dit_inj.ok(); ++dit_inj) {
         m_phi_injected[dit_inj].setVal(phi_val);
     }
   
   
   /*
   DataIterator dit0 = phi_injected.dataIterator();
   for (dit0.begin(); dit0.ok(); ++dit0) {
     const Box& this_box = phi_injected[dit0].box();
     std::string str_out;
     int ind = dit0.i().intCode();
     string f_name = "output" + std::to_string(ind) + ".txt";
     ofstream fout(f_name);
     IntVect iv_lo = this_box.smallEnd();
     IntVect iv_hi = this_box.bigEnd();
     fout<<"phi_injected"<<std::endl;
     str_out = "lo: "+std::to_string(iv_lo[0]) + " " + std::to_string(iv_lo[1]) + " " + std::to_string(iv_lo[2]) + " " + std::to_string(iv_lo[3]);
     fout<<str_out<<std::endl;
     str_out = "hi: "+std::to_string(iv_hi[0]) + " " + std::to_string(iv_hi[1]) + " " + std::to_string(iv_hi[2]) + " " + std::to_string(iv_hi[3]);
     fout<<str_out<<std::endl;
     fout.close();
   }
      
   LevelData<FArrayBox> phi_injected_tmp;
   phi_injected_tmp.define(m_grids_inj, 1, IntVect::Zero);
   DataIterator dit1 = phi_injected_tmp.dataIterator();
   for (dit1.begin(); dit1.ok(); ++dit1) {
     phi_injected_tmp[dit1].setVal(0.0);
     const Box& this_box = phi_injected_tmp[dit1].box();
     std::string str_out;
     int ind = dit1.i().intCode();
     string f_name = "output" + std::to_string(ind) + ".txt";
     ofstream fout(f_name, std::ios_base::app);
     IntVect iv_lo = this_box.smallEnd();
     IntVect iv_hi = this_box.bigEnd();
     fout<<"phi_injected_tmp"<<std::endl;
     str_out = "lo: "+std::to_string(iv_lo[0]) + " " + std::to_string(iv_lo[1]) + " " + std::to_string(iv_lo[2]) + " " + std::to_string(iv_lo[3]);
     fout<<str_out<<std::endl;
     str_out = "hi: "+std::to_string(iv_hi[0]) + " " + std::to_string(iv_hi[1]) + " " + std::to_string(iv_hi[2]) + " " + std::to_string(iv_hi[3]);
     fout<<str_out<<std::endl;
     fout.close();     
     
   }   
     MayDay::Error( "DONE" );         
     */         
              
  


         
   // Create valid-cell box that extends ghost_vect number
   // of cells away from the boundary
   const Box& domain_box = (geometry.domain()).domainBox();
   const IntVect& ghost_vect( soln.ghostVect() );
   const IntVect& ghost_vect_fr( 4*IntVect::Unit ); //ghost vect for reflected BC, needed 4 ghosts 
   const DisjointBoxLayout& dbl = soln.getBoxes();

   
   IntVect lo_end(domain_box.smallEnd());
   IntVect hi_end(domain_box.bigEnd());
   

   int N_ghost = 4;

   // Get dv_par
   const VEL::VelCoordSys& vel_coords = geometry.velSpaceCoordSys();
   const VEL::RealVect& vel_dx = vel_coords.dx();
   const Real dv_par = vel_dx[0];

   int v_min = lo_end[VPARALLEL_DIR];
   int v_max = hi_end[VPARALLEL_DIR];
   int N_v = v_max - v_min + 1 + 2*N_ghost;
   Real* arr_v = new Real[N_v];
   for (int i=0; i<N_v; ++i) {
     arr_v[i] = (i + 0.5 - N_v/2)*dv_par;
   }
   int y_min = lo_end[1];
   int y_max = hi_end[1];
   
   ////// added -1 as the last index is included if (side == Side::Lo) hi_end[dir] = lo_end[dir] + ghost_vect[dir];
   if (side == Side::Lo) {hi_end[dir] = lo_end[dir] + ghost_vect_fr[dir] - 1;} //VG changed ghost_vect to ghost_vect_fr to get 4 cells
   ////// added +1 as the last index is included if (side == Side::Hi) lo_end[dir] = hi_end[dir] - ghost_vect[dir];
   if (side == Side::Hi) {lo_end[dir] = hi_end[dir] - ghost_vect_fr[dir] + 1;}
         
   Box refl_bnd_box(lo_end, hi_end);
         
   // Create the flipped data object we will need for reflecting particles (nominally, electrons)
   //  below the potential barrier
   int reflectDir = VPARALLEL_DIR;
   int reflectCoord = 0;
   Vector<Tuple<DataIndex, 2> > boxCorrelation;
   
   DisjointBoxLayout flippedGrids;
   getFlippedGrids(flippedGrids, boxCorrelation, dbl, refl_bnd_box, reflectDir, reflectCoord);
   LevelData<FArrayBox> flippedData(flippedGrids, 1);
   soln.copyTo(flippedData);

// tmp
LevelData<FArrayBox> flippedData_gc(flippedGrids, 1, 2*IntVect::Unit); // added 4 ghosts
DataIterator fdit2 = flippedData.dataIterator();
   for (fdit2.begin(); fdit2.ok(); ++fdit2) {
	   flippedData_gc[fdit2].copy(flippedData[fdit2]);
   }
   flippedData_gc.exchange();
      
   // Iterate over patches of flipped data, a small subset of the full data
   DataIterator fdit = flippedData.dataIterator();
   for (fdit.begin(); fdit.ok(); ++fdit) {
            
      // find the iterator value for the UNFLIPPED data corresponding to this flipped data
      DataIndex regDataIndex;
      for (int n=0; n<boxCorrelation.size(); n++)
      {
         if (boxCorrelation[n][1] == fdit() )
         {
            regDataIndex = boxCorrelation[n][0];
         }
      }
            
      FArrayBox& this_dfn = soln[regDataIndex];
      const Box& this_box = dbl[regDataIndex];
      FArrayBox& this_phi = m_phi_injected[regDataIndex];
      
      // Because m_grids_inj corresponds to ghost cell layer,
      // need to shift phi_bc to the last valid cell
      if (m_compute_potential_BC) {
         this_phi.shift(dir, -sign(side));
      }

      const FluxBox& this_vel = a_velocity[regDataIndex];
                  
      Box boundaryBox = adjCellBox(this_box, dir, side, ghost_vect[dir]);
            
      FArrayBox velocityRealCoords(boundaryBox, VEL_DIM);
      const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(this_box);
      block_coord_sys.getVelocityRealCoords(velocityRealCoords);
            
      const int SIDE(side);
      if(m_advect_scheme=="uw1") {
        FORT_SET_LOGICAL_SHEATH_BC_SPECIAL2(CHF_FRA1(this_dfn,0),
                                            CHF_BOX(boundaryBox),
                                            CHF_CONST_FRA1(flippedData[fdit],0),
                                            CHF_CONST_R1D(arr_v, N_v),
                                            CHF_CONST_FRA1(this_vel[dir],dir),
                                            CHF_CONST_FRA1(this_phi,0),
                                            CHF_CONST_REAL(mass),
                                            CHF_CONST_REAL(charge),
                                            CHF_CONST_INT(SIDE) );
      }
      //else if(m_advect_scheme=="uw3" && a_species_index==0) {
          else if(m_advect_scheme=="uw3" && a_species_index>-1) {
        FORT_SET_LOGICAL_SHEATH_BC_UW3(CHF_FRA1(this_dfn,0),
                                       CHF_BOX(boundaryBox),
                                       CHF_CONST_FRA1(flippedData_gc[fdit],0),//CHF_CONST_FRA1(flippedData[fdit],0),
                                       CHF_CONST_R1D(arr_v, N_v),
                                       CHF_CONST_FRA1(this_vel[dir],dir),
                                       CHF_CONST_FRA1(this_phi,0),
                                       CHF_CONST_REAL(mass),
                                       CHF_CONST_REAL(charge),
                                       CHF_CONST_INT(SIDE) );
      }
      else if(m_advect_scheme=="weno" && a_species_index>-1) {
        FORT_SET_LOGICAL_SHEATH_BC_WENO(CHF_FRA1(this_dfn,0),
                                       CHF_BOX(boundaryBox),
                                       CHF_CONST_FRA1(flippedData_gc[fdit],0),//CHF_CONST_FRA1(flippedData[fdit],0),
                                       CHF_CONST_R1D(arr_v, N_v),
                                       CHF_CONST_FRA1(this_vel[dir],dir),
                                       CHF_CONST_FRA1(this_phi,0),
                                       CHF_CONST_REAL(mass),
                                       CHF_CONST_REAL(charge),
                                       CHF_CONST_INT(SIDE) );
      }
        else if(0 && m_advect_scheme=="uw3" && a_species_index==1) {
        FORT_SET_LOGICAL_SHEATH_BC_UW3_ELE(CHF_FRA1(this_dfn,0),
                                       CHF_BOX(boundaryBox),
                                       CHF_CONST_FRA1(flippedData_gc[fdit],0),//CHF_CONST_FRA1(flippedData[fdit],0),
                                       CHF_CONST_R1D(arr_v, N_v),
                                       CHF_CONST_FRA1(this_vel[dir],dir),
                                       CHF_CONST_FRA1(this_phi,0),
                                       CHF_CONST_REAL(mass),
                                       CHF_CONST_REAL(charge),
                                       CHF_CONST_INT(SIDE) );
      }
      else {MayDay::Error( "Advection scheme is not implemented in the BC class" );}
         
   }

   delete[] arr_v;
   // TMP test
   Real val_out = computeAdvection(a_species, a_species_index, a_velocity);
   p_res[0] = val_out;
   
   /*  
   if(a_species_index==0){
     
     //if (procID()==0) {cout<<"Ions are done.\n";}
   
     this->applyBC_iter( p_res, a_species, 1, a_velocity, a_phi);
   
     Real val_out = computeAdvection(a_species, 0, a_velocity);
     p_res[0] = val_out;
     if (procID()==0) {
       //cout<<"advection scheme:  "<<m_advect_scheme<<endl;
       //cout<<"flux_0: "<<val_out<<endl;
     }
     val_out = computeAdvection(a_species, 1, a_velocity);
     p_res[1] = val_out;
     if (procID()==0) {
       //cout<<"advection scheme:  "<<m_advect_scheme<<endl;
       //cout<<"flux_1: "<<val_out<<endl;
     }
  }
  */
  //if(a_species_index==1) {
  //if (procID()==0) {cout<<"Electrons are done.\n";}
  //} 
}

//----------------------------------------------------------------------------------------------------------

void
InsulatingSheathBC_new::parseParameters(const ParmParse& a_pp)
{
   a_pp.query("compute_potential_BC", m_compute_potential_BC);
   a_pp.query("sheath_bc_type", m_sheath_bc_type);
   a_pp.query("advection_scheme", m_advect_scheme);
   a_pp.query("iter_number", m_iter_number);
   a_pp.query("tolerance", m_tolerance);
   a_pp.query("debug", m_debug);
   a_pp.query("newton_solver", m_newton);

}

//----------------------------------------------------------------------------------------------------------

Real
InsulatingSheathBC_new::computeAdvection(const KineticSpeciesPtrVect& a_species,
                                         const int& a_species_index,
                                         const LevelData<FluxBox>& a_velocity) const
{
//const PhaseGeom& geometry( a_species.phaseSpaceGeometry() );



   const int& dir( m_bdry_layout->dir() );
   const Side::LoHiSide& side( m_bdry_layout->side() );
   
   //if (procID()==0)  {cout<<"dir = "<<dir<<"   "<<"side = "<<sign(side)<<endl;}
   
   KineticSpecies& species_physical( *(a_species[a_species_index]) );
   LevelData<FArrayBox>& soln( species_physical.distributionFunction() );
   const PhaseGeom& geometry( species_physical.phaseSpaceGeometry() );
   const double mass = species_physical.mass();
   const double charge = species_physical.charge();
   //const double factor = 1.0/mass;
   
      // Testing velocity
//   CFG::IntVect iv_plot;
//   iv_plot[0] = 0;
//   iv_plot[1] = 1;
//   geometry.plotAtConfigurationIndex("velocty_data", iv_plot, a_velocity, 0.0);
//   MayDay::Error( "DONE" );
   // End testing 
   
   /* DOES NO work either
   const double mass = species_physical.mass();
   const double charge = species_physical.charge();

   DisjointBoxLayout tmp_layout;
   tmp_layout.deepCopy(m_grids_full);
    for (int ivec p= 0; ivec < m_grids_full.rawPtr()->size(); ++ivec)
    { Box& base_box = *const_cast<Box*>(&((*m_grids_full.rawPtr())[ivec].box));
       base_box.shift(dir,-sign(side));
      }
   tmp_layout.closeNoSort();
   */

   DisjointBoxLayout tmp_layout;
   if (side == Side::Lo) {adjCellHi(tmp_layout, m_grids_full, dir, 1);}
   if (side == Side::Hi) {adjCellLo(tmp_layout, m_grids_full, dir, 1);}// Somehow it adds a cell to the direction of the interior of the domain
   // addCellHi add +1 to top, addCellLo adds +1 to bottom

   
  /*
  HAND BUILT BOXLAYOUT DOES NOT WORK DUE TO M_INDEX
   Vector<Box> v_box = m_grids_full.boxArray();
   Vector<int> v_int = m_grids_full.procIDs();
   //for (auto & element : v_box) {
     //// Rely on the fact there is only one side used at moment
     //element.shift(dir,-sign(side));
     //}
   //for (auto it = begin (v_box); it != end (v_box); ++it) {
   //  it->shift(dir,-sign(side));
   //}
   for(int i=0; i < v_box.size(); i++){
     v_box[i].shift(dir,-sign(side));;
   }
   DisjointBoxLayout tmp_layout;
   tmp_layout.define(v_box, v_int);
   */
   
   
   // Start tmp test of layout
   /*
   DataIterator dit1 = m_grids_full.dataIterator();
   for (dit1.begin(); dit1.ok(); ++dit1) {
     const Box& this_box = m_grids_full[dit1];
     std::string str_out;
     int ind = dit1.i().intCode();
     string f_name = "output" + std::to_string(ind) + ".txt";
     ofstream fout(f_name);
     //IntVect iv_lo = this_box.loVect();
     IntVect iv_lo = this_box.smallEnd();
     IntVect iv_hi = this_box.bigEnd();
     //IntVect iv_hi = this_box.hiVect();
     int pID = m_grids_full.procID(dit1.i());
     fout<<"box_m_grids   "<<pID<<std::endl;
     str_out = "lo: "+std::to_string(iv_lo[0]) + " " + std::to_string(iv_lo[1]) + " " + std::to_string(iv_lo[2]) + " " + std::to_string(iv_lo[3]);
     fout<<str_out<<std::endl;
     str_out = "hi: "+std::to_string(iv_hi[0]) + " " + std::to_string(iv_hi[1]) + " " + std::to_string(iv_hi[2]) + " " + std::to_string(iv_hi[3]);
     fout<<str_out<<std::endl;
     fout<<std::endl;
     fout.close();
   } 
    
   DataIterator dit2 = tmp_layout.dataIterator();
   for (dit2.begin(); dit2.ok(); ++dit2) {
     const Box& this_box = tmp_layout[dit2];
     std::string str_out;
     int ind = dit2.i().intCode();
     string f_name = "output" + std::to_string(ind) + ".txt";
     ofstream fout(f_name, std::ofstream::out | std::ofstream::app);
     //IntVect iv_lo = this_box.loVect();
     IntVect iv_lo = this_box.smallEnd();
     IntVect iv_hi = this_box.bigEnd();
     //IntVect iv_hi = this_box.hiVect();
     int pID = tmp_layout.procID(dit2.i());
     fout<<"box_tmp   "<<pID<<std::endl;
     str_out = "lo: "+std::to_string(iv_lo[0]) + " " + std::to_string(iv_lo[1]) + " " + std::to_string(iv_lo[2]) + " " + std::to_string(iv_lo[3]);
     fout<<str_out<<std::endl;
     str_out = "hi: "+std::to_string(iv_hi[0]) + " " + std::to_string(iv_hi[1]) + " " + std::to_string(iv_hi[2]) + " " + std::to_string(iv_hi[3]);
     fout<<str_out<<std::endl;
     fout<<std::endl;
     fout.close();
   }
   
   MayDay::Error( "Test done." );
   */
   // Finish test
   
   //LevelData<FArrayBox> face_dfn(m_grids_full, 1, IntVect::Zero);
   LevelData<FArrayBox> face_dfn(tmp_layout, 1, IntVect::Zero);
   //DataIterator dit = face_dfn.dataIterator();      
   DataIterator dit = m_grids_full.dataIterator();   
   for (dit.begin(); dit.ok(); ++dit) {
     const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
     //const DataIndex& interior_dit(dit.i());
     //const DataIndex& interior_dit (DataIndex(dit));
     const FArrayBox& this_dfn( soln[interior_dit] );
     const FluxBox& this_vel = a_velocity[interior_dit];
     const Box& this_box = this_dfn.box();
     FArrayBox& this_face_dfn_dir = face_dfn[dit];
     
     //FArrayBox& this_face_phi_dir( this_dfn[dir] );
      const FArrayBox& this_norm_vel_dir( this_vel[dir] );
            
      
  const Box& vel_box = this_norm_vel_dir.box();
  Box res_box(this_face_dfn_dir.box());
      
/*
            /// TMP output to check boxes
              std::string str_out;
              //int ind = interior_dit.intCode();
              int ind = dit.i().intCode();
  
              string f_name = "output" + std::to_string(ind) + ".txt";
              ofstream fout(f_name);
              //IntVect iv_lo = this_box.loVect();
              IntVect iv_lo = this_box.smallEnd();
              IntVect iv_hi = this_box.bigEnd();
              //IntVect iv_hi = this_box.hiVect();
              fout<<"dfn_box"<<std::endl;
              str_out = "lo: "+std::to_string(iv_lo[0]) + " " + std::to_string(iv_lo[1]) + " " + std::to_string(iv_lo[2]) + " " + std::to_string(iv_lo[3]);
              fout<<str_out<<std::endl;
              str_out = "hi: "+std::to_string(iv_hi[0]) + " " + std::to_string(iv_hi[1]) + " " + std::to_string(iv_hi[2]) + " " + std::to_string(iv_hi[3]);
              fout<<str_out<<std::endl;
              fout<<std::endl;
              
              const Box& vel_box = this_norm_vel_dir.box();

              //iv_lo = vel_box.loVect();
              //iv_hi = vel_box.hiVect();
              iv_lo = vel_box.smallEnd();
              iv_hi = vel_box.bigEnd();
              fout<<"vel_box"<<std::endl;
              str_out = "lo: "+std::to_string(iv_lo[0]) + " " + std::to_string(iv_lo[1]) + " " + std::to_string(iv_lo[2]) + " " + std::to_string(iv_lo[3]);
              fout<<str_out<<std::endl;
              str_out = "hi: "+std::to_string(iv_hi[0]) + " " + std::to_string(iv_hi[1]) + " " + std::to_string(iv_hi[2]) + " " + std::to_string(iv_hi[3]);
              fout<<str_out<<std::endl;
              fout<<std::endl;
             
	      Box res_box(this_face_dfn_dir.box());
              //iv_lo = res_box.loVect();
              //iv_hi = res_box.hiVect();
              iv_lo = res_box.smallEnd();
              iv_hi = res_box.bigEnd();
              
              fout<<"res_box"<<std::endl;
              str_out = "lo: "+std::to_string(iv_lo[0]) + " " + std::to_string(iv_lo[1]) + " " + std::to_string(iv_lo[2]) + " " + std::to_string(iv_lo[3]);
              fout<<str_out<<std::endl;
              str_out = "hi: "+std::to_string(iv_hi[0]) + " " + std::to_string(iv_hi[1]) + " " + std::to_string(iv_hi[2]) + " " + std::to_string(iv_hi[3]);
              fout<<str_out<<std::endl;
              fout<<std::endl;                            
              */ 
            
            
      //res_box.surroundingNodes( dir );     
      // Continue
      //res_box.shift(dir, -sign(side));
          if (m_advect_scheme=="uw1") {
      for (BoxIterator bit(res_box); bit.ok(); ++bit) {
        IntVect iv = bit();
        IntVect iv_v = iv;
        iv_v[dir] += (sign(side)+1)/2;
        Real this_v_norm = this_norm_vel_dir(iv_v,dir);
        // Outgoing flux
        if ( sign(side)*this_v_norm>0 ) { this_face_dfn_dir(iv,0) = this_dfn(iv,0); }
        else {
          IntVect iv_p = iv;
          iv_p[dir] += sign(side);
          this_face_dfn_dir(iv,0) = this_dfn(iv_p,0);
        }
      }
    }
    
    else if (m_advect_scheme=="uw3"){
          for (BoxIterator bit(res_box); bit.ok(); ++bit) {
            IntVect iv = bit();
            IntVect iv_v = iv;
            iv_v[dir] += (sign(side)+1)/2;
            Real this_v_norm = this_norm_vel_dir(iv_v,dir);
            if(sign(side)*this_v_norm>0)//outgoing flux
            {
                IntVect iv_p = iv;
                IntVect iv_m = iv;
                iv_p[dir] += sign(side);
                iv_m[dir] -= sign(side);
                this_face_dfn_dir(iv,0) = (2*this_dfn(iv_p,0) + 5*this_dfn(iv,0) - this_dfn(iv_m,0)) / 6.0;
            }
            else{
                IntVect iv_p = iv;
                IntVect iv_pp = iv;
                iv_p[dir] += sign(side);
                iv_pp[dir] += 2*sign(side);
                this_face_dfn_dir(iv,0) = (2*this_dfn(iv,0) + 5*this_dfn(iv_p,0) - this_dfn(iv_pp,0)) / 6.0;
            }
          }
        }
        else if (m_advect_scheme=="weno"){
          for (BoxIterator bit(res_box); bit.ok(); ++bit) {
            IntVect iv = bit();
            IntVect iv_v = iv;
            iv_v[dir] += (sign(side)+1)/2;
            Real this_v_norm = this_norm_vel_dir(iv_v,dir);
            if(sign(side)*this_v_norm>0)//outgoing flux
            {
                double ptr[5];
                IntVect iv_p = iv;
                IntVect iv_m = iv;
                iv_p[dir] += sign(side);
                iv_m[dir] -= sign(side);
                ptr[2] = this_dfn(iv,0);
                ptr[1] = this_dfn(iv_m,0);
                ptr[3] = this_dfn(iv_p,0);
                iv_p[dir] += sign(side);
                iv_m[dir] -= sign(side);
                ptr[0] = this_dfn(iv_m,0);
                ptr[4] = this_dfn(iv_p,0);
                this_face_dfn_dir(iv,0) = computeFaceValueWENO(ptr);
            }
            else{
              double ptr[5];
                IntVect iv_p = iv;
                IntVect iv_pp = iv;
                iv_p[dir] += sign(side);
                iv_pp[dir] += 2*sign(side);
                ptr[2] = this_dfn(iv_p,0);
                ptr[1] = this_dfn(iv_pp,0);
                ptr[3] = this_dfn(iv,0);
                iv_pp[dir] += sign(side);
                iv_p[dir] -= 2*sign(side);
                ptr[0] = this_dfn(iv_pp,0);
                ptr[4] = this_dfn(iv_p,0);
                this_face_dfn_dir(iv,0) = computeFaceValueWENO(ptr);
                //(2*this_dfn(iv,0) + 5*this_dfn(iv_p,0) - this_dfn(iv_pp,0)) / 6.0;
            }
          }
        }
      else
      if(0 && m_advect_scheme=="uw1") {
         FORT_UW1FACE( CHF_FRA( this_face_dfn_dir ),
                       CHF_CONST_FRA( this_dfn ),
                       CHF_CONST_FRA1( this_norm_vel_dir, dir ),
                       CHF_BOX( res_box ),
                       CHF_CONST_INT( dir ) );
      } else
      if(0 && m_advect_scheme=="uw3") {
         FORT_UW3FACE( CHF_FRA( this_face_dfn_dir ),
                       CHF_CONST_FRA( this_dfn ),
                       CHF_CONST_FRA1( this_norm_vel_dir, dir ),
                       CHF_BOX( res_box ),
                       CHF_CONST_INT( dir ) );
	 /*} else
     if(0 && m_advect_scheme=="uw3f") {
         FORT_UW3FFACE( CHF_FRA( this_face_dfn_dir ),
                       CHF_CONST_FRA( this_dfn ),
                       CHF_CONST_FRA1( this_norm_vel_dir, dir ),
                       CHF_BOX( res_box ),
                       CHF_CONST_INT( dir ) );
	 */
      } else
     if(0 && m_advect_scheme=="uw5") {
         FORT_UW5FACE( CHF_FRA( this_face_dfn_dir ),
                       CHF_CONST_FRA( this_dfn ),
                       CHF_CONST_FRA1( this_norm_vel_dir, dir ),
                       CHF_BOX( res_box ),
                       CHF_CONST_INT( dir ) );

     } else
     if(0 && m_advect_scheme=="bweno") {
         FORT_BWENOFACE( CHF_FRA( this_face_dfn_dir ),
                       CHF_CONST_FRA( this_dfn ),
                       CHF_CONST_FRA1( this_norm_vel_dir, dir ),
                       CHF_BOX( res_box ),
                       CHF_CONST_INT( dir ) );

     }
     else {MayDay::Error( "No advection scheme specified in the BC class!" );}
     //this_face_dfn_dir.mult(this_norm_vel_dir);
     
     //res_box.shift(dir, sign(side));
    
    // Debugging again
    /*for (BoxIterator bit(res_box); bit.ok(); ++bit) {
           IntVect iv = bit();
           IntVect iv_v = iv;
           IntVect iv_m = iv;
           iv_v[dir] -= sign(side);
           iv_m[dir] += sign(side);
           str_out = "iv: "+std::to_string(iv[0]) + " " + std::to_string(iv[1]) + " " + std::to_string(iv[2]) + " " + std::to_string(iv[3])+":   ";
           str_out += std::to_string(this_dfn(iv_v,0)) + "    " + std::to_string(this_dfn(iv,0)) + "    " + std::to_string(this_dfn(iv_m,0)) + "    " + std::to_string(this_face_dfn_dir(iv,0)) + "    "+ std::to_string(this_norm_vel_dir(iv,dir));
           fout<<str_out<<std::endl;
    }
    
    fout.close();
    // Finished debuging
    */
    Real tmp;
    for (BoxIterator bit(res_box); bit.ok(); ++bit) {
           IntVect iv = bit();
                 IntVect iv_v = iv;
                  // Shift to one cell in y space. Non geometry invariant, must be fixed!!!
                 iv_v[dir] += (sign(side)+1)/2; // This shift is for the big end of the velbox that is +1 due to flux box
           tmp = this_norm_vel_dir(iv,dir);
           tmp = this_face_dfn_dir(iv,0);
           this_face_dfn_dir(iv,0) *= this_norm_vel_dir(iv_v,dir);
    }
    
   }
   Real res = simple_sum( face_dfn );
   return charge*res/mass/mass;
}

void
InsulatingSheathBC_new::computeHeatFluxes( Real * p_res,
                                           const KineticSpeciesPtrVect& a_species,
                                           const int& a_species_index,
                                           const LevelData<FluxBox>& a_velocity) const
{
  // Get boundary layout data
  const int& dir( m_bdry_layout->dir() );
  const Side::LoHiSide& side( m_bdry_layout->side() );

  // Get kinetic species data
  KineticSpecies& species_physical( *(a_species[a_species_index]) );
  LevelData<FArrayBox>& soln( species_physical.distributionFunction() );
  const PhaseGeom& geometry( species_physical.phaseSpaceGeometry() );
  const double mass = species_physical.mass();
  const double charge = species_physical.charge();
  //const double factor = 1.0/mass;
  const VEL::VelCoordSys& vel_coords = geometry.velSpaceCoordSys();
  const LevelData<FArrayBox>& B_injected( geometry.getBFieldMagnitude() );

  // Create an auxiliary layout that actually takes valid cells in the computational domain
  DisjointBoxLayout tmp_layout;
  if (side == Side::Lo) {adjCellHi(tmp_layout, m_grids_full, dir, 1);}
  if (side == Side::Hi) {adjCellLo(tmp_layout, m_grids_full, dir, 1);}

  // This is LevelData on the boundary
  LevelData<FArrayBox> face_dfn(tmp_layout, 1, IntVect::Zero); // This LevelData is for the flux
  LevelData<FArrayBox> face_dfn_mv2(tmp_layout, 1, IntVect::Zero); // This LevelData is for the mv^2/2 flux
  LevelData<FArrayBox> face_dfn_muB(tmp_layout, 1, IntVect::Zero); // This LevelData is for the mu B flux
  DataIterator dit = m_grids_full.dataIterator();   
  for (dit.begin(); dit.ok(); ++dit) {
    // We need interior dit, since solution and a_velocity are located in the computational domain, where boundary laypit dit is invalid
    const DataIndex& interior_dit( m_bdry_layout->dataIndex(dit) );
    const FArrayBox& this_dfn( soln[interior_dit] );
    const FluxBox& this_vel = a_velocity[interior_dit]; // this velocity is mapped but we will need it for the projection
    const Box& this_box = this_dfn.box();
    FArrayBox& this_face_dfn_dir = face_dfn[dit];


    const FArrayBox& this_norm_vel_dir( this_vel[dir] );
    const Box& vel_box = this_norm_vel_dir.box();
    Box res_box(this_face_dfn_dir.box());
    // Ideally, it should be done using FORT_UW1FACE and FORT_UW3FACE functions, but we update it later?????
    if (m_advect_scheme=="uw1") {
      for (BoxIterator bit(res_box); bit.ok(); ++bit) {
        IntVect iv = bit();
        IntVect iv_v = iv;
        iv_v[dir] += (sign(side)+1)/2;
        Real this_v_norm = this_norm_vel_dir(iv_v,dir);
        // Outgoing flux
        if ( sign(side)*this_v_norm>0 ) { this_face_dfn_dir(iv,0) = this_dfn(iv,0); }
        else {
          IntVect iv_p = iv;
          iv_p[dir] += sign(side);
          this_face_dfn_dir(iv,0) = this_dfn(iv_p,0);
        }
      }
    }
    else if (m_advect_scheme=="uw3") {
      for (BoxIterator bit(res_box); bit.ok(); ++bit) {
        IntVect iv = bit();
        IntVect iv_v = iv;
        iv_v[dir] += (sign(side)+1)/2;
        Real this_v_norm = this_norm_vel_dir(iv_v,dir);
        // Explicitly using upwind3 stencil u_face = (2*u[+1] + 5*u[0] - u[-1])/6, where u is cell-AVERAGED value of the PDF
        if(sign(side)*this_v_norm>0) //outgoing flux
        {
          IntVect iv_p = iv;
          IntVect iv_m = iv;
          iv_p[dir] += sign(side);
          iv_m[dir] -= sign(side);
          this_face_dfn_dir(iv,0) = (2*this_dfn(iv_p,0) + 5*this_dfn(iv,0) - this_dfn(iv_m,0)) / 6.0;
        }
        else{
          IntVect iv_p = iv;
          IntVect iv_pp = iv;
          iv_p[dir] += sign(side);
          iv_pp[dir] += 2*sign(side);
          this_face_dfn_dir(iv,0) = (2*this_dfn(iv,0) + 5*this_dfn(iv_p,0) - this_dfn(iv_pp,0)) / 6.0;
        }
      }
    }
            else if (m_advect_scheme=="weno"){
          for (BoxIterator bit(res_box); bit.ok(); ++bit) {
            IntVect iv = bit();
            IntVect iv_v = iv;
            iv_v[dir] += (sign(side)+1)/2;
            Real this_v_norm = this_norm_vel_dir(iv_v,dir);
            if(sign(side)*this_v_norm>0)//outgoing flux
            {
                double ptr[5];
                IntVect iv_p = iv;
                IntVect iv_m = iv;
                iv_p[dir] += sign(side);
                iv_m[dir] -= sign(side);
                ptr[2] = this_dfn(iv,0);
                ptr[1] = this_dfn(iv_m,0);
                ptr[3] = this_dfn(iv_p,0);
                iv_p[dir] += sign(side);
                iv_m[dir] -= sign(side);
                ptr[0] = this_dfn(iv_m,0);
                ptr[4] = this_dfn(iv_p,0);
                this_face_dfn_dir(iv,0) = computeFaceValueWENO(ptr);
            }
            else{
              double ptr[5];
                IntVect iv_p = iv;
                IntVect iv_pp = iv;
                iv_p[dir] += sign(side);
                iv_pp[dir] += 2*sign(side);
                ptr[2] = this_dfn(iv_p,0);
                ptr[1] = this_dfn(iv_pp,0);
                ptr[3] = this_dfn(iv,0);
                iv_pp[dir] += sign(side);
                iv_p[dir] -= 2*sign(side);
                ptr[0] = this_dfn(iv_pp,0);
                ptr[4] = this_dfn(iv_p,0);
                this_face_dfn_dir(iv,0) = computeFaceValueWENO(ptr);
                //(2*this_dfn(iv,0) + 5*this_dfn(iv_p,0) - this_dfn(iv_pp,0)) / 6.0;
            }
          }
        }
      else
         if(0 && m_advect_scheme=="uw1") {
      FORT_UW1FACE( CHF_FRA( this_face_dfn_dir ),
                    CHF_CONST_FRA( this_dfn ),
                    CHF_CONST_FRA1( this_norm_vel_dir, dir ),
                    CHF_BOX( res_box ),
                    CHF_CONST_INT( dir ) );
    } else
    if(0 && m_advect_scheme=="uw3") {
      FORT_UW3FACE( CHF_FRA( this_face_dfn_dir ),
                    CHF_CONST_FRA( this_dfn ),
                    CHF_CONST_FRA1( this_norm_vel_dir, dir ),
                    CHF_BOX( res_box ),
                    CHF_CONST_INT( dir ) );
    }
    else {MayDay::Error( "No advection scheme specified in the BC class!" );}

    // Get the velocity
    FArrayBox velocityRealCoords(res_box, VEL_DIM);
    const PhaseBlockCoordSys& block_coord_sys = geometry.getBlockCoordSys(res_box);
    block_coord_sys.getVelocityRealCoords(velocityRealCoords);
      
    
    Real B_field = 3.0; // HARDWIRED HERE, we need to get actual B
    const FArrayBox& this_Bfield( B_injected[interior_dit] );
    const Box& b_box = this_Bfield.box();
    IntVect s_end =  b_box.smallEnd();
    
    Real v_par, mu;
    for (BoxIterator bit(res_box); bit.ok(); ++bit) {
      IntVect iv = bit();
      IntVect iv_v = iv;
      IntVect iv_b = iv;
      // Make v_par and mu components legit
      iv_b[CFG_DIM] = s_end[CFG_DIM];
      iv_b[CFG_DIM+1] = s_end[CFG_DIM+1];
      // Shift to one cell in y space. Non geometry invariant, must be fixed!!!
      iv_v[dir] += (sign(side)+1)/2; // This shift is for the big end of the velbox that is +1 due to flux box
      // it should be iv_v[1] -= sign(side);
      // Assuming that cell centered values of v_=parr and mu are the same that on the boundary
      v_par = velocityRealCoords(iv,0);
      mu = velocityRealCoords(iv,1);
      this_face_dfn_dir(iv,0) *= this_norm_vel_dir(iv_v,dir);
      face_dfn_mv2[dit](iv,0) = this_face_dfn_dir(iv,0)*mass*v_par*v_par/2;
      face_dfn_muB[dit](iv,0) = this_face_dfn_dir(iv,0)*mu*this_Bfield(iv_b,0)/2;
    }
    

  }
  Real factor_norm = 1.0/mass/mass/4.0; // Here, 4 = Nx, since simple_sum will take a sum over all indexes and we have Nx=4 in the input file
  p_res[0] = factor_norm*simple_sum( face_dfn, false);
  p_res[1] = factor_norm*simple_sum( face_dfn_mv2, false );
  p_res[2] = factor_norm*simple_sum( face_dfn_muB, false );

}


#include "NamespaceFooter.H"


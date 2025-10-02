#include "ImportData.H"

#include "MagGeom.H"
#include "LogRectCoordSys.H"
#include "SNCoreCoordSys.H"
#include "SingleNullCoordSys.H"

#include "newMappedGridIO.H"

#include "NamespaceHeader.H"

ImportData::ImportData(ParmParse& a_pp,
                       const int& a_verbosity )
  : GridFunction(a_verbosity),
    m_scale_factor(1.0),
    m_floor(0.)
{
   parseParameters( a_pp );
}


inline
void ImportData::parseParameters( ParmParse& a_pp )
{
   if ( a_pp.contains("data_file")) {
        a_pp.get("data_file", m_data_file);
   }
   else {
       MayDay::Error("OneDimData:: No data file specified");
   }

   a_pp.query("scale_factor", m_scale_factor);
   a_pp.query("floor", m_floor); 

   if (m_verbosity) {
      printParameters();
   }
}

void ImportData::assign(LevelData<FArrayBox>&       a_data,
                        const MultiBlockLevelGeom&  a_geometry,
                        const Real&                 a_time,
                        const bool&                 a_cell_averages ) const
{
   checkGeometryValidity( a_geometry );
   
   //Read-in hdf5 data file
   
   Vector<DisjointBoxLayout> vectGrids;
   Vector<LevelData<FArrayBox>* >  vectData;
   Vector<string> vectNames;
   Box domain;
   RealVect dx;
   Real dt;
   Real time;
   Vector<IntVect> refRatio;
   int numLevels;
   
   ReadAnisotropicAMRHierarchyHDF5(m_data_file,
                                   vectGrids,
                                   vectData,
                                   vectNames,
                                   domain,
                                   dx,
                                   dt,
                                   time,
                                   refRatio,
                                   numLevels);
   
   LevelData<FArrayBox>& imported_data = *vectData[0];

   imported_data.copyTo(a_data);

   const DisjointBoxLayout& grids = a_data.disjointBoxLayout();
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
     a_data[dit].mult(m_scale_factor);
   }
   
   
   //Handle ghost cells
   if (a_data.ghostVect() >= 2*IntVect::Unit) {
      ((MagGeom&)a_geometry).extrapolateToPhysicalGhosts(a_data, true);
   }

   if (a_data.ghostVect() == IntVect::Unit) {
      ((MagGeom&)a_geometry).extrapolateToPhysicalGhosts(a_data, false);
   }

   if (a_data.ghostVect() >= 3*IntVect::Unit) {
      if ( procID()==0 ) MayDay::Warning("ImportData: the object has more than two layers of ghost cells. Only two layers are filled by extrapolation");
   }

   //Imposing floor
   for (DataIterator dit( grids.dataIterator() ); dit.ok(); ++dit) {
     BoxIterator bit(a_data[dit].box());
     for (bit.begin(); bit.ok(); ++bit) {   
       IntVect iv = bit();
       for (int n_comp=0; n_comp<a_data.nComp(); ++n_comp) {
	 if (abs(a_data[dit](iv,n_comp)) < m_floor) {
	   a_data[dit](iv,n_comp) = m_floor;
	 }
       }
     }
   }
}

void ImportData::checkGeometryValidity( const MultiBlockLevelGeom& a_geometry ) const
{
   const MultiBlockCoordSys& coord_sys( *(a_geometry.coordSysPtr()) );
   bool unknown_geometry( typeid(coord_sys) != typeid(LogRectCoordSys) );
   unknown_geometry &= (typeid(coord_sys) != typeid(SingleNullCoordSys));
   unknown_geometry &= (typeid(coord_sys) != typeid(SNCoreCoordSys));
   
   if ( unknown_geometry ) {
      const std::string msg( "ImportData: Attempt to use unknown geometry. ");
      MayDay::Error( msg.c_str() );
   }
}

void ImportData::setPointwise(FArrayBox&                 a_data,
			      const MultiBlockLevelGeom& a_geometry,
			      const FArrayBox&           a_real_coords,
			      const FArrayBox&           a_normalized_flux,
			      const int                  a_block_number) const
{
  const std::string msg( "ImportData::setPointwise: not implemented !!! ImportData-type functions cannot be used for BCs ");
  MayDay::Error( msg.c_str() );
}


void ImportData::printParameters() const
{
   if (procID()==0) {
      std::cout << "ImportData grid function parameters:" << std::endl;
      std::cout << "  data file name: "   << m_data_file   << std::endl;
      std::cout << std::endl;
   }
}


#include "NamespaceFooter.H"

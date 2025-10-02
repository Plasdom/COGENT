#ifndef _IMPORTDATA_H_
#define _IMPORTDATA_H_

#include "GridFunction.H"
#include "NamespaceHeader.H"

/**
 * ImportData function class.
 *
 * Presently, assigns data from and hdf5 file.
 * Assumes that the same DBL was used to create the hdf5 file
 * and for the LevelData<FArrayBox> object we import data into.
 * Fills the ghost cells by extrapolation.
 *
 * The following represents a sample input entry for this function choice.
 *
 * \verbatim
 * type = "ImportData"
 * data_file = "path_to_file"
 *
*/
class ImportData : public GridFunction
{
   public:

      /// Constructor using a database.
      /**
       * Constructs the object by getting parameter from a ParmParse database.
       *
       * @param[in] pp ParmParse database that contains the parameters for
       *               the new GridFunction object.
       * @param[in] verbosity Flag to turn on additional levels of output
       */
      ImportData( ParmParse& pp, const int& verbosity );

      /// Destructor.
      /**
       */
      virtual ~ImportData() {;}

      /// Print object parameters.
      /**
       */
      virtual void printParameters() const;

      virtual void assign(LevelData<FArrayBox>& data,
                  const MultiBlockLevelGeom& geometry,
                  const Real& time,
                  const bool& as_cell_averages = true ) const;

   private:

      // prohibit copying
      ImportData( const ImportData& );
      ImportData& operator=( const ImportData& );

   
   
      /// Set pointwise values at cell centers.
      /**
       */
      virtual void setPointwise(FArrayBox&                  dfn,
                                const MultiBlockLevelGeom&  geometry,
                                const FArrayBox&            real_coords,
                                const FArrayBox&            normalized_flux,
                                const int                   block_number ) const;
   
      /// Parse the input database for parameters.
      /**
       */
      inline void parseParameters( ParmParse& pp );
   
      /// Check geometry validity.
      /**
       */
      virtual void checkGeometryValidity( const MultiBlockLevelGeom& geometry ) const;

      std::string m_data_file;

      Real m_scale_factor;
      Real m_floor;
};

#include "NamespaceFooter.H"

#endif

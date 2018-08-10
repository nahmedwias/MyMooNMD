/** ************************************************************************ 
*
* @name       POD_TCDR2D
* @brief      POD for time-dep. convection-diffusion scalar problems in 2D
*
* @author     Swetlana Giere
* @date       08.03.2017 (start of implementation)
*
****************************************************************************/

#ifndef __POD_TCDR2D__
#define __POD_TCDR2D__

// boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <MooNMD_Io.h>
#include <Matrix.h>
#include <POD.h>
#include <BlockFEMatrix.h>
#include <FEFunction2D.h>
#include <BlockVector.h>
#include <PostProcessing2D.h>
#include <ParameterDatabase.h>
#include <Example_TimeCD2D.h>


#include <fstream>
#include <memory>
#include <vector>

using namespace std;
namespace ublas = boost::numeric::ublas;


class POD_TCDR2D : public POD
{
  protected:
	/** @brief Definition of the used example */
	const Example_TimeCD2D example;
	/** @brief Finite Element space */
	TFESpace2D fe_space;
	/** @brief gramian matrix (needed for POD computation) */
	BlockFEMatrix gramian_matrix;
	BlockVector pod_mode;
	/** @brief Finite Element function */
	TFEFunction2D fe_function;

	/** @brief class for output handling */
	PostProcessing2D outputWriter;

	/** @brief set parameters in database
	*
	* This functions checks if the parameters in the database are meaningful
	* and resets them otherwise. The hope is that after calling this function
	* this class is fully functional.
	*
	* If some parameters are set to unsupported values, an error occurs and
	* throws an exception.
	*/
	void set_parameters();

	/** @brief Assemble gramian matrix
	 *
	 * Assemble gramian matrix, which describes the inner product, with respect
	 * to which the POD basis will be computed.
	*/
	void assemble_gramian();

  public:
      /** @brief constructor
       *
       */
      POD_TCDR2D(TCollection& coll, const ParameterDatabase& param_db,
      		const Example_TimeCD2D& ex);

      /** @brief default destructor
      */
      ~POD_TCDR2D();

      /**
      * @brief Compute POD basis
      *
      * Compute POD basis from snapshots with the compute_basis() routine
      * from the class POD. Within the function, appropriate gramian matrix
      * according to the 'pod_inner_product' parameter is assembled and
      * incorporated into computation of the POD basis.
      */
      void compute_pod_basis();

      /**
      * @brief Read POD basis from file
      *
      * Read POD basis from file. The functionality could be interesting if POD
      * basis is already computed and stored into file and one
      * only wants to write vtk files of the POD basis functions.
      */
      void read_basis();

      /// @brief write POD basis to file and write vtk files
      void output();

       // getters and setters
      const Example_TimeCD2D& get_example() const
      { return example; }
      const TFEFunction2D & get_function() const
      { return fe_function; }
      TFEFunction2D & get_function()
      { return fe_function; }
      const TFESpace2D & get_space() const
      { return fe_space; }
};

#endif

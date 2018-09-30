/** ************************************************************************
*
* @name       ROM_TCDR2D
* @brief      All about ROM for time-dep. convection-diffusion scalar problems
* 
* NOTE: Current implementation is only valid for problems with time-independent
* coefficients (diffusion, convection, reaction) and source term. In case of
* time-dependent coefficients, corresponding matrices have to be assembled in every
* time step, which is very inefficient. Alternatively, if coefficients can be
* splitted into temporal and space parts, then the corresponding matrices involving
* only the space part can be precomputed offline and within the time loop they only
* must be updated with the corresponding temporal part of the coefficients.
*
* @author     Swetlana Giere
* @date       15.03.2017 (start of implementation)
*
****************************************************************************/


//TODO Implement regularized ROM initial condition


#ifndef __ROM_TCDR2D__
#define __ROM_TCDR2D__

// boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <MooNMD_Io.h>
#include <Matrix.h>
#include <ROM.h>
#include <BlockFEMatrix.h>
#include <FEFunction2D.h>
#include <BlockVector.h>
#include <DataWriter.h>
#include <ParameterDatabase.h>
#include <Example_TimeCD2D.h>
#include <TCD2D_POD_LocalAssemble.h>


#include <fstream>
#include <memory>
#include <vector>

using namespace std;
namespace ublas = boost::numeric::ublas;

class ROM_TCDR2D : public ROM
{
  private:

	/** @brief Definition of the used example */
	const Example_TimeCD2D example;
	/** @brief Finite Element space */
	TFESpace2D fe_space;
	/** @brief Gramian matrix (needed for reduction of solution) */
	BlockFEMatrix gramian_matrix;
	/** @brief full-order solution */
	BlockVector full_sol;
	/** @brief Finite Element function */
	TFEFunction2D fe_function;
	/** @brief Class for output handling */
	DataWriter2D outputWriter;
	/** @brief Store errors to compute accumulated error norms */
	std::vector<double> errors;

	/** @brief Reduced system matrix */
	ublas::matrix<double> sys_mat;
	/** @brief Reduced system rhs */
	ublas::vector<double> sys_rhs;
	/** @brief ROM solution */
    ublas::vector<double> red_sol;

    /** @brief Reduced mass matrix */
    ublas::matrix<double> mass_mat;
    /** @brief Reduced convection-diffusion-reaction matrix */
    ublas::matrix<double> cdr_mat;
    /** @brief Reduced convection-diffusion-reaction matrix times snaps mean */
    ublas::vector<double> cdr_mat_mean;
    /** @brief Reduced source term (of its variational formulation)*/
    ublas::vector<double> source;

	/** @brief assemble FE matrix by providing a local assembling routine 'local_assemble_param'
	*/
	void assemble(BlockFEMatrix &mat, AssembleFctParam *local_assemble_param);

	/** @brief assemble FE rhs by providing a local assembling routine 'local_assemble_param'
	*/
	void assemble(BlockVector &rhs, AssembleFctParam *local_assemble_param);

	/**
	 * @brief Assemble gramian matrix for reduction of finite element solution
	 *
	 * Assemble gramian matrix which will be used for the reduction of the
	 * full-order finite element solution. Note that the matrix should be the
	 * same as used to compute the underlying POD basis. Value of parameter
	 * 'pod_inner_product' is used to decide which matrix has to be assembled.
	 */
	void assemble_gramian();

  public:

	/** @brief Constructor
	*/
	ROM_TCDR2D(TCollection& coll, const ParameterDatabase& param_db,
	      	   const Example_TimeCD2D& ex);

	/** @brief Default destructor
	*/
	~ROM_TCDR2D();

	/** @brief Check and set parameters in database
	*
	* This functions checks if the parameters in the database are meaningful
	* and resets them otherwise. The hope is that after calling this function
	* this class is fully functional.
	*
	* If some parameters are set to unsupported values, an error occurs and
	* throws an exception.
	*/
	void set_parameters();

	/** @brief Assemble required full-order matrices/vectors and reduce them
	 *
	 * Assemble mass, cdr(convection-diffusion-reaction), gramian matrices. If
	 * parameter DISCTYPE is set to 2 then the corresponding matrices/vectors with
	 * including SUPG parts are assembled.
	 * Reduce mass and cdr matrices -> mass_mat, cdr_mat, cdr_mat_mean
	 * Assemble source term and reduce it -> source.
	 * If only_rhs is set to true then only the source term is assmelbed and reduced
	 * (coud be useful within the time loop). Otherwise, both matrices and rhs are
	 * processed.
	 *
	 * @param only_rhs Flag indicates if both matrices and rhs have to be processed or just rhs
	*/
    void assemble_reduce(bool only_rhs=false);

    /** @brief Compute ROM initial condition
     *
     * With parameter 'rom_init_regularized' set to false, ROM initial condition
     * is computed in the standard way (L2 projection into POD space). Otherwise,
     * the regularized ROM initial condition (see PhD Thesis of S.Giere, Sec. 3.2.2.) is
     * computed. The filter width in the Helmholtz equation is controlled by database
     * parameter 'differential_filter_width'.
    */
    void compute_initial_solution();

    /** @brief Compute system matrix
     *
     * System matrix sys_mat is computed and its inverse is pre-computed for a
     * more efficient solution of the linear system (only valid for
     * time-independent problem coefficients).
    */
    void set_system_matrix();

    /** @brief Compute system rhs
     *
     * @note Note that assembling in every time step in very inefficient in the
     * ROM context. Try to avoid it, e.g., by storing all FE source terms at all
     * times and reduce them offline and online just use the reduced vectors.
     * Alternatively, for source terms given in the separated time-space form,
     * one could pre-assemble the space part and reduce it offline and within the
     * time loop just multiply it with the corresponding temporal coefficient
     * (e.g., see PhD Thesis of S.Giere, Sec. 4.2.).
     *
     * @param reassemble Flag indicating if the source term has to be reassembled
    */
    void set_system_rhs(bool reassemble=true);

    /** @brief solve the system
    */
    void solve();

    /// @brief measure errors and write solution to vtk
    void output(int time_step);
};

#endif

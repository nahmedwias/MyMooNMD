/** ***********************************************************************************
* @(#)ROM.h   14/03/2017
* 
* Class:      ROM (independent of problem's dimension)
*
* Purpose:    This class provides different routines for basis transformation (between
* 			  the underlying finite element and POD space), solving the ROM system.
*
*
***************************************************************************************/

#ifndef __ROM__
#define __ROM__

// boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>

#include <MooNMD_Io.h>
#include <Matrix.h>
#include <ParameterDatabase.h>
#include <POD.h>
#include <fstream>
#include <memory>


using namespace std;
namespace ublas = boost::numeric::ublas;
typedef ublas::permutation_matrix<std::size_t> pmatrix;
typedef ublas::identity_matrix<double> idmatrix;

class ROM : public POD
{
  private:

	/* inverse matrix for solving linear system
	 * (make sense for time-independent system matrix)
	 */
	ublas::matrix<double> inv_mat;

  public:
	/**
	* @brief default constructor
	*/
	ROM();

	/**
	* @brief constructor
	*/
	ROM(const ParameterDatabase& param_db);


	/**
	* @brief default constructor
	*/
	~ROM();

	/**
	* @brief Reduce full sparse matrix using POD basis
	*
	* red_mat = pod_basis^T * full_mat * pod_basis
	*
	* @param full_mat Sparse ParMooN matrix to reduce
	* @param red_mat  Reduced matrix to be computed
	*/
	void reduce(std::shared_ptr<TMatrix> full_mat, ublas::matrix<double> &red_mat);
	ublas::matrix<double> reduce(std::shared_ptr<TMatrix> full_mat);

	/**
	* @brief Reduce full-order sparse matrix multiplied with snapshots mean using POD basis
	*
	* red_vec = pod_basis^T * full_mat * snaps_mean
	*
	* @param full_mat Sparse ParMooN matrix
	* @param red_vec  Reduced array to be computed
	*
	*/
	void reduce_mat_mean(std::shared_ptr<TMatrix> full_mat, ublas::vector<double> &red_vec);
	ublas::vector<double> reduce_mat_mean(std::shared_ptr<TMatrix> full_mat);
	
	/**
	* @brief Reduce full-order vector using POD basis
	*
	* red_vec = pod_basis^T * full_vec
	*
	* @param full_vec Full-order ublas vector
	* @param red_vec  Reduced-order ublas vector
	*
	*/
	void reduce(const ublas::vector<double> &full_vec, ublas::vector<double> &red_vec);
	ublas::vector<double> reduce(const ublas::vector<double> &full_vec);
	ublas::vector<double> reduce(double *full_vec);

	/**
	* @brief Reduce full-order solution
	*
	* Reduce full-order solution by first subtracting from it the snapshots'
	* mean (rom_db["pod_fluct"] == true) and then projecting it onto POD basis
	*
	* red_sol = pod_basis^T * Gram_mat * (full_sol - mean_snaps),
	*
	* where Gram_mat represents the inner product used for the computation of POD basis.
	* NOTE: For the pointer full_soll enough memory must be allocated.
	*
	* @param full_sol Finite element solution
	* @param red_sol  ROM solution
	*
	*/
	void reduce_solution(double *full_sol, ublas::vector<double> &red_sol);
	
	/**
	* @brief Compute full-order solution from the ROM solution
	*
	* Compute full-order solution form the ROM solution by first
	* multiplying the ROM solution with the POD matrix and then
	* adding to it the snapshots' mean if rom_db["pod_fluct"] == true.
	*
	* full_sol = mean_snaps + pod_basis * red_sol
	*
	* NOTE: For the pointer full_soll enough memory must be allocated.
	*
	* @param red_sol  ROM solution
	* @param full_sol Finite element solution
	*
	*/
	void get_full_solution(const ublas::vector<double> &red_sol, double *full_sol);

	/**
	* @brief Solve reduced (small) linear system using UBLAS routine
	*
	* mat * sol = rhs
	*
	* @param mat Sytem matrix
	* @param rhs System right-hand side
	*
	* @return sol Solution
	*/
	ublas::vector<double> solve(const ublas::matrix<double> &mat, 
			                    const ublas::vector<double> &rhs);

	/**
	 * @brief Precompute inverse matrix of the system matrix
	 * It makes sense only for problems with a time-independent system matrix
	 */
	void prepare_solve(const ublas::matrix<double> &mat);
	/**
     * @brief Solve linear system
     *
     * It can be used only in combination with prepare_solve()
     *
     * @param  rhs System right-hand side
     *
     * @return Solution
     */
   	ublas::vector<double> solve(const ublas::vector<double> &rhs);

	/**
	* @brief Set gramian matrix (needed for reduction of full-order solution)
	*
	*/
	void set_gramian(std::shared_ptr<TMatrix> gram_mat);
};

#endif

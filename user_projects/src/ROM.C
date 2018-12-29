/** ***********************************************************************************
*
* Class:      ROM (independent of problem's dimension)
*
* Purpose:    This class provides different routines for basis transformation (between
* 			  the underlying finite element and POD space), solving the ROM system.
*
*
***************************************************************************************/

#include <ROM.h>
#include <ParameterDatabase.h>


/**************************************************************************** */
ROM::ROM(const ParameterDatabase& param_db)
: POD(param_db)
{
  this->inv_mat.clear();
  POD::read_basis();
}

/**************************************************************************** */
ROM::~ROM()
{
}

/**************************************************************************** */
void ROM::reduce(std::shared_ptr<TMatrix> full_mat, ublas::matrix<double> &red_mat)
{ 
  ublas::compressed_matrix<double> ublas_mat;
  convert_to_ublas(full_mat, ublas_mat);
  
  red_mat.resize(POD::rank,POD::rank);
  
  ublas::matrix<double> tmp_mat (ublas_mat.size1(), POD::rank);
  
  noalias(tmp_mat) = prod(ublas_mat, POD::pod_basis);
  noalias(red_mat) = prod(trans(POD::pod_basis), tmp_mat);
}

/**************************************************************************** */
ublas::matrix<double> ROM::reduce(std::shared_ptr<TMatrix> full_mat)
{ 
  ublas::matrix<double> red_mat; 
  
  reduce(full_mat, red_mat);
  return red_mat;
}

/**************************************************************************** */
void ROM::reduce_mat_mean(std::shared_ptr<TMatrix> full_mat, ublas::vector<double> &red_vec)
{
  if(!POD::rom_db["pod_fluct"])
  {
	ErrThrow("ROM::reduce_mat_mean(..): Function is only available if"
			 "parameter 'pod_fluct' is set to true.");
  }
  
  ublas::compressed_matrix<double> ublas_mat;
  convert_to_ublas(full_mat, ublas_mat);
  
  red_vec.resize(POD::rank);
  
  ublas::vector<double> tmp_vec (ublas_mat.size1());
  tmp_vec = prod(ublas_mat, POD::snaps_mean);
  red_vec = prod(trans(POD::pod_basis), tmp_vec);
}

/**************************************************************************** */
ublas::vector<double> ROM::reduce_mat_mean(std::shared_ptr<TMatrix> full_mat)
{
  ublas::vector<double> red_vec;
  
  reduce_mat_mean(full_mat, red_vec);
  return red_vec;
}

/**************************************************************************** */
void ROM::reduce(const ublas::vector<double> &full_vec, ublas::vector<double> &red_vec)
{
  if(full_vec.size() != POD::length)
  {
    ErrThrow("ROM::reduce(const ublas::vector<double>&, const ublas::vector<double>&): "
    		 "Length of input full-order vector must coincide with the length (dof) of "
    		 "POD basis. Length of input full-order vector: ", full_vec.size(), 
			 ", length of POD basis: ", POD::length);
  }
  red_vec.resize(POD::rank);
  red_vec.clear();
  red_vec = prod(trans(POD::pod_basis),full_vec);
}

/**************************************************************************** */
ublas::vector<double> ROM::reduce(const ublas::vector<double> &full_vec)
{
  ublas::vector<double> red_vec;
  reduce(full_vec, red_vec);
  return red_vec;
}

/**************************************************************************** */
ublas::vector<double> ROM::reduce(double *full_vec)
{
  ublas::vector<double> red_vec;
  ublas::vector<double> tmp_vec(POD::length);

  for (int i=0; i<POD::length; ++i)
  {
	tmp_vec(i) = full_vec[i];
  }
  reduce(tmp_vec, red_vec);
  return red_vec;
}

/**************************************************************************** */
void ROM::reduce_solution(double *full_sol, ublas::vector<double> &red_sol)
{
  ublas::vector<double> tmp_vec(POD::length);

  for (int i=0; i<POD::length; ++i)
  {
	if(POD::rom_db["pod_fluct"])
	  tmp_vec(i) = full_sol[i] - POD::snaps_mean(i);
	else
      tmp_vec(i) = full_sol[i];
  }
  if(POD::rom_db["pod_inner_product"].get<std::string>() != "eucl")
  {
	tmp_vec = prod(POD::gramian_mat, tmp_vec);
  }
  red_sol = prod(trans(POD::pod_basis), tmp_vec);
}

/**************************************************************************** */
void ROM::get_full_solution(const ublas::vector<double> &red_sol, double *full_sol)
{
  if(red_sol.size() != POD::rank)
  {
	ErrThrow("ROM::get_full_solution(..): Length of input reduced-order "
			"vector must conicide with the dimension of the POD basis. "
			 "Length of input vector: ", red_sol.size(), 
			 ", dimension of POD basis: ", POD::rank);
  }
  ublas::vector<double> tmp_vec(POD::length);
  tmp_vec = prod(POD::pod_basis, red_sol);
  
  for (int i=0; i<POD::length; ++i)
  {
  	if(POD::rom_db["pod_fluct"])
  	  full_sol[i] = tmp_vec(i) + POD::snaps_mean(i);
  	else
      full_sol[i] = tmp_vec(i);
  }
}

/**************************************************************************** */
void ROM::prepare_solve(const ublas::matrix<double> &mat)
{
  ublas::matrix<double> tmp_mat = mat;
  ublas::matrix<double> _inv_mat( tmp_mat.size1(),tmp_mat.size2());
  pmatrix pm( tmp_mat.size1() );
  int res = lu_factorize( tmp_mat, pm );
  _inv_mat.assign( idmatrix(tmp_mat.size1()) );
  lu_substitute(tmp_mat, pm, _inv_mat);
  this->inv_mat = _inv_mat;
}

/**************************************************************************** */
ublas::vector<double> ROM::solve(const ublas::vector<double> &rhs)
{
  if(this->inv_mat.size1() != POD::rank)
  {
    ErrThrow("ROM::solve(const ublas::vector<double>, ublas::vector<double>): "
    		 "Call first prepare_solve(const ublas::matrix<double>) to obtain "
    		 "the inverse matrix of the system matrix.");
  }
  ublas::vector<double> sol(POD::rank);
  sol = prod(this->inv_mat, rhs);
  return sol;
}

/**************************************************************************** */
ublas::vector<double> ROM::solve(const ublas::matrix<double> &mat, 
		                         const ublas::vector<double> &rhs)
{
  prepare_solve(mat);
  ublas::vector<double> sol = solve(rhs);
  return sol;
}

/**************************************************************************** */
void ROM::set_gramian(std::shared_ptr<TMatrix> gram_mat)
{
  POD::set_gramian(gram_mat);
}

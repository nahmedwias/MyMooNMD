#include <Iteration_jacobi.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>
#include <algorithm>

template <class L>
void extract_diagonal_entries(std::vector<double> & diag_entries, const L & mat)
{
  diag_entries = mat.get_diagonal();
  // make sure there are no zeros on the diagonal
  if(std::find_if(diag_entries.begin(), diag_entries.end(), 
                 [](const double& d){return d == 0.;}) != diag_entries.end())
    ErrThrow("zero entry on the diagonal using Jacobi solver/preconditioner");
}

// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_jacobi<L, V>::Iteration_jacobi(const L& mat)
 : IterativeMethod<L, V>(std::make_shared<NoPreconditioner<V>>(), "Jacobi"),
   Preconditioner<V>(),
   diagonal_entries(mat.get_n_total_rows(), 0.), linear_operator(mat)
{
  // extract the diagonal entries of the matrix
  extract_diagonal_entries(this->diagonal_entries, mat);
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class Vector>
std::pair<unsigned int, double> Iteration_jacobi<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{
  // x^{k+1} = w D^{-1}(b - Rx^k) + (1-w)x^k
  // D is the diagonal, R is A-D, w is the damping parameter
  // Since we can only multiply with all of A, we do the following:
  // -> x^{k+1} = w D^{-1}(b-Ax^k) + x^k
  
  // extract the diagonal entries of the matrix
  extract_diagonal_entries(this->diagonal_entries, A);
  size_t length = this->diagonal_entries.size();
  
  Vector r(rhs);
  A.apply_scaled_add(solution, r, -1.);
  // now: r = rhs - A*solution
  
  double normb = rhs.norm();
  if (normb == 0.0)
    normb = 1;
  // check for convergence
  double resid = r.norm() / normb;
  // safe initial residual, used to check stopping criteria later
  if(this->converged(resid, 0)) 
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  
  for(unsigned int i = 1; i <= this->max_n_iterations; ++i) 
  {
    // now: r = rhs - A*solution
    for(size_t j = 0; j < length; ++j)
      r[j] *= this->damping / this->diagonal_entries[j];
    // now: r = w D^{-1}(rhs - A*solution)
    r += solution;
    // now: r = w D^{-1}(rhs - A*solution) + solution
    std::swap(solution, r); // solution is now the new solution
    
    // compute residual
    r = rhs;
    A.apply_scaled_add(solution, r, -1.);
    
    // check for convergence
    resid = r.norm() / normb;
    if(this->converged(resid, i))
    {
      return std::pair<unsigned int, double>(i, resid);
    }
  }
  
  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Iteration_jacobi<L, V>::update(const L& A)
{
  if(&this->linear_operator != &A)
    ErrThrow("You are trying to update a Iteration_jacobi object with a matrix "
             "which is not the one you constructed it with. Please create a "
             "new Iteration_jacobi object");
  // update the diagonal entries
  extract_diagonal_entries(this->diagonal_entries, A);
  this->IterativeMethod<L, V>::update(A);
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Iteration_jacobi<L, V>::apply(const V & z, V & r) const
{
  // as if no preconditioner is used (this initialized r, in case it has not 
  // been done already)
  r = z;
  for(size_t j = 0, length = this->diagonal_entries.size(); j < length; ++j)
    r[j] /= this->diagonal_entries[j];
}

/* ************************************************************************** */
// L - LinearOperator
template <class L, class V>
const L& Iteration_jacobi<L, V>::get_operator() const
{
    return linear_operator;
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_jacobi<BlockFEMatrix, BlockVector>;
template class Iteration_jacobi<BlockMatrix,   BlockVector>;

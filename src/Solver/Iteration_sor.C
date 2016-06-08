#include <Iteration_sor.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>
#include <algorithm>

// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_sor<L, V>::Iteration_sor(const L& mat, int flag, double w)
 : IterativeMethod<L, V>(std::make_shared<NoPreconditioner<V>>(),
                         flag == 2 ? "ssor" : "sor"),
   Preconditioner<V>(),
   sor_type(flag), omega(w), linear_operator(mat)
{
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class Vector>
std::pair<unsigned int, double> Iteration_sor<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{
  // (D + wL) x^{k+1} = w b - (wU + (w-1)D)x^k
  // D is the diagonal, L is the lower triangular part of A, U is the upper 
  // triangular part of A, w is the relaxation parameter
  // Since we can only multiply with all of A, we do the following:
  // -> x^{k+1} = w (D+wL)^{-1}(b-Ax^k) + x^k
  
  Vector r(rhs);
  A.apply_scaled_add(solution, r, -1.);
  // now: r = rhs - A*solution
  
  double normb = norm(rhs);
  if (normb == 0.0)
    normb = 1;
  // check for convergence
  double resid = norm(r) / normb;
  if(this->converged(resid, resid, 0)) 
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  // safe initial residual, used to check stopping criteria later
  double resid0 = resid;
  
  for(unsigned int i = 1; i <= this->max_n_iterations; ++i) 
  {
    // one (s)sor step
    A.sor_sweep(rhs, solution, this->omega, this->sor_type);
    //solution.print(std::string("sol_") + std::to_string(i));
    // compute new residual
    r = rhs;
    A.apply_scaled_add(solution, r, -1.);
    
    // check for convergence
    resid = norm(r) / normb;
    Output::print<4>(this->name, " iteration ", i, " ", resid, "\t",
                     solution.norm());
    if(this->converged(resid, resid0, i))
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
void Iteration_sor<L, V>::apply(const V & z, V & r) const
{
  // initialize r, in case it has not been done already
  r.copy_structure(z);
  r.reset();
  this->linear_operator.sor_sweep(z, r, this->omega, this->sor_type);
}

/* ************************************************************************** */
// L - LinearOperator
template <class L, class V>
const L& Iteration_sor<L, V>::get_operator() const
{
    return linear_operator;
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_sor<BlockFEMatrix, BlockVector>;
template class Iteration_sor<BlockMatrix,   BlockVector>;
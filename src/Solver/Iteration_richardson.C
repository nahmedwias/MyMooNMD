#include <Iteration_richardson.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Iteration_jacobi.h>
#include <MooNMD_Io.h>

// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_richardson<L, V>::Iteration_richardson(
  std::shared_ptr<Preconditioner<V>> prec)
 : IterativeMethod<L, V>(prec, "Richardson")
{
  if(!prec)
  {
    ErrThrow("No preconditioner specified. Choose NoPreconditioner<Vector> if "
             "you don't need one");
  }
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class Vector>
std::pair<unsigned int, double> Iteration_richardson<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{
  Vector z;
  
  double normb = norm(rhs);
  if (normb == 0.0)
    normb = 1;
  
  //Vector r = rhs - A*solution;
  Vector r(rhs); // copy values
  A.apply_scaled_add(solution, r, -1.0);
  double resid = norm(r) / normb;
  if(this->converged(resid, resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  // safe initial residual, used to check stopping criteria later
  double resid0 = resid;
  
  for (unsigned int i = 1; i <= this->max_n_iterations; i++)
  {
    //z = M.solve(r);
    this->prec->apply(r, z);
    // update solution by w*z (w = damping factor)
    solution.add_scaled(z, this->damping);
    
    //r = b - A * x;
    r = rhs;
    A.apply_scaled_add(solution, r, -1.0);
    
    resid = norm(r) / normb;
    Output::print<4>("richardson iteration ", i, " ", resid);
    if(this->converged(resid, resid0, i))
    {
      return std::pair<unsigned int, double>(i, resid);
    }
  }
  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_richardson<BlockFEMatrix, BlockVector>;
template class Iteration_richardson<BlockMatrix,   BlockVector>;
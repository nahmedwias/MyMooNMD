#include <Iteration_cg.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>


// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_cg<L, V>::Iteration_cg(std::shared_ptr<Preconditioner<V>> prec)
 : IterativeMethod<L, V>(prec, "cg")
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
std::pair<unsigned int, double> Iteration_cg<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{
  Vector p, z;
  Vector q(rhs); // this does a copy. We don't need the entries though.
  double alpha, beta, rho;
  double rho_1=0.; // pour le compilateur.

  //Vector r = rhs - A*solution;
  Vector r(rhs); // copy values
  A.apply_scaled_add(solution, r, -1.0);
  
  double normb = rhs.norm();
  if (normb == 0.0)
    normb = 1;
  double resid = r.norm() / normb;
  // safe initial residual, used to check stopping criteria later
  if(this->converged(resid, 0)) 
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  
  for (unsigned int i = 1; i <= this->max_n_iterations; ++i) 
  {
    //z = M.solve(r);
    this->prec->apply(r, z);
    rho = dot(r, z);

    if (i == 1)
      p = z;
    else
    {
      beta = rho / rho_1;
      //p = z + beta * p;
      p *= beta;
      p += z;
    }

    //q = A*p;
    A.apply(p, q);
    alpha = rho / dot(p, q);

    //solution += alpha * p;
    solution.add_scaled(p, alpha);
    
    //r -= alpha * q;
    r.add_scaled(q, -alpha);
    
    resid = r.norm() / normb;
    if(this->converged(resid, i))
    {
      return std::pair<unsigned int, double>(i, resid);
    }
    rho_1 = rho;
  }
  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_cg<BlockMatrix,   BlockVector>;
template class Iteration_cg<BlockFEMatrix, BlockVector>;

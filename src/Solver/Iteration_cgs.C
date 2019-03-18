#include <Iteration_cgs.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Iteration_jacobi.h>
#include <MooNMD_Io.h>


// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_cgs<L, V>::Iteration_cgs(std::shared_ptr<Preconditioner<V>> prec)
 : IterativeMethod<L, V>(prec, "cgs")
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
std::pair<unsigned int, double> Iteration_cgs<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{
  double alpha, beta, rho_1, rho_2;
  Vector p, phat(rhs), q, u, uhat(rhs);
  Vector qhat(rhs), vhat(rhs); // this does a copy. We don't need the entries.

  //Vector r = rhs - A*solution;
  Vector r(rhs); // copy values
  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - A*solution
  
  //Vector rtilde = r;
  Vector rtilde(r); // copy values
  
  double normb = rhs.norm();
  if (normb == 0.0)
    normb = 1;
  
  double resid = r.norm() / normb;
  // safe initial residual, used to check stopping criteria later
  if(this->converged(resid, 0)) 
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  
  for (unsigned int i = 1; i <= this->max_n_iterations; i++)
  {
    rho_1 = dot(rtilde, r);
    if (rho_1 == 0) // this should not ever happen 
    {
      return std::pair<unsigned int, double>(i, r.norm() / normb);
    }
    if (i == 1)
    {
      u = r;
      p = u;
    } 
    else
    {
      beta = rho_1 / rho_2;
      //u = r + beta * q;
      u = q;
      u *= beta;
      u += r;
      //p = u + beta * (q + beta * p);
      p *= beta;
      p += q;
      p *= beta;
      p += u;
    }
    //phat = M.solve(p);
    this->prec->apply(p, phat);
    //vhat = A*phat;
    A.apply(phat, vhat);
    alpha = rho_1 / dot(rtilde, vhat);
    //q = u - alpha * vhat;
    q = u;
    q.add_scaled(vhat, -alpha);
    //uhat = M.solve(u + q);
    u += q;
    this->prec->apply(u, uhat);
    //solution += alpha * uhat;
    solution.add_scaled(uhat, alpha);
    //qhat = A * uhat;
    A.apply(uhat, qhat);
    //r -= alpha * qhat;
    qhat *= -alpha;
    r += qhat;
    rho_2 = rho_1;
    
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
// explicit instantiations
template class Iteration_cgs<BlockMatrix,   BlockVector>;
template class Iteration_cgs<BlockFEMatrix, BlockVector>;

#include <Iteration_bicgstab.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>
#include <typeinfo>

// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_bicgstab<L, V>::Iteration_bicgstab(
  std::shared_ptr<Preconditioner<V>> prec)
 : IterativeMethod<L, V>(prec, "bi-cgstab")
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
std::pair<unsigned int, double> Iteration_bicgstab<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{
  double rho_1, rho_2 = 1., alpha = 1., beta, omega = 1.;
  Vector p, phat(rhs), s, shat(rhs);
  Vector t(rhs), v(rhs); // this does a copy. We don't need the entries though.
  
  //Vector r = rhs - A * solution;
  Vector r(rhs); // copy values
  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - A*solution
  //Vector rtilde = r;
  Vector rtilde(r);
  
  double normb = norm(rhs);
  if (normb == 0.0)
    normb = 1;
  double resid = norm(r) / normb;
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
      return std::pair<unsigned int, double>(i, norm(r) / normb);
    }
    if (i == 1)
      p = r;
    else
    {
      beta = (rho_1/rho_2) * (alpha/omega);
      //p = r + beta * (p - omega * v);
      p.add_scaled(v, -omega);
      p *= beta;
      p += r;
    }
    //phat = M.solve(p);
    this->prec->apply(p, phat);
    //v = A * phat;
    A.apply(phat, v);
    alpha = rho_1 / dot(rtilde, v);
    //s = r - alpha * v;
    s = r;
    s.add_scaled(v, -alpha);
    resid = norm(s) / normb;
    if(this->converged(resid, i))
    {
      //solution += alpha * phat;
      solution.add_scaled(phat, alpha);
      return std::pair<unsigned int, double>(i, resid);
    }
    //shat = M.solve(s);
    this->prec->apply(s, shat);
    //t = A * shat;
    A.apply(shat, t);
    omega = dot(t,s) / dot(t,t);
    //solution += alpha * phat + omega * shat;
    solution.add_scaled(phat, alpha);
    solution.add_scaled(shat, omega);
    //r = s - omega * t;
    r = s;
    r.add_scaled(t, -omega);
    
    rho_2 = rho_1;
    resid = norm(r) / normb;
    Output::print<4>("bi-cgstab iteration ", i, " ", resid);
    if(this->converged(resid, i))
    {
      return std::pair<unsigned int, double>(i, resid);
    }
    if(omega == 0)
    {
      ErrThrow("unhappy breakdown in Bi-CGStab. iteration ", i, ", residual ",
               resid);
    }
  }
  
  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_bicgstab<BlockMatrix, BlockVector>;
template class Iteration_bicgstab<BlockFEMatrix, BlockVector>;
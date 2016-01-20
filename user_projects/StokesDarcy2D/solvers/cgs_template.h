//*****************************************************************
// Iterative template routine -- CGS
//
// CGS solves the unsymmetric linear system Ax = b 
// using the Conjugate Gradient Squared method
//
// CGS follows the algorithm described on p. 26 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************

template < class Matrix, class Vector, class Preconditioner, class Real >
int 
CGS(const Matrix &A, Vector &x, const Vector &b,
    const Preconditioner &M, unsigned int &max_iter, Real &tol)
{
  Real resid;
  Real rho_1, rho_2, alpha, beta;
  Vector p, phat(b), q, qhat, vhat, u, uhat(b);

  Real normb = norm(b);
  //Vector r = b - A*x;
  Vector r(x); // r is only a structural copy of x, all values are set to 0
  r = b;       // copy values
  A.apply_scaled_add(x, r, -1.0); // now r = b - Ax
  
  //Vector rtilde = r;
  Vector rtilde(r); // only a structural copy of r, all values are set to 0
  rtilde = r; // copy values

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = norm(r) / normb) <= tol)
  {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++)
  {
    rho_1 = dot(rtilde, r);
    if (rho_1 == 0)
    {
      tol = norm(r) / normb;
      max_iter = i;
      return 2;
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
    M.solve(p, phat);
    //vhat = A*phat;
    A.apply(phat, vhat);
    alpha = rho_1 / dot(rtilde, vhat);
    //q = u - alpha * vhat;
    q = vhat;
    q *= -alpha;
    q += u;
    //uhat = M.solve(u + q);
    u += q;
    M.solve(u, uhat);
    //x += alpha * uhat;
    uhat *= alpha;
    x += uhat;
    uhat *= 1.0/alpha;
    //qhat = A * uhat;
    A.apply(uhat, qhat);
    //r -= alpha * qhat;
    qhat *= -alpha;
    r += qhat;
    rho_2 = rho_1;
    if ((resid = norm(r) / normb) < tol)
    {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }
  
  tol = resid;
  return 1;
}


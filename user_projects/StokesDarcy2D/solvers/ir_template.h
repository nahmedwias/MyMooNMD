//*****************************************************************
// Iterative template routine -- Preconditioned Richardson
//
// IR solves the unsymmetric linear system Ax = b using 
// Iterative Refinement (preconditioned Richardson iteration).
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
IR(const Matrix &A, Vector &x, const Vector &b,
   const Preconditioner &M, unsigned int &max_iter, Real &tol)
{
  Real resid;
  Vector z;

  Real normb = norm(b);
  //Vector r = b - A*x;
  Vector r(x); // r is only a structural copy of x, all values are set to 0
  r = b;       // copy values
  A.apply_scaled_add(x, r, -1.0);
  
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
    //z = M.solve(r);
    M.solve(r,z);
    x += z;
    
    //r = b - A * x;
    r = b;
    A.apply_scaled_add(x, r, -1.0);
    
    resid = norm(r) / normb;
    if (resid <= tol)
    {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;
}




#ifndef _GMRES_H_
#define _GMRES_H_
//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
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

#include <Matrix.h>

template<class Real> 
void GeneratePlaneRotation(const Real &dx, const Real &dy, Real &cs, Real &sn)
{
  if (dy == 0.0) 
  {
    cs = 1.0;
    sn = 0.0;
  }
  else if (abs(dy) > abs(dx)) 
  {
    Real temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  }
  else 
  {
    Real temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}

template<class Real> 
void ApplyPlaneRotation(Real &dx, Real &dy, const Real &cs, const Real &sn)
{
  Real temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}

template < class Matrix, class Vector >
void 
Update(Vector &x, const int k, const Matrix &h, const Vector &s, Vector v[])
{
  Vector y(s); // y is only a structural copy of s, all values are set to 0
  y = s; // copy values

  // Backsolve:  
  for (int i = k; i >= 0; i--) 
  {
    y(i) /= h(i,i);
    for (int j = i - 1; j >= 0; j--)
      y(j) -= h(j,i) * y(i);
  }

  for (int j = 0; j <= k; j++)
  {
    //x += v[j] * y(j);
    v[j] *= y(j);
    x += v[j];
    v[j] *= 1.0 / y(j);
  }
}


template < class Real >
Real 
inline abs(Real x)
{
  return (x > 0 ? x : -x);
}


template < class Operator, class Vector, class Preconditioner,
           class Matrix, class Real >
int 
GMRES_left(const Operator &A, Vector &x, const Vector &b,
           const Preconditioner &M, Matrix &H, unsigned int &restart,
           unsigned int &max_iter, Real &tol)
{
  if(restart > max_iter)
    restart = max_iter;
  Real red_factor = std::max(TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR,
                             TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE);
  Real resid, resid0; // current and initial residual
  unsigned int i, j = 1, k;
  Vector s(restart+1), cs(restart+1), sn(restart+1), w;
  
  
  //Real normb = norm(M.solve(b));
  Vector a; // intermediate vector
  M.solve(b, a);
  //Vector r = M.solve(b - A * x);
  Vector r;
  a = b; // reuse Vector 'a' to compute residual
  A.apply_scaled_add(x, a, -1.0); // now a = b - Ax
  M.solve(a, r);
  
  Real beta = norm(r);
  
  if ((resid = norm(r)) <= tol) 
  {
    tol = resid;
    max_iter = 0;
    return 0;
  }
  resid0 = resid;
  
  Vector *v = new Vector[restart+1];
  
  while (j <= max_iter) 
  {
    //v[0] = r * (1.0 / beta);    // ??? r / beta
    v[0] = r;
    v[0] *= 1.0 / beta;
    s = 0.0;
    s(0) = beta;
    
    for (i = 0; i < restart && j <= max_iter; i++, j++) 
    {
      //w = M.solve(A * v[i]);
      A.apply(v[i], a); // a = A*v[i] // reuse Vector 'a'
      M.solve(a, w);
      
      for (k = 0; k <= i; k++) 
      {
        H(k, i) = dot(w, v[k]);
        //w -= H(k, i) * v[k];
        if(H(k, i) != 0.0)
        {
          v[k] *= H(k, i);
          w -= v[k];
          v[k] *= 1.0 / H(k, i);
        }
      }
      H(i+1, i) = norm(w);
      //v[i+1] = w * (1.0 / H(i+1, i)); // ??? w / H(i+1, i)
      v[i+1] = w;
      v[i+1] *= 1.0 / H(i+1, i);

      for (k = 0; k < i; k++)
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));
      
      if((resid = abs(s(i+1))) < tol || resid < red_factor*resid0) 
      {
        Update(x, i, H, s, v);
        tol = resid;
        max_iter = j;
        delete [] v;
        return 0;
      }
      if(TDatabase::ParamDB->SC_VERBOSE > 2)
      {
        OutPut("gmres iteration " << j << " " << resid << endl);
      }
    }
    Update(x, restart - 1, H, s, v);
    //r = M.solve(b - A * x);
    a = b;
    A.apply_scaled_add(x, a, -1.0);
    M.solve(a, r);
    
    beta = norm(r);
    if(fabs(beta-resid)>0.01*beta)
    {
      OutPut("restart residual changed " << beta << "  " << resid << endl);
    }
    if ((resid = beta) < tol || resid < red_factor*resid0) 
    {
      tol = resid;
      max_iter = j;
      delete [] v;
      return 0;
    }
    // else restart
  }
  
  tol = resid;
  delete [] v;
  return 1;
}




template < class Operator, class Vector, class Preconditioner,
           class Matrix, class Real >
int 
GMRES_right(const Operator &A, Vector &x, const Vector &b,
            const Preconditioner &M, Matrix &H, unsigned int &restart,
            unsigned int &max_iter, Real &tol)
{
  if(restart > max_iter)
    restart = max_iter;
  Real red_factor = std::max(TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR,
                             TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE);
  Real resid, resid0; // current and initial residual
  unsigned int i, j = 1, k;
  Vector s(restart+1), cs(restart+1), sn(restart+1), w(b);
  
  
  //Vector r = b - A * x;
  Vector r;
  r = b; // reuse Vector 'a' to compute residual
  A.apply_scaled_add(x, r, -1.0); // now r = b - Ax
  
  Real beta = norm(r);
  
  if ((resid = beta) <= tol) 
  {
    tol = resid;
    max_iter = 0;
    return 0;
  }
  resid0 = resid;
  
  Vector *v = new Vector[restart+1];
  
  while (j <= max_iter) 
  {
    //v[0] = r * (1.0 / beta);
    v[0] = r;
    v[0] *= 1.0 / beta;
    s = 0.0;
    s(0) = beta;
    
    for (i = 0; i < restart && j <= max_iter; i++, j++) 
    {
      //r = A * M.solve(v[i]);
      w = 0.0;
      M.solve(v[i], w); // w = M.solve(v[i])
      A.apply(w, r); // r = A * w
      
      
      for (k = 0; k <= i; k++) 
      {
        H(k, i) = dot(r, v[k]);
        //r -= H(k, i) * v[k];
        if(H(k, i) != 0.0)
        {
          v[k] *= H(k, i);
          r -= v[k];
          v[k] *= 1.0 / H(k,i);
        }
      }
      H(i+1, i) = norm(r);
      //v[i+1] = r * (1.0 / H(i+1, i));
      v[i+1] = r;
      v[i+1] *= 1.0 / H(i+1, i);

      for (k = 0; k < i; k++)
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));
      
      if((resid = abs(s(i+1))) < tol || resid < red_factor*resid0) 
      {
        w = 0.0; // reuse w for update of x
        Update(w, i, H, s, v);
        M.solve(w, r); // r = M.solve(w);
        x += r;
        
        tol = resid;
        max_iter = j;
        delete [] v;
        return 0;
      }
      if(TDatabase::ParamDB->SC_VERBOSE > 2)
      {
        OutPut("gmres iteration " << j << " " << resid << endl);
      }
    }
    w = 0.0;
    Update(w, restart - 1, H, s, v);
    M.solve(w, r); // r = M.solve(w)
    x += r; // update solution
    // compute new residual
    r = b;
    A.apply_scaled_add(x, r, -1.0); // r = b - A * x;
    
    beta = norm(r);
    if(fabs(beta-resid)>0.01*beta)
    {
      OutPut("restart residual changed " << beta << "  " << resid << endl);
    }
    if ((resid = beta) < tol || resid < red_factor*resid0) 
    {
      tol = resid;
      max_iter = j;
      delete [] v;
      return 0;
    }
    // else restart
  }
  
  tol = resid;
  delete [] v;
  return 1;
}

template < class Operator, class Vector, class PreconditioningStrategy,
           class Matrix, class Real >
int
FGMRES(const Operator &A, Vector &x, const Vector &b,
       const PreconditioningStrategy &P, Matrix &H, unsigned int &restart,
       unsigned int &max_iter, Real &tol)
{
  if(restart > max_iter)
    restart = max_iter;
  Real red_factor = std::max(TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR,
                             TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE);
  Real resid, resid0; // current and initial residual
  unsigned int i, j = 1, k;
  Vector s(restart+1), cs(restart+1), sn(restart+1);


  //Vector r = b - A * x;
  Vector r;
  r = b; // reuse Vector 'r' to compute residual
  A.apply_scaled_add(x, r, -1.0); // now r = b - Ax

  Real beta = norm(r);

  if ((resid = beta) <= tol)
  {
    tol = resid;
    max_iter = 0;
    return 0;
  }
  resid0 = resid;

  Vector *v = new Vector[restart+1];
  Vector *z = new Vector[restart+1]; //array to store the outputs of the preconditioning processes

  while (j <= max_iter)
  {
    //v[0] = r * (1.0 / beta);
    v[0] = r;
    v[0] *= 1.0 / beta;
    s = 0.0;
    s(0) = beta;

    for (i = 0; i < restart && j <= max_iter; i++, j++)
    {
      //z[i] = A * M.solve(v[i]); //where M.solve(v[i]) is replaced by application of a preconditioning strategy
      z[i]=r; //copy structure of r
      z[i] = 0.0;
      P.solve(i,j,v[i],z[i]); // apply a preconditioning strategy with rhs v[i] to obtain z[i]
      A.apply(z[i], r); // r = A * z[i]


      for (k = 0; k <= i; k++)
      {
        H(k, i) = dot(r, v[k]);
        //r -= H(k, i) * v[k];
        if(H(k, i) != 0.0)
        {
          v[k] *= H(k, i);
          r -= v[k];
          v[k] *= 1.0 / H(k,i);
        }
      }
      H(i+1, i) = norm(r);

      //v[i+1] = r * (1.0 / H(i+1, i)); //can lead to "unhappy" crash in FGMRES if H(i+1, i) == 0
      v[i+1] = r;
      if(H(i+1, i) != 0){
    	  v[i+1] *= 1.0 / H(i+1, i);
      } else {
    	  Error("Unhappy breakdown in FGMRES Iteration " << j << ". " << endl);
    	  exit(-1);
      }

      for (k = 0; k < i; k++)
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));

      GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));

      if((resid = abs(s(i+1))) < tol || resid < red_factor*resid0)
      {
    	//use r as auxiliary here
        r = 0.0; //set entries to zero
        Update(r, i, H, s, z); // solve the least squares problem by backsolve
        x += r; //and update the solution

        tol = resid;
        max_iter = j;
        delete [] v;
        delete [] z;
        return 0;
      }
      if(TDatabase::ParamDB->SC_VERBOSE > 2)
      {
        OutPut("FGMRES iteration " << j << " " << resid << endl);
      }
    }
	//use r as auxiliary here
    r = 0.0; //set entries to zero - just in case
    Update(r, restart - 1, H, s, z); // solve the least squares problem by backsolve
    x += r; //and update the solution

    // compute new residual r
    r = b;
    A.apply_scaled_add(x, r, -1.0); // r = b - A * x;
    // and new defect beta
    beta = norm(r);

    if(fabs(beta-resid)>0.01*beta)
    {
      OutPut("restart residual changed " << beta << "  " << resid << endl);
    }
    if ((resid = beta) < tol || resid < red_factor*resid0)
    {
      tol = resid;
      max_iter = j;
      delete [] v;
      delete [] z;
      return 0;
    }
    // else restart
  }

  tol = resid;
  delete [] v;
  delete [] z;
  return 1;
}


class TriangularMatrix // upper triangular matrix and one off diagonal
{
public:
  TMatrix * M;
  
  TriangularMatrix(const unsigned int size)
  {
    int n_elements = (size*size+size)/2 + size-1;
    
    int * RowPtr = new int[size+1];
    int * ColPtr = new int[n_elements];
    RowPtr[0] = 0;
    for(unsigned int row = 1; row <= size; row++)
    {
      int n = size-row+2; // number of entries in this row
      if(row==1) n  -= 1;
      RowPtr[row] = RowPtr[row-1] + n;
      for(int col = 0; col < n; col++)
      {
        ColPtr[RowPtr[row-1]+col] = col + row - 2;
        if(row==1) ColPtr[RowPtr[row-1]+col] += 1; 
      }
    }
    
    
    double * elements = new double[n_elements];
    std::fill(elements, elements+n_elements, 0.0);
    TStructure * structure = new TStructure(size, size, n_elements, ColPtr, 
                                            RowPtr);
    M = new TMatrix(structure, elements);
  };
  
  double& operator()(int i, int j)
  {
    if(i > j+1)
    {
      cout << "ERROR in TriangularMatrix : index (" << i << "," << j << ") is "
           << "not upper triangular\n";
      exit(0);
    }
    return (*M)(i,j);
  }
  
  const double& operator()(int i, int j) const
  {
    if(i > j+1)
    {
      cout << "ERROR in TriangularMatrix : index (" << i << "," << j << ") is "
           << "not upper triangular\n";
      exit(0);
    }
    return (*M)(i,j);
  }
  
  ~TriangularMatrix()
  {
    delete M->GetStructure();
    delete M;
  }
  
  void print_full() const
  {
    this->M->PrintFull("H");
  }
};


#endif

#include <Iteration_gmres.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <cmath> // std::abs
#include <MooNMD_Io.h>

template <class LinearOperator, class Vector>
Iteration_gmres<LinearOperator, Vector>::Iteration_gmres(
  std::shared_ptr<Preconditioner<Vector>> prec, gmres_type t)
 : IterativeMethod<LinearOperator, Vector>(prec, "gmres"), type(t)
{
  if(!prec)
  {
    ErrThrow("No preconditioner specified. Choose NoPreconditioner<Vector> if "
             "you don't need one");
  }
  // change name according to the type and check if type is known
  switch(this->type)
  {
    case gmres_type::left:
      this->name = "left gmres";
      break;
    case gmres_type::right:
      this->name = "right gmres";
      break;
    case gmres_type::flexible:
      this->name = "flexible gmres";
      break;
    default:
      ErrThrow("unknown gmres type ");
      break;
  }
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::iterate(const LinearOperator & A,
                                                 const Vector & rhs,
                                                 Vector & solution)
{
  switch(this->type)
  {
    case gmres_type::left:
      return left_gmres(A, rhs, solution);
      break;
    case gmres_type::right:
      return right_gmres(A, rhs, solution);
      break;
    case gmres_type::flexible:
      return flexible_gmres(A, rhs, solution);
      break;
    default:
      ErrThrow("unknown gmres type ");
      break;
  }
}

/* ************************************************************************** */
/// @brief a simple class to store an upper triangular matrix and one off 
/// diagonal
///
/// Example:
///     ( * * * * * )
///     ( * * * * * )
///     (   * * * * )
///     (     * * * )
///     (       * * )
/// This is needed for gmres.
class TriangularMatrix
{
  public:
    // one vector for each row
    std::vector<std::vector<double>> entries;
    
    TriangularMatrix(const unsigned int size) : entries(size)
    {
      Output::print<5>("TriangularMatrix constructor ", size);
      //unsigned int n_elements = (size*size+size)/2 + size-1;
      // the first two rows have `size` many columns
      entries[0] = std::vector<double>(size, 0.);
      entries[1] = std::vector<double>(size, 0.);
      // the following rows have one column less the their previous row, the last
      // row has two entries
      for(unsigned int row = 2; row < size; ++row)
      {
        size_t n = size-row+1; // number of entries in this row
        entries[row] = std::vector<double>(n);
      }
    };
    
    double& operator()(int i, int j)
    {
      if(i > j+1)
        ErrThrow("TriangularMatrix: index (", i, ",", j, ") is not upper "
                 "triangular");
      // notice that index 0 in  entries[i] corresponds to column i-1
      j = 1 + j - i;
      return entries.at(i).at(j);
    }
    
    const double& operator()(int i, int j) const
    {
      if(i > j+1)
        ErrThrow("TriangularMatrix: index (", i, ",", j, ") is not upper "
                 "triangular");
      // notice that index 0 in  entries[i] corresponds to column i-1
      j = 1 + j - i;
      return entries.at(i).at(j);
    }
};

/* ************************************************************************** */
void GeneratePlaneRotation(const double &dx, const double &dy, double &cs, 
                           double &sn)
{
  if (dy == 0.0) 
  {
    cs = 1.0;
    sn = 0.0;
  }
  else if (std::abs(dy) > std::abs(dx)) 
  {
    double temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  }
  else 
  {
    double temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}
void ApplyPlaneRotation(double &dx, double &dy, const double &cs,
                        const double &sn)
{
  double temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}
template <class Matrix, class Vector>
void Update(Vector &x, const int k, const Matrix &h, std::vector<double> y,
            Vector v[])
{
  // Backsolve:  
  for (int i = k; i >= 0; i--) 
  {
    y.at(i) /= h(i,i);
    for (int j = i - 1; j >= 0; j--)
      y[j] -= h(j,i) * y[i];
  }

  for (int j = 0; j <= k; j++)
  {
    //x += v[j] * y[j];
    v[j] *= y[j];
    x += v[j];
    v[j] *= 1.0 / y[j];
  }
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::left_gmres(const LinearOperator & A,
                                                    const Vector & rhs,
                                                    Vector & solution)
{
  std::vector<double> s(this->restart+1);
  std::vector<double> cs(this->restart+1);
  std::vector<double> sn(this->restart+1);
  Vector w;
  
  TriangularMatrix H(this->restart+1);
  
  //Vector r = prec.apply(rhs - A * solution);
  Vector r;
  Vector a = rhs; // intermediate vector
  A.apply_scaled_add(solution, a, -1.0); // now a = rhs - A*solution
  this->prec->apply(a, r);
  
  double resid = norm(r); // compute initial residual
  double beta = resid; // initialize beta as initial residual
  // safe initial residual, used to check stopping criteria later
  this->initial_residual = resid;
  if(this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  
  Vector *v = new Vector[this->restart+1];
  
  unsigned int j = 1; // iteration index
  while (j <= this->max_n_iterations) 
  {
    //v[0] = r * (1.0 / beta);    // ??? r / beta
    v[0] = r;
    v[0] *= 1.0 / beta;
    std::fill(s.begin(), s.end(), 0.0); // set all entries in s to zero
    s[0] = beta;
    
    for(unsigned int i = 0; i < this->restart && j <= this->max_n_iterations;
        i++, j++) 
    {
      //w = M.solve(A * v[i]);
      A.apply(v[i], a); // a = A*v[i] // reuse Vector 'a'
      this->prec->apply(a, w);
      
      for (unsigned int k = 0; k <= i; k++) 
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

      for (unsigned int k = 0; k < i; k++)
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
      
      resid = std::abs(s[i+1]);
      if(this->converged(resid, j))
      {
        Update(solution, i, H, s, v);
        delete [] v;
        return std::pair<unsigned int, double>(j, resid);
      }
    }
    Update(solution, this->restart - 1, H, s, v);
    //r = M.solve(b - A * x);
    a = rhs;
    A.apply_scaled_add(solution, a, -1.0);
    this->prec->apply(a, r);
    
    beta = norm(r);
    if(fabs(beta-resid)>0.01*beta)
    {
      Output::print("restart residual changed ", beta, "  ", resid);
    }
    if(this->converged(resid, j))
    {
      delete [] v;
      return std::pair<unsigned int, double>(j, resid);
    }
    // else restart
  }
  
  delete [] v;
  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::right_gmres(const LinearOperator & A,
                                                     const Vector & rhs,
                                                     Vector & solution)
{
  std::vector<double> s(this->restart+1);
  std::vector<double> cs(this->restart+1);
  std::vector<double> sn(this->restart+1);
  Vector w(rhs);
  
  TriangularMatrix H(this->restart+1);
  
  //Vector r = b - A * solution;
  Vector r;
  r = rhs; // reuse Vector 'a' to compute residual
  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - Ax
  
  double resid = norm(r); // compute initial residual
  double beta = norm(r); // initialize beta as initial residual
  // safe initial residual, used to check stopping criteria later
  this->initial_residual = resid;
  if(this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  
  Vector *v = new Vector[this->restart+1];
  
  unsigned int j = 1; // iteration index
  while (j <= this->max_n_iterations) 
  {
    //v[0] = r * (1.0 / beta);
    v[0] = r;
    v[0] *= 1.0 / beta;
    std::fill(s.begin(), s.end(), 0.0); // set all entries in s to zero
    s[0] = beta;
    
    for(unsigned int i = 0; i < this->restart && j <= this->max_n_iterations;
        i++, j++) 
    {
      //r = A * M.solve(v[i]);
      w = 0.0;
      this->prec->apply(v[i], w); // w = M.solve(v[i])
      A.apply(w, r); // r = A * w
      
      
      for (unsigned int k = 0; k <= i; k++) 
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

      for (unsigned int k = 0; k < i; k++)
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
      
      resid = std::abs(s[i+1]);
      if(this->converged(resid, j)) 
      {
        w = 0.0; // reuse w for update of solution
        Update(w, i, H, s, v);
        this->prec->apply(w, r); // r = M.solve(w);
        solution += r;
        
        delete [] v;
        return std::pair<unsigned int, double>(j, resid);
      }
    }
    w = 0.0;
    Update(w, this->restart - 1, H, s, v);
    this->prec->apply(w, r); // r = M.solve(w)
    solution += r; // update solution
    // compute new residual
    r = rhs;
    A.apply_scaled_add(solution, r, -1.0); // r = rhs - A * solution;
    
    beta = norm(r);
    if(fabs(beta-resid)>0.01*beta)
    {
      Output::print<1>("restart residual changed ", beta, "  ", resid);
    }
    if(this->converged(resid, j)) 
    {
      delete [] v;
      return std::pair<unsigned int, double>(j, resid);
    }
    // else restart
  }
  
  // not converged
  delete [] v;
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::flexible_gmres(const LinearOperator& A,
                                                        const Vector & rhs,
                                                        Vector & solution)
{
  std::vector<double> s(this->restart+1);
  std::vector<double> cs(this->restart+1);
  std::vector<double> sn(this->restart+1);
  
  TriangularMatrix H(this->restart+1);
  
  //Vector r = rhs - A * solution;
  Vector r(rhs); // copy values
  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - A*solution
  
  double resid = norm(r); // compute initial residual
  double beta = resid; // initialize beta as initial residual
  // safe initial residual, used to check stopping criteria later
  this->initial_residual = resid;
  if(this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  
  Vector *v = new Vector[this->restart+1];
  //array to store the outputs of the preconditioning processes
  Vector *z = new Vector[this->restart+1];
  
  unsigned int j = 1;
  while (j <= this->max_n_iterations)
  {
    //v[0] = r * (1.0 / beta);
    v[0] = r;
    v[0] *= 1.0 / beta;
    std::fill(s.begin(), s.end(), 0.0); // set all entries in s to zero
    s[0] = beta;
    
    for(unsigned int i = 0; i < this->restart && j <= this->max_n_iterations;
        i++, j++)
    {
      //z[i] = A * M.solve(v[i]); //where M.solve(v[i]) is replaced by application of a preconditioning strategy
      z[i]=r; //copy structure of r
      z[i] = 0.0;
      // apply a preconditioning strategy with rhs v[i] to obtain z[i]
      this->prec->apply(i, j, v[i], z[i]);
      A.apply(z[i], r); // r = A * z[i]
      
      for (unsigned int k = 0; k <= i; k++)
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
      
      //can lead to "unhappy" crash in FGMRES if H(i+1, i) == 0
      //v[i+1] = r * (1.0 / H(i+1, i)); 
      v[i+1] = r;
      if(H(i+1, i) != 0)
      {
        v[i+1] *= 1.0 / H(i+1, i);
      }
      else
      {
        ErrThrow("Unhappy breakdown in flexible gmres Iteration ", j);
      }
      
      for (unsigned int k = 0; k < i; k++)
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
      
      resid = std::abs(s[i+1]);
      if(this->converged(resid, j))
      {
        //use r as auxiliary here
        r = 0.0; //set entries to zero
        Update(r, i, H, s, z); // solve the least squares problem by backsolve
        solution += r; //and update the solution
        
        delete [] v;
        delete [] z;
        return std::pair<unsigned int, double>(j, resid);
      }
    }
    //use r as auxiliary here
    r = 0.0; //set entries to zero - just in case
    // solve the least squares problem by backsolve
    Update(r, this->restart - 1, H, s, z); 
    solution += r; //and update the solution
    
    // compute new residual r
    r = rhs;
    A.apply_scaled_add(solution, r, -1.0); // r = rhs - A * solution;
    // and new defect beta
    beta = norm(r);
    
    if(std::abs(beta - resid) > 0.01*beta)
    {
      Output::print<1>("restart residual changed ", beta, "  ", resid);
    }
    if(this->converged(resid, j))
    {
      delete [] v;
      delete [] z;
      return std::pair<unsigned int, double>(j, resid);
    }
    // else restart
  }
  
  // not converged
  delete [] v;
  delete [] z;
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_gmres<BlockMatrix,   BlockVector>;
template class Iteration_gmres<BlockFEMatrix, BlockVector>;

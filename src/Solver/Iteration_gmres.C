#include <Iteration_gmres.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <cmath> // std::abs
#include <MooNMD_Io.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#include <mpi.h>
#endif

std::string type_name(gmres_type t)
{
  switch(t)
  {
    case gmres_type::left:
      return "left gmres";
      break;
    case gmres_type::right:
      return "right gmres";
      break;
    case gmres_type::flexible:
      return "flexible gmres";
      break;
    default:
      ErrThrow("unknown gmres type ");
      break;
  }
}

template <class LinearOperator, class Vector>
Iteration_gmres<LinearOperator, Vector>::Iteration_gmres(
  std::shared_ptr<Preconditioner<Vector>> prec, gmres_type t)
 : IterativeMethod<LinearOperator, Vector>(prec, type_name(t)), type(t), s(),
   cs(), sn(), v(), z()
{
  if(!prec)
  {
    ErrThrow("No preconditioner specified. Choose NoPreconditioner<Vector> if "
             "you don't need one");
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
      // the following rows have one column less the their previous row, the
      // last row has two entries
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
            std::vector<Vector>& v)
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
    x.add_scaled(v[j], y[j]);
  }
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::left_gmres(const LinearOperator & A,
                                                    const Vector & rhs,
                                                    Vector & solution)
{
  this->s.resize(this->restart+1);
  this->cs.resize(this->restart+1);
  this->sn.resize(this->restart+1);
  Vector w;
  
  TriangularMatrix H(this->restart+1);
  
  //Vector r = prec.apply(rhs - A * solution);
  Vector r;
  Vector a = rhs; // intermediate vector
  A.apply_scaled_add(solution, a, -1.0); // now a = rhs - A*solution
  this->prec->apply(a, r);
  
  double resid = r.norm(); // compute initial residual
  double beta = resid; // initialize beta as initial residual
  // safe initial residual, used to check stopping criteria later
  if(this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  
  this->v.resize(this->restart+1);
  
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
      H(i+1, i) = w.norm();
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
        return std::pair<unsigned int, double>(j, resid);
      }
    }
    Update(solution, this->restart - 1, H, s, v);
    //r = M.solve(b - A * x);
    a = rhs;
    A.apply_scaled_add(solution, a, -1.0);
    this->prec->apply(a, r);
    
    beta = r.norm();
    if(fabs(beta-resid)>0.01*beta)
    {
      Output::print("restart residual changed ", beta, "  ", resid);
    }
    if(this->converged(resid, j))
    {
      return std::pair<unsigned int, double>(j, resid);
    }
    // else restart
  }

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
  this->s.resize(this->restart+1);
  this->cs.resize(this->restart+1);
  this->sn.resize(this->restart+1);
  Vector w(rhs);
  
  TriangularMatrix H(this->restart+1);
  
  //Vector r = b - A * solution;
  Vector r;
  r = rhs; // reuse Vector 'a' to compute residual
  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - Ax
  
  double resid =r.norm(); // compute initial residual
  double beta = r.norm(); // initialize beta as initial residual
  // safe initial residual, used to check stopping criteria later
  if(this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  
  this->v.resize(this->restart+1);
  
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
      H(i+1, i) = r.norm();
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
    
    beta = r.norm();
    if(fabs(beta-resid)>0.01*beta)
    {
      Output::print<1>("restart residual changed ", beta, "  ", resid);
    }
    if(this->converged(resid, j)) 
    {
      return std::pair<unsigned int, double>(j, resid);
    }
    // else restart
  }
  
  // not converged
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::flexible_gmres(const LinearOperator& A,
                                                        const Vector & rhs,
                                                        Vector & solution)
{
  //MPI: rhs and solution in consistency level 0 for computation of global norm
#ifdef _MPI
  int size, my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::vector<const TParFECommunicator3D*> comms = A.get_communicators();
#endif

  this->s.resize(this->restart+1);
  this->cs.resize(this->restart+1);
  this->sn.resize(this->restart+1);
  
  TriangularMatrix H(this->restart+1);

  //Vector r = rhs - A * solution;
  Vector r(rhs); // copy values
#ifdef _MPI
    //MPI: solution in consistency level 3 for vector.matrix multiplication
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(solution.block(bl), 3);
    }
#endif
  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - A*solution

  //compute initial residual
#ifndef _MPI
  double resid = r.norm();
#elif _MPI
  double resid = r.norm(comms);
#endif

  double beta = resid; // initialize beta as initial residual

  // safe initial residual, used to check stopping criteria later
  if(this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }
  
  this->v.resize(this->restart+1);
  //array to store the outputs of the preconditioning processes
  this->z.resize(this->restart+1);
  
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
#ifdef _MPI
    //MPI: solution in consistency level 3 for computation of global norm
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(z[i].block(bl), 3);
    }
#endif
      A.apply(z[i], r); // r = A * z[i]
      
      for (unsigned int k = 0; k <= i; k++)
      {
#ifndef _MPI
        H(k, i) = dot(r, v[k]);
#elif _MPI
        H(k, i) = dot(r, v[k], comms);
#endif
        //r -= H(k, i) * v[k];
        r.add_scaled(v[k], -H(k,i));
      }
#ifndef _MPI
      H(i+1, i) = r.norm();
#elif _MPI
      H(i+1, i) = r.norm(comms);
#endif
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
        
        return std::pair<unsigned int, double>(j, resid);
      }
    }
    //use r as auxiliary here
    r = 0.0; //set entries to zero - just in case
    // solve the least squares problem by backsolve
    Update(r, this->restart - 1, H, s, z); 
    solution += r; //and update the solution
    
#ifdef _MPI
    //MPI: solution in consistency level 3 for computation of global norm
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(solution.block(bl), 3);
    }
#endif
    // compute new residual r
    r = rhs;
    A.apply_scaled_add(solution, r, -1.0); // r = rhs - A * solution;
    // and new defect beta
#ifndef _MPI
    beta = r.norm();
#elif _MPI
    beta = r.norm(comms);
#endif
    
    if(std::abs(beta - resid) > 0.01*beta)
    {
      Output::print<1>("restart residual changed ", beta, "  ", resid);
    }
    if(this->converged(resid, j))
    {
      return std::pair<unsigned int, double>(j, resid);
    }
    // else restart
  }
  // not converged
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_gmres<BlockFEMatrix,   BlockVector>;
// In MPI case we are so dependent on the connection of Matrix and FESpace, that
// it does not make sense to instantiate the function for BlockMatrix.
#ifndef _MPI
template class Iteration_gmres<BlockMatrix,   BlockVector>;
#endif

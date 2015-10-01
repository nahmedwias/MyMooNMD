#ifndef __TEMPLATE_SOLVER__
#define __TEMPLATE_SOLVER__

/** *********************************************************************** 
* @(#)solver.h        
* 
* Class:   solver
*
* Purpose: class to solve A*x = b, where A is a linear operator (e.g. a matrix), 
*          b is a vector representing the right hand side, and x is a vector 
*          representing the solution. 
*
*          The vector-class 'Vector' and the linear operator-class 
*          'LinearOperator' have to implemenent the following 
*          operators / functions:
*          - friend double norm(Vector b); // the (eucledian) norm of the vector
*          - Vector(const Vector& r);  // constructor, does not copy values
*          - Vector(); // empty constructor
*          - Vector& operator=(const Vector& r); // copy values
*          - Vector& operator*=(double a);   // scale all values by 'a'
*          - Vector& operator+=(const Vector& b); // add other vector
*          - Vector& operator-=(const Vector& b); // substract other vector
*          - friend double dot(const Vector& a, const Vector& b);//inner product
*          - void LinearOperator::apply(const Vector & x, Vector & y) const;
*            // set y = A*x
*          - void LinearOperator::apply_scaled_add(const Vector & x, Vector & y,
*                                                  double a)
*            // set y = y + a * A*x
*
* @author     Ulrich Wilbrandt
* @date       20.11.13
*
***************************************************************************/
#include <Database.h>
#include <MooNMD_Io.h>

#include <cg_template.h>
#include <cgs_template.h>
#include <ir_template.h>
#include <gmres_template.h>
#include <UmfpackSolverWrapper.h>
#include <memory>

enum solver_type {direct_umfpack, direct_pardiso, richardson, cg, cgs, gmres_left, gmres_right, fgmres};

template < class LinearOperator, class Vector, class Preconditioner>
class solver
{
  protected:
    LinearOperator* A;
    Vector* x; // solution
    Vector* b; // right hand side
    Preconditioner *M;
    solver_type type;
    
    // store the number of used iterations
    unsigned int n_iterations;
    
    // for direct solvers (usually nullptr)
    UmfpackSolverWrapper* umfpack_solver;
    
    std::shared_ptr<TMatrix> Mat;

  public:
    // constructor
    // TODO: Implement factorisation for direct solvers
    solver(LinearOperator* AA, Vector* xx, Vector* bb, Preconditioner* MM);
 
    ~solver();

    // apply solver, default values are replaced by values from database (i.e. 
    // from input file)
    void solve(unsigned int max_it = 0, double tol = -1.0, double red = -1.0,
               unsigned int restart = 0);
    
    void set_type(solver_type t) { type = t; }
    
    unsigned int get_n_iterations() const { return n_iterations; }
};

/** ************************************************************************ */
template < class LinearOperator, class Vector, class Preconditioner >
solver < LinearOperator, Vector, Preconditioner >::solver(LinearOperator* AA,
                                                          Vector* xx,
                                                          Vector* bb,
                                                          Preconditioner* MM)
 : A(AA), x(xx), b(bb), M(MM), n_iterations(0), umfpack_solver(nullptr),
   Mat(nullptr)
{
  //OutPut("constructor of template solver class\n");
  // set initial iterate (the contructor above only copied the structure, not
  // the values)
  xx = x;

  if(TDatabase::ParamDB->SOLVER_TYPE == 2)
  {
    Mat = A->get_combined_matrix();
    type = direct_umfpack; // pardiso not supported currently
    umfpack_solver = new UmfpackSolverWrapper(Mat.get());
  }
  else if(TDatabase::ParamDB->SC_SOLVER_SCALAR == 11)
    type = richardson;
  else if(TDatabase::ParamDB->SC_SOLVER_SCALAR == 12)
    type = cg;
  else if(TDatabase::ParamDB->SC_SOLVER_SCALAR == 13)
    type = cgs;
  else if(TDatabase::ParamDB->SC_SOLVER_SCALAR == 14)
	  type = gmres_left;// left gmres
  else if ( TDatabase::ParamDB->SC_SOLVER_SCALAR == 15)
	  type = gmres_right;// right gmres
  else if(TDatabase::ParamDB->SC_SOLVER_SCALAR == 16)
      type = fgmres;    // flexible gmres

  else
    ErrThrow("unknown solver settings");
};

/** ************************************************************************ */
template < class LinearOperator, class Vector, class Preconditioner>
void solver<LinearOperator, Vector, Preconditioner>::solve(unsigned int max_it, 
                                                           double tol,
                                                           double red, 
                                                           unsigned int restart)
{
  //OutPut("solve using template solver class\n");
  if (max_it == 0)
    max_it = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
  if(tol < 0.0)
    tol = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
  if(red < 0.0)
    red = TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
  if(restart == 0)
    restart = TDatabase::ParamDB->SC_GMRES_RESTART;
  
  switch(type)
  {
    case direct_umfpack:
    case direct_pardiso:
    {
      umfpack_solver->solve(b->get_entries(), x->get_entries());
      break;
    }
    case richardson: // aka iterative refinement
    {
      int ret = IR( *A, *x, *b, *M, max_it, tol);
      this->n_iterations += max_it;
      if(TDatabase::ParamDB->SC_VERBOSE && ret == 0)
        OutPut(" richardson iterations " << max_it << "\ttolerance " << tol 
               << endl);
      if(ret == 1)
        OutPut(" richardson did not converge within " << max_it << 
               " iterations, current residual " << tol << endl);
      break;
    }
    case cg:
    {
      int ret = CG( *A, *x, *b, *M, max_it, tol);
      this->n_iterations += max_it;
      if(TDatabase::ParamDB->SC_VERBOSE && ret == 0)
      {
        OutPut(" cg iterations " << max_it << "\ttolerance " << tol << endl);
      }
      if(ret == 1)
      {
        OutPut(" cg did not converge within " << max_it << 
               " iterations, current residual " << tol << endl);
      }
      break;
    }
    case cgs:
    {
      int ret = CGS( *A, *x, *b, *M, max_it, tol);
      this->n_iterations += max_it;
      if(TDatabase::ParamDB->SC_VERBOSE && ret == 0)
      {
        OutPut(" cgs iterations " << max_it << "\ttolerance " << tol << endl);
      }
      else if(ret == 1)
      {
        OutPut(" cgs did not converge within " << max_it << 
               " iterations, current residual " << tol << endl);
      }
      else
        OutPut("WARNING! something wrong with cgs. " << max_it << "  " << tol 
               << endl);
      break;
    }
    case gmres_left:
    {
      if(restart > max_it) 
        restart = max_it;
      TriangularMatrix H(restart + 1); // triangular matrix
      int ret;
      //apply GMRES_left
      ret = GMRES_left(*A, *x, *b, *M, H, restart, max_it, tol);

      this->n_iterations += max_it;
      if(TDatabase::ParamDB->SC_VERBOSE && ret == 0)
      {
        OutPut(" gmres_left iterations " << max_it << "\ttolerance " << tol 
<<endl);
      }
      if(ret == 1)
      {
        OutPut(" gmres_left did not converge within " << max_it <<
               " iterations, current residual " << tol << endl);
      }
      break;
    }
    case gmres_right:
    {
      if(restart > max_it)
        restart = max_it;
      TriangularMatrix H(restart + 1); // triangular matrix
      int ret;
      //apply GMRES_right
      ret = GMRES_right(*A, *x, *b, *M, H, restart, max_it, tol);

      this->n_iterations += max_it;
      if(TDatabase::ParamDB->SC_VERBOSE && ret == 0)
      {
        OutPut(" gmres_right iterations " << max_it << "\ttolerance " << tol 
               << endl);
      }
      if(ret == 1)
      {
        OutPut(" gmres_right did not converge within " << max_it <<
               " iterations, current residual " << tol << endl);
      }
      break;
    }
    case fgmres:
    {
      if(restart > max_it)
        restart = max_it;
      TriangularMatrix H(restart + 1); // triangular matrix
      int ret;
      ret = FGMRES(*A, *x, *b, *M, H, restart, max_it, tol); //apply FGMRES

      this->n_iterations += max_it;
      if(TDatabase::ParamDB->SC_VERBOSE && ret == 0)
      {
        OutPut(" fgmres iterations " << max_it << "\ttolerance " << tol <<endl);
      }
      if(ret == 1)
      {
        OutPut(" fgmres did not converge within " << max_it <<
               " iterations, current residual " << tol << endl);
      }
      break;
    }
  }
}

/** ************************************************************************ */
template < class LinearOperator, class Vector, class Preconditioner>
solver<LinearOperator, Vector, Preconditioner>::~solver()
{
  
}

#endif

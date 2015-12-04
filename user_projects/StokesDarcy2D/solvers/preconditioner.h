#ifndef __TEMPLATE_PRECONDITIONER__
#define __TEMPLATE_PRECONDITIONER__

/** *********************************************************************** 
* @(#)preconditioner.h        
* 
* Class:   preconditioner
*
* Purpose: class to solve P*x = z, where P is preconditioner (e.g. a matrix), 
*          z is a vector representing the right hand side, and x is a vector 
*          representing the solution. 
*
*          The vector-class 'Vector' and the linear operator-class 
*          'LinearOperator' have to implemenent the following 
*          operators / functions:
*          - Vector& operator=(const Vector& r); // copy values
*          - const double & LinearOperator::operator()(const int i, const int j) const;
*
*          The last function is needed to directly access values which is needed
*          if e.g. a Jacobi preconditioner is used. This is of course only 
*          possible, if the linear operator is a matrix. Otherwise simply 
*          implement a dummy function and don't use the Jacobi preconditioner
*          which would not make sence in this case anyway.
*
* @author     Ulrich Wilbrandt
* @date       20.11.13
*
***************************************************************************/

#include <Database.h>

enum preconditioner_type {no_preconditioner, Jacobi};

template < class LinearOperator, class Vector >
class preconditioner
{
  LinearOperator * M;
  preconditioner_type type;

public:
  preconditioner(LinearOperator * m);
  void solve(const Vector &z, Vector &r) const;

  /** @brief Method to use the preconditioner in FGMRES.
   *
   * So far i and j get ignored and simply solve(z,r) is called.
   *
   * @param i Number of current iteration since last restart in FGMRES.
   * @param j Number of current iteration in FGMRES.
   * @param z The right hand side  of the preconditioning.
   * @param r The obtained vector.
   */
  void solve(int i, int j, const Vector &z, Vector &r) const{
	  solve(z,r);
  }
};


template < class LinearOperator, class Vector >
preconditioner<LinearOperator, Vector>::preconditioner(LinearOperator * m): M(m)
{
  if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 1)
    type = Jacobi;
  else if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == -4711)
    type = no_preconditioner;
  else
  {
    type = no_preconditioner;
    if(TDatabase::ParamDB->SOLVER_TYPE != 2) // not direct solver
      Output::print<1>("no preconditioner secified");
  }
}

template < class LinearOperator, class Vector >
void preconditioner<LinearOperator, Vector>::solve(const Vector &z, Vector &r) const
{
  switch(type)
  {
    case no_preconditioner:
      r = z;
      break;
    case Jacobi:
    {
      r = z;
      unsigned int l = z.length();
      for(unsigned int i = 0; i < l; i++)
      {
        r[i] /= (*M)(i,i);
      }
      break;
    }
  }
}



#endif

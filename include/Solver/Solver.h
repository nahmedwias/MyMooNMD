#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <ParameterDatabase.h>
#include <DirectSolver.h>
#include <memory>

template <class LinearOperator, class Vector>
class Solver
{
  public:
    
    Solver(const ParameterDatabase& param_db);
    
    void update_matrix(const LinearOperator& matrix);
    
    /// @brief solve after calling `Solver::update_matrix`
    ///
    /// This only makes sense for direct solvers, where the factorization is 
    /// stored. This way you can solve many times for different right hand 
    /// sides.
    ///
    /// Using iterative solvers one should use the other Solver::solve method
    /// which also has the matrix as an argument.
    void solve(const Vector& rhs, Vector& solution);
    
    void solve(const LinearOperator& matrix, const Vector& rhs,
               Vector& solution);
    
    const ParameterDatabase& get_db();
    
  protected:
    
    ParameterDatabase db;
    
    std::unique_ptr<DirectSolver> direct_solver;
};

#endif // __SOLVER_H__
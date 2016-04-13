#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <ParameterDatabase.h>
#include <DirectSolver.h>
#include <BlockVector.h>
#include <BlockMatrix.h>
#include <memory>

class Solver
{
  public:
    
    Solver(const ParameterDatabase& param_db);
    
    void update_matrix(const BlockMatrix& matrix);
    
    void solve(const BlockVector& rhs, BlockVector& solution);
    
    const ParameterDatabase& get_db();
    
  protected:
    
    ParameterDatabase db;
    
    std::unique_ptr<DirectSolver> direct_solver;
};

#endif // __SOLVER_H__
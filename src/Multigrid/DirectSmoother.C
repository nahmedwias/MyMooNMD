/**
 * @file DirectSmoother.C
 *
 * @date 2016/10/11
 * @author Clemens Bartsch
 */

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <DirectSmoother.h>

DirectSmoother::DirectSmoother()
: solver_()
{
  ;
}

void DirectSmoother::smooth(const BlockVector& rhs, BlockVector& solution)
{
  solver_->solve(rhs, solution);
}

void DirectSmoother::update(const BlockFEMatrix& matrix)
{
  DirectSolver::DirectSolverTypes type = DirectSolver::DirectSolverTypes::umfpack;
  solver_.reset(new DirectSolver(matrix, type));
}

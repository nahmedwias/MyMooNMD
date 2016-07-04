/**
 * @file DirectSmoother.C
 *
 * @date 2016/10/11
 * @author Clemens Bartsch
 */

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <DirectSmoother.h>

#ifdef _MPI
#include <MumpsWrapper.h>
#endif

DirectSmoother::DirectSmoother()
: solver_()
#ifdef _MPI
  , comms_()
#endif
{

}

void DirectSmoother::smooth(const BlockVector& rhs, BlockVector& solution)
{
#ifdef _MPI
  solver_->solve(rhs,solution, comms_);
#else
  solver_->solve(rhs, solution);
#endif
}

void DirectSmoother::update(const BlockFEMatrix& matrix)
{
#ifdef _MPI
  solver_.reset(new MumpsWrapper(matrix, matrix.get_communicators()));
  comms_ = matrix.get_communicators();
#else
  DirectSolver::DirectSolverTypes type = DirectSolver::DirectSolverTypes::umfpack;
  solver_.reset(new DirectSolver(matrix, type));
#endif
}

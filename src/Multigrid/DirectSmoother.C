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
//  Output::dash("DirectSmoother::smooth");
//  std::shared_ptr<TMatrix> mat_solver = solver_->get_matrix();
//  Output::print("Norm mat solver: ", mat_solver->GetNorm());
  solver_->solve(rhs, solution);

}

void DirectSmoother::update(const BlockFEMatrix& matrix)
{
//  Output::dash("DirectSmoother::update");
  DirectSolver::DirectSolverTypes type = DirectSolver::DirectSolverTypes::umfpack;
  solver_.reset(new DirectSolver(matrix, type));

//  //CB DEBUG
//  std::shared_ptr<TMatrix> mat_solver = solver_->get_matrix();
//  Output::print("Norm mat solver: ", mat_solver->GetNorm());
//  std::shared_ptr<TMatrix> mat_input = matrix.get_combined_matrix();
//  Output::print("Norm mat input: ", mat_input->GetNorm());
//
//  mat_input->add_scaled(*mat_input.get(), -1.0);
//  Output::print("Norm mat diff: ", mat_input->GetNorm());
//
//  //END DEBUG
}

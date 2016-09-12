/**
 * @file SORSmoother.C
 * Implementation of class SORSmoother.
 *
 * @date 2016/09/09
 * @author Clemens Bartsch
 */

#include <SORSmoother.h>
#include <BlockVector.h>
#include <BlockFEMatrix.h>
#include <Iteration_sor.h>
#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

SORSmoother::SORSmoother(double omega, size_t sor_strat)
: sor_(nullptr), omega_(omega), sor_strat_(sor_strat)
{
}

void SORSmoother::smooth(const BlockVector& rhs, BlockVector& solution)
{
  //MPI: solution and rhs enter as level 3 consistent

  // Calculate current defect.
  BlockVector defect = rhs;
  sor_->get_operator().apply_scaled_add(solution, defect, -1.0);

  //MPI: defect now level 1 consistent

  sor_->apply_smoother(defect, defect); //one sor sweep - like in preconditioning

  solution.add_scaled(defect, 1.0);
  ///MPI: solution now in level 0 consistency
}

void SORSmoother::update(const BlockFEMatrix& matrix)
{
  //Reset the sor object, w.
  sor_.reset(new Iteration_sor<BlockFEMatrix, BlockVector>(matrix, sor_strat_, omega_));
}



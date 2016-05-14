/**
 * @file JacobiSmoother.C
 * Implementation of class JacobiSmoother.
 *
 * @date 2016/05/13
 * @author Clemens Bartsch
 */

#include <JacobiSmoother.h>
#include <BlockVector.h>
#include <BlockFEMatrix.h>
#include <Iteration_jacobi.h>

JacobiSmoother::JacobiSmoother()
: jacobi(nullptr)
{
  ;
}

void JacobiSmoother::smooth(const BlockVector& rhs, BlockVector& solution)
{
  //Apply one step of a jacobi iteration with the matrix stored therein.
  BlockVector defect = rhs; //copy
  jacobi->get_operator().apply_scaled_add(solution, defect, -1.0);


  jacobi->apply(defect, defect);
  solution.add_scaled(defect, 1.0);
}

void JacobiSmoother::update(const BlockFEMatrix& matrix)
{
  //Reset the diagonal matrix stored in the jacobi object
  jacobi.reset(new Iteration_jacobi<BlockFEMatrix, BlockVector>(matrix));
}



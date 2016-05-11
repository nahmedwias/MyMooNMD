/**
 * @file Implementation of class MultigridLevel.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#include <BlockFEMatrix.h>
#include <MooNMD_Io.h>
#include <MultigridLevel.h>

MultigridLevel::MultigridLevel(BlockFEMatrix* matrix,
                               double beta, SmootherCode sm)
:  matrix_(matrix),
   defect_(*matrix, true), rhs_(*matrix,true), solution_(*matrix, false),
   smoother_(nullptr), beta_(beta)
{
  Output::info("MultigridLevel", "Constructed a MultigridLevel object.");

  //Determine which smoother to use and construct the object.
  switch(sm)
  {
    case SmootherCode::DIRECT_SOLVE:
      Output::info("MultigridLevel", "Dummy: DIRECT_SOLVE smoother");
      break;
    case SmootherCode::NODAL_VANKA:
      Output::info("MultigridLevel", "Dummy: NODAL_VANKA smoother");
      break;
    case SmootherCode::CELL_VANKA:
      Output::info("MultigridLevel", "Dummy: CELL_VANKA smoother");
      break;
    case SmootherCode::BATCH_VANKA:
      Output::info("MultigridLevel", "Dummy: BATCH_VANKA smoother");
      break;
    default:
      ErrThrow("Unknown SmootherCode!");
  }
}


void MultigridLevel::apply_smoother(BlockVector& sol) const
{
  smoother_->smooth(rhs_, sol);
}

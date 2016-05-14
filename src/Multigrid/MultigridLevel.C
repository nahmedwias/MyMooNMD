/**
 * @file Implementation of class MultigridLevel.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#include <BlockFEMatrix.h>
#include <DirectSmoother.h>
#include <JacobiSmoother.h>
#include <MooNMD_Io.h>
#include <MultigridLevel.h>
#include <NoSmoother.h>

MultigridLevel::MultigridLevel(BlockFEMatrix* matrix,
                               SmootherCode sm)
:  matrix_(matrix),
   defect_(*matrix, true), residual_(1e10),
   rhs_(*matrix,true), solution_(*matrix, false),
   smoother_(nullptr)
{
  Output::info("MultigridLevel", "Constructed a MultigridLevel object.");

  //Determine which smoother to use and construct the object.
  switch(sm)
  {
    case SmootherCode::DIRECT_SOLVE:
      smoother_ = std::make_shared<DirectSmoother>();
      break;
    case SmootherCode::JACOBI:
      smoother_ = std::make_shared<JacobiSmoother>();
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
    case SmootherCode::NO_SMOOTHER:
      smoother_ = std::make_shared<NoSmoother>();
      break;
    default:
      ErrThrow("Unknown SmootherCode!");
  }
}


void MultigridLevel::apply_smoother()
{
  smoother_->smooth(rhs_, solution_);
}

void MultigridLevel::calculate_defect()
{
  defect_ = rhs_; //copy
  matrix_->apply_scaled_add(solution_, defect_, -1.0);
  residual_ = sqrt(dot(defect_,defect_));
}

void MultigridLevel::update_smoother()
{
  smoother_->update(*matrix_);
}

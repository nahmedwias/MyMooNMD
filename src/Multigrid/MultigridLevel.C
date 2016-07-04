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
#include <ParameterDatabase.h>
#include <VankaSmootherNew.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

MultigridLevel::MultigridLevel(BlockFEMatrix* matrix,
                               SmootherCode sm,
                               const ParameterDatabase& db)
:  matrix_(matrix),
   defect_(*matrix, true), residual_(1e10),
   rhs_(*matrix,true), solution_(*matrix, false),
   smoother_(nullptr)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank =0;
#endif
  Output::info<4>("MultigridLevel", "Constructed a MultigridLevel object. matrix "
               "dimensions: (", matrix->get_n_total_rows(), ",",
               matrix->get_n_total_columns(), "), n_cells ",
               matrix->get_ansatz_space(0,0).GetCollection()->GetN_Cells());

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
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmootherNew>(VankaType::NODAL, damp);
      break;
    }
    case SmootherCode::CELL_VANKA:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmootherNew>(VankaType::CELL, damp);
      break;
    }
    case SmootherCode::BATCH_VANKA:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmootherNew>(VankaType::BATCH, damp);
      break;
    }
    //Vanka smoothers with storage of local systems
    case SmootherCode::NODAL_VANKA_STORE:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmootherNew>(VankaType::NODAL, damp, true);
      break;
    }
    case SmootherCode::CELL_VANKA_STORE:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmootherNew>(VankaType::CELL, damp, true);
      break;
    }
    case SmootherCode::BATCH_VANKA_STORE:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmootherNew>(VankaType::BATCH, damp, true);
      break;
    }
    default:
      ErrThrow("Unknown SmootherCode!");
  }
}


void MultigridLevel::apply_smoother()
{
  smoother_->smooth(rhs_, solution_);

#ifdef _MPI
  std::vector<const TParFECommunicator3D*> comms = matrix_->get_communicators();
  for(size_t bl =0; bl < comms.size(); ++bl)
  {
    comms[bl]->consistency_update(solution_.block(bl), 3); //restore level 3 consistency of solution_
  }
#endif
}

void MultigridLevel::calculate_defect()
{
  defect_ = rhs_;

#ifdef _MPI
  std::vector<const TParFECommunicator3D*> comms = matrix_->get_communicators();
  for(size_t bl =0; bl < comms.size(); ++bl)
  {
    comms[bl]->consistency_update(solution_.block(bl), 2); //restore level 2 consistency of solution_
  }
#endif

  matrix_->apply_scaled_add(solution_, defect_, -1.0);

#ifdef _MPI
  residual_ = defect_.norm_global(comms);
#else
  residual_ = sqrt(dot(defect_,defect_));
#endif
}

void MultigridLevel::update_smoother()
{
  smoother_->update(*matrix_);
}

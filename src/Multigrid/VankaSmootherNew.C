/*
 * VankaSmootherNew.C
 *
 *  Created on: May 16, 2016
 *      Author: bartsch
 */

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <DirectSolver.h>
#include <MooNMD_Io.h>
#include <VankaSmootherNew.h>

#include <memory>

// Experimental Macro to avoid double code.
#ifdef __2D__
#define TFESpaceXD TFESpace2D
#elif __3D__
#define TFESpaceXD TFESpace3D
#endif

//! Default constructor.
VankaSmootherNew::VankaSmootherNew(VankaType type, double damp_factor)
: type_(type), dimension_(0), damp_factor_(damp_factor),
  matrix_global_(nullptr), press_dofs_local_(0), velo_dofs_local_(0)
{

}

void VankaSmootherNew::update(const BlockFEMatrix& matrix)
{
  //Check if matrix looks like saddle point problem (i.e.: all but last block
  // row space are the same)
  size_t n_blocks = matrix.get_n_cell_rows();
  if(n_blocks <= 1)
    ErrThrow("Matrix has too few cell rows for VankaSmoother!");
  const TFESpaceXD* first_space = &matrix.get_test_space(0,0);
  const TFESpaceXD* last_space = &matrix.get_test_space(n_blocks - 1 , 0);
  for(size_t i = 0; i < n_blocks - 1 ; ++i)
  {
    if(&matrix.get_test_space(i,0) != first_space)
    {
      ErrThrow("So far VankaSmoother will only operate for saddle point matrices."
          " This matrix does not look like one! ", i);
    }
  }
  //Find out if spaces changed, and if so: reassort dof batches
  //pressure (when pressure space has changed)
  if(last_space != pressure_space_)
  {
    pressure_space_ = last_space;
    set_up_pressure_batches(*pressure_space_);
  }
  //velocity ( when either space changed)
  if(first_space != velocity_space_ || last_space != pressure_space_)
  {
    velocity_space_ = first_space;
    bool is_transposed{true}; //out-variable for get_block
    std::shared_ptr<const TMatrix> coupling_block
      = matrix.get_block(n_blocks - 1, 0, is_transposed); //the block B1
    if(is_transposed)
      ErrThrow("Oy vey! That coupling block is transposed, "
          "that's not what I expected.")
    set_up_velocity_batches(*coupling_block.get(), *velocity_space_);
  }

  //Reset the stored global matrix on which all the work is done
  dimension_ = n_blocks - 1;
  matrix_global_=matrix.get_combined_matrix();
}

//The implementation of the smoothing step is very procedural in nature.
void VankaSmootherNew::smooth(const BlockVector& rhs, BlockVector& solution )
{
  if(rhs.n_blocks() != dimension_ + 1)
    ErrThrow("VankaSmoother: rhs dimension does not fit!");

  if(solution.n_blocks() != dimension_ + 1)
    ErrThrow("VankaSmoother: solution dimension does not fit!");

  //loop over all local systems
  for(size_t i = 0; i < press_dofs_local_.size() ; ++i)
  {
    const DofBatch& velo_dofs = velo_dofs_local_.at(i);
    const DofBatch& press_dofs = press_dofs_local_.at(i);
    size_t n_velo_dofs_local = velo_dofs.getSize();
    size_t n_press_dofs_local = press_dofs.getSize();

    size_t n_velo_dofs_global = this->velocity_space_->GetN_DegreesOfFreedom();

    size_t size_local = dimension_* n_velo_dofs_local + n_press_dofs_local;

    //read: all_dofs[i] = global dof corresponding to local dof i (added over all block rows)
    std::vector<int> dof_map(size_local, 0);

    //local solution and rhs
    std::vector<double> rhs_local(size_local, 0.0);
    std::vector<double> solution_local(size_local, 0.0);

    //These vectors will form the local matrix.
    std::vector<int> rowptr_local;
    std::vector<int> kcol_local;
    std::vector<double> entries_local;

    //reserve space for local matrices
    rowptr_local.reserve(size_local + 1);
    kcol_local.reserve(size_local*size_local); //this is too much space! local matrix will be sparse, too
    entries_local.reserve(size_local*size_local);

    /* ******** Fill dof map. *********** */
    for(size_t k = 0; k < size_local; ++k)
    {
      if(k < dimension_* n_velo_dofs_local) //is a velo dof
      {
        int k_modulo_n_velo_dofs_local = k%n_velo_dofs_local;
        int block = k/n_velo_dofs_local; //integer division
        dof_map.at(k) = block * n_velo_dofs_global
            + velo_dofs.getDof(k_modulo_n_velo_dofs_local);
      }
      else //is a pressure dof
      {
        int block = dimension_;
        int k_minus_dim_times_velo = k - dimension_*n_velo_dofs_local;
        dof_map.at(k) = block * n_velo_dofs_global
            + press_dofs.getDof(k_minus_dim_times_velo);
      }
    }

    /* ******** Set up local matrix *********** */
    rowptr_local.push_back(0);
    for(size_t dof_loc = 0; dof_loc < size_local; ++dof_loc) //loop over all local rows
    {
      size_t dof_glo = dof_map.at(dof_loc); //the corresponding global dof

      int begin_r_glo = matrix_global_->GetRowPtr()[dof_glo];
      int end_r_glo = matrix_global_->GetRowPtr()[dof_glo + 1];
      size_t n_entries_in_row_local = 0;

      size_t c_loc = 0;// start with the 0th local column - exploit that
                       // the global KCol Array and the dof_map array are sorted

      for( int i = begin_r_glo ; i < end_r_glo ; ++i )
      {
        double entry = matrix_global_->GetEntries()[i];
        if(entry != 0) //don't copy zeroes.
        {
          int c_glo = matrix_global_->GetKCol()[i];
          //find out what the local column is
          while(dof_map.at(c_loc) < c_glo)
          {
            ++c_loc;
            if(c_loc == size_local)
              break;
          }
          if(c_loc == size_local)
            break; //break loop, we're behind the end!
          if(dof_map.at(c_loc) == c_glo) //the dof c_glo is of interest for the local system
          {
            kcol_local.push_back(c_loc);
            entries_local.push_back(entry);
            //count up the number of entries in the row
            ++n_entries_in_row_local;
          }
          //else just go on
        }
      }
      rowptr_local.push_back(rowptr_local.back() + n_entries_in_row_local);

    }
    std::shared_ptr<TStructure> structure_local
    = std::make_shared<TStructure>(size_local, entries_local.size(),
                                   &kcol_local.at(0), &rowptr_local.at(0));
    std::shared_ptr<TMatrix> matrix_local
    = std::make_shared<TMatrix>(structure_local);
    matrix_local->setEntries(entries_local);


    /* ******** Set up local right hand side. *********** */
    /* (Local right hand side is global defect in local rows) */
    for(size_t dof_loc = 0; dof_loc < size_local; ++dof_loc) //loop over all local rows
    {

      size_t dof_glo = dof_map.at(dof_loc); //the corresponding global dof
      int begin_r_glo = matrix_global_->GetRowPtr()[dof_glo];
      int end_r_glo = matrix_global_->GetRowPtr()[dof_glo + 1];

      double temp = rhs.get_entries()[dof_glo];
      for( int i = begin_r_glo ; i < end_r_glo ; ++i )
      {
        int c_glo = matrix_global_->GetKCol()[i];
        double matrix_glo_entry = matrix_global_->GetEntries()[i];
        temp -= matrix_glo_entry * solution.get_entries()[c_glo];
      }

      rhs_local.at(dof_loc) = temp;
    }

    /* ******** Solve the local system with UMFPACK. *********** */
    DirectSolver ds(matrix_local, DirectSolver::DirectSolverTypes::umfpack);
    ds.solve(&rhs_local.at(0), &solution_local.at(0));

    /* ******** Add damped local solution to global solution. *********** */
    for(size_t dof_loc = 0; dof_loc < size_local; ++dof_loc) //loop over all local rows
    {
      double damp = damp_factor_;
      size_t dof_glo = dof_map.at(dof_loc); //the corresponding global dof
      solution.get_entries()[dof_glo] += damp*solution_local.at(dof_loc);
    }
  }
}

/*!
 * @brief Assembling of the Pressure Dofs.
 *
 * @param pressureSpacePtr A pointer to the pressure space the dofs refer to.
 *
 * This method determines the size, the setup and the order of the pressure dof batches.
 * Playing with it renders different Vanka smoothers!
 */
void VankaSmootherNew::set_up_pressure_batches(const TFESpace& pressureSpace){
  switch (type_) {
    case (VankaType::NODAL):{ // pressure node oriented Vanka
      //There are as many batches as cells in this case.
      int nBatches = pressureSpace.GetN_DegreesOfFreedom();
      //Loop over cells.
      for (int i=0;i<nBatches;i++){
        //To every cell belongs a pressure batch. Create that batch here.
        DofBatch currentPressureBatch{};
        currentPressureBatch.addDof(i);
        //Make the pressure batch nice and clean and copy it into the list.
        currentPressureBatch.tidyUp();
        press_dofs_local_.push_back(currentPressureBatch);
      }
      //End loop over cells
    } break;

    case (VankaType::CELL):
    case (VankaType::BATCH):{ //cell and cellbatch oriented Vanka treat the pressure dofs the same.

      //There are as many batches as cells in this case.
      int nBatches = pressureSpace.GetN_Cells();

      // Store a pointer to the array of globl dofs and the array of begin indices (for each cell).
      int* allDofArray = pressureSpace.GetGlobalNumbers();
      int* beginIndexArray = pressureSpace.GetBeginIndex();

      //Loop over cells.
      for (int i=0;i<nBatches;i++){
        //To every cell belongs a pressure batch. Create that batch here.
        DofBatch currentPressureBatch{};
        //Loop over cell dofs.
        for (int j=beginIndexArray[i];j<beginIndexArray[i+1];j++){
          //Put the current dof into current pressure batch
          currentPressureBatch.addDof(allDofArray[j]);
        }
        // End loop over cell dofs.

        //Make the pressure batch nice and clean and copy it into the list.
        currentPressureBatch.tidyUp();
        press_dofs_local_.push_back(currentPressureBatch);
      }
      //End loop over cells
    }
    break;

    default: {
      Error("Unknown or unimplemented Vanka smoother type! " << endl);
      exit(-1);
    } break;


  }
}

/*!
 * @brief Assembling of the Velocity Dofs.
 *
 * @param[in]   pressureVelocityMatrix A matrix block from which the velo-pressure coupling
 *      can be determined (via non-zero entries).
 * @param[in] velocitySpace the velocity space which is "needed" for by-the-book assembling in CELL case.
 *
 * Look for non-zero entries in the matrix pressureVelocityMatrix.
 * To each row (crsp. pressure dof) put the columns (crsp. velo dofs)
 * into the corresponding velocity dof batch.
 *
 * This method is the same for all nodal vankas, but think carefully about which matrix block to pass in
 * the different NSTYPE cases.
 */
void VankaSmootherNew::set_up_velocity_batches(const TMatrix& pressureVelocityMatrix,
                                            const TFESpace& velocitySpace){
  switch (type_) {
    case VankaType::NODAL:
    case VankaType::BATCH: //nodal and cellbatch Vanka treat assorting of velocity batches the same.
    {
      // Hold pointers for lookup in the coupling matrix.
      const int* columnsPointer = pressureVelocityMatrix.GetStructure().GetKCol();
      const int* rowPointer = pressureVelocityMatrix.GetStructure().GetRowPtr();

      // Loop over pressure batches.
      for(auto pBatchesIterator : press_dofs_local_){
        // Construct the new corresponding velocity batch.
        DofBatch currentVeloBatch{};
        //Loop over pressure dofs in current pressure batch
        for (auto pDof : pBatchesIterator){
          //Loop over all the places in columnsPointer, where entries are connected to row pDof.
          for (int i = rowPointer[pDof]; i != rowPointer[pDof+1]; i++){
            //Every column index that appears in the sparsity structure
            //is a velo dof connected to the pressure dof.
            currentVeloBatch.addDof(columnsPointer[i]);
          }

        }
        //End loop over pressure dofs in current pressure batch.

        //Make the velocity batch nice and clean and add a copy to the list.
        currentVeloBatch.tidyUp();
        velo_dofs_local_.push_back(currentVeloBatch);
      }
      // End loop over pressure batches
      break;
    }
    case VankaType::CELL:
    { //for cell vanka gather all velo dofs in one cell into one batch
      //There are as many batches as cells in this case.
      int nBatches = velocitySpace.GetN_Cells();

      // Store a pointer to the array of globl dofs and the array of begin indices (for each cell).
      int* allDofArray = velocitySpace.GetGlobalNumbers();
      int* beginIndexArray = velocitySpace.GetBeginIndex();

      //Loop over cells.
      for (int i=0;i<nBatches;i++){
        //To every cell belongs a pressure batch. Create that batch here.
        DofBatch currentVeloBatch{};
        //Loop over cell dofs.
        for (int j=beginIndexArray[i];j<beginIndexArray[i+1];j++){
          //Put the current dof into current pressure batch
          currentVeloBatch.addDof(allDofArray[j]);
        }
        // End loop over cell dofs.

        //Make the pressure batch nice and clean and copy it into the list.
        currentVeloBatch.tidyUp();
        velo_dofs_local_.push_back(currentVeloBatch);
      }
      //End loop over cells

      break;
    }
    default:
    {
      ErrThrow("Unknown or unimplemented Vanka smoother type! ");
      break;
    }

  }
}

#undef TFESpaceXD


/*!
 * VankaSmoother.C
 *
 * @date Mar 31, 2015
 * @author Clemens Bartsch
 *
 * Implementation of class VankaSmoother declared in VankaSmoother.h.
 */
#include <VankaSmoother.h>

#include <MooNMD_Io.h>
#include <Matrix.h>
#include <FESpace.h>
#include <Database.h>

#include <DofBatch.h>
#include <FortranStyleMatrix.h>

#include <cmath>
/*!
 * @brief Constructor.
 *
 * @param[in] nBatches The number of dof batches.
 * @param[in] level The number of the mg level this acts on.
 *
 * Determines the type of the smoother to be used, depending on TDatabase::ParamDB->SC_SMOOTHER_SADDLE
 * (or TDatabase::ParamDB->SC_COARSE_SMOOTHER_SADDLE if level equals 0).
 * 20, 30 cell Vanka
 * 21, 31 nodal Vanka
 * 22, 32 cell batch Vanka
 *
 * Takes the number of dof batches (which equals the number of local systems) and reserves space
 * for the data lists.
 */
VankaSmoother::VankaSmoother(size_t nBatches, size_t level){
	//determine the type of the Vanka object
	if (level==0){ //this is the coarsest level
		switch (TDatabase::ParamDB->SC_COARSE_SMOOTHER_SADDLE){
		//without storing local matrices
		case 20:
			storesMatrices_=false;
			type_ = VankaType::CELL;
			break;
		case 21:
			type_ = VankaType::NODAL;
			storesMatrices_=false;
			break;
		case 22:
			type_ = VankaType::CELLBATCH;
			storesMatrices_=false;
			break;
			//with storing local matrices
		case 30:
			storesMatrices_=true;
			type_ = VankaType::CELL;
			break;
		case 31:
			storesMatrices_=true;
			type_ = VankaType::NODAL;
			break;
		case 32:
			storesMatrices_=true;
			type_ = VankaType::CELLBATCH;
			break;
		default:
			storesMatrices_=false;
			type_ = VankaType::NONE;
		}
	} else { //this is one of the finer grids
		switch (TDatabase::ParamDB->SC_SMOOTHER_SADDLE){
		//without storing local matrices
		case 20:
			storesMatrices_=false;
			type_ = VankaType::CELL;
			break;
		case 21:
			type_ = VankaType::NODAL;
			storesMatrices_=false;
			break;
		case 22:
			type_ = VankaType::CELLBATCH;
			storesMatrices_=false;
			break;
			//with storing local matrices
		case 30:
			storesMatrices_=true;
			type_ = VankaType::CELL;
			break;
		case 31:
			storesMatrices_=true;
			type_ = VankaType::NODAL;
			break;
		case 32:
			storesMatrices_=true;
			type_ = VankaType::CELLBATCH;
			break;
		default:
			storesMatrices_=false;
			type_ = VankaType::NONE;
		}
	}
	if (type_ != VankaType::NONE){
		// Reserve space for the lists. TODO This requires push_back instead of random access when inserting elements!
		pressureBatches_.reserve(nBatches);
		velocityBatches_.reserve(nBatches);
		if(storesMatrices_){
			localMatrices_.reserve(nBatches);
		} else { // we store at most 1 system a time.
			localMatrices_.reserve(1);
		}
	}
}

/*!
 * @brief Add batch of pressure dofs to the end of the list of pressure batches.
 *
 * @param[in] batch The dof batch to be added.
 */
void VankaSmoother::addPressureBatch(const DofBatch& batch){
	pressureBatches_.push_back(batch);
}

/*!
 * @brief Add batch of velocity dofs to the end of the list of velocity batches.
 *
 * @param[in] batch The dof batch to be added.
 */
void VankaSmoother::addVelocityBatch(const DofBatch& batch){
	velocityBatches_.push_back(batch);
}


/*!
 * @brief Assembling of the Pressure Dofs.
 *
 * @param pressureSpacePtr A pointer to the pressure space the dofs refer to.
 *
 * This method determines the size, the setup and the order of the pressure dof batches.
 * Playing with it renders different Vanka smoothers!
 */
void VankaSmoother::assortPressureBatches(const TFESpace& pressureSpace){
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
			addPressureBatch(currentPressureBatch);
		}
		//End loop over cells
	} break;

	case (VankaType::CELL):
	case (VankaType::CELLBATCH):{ //cell and cellbatch oriented Vanka treat the pressure dofs the same.

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
			addPressureBatch(currentPressureBatch);
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
 * @param[in] 	pressureVelocityMatrix A matrix block from which the velo-pressure coupling
 * 			can be determined (via non-zero entries).
 * @param[in] velocitySpace the velocity space which is "needed" for by-the-book assembling in CELL case.
 *
 * Look for non-zero entries in the matrix pressureVelocityMatrix.
 * To each row (crsp. pressure dof) put the columns (crsp. velo dofs)
 * into the corresponding velocity dof batch.
 *
 * This method is the same for all nodal vankas, but think carefully about which matrix block to pass in
 * the different NSTYPE cases.
 */
void VankaSmoother::assortVelocityBatches(const TMatrix& pressureVelocityMatrix,
		const TFESpace& velocitySpace){
	switch (type_) {
	case VankaType::NODAL:
	case VankaType::CELLBATCH: //nodal and cellbatch Vanka treat assorting of velocity batches the same.
	{
	// Hold pointers for lookup in the coupling matrix.
	const int* columnsPointer = pressureVelocityMatrix.GetStructure().GetKCol();
	const int* rowPointer = pressureVelocityMatrix.GetStructure().GetRowPtr();

	// Loop over pressure batches.
	for(auto pBatchesIterator : pressureBatches_){
		// Construct the new corresponding velocity batch.
		DofBatch currentVelocityBatch{};
		//Loop over pressure dofs in current pressure batch
		for (auto pDof : pBatchesIterator){
			//Loop over all the places in columnsPointer, where entries are connected to row pDof.
			for (int i = rowPointer[pDof]; i != rowPointer[pDof+1]; i++){
				//Every column index that appears in the sparsity structure
				//is a velo dof connected to the pressure dof.
				currentVelocityBatch.addDof(columnsPointer[i]);
			}

		}
		//End loop over pressure dofs in current pressure batch.

		//Make the velocity batch nice and clean and add a copy to the list.
		currentVelocityBatch.tidyUp();
		addVelocityBatch(currentVelocityBatch);
	}
	// End loop over pressure batches
	break;
	}
	case VankaType::CELL:
	{	//for cell vanka gather all velo dofs in one cell into one batch
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
			addVelocityBatch(currentVeloBatch);
		}
		//End loop over cells

	break;
	}
	default:
	{
		Error("Unknown or unimplemented Vanka smoother type! " << endl);
		exit(-1);
		break;
	}

	}
}

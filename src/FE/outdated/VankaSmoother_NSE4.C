/*
 * VankaSmoother_NSE4.C
 *
 *  @date Mar 31, 2015
 *  @author Clemens Bartsch
 *
 *  Implementation of class VankaSmoother_NSE4 declared in VankaSmoother_NSE4.h.
 */

#include <VankaSmoother_NSE4.h>

#include <NSE_MGLevel4.h>
#ifdef __3D__
#include <Matrix3D.h>
#else
#include <Matrix2D.h>
#endif
#include <FortranStyleMatrix.h>
#include <Database.h>

	/*!
	 * @brief Constructor.
	 *
	 * @param[in] multiGridLevel The MG level the smoother acts on.
	 * @param[in] nBatches The expected or known number of batches.
	 * @param[in] levelNo The number of the MG level. Needed
	 * to determine which database parameter to use.
	 */
	VankaSmoother_NSE4::VankaSmoother_NSE4(const TNSE_MGLevel4& multiGridLevel,
			size_t nBatches, size_t levelNo) : VankaSmoother(nBatches, levelNo),
			multiGridLevel_(multiGridLevel)
	{
	}

	/*!
	 * @brief For a constructed Vanka Object set up the dof batches and
	 * (if storage is chosen) the local matrices.
	 */
	void VankaSmoother_NSE4::initialize(){
		// Clear the lists.
		pressureBatches_.clear();
		velocityBatches_.clear();
		localMatrices_.clear();

		// Assort the pressure batches.
		assortPressureBatches(*multiGridLevel_.GetPSpace());

		//Assort the velo batches. This is done by searching through B1 in the NSTYPE 4 case.
		assortVelocityBatches(*(TMatrix*)multiGridLevel_.getB1(),*multiGridLevel_.GetUSpace());

		// Set up the local Matrices for a full Vanka if supposed to store them.
		if(storesMatrices_){
			// Loop over number of batch pairs
			for (size_t j = 0; j != pressureBatches_.size(); ++j){
				setUpLocalMatrix(j);
			}
		}
	}

	/*! @brief This is the method to solve all local systems.
	 *  The solution is directly written into "currentSolution".
	 *
	 *  @param[in,out] currentSolution The current global solution array
	 *  which gets modified by this method.
	 *  @param[in] currentRHS The current global right hand side
	 */
	void VankaSmoother_NSE4::solveLocalSystems(double *currentSolution,
			const double* const currentRhs){
		// loop over local systems
		for(size_t systemIndex = 0; systemIndex < pressureBatches_.size() ;++systemIndex){

			//Set up vector for local defect/local solution.
			#ifdef __3D__
			size_t nDofs = 3*velocityBatches_.at(systemIndex).getSize() + pressureBatches_.at(systemIndex).getSize();
			#else
			size_t nDofs = 2*velocityBatches_.at(systemIndex).getSize() + pressureBatches_.at(systemIndex).getSize();
			#endif
			double* localDefect = new double[nDofs];


			if(storesMatrices_){ //version with stored local matrices
				// put up current local rhs
				setUpLocalRhs(systemIndex, currentSolution,
						currentRhs, localDefect);
				//solve the prepared system, the matrix is already stored
				solveOneLocalSystem(systemIndex,
						currentSolution, localDefect);

			} else { // version where local matrices have to be put up especially
				//Set up matrix and right hand side at once, store matrix in localMatrices_[0] and solve.
				DofBatch velo0(velocityBatches_.at(0));
				DofBatch press0(pressureBatches_.at(0));
				velocityBatches_.at(0)=velocityBatches_.at(systemIndex);
				pressureBatches_.at(0)=pressureBatches_.at(systemIndex);

				//set up local matrix and store it at place "0" in the list
				setUpLocalMatrixAndRhs(0, currentSolution,
						currentRhs, localDefect);
				//solve the prepared system.
				solveOneLocalSystem(0, currentSolution, localDefect);

				//restore original state of the lists
				velocityBatches_.at(0)=velo0;
				pressureBatches_.at(0)=press0;
				localMatrices_.pop_back();
			}

			delete[] localDefect; localDefect=nullptr;
		} //end loop over local systems
	}

	/*! *  @brief Put up the current local system from the current defect and solution and solve it.
	 *
	 *  @param[in] systemIndex The index of the system to solve.
	 *  @param[in, out] globalSolution The current global solution.
	 *  @param[in, out] localDefect Local Defect. Gets modified in the process
	 *  and is not needed anymore afterwards.
	 *
	 *  Performs one local Vanka step, uses the precomputed dof batches
	 *  and local systems. Solves the systems with the LAPACK solver.
	 */
	void VankaSmoother_NSE4::solveOneLocalSystem(size_t systemIndex,
			double* globalSolution, double* localDefect
			) const{

		// Hold some references for convenience.
		const DofBatch& pBatch = pressureBatches_.at(systemIndex); //current pressure batch
		const DofBatch& vBatch = velocityBatches_.at(systemIndex); //to current velocity batch
		size_t nPressureDofs = pBatch.getSize(); //number of dofs in the pressure batch
		size_t nVelocityDofs = vBatch.getSize(); //number of dofs in the velocity batch

		//Pluck apart the solution array
		double* globalSolutionUx = globalSolution;
		double* globalSolutionUy = &globalSolution[multiGridLevel_.GetN_UDOF()];
		#ifdef __3D__
		double* globalSolutionUz = &globalSolution[2*multiGridLevel_.GetN_UDOF()];
		double* globalSolutionP = &globalSolution[3*multiGridLevel_.GetN_UDOF()];
		#else
		double* globalSolutionP = &globalSolution[2*multiGridLevel_.GetN_UDOF()];
		#endif
		//Solve the local system.
		localMatrices_.at(systemIndex).solve(localDefect);

		// Add the local error vector to the global solution at the right places
		double damp = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
		for(size_t i=0; i < nVelocityDofs ;++i){
			//write to global defect ux
			globalSolutionUx[vBatch.getDof(i)] += damp*localDefect[i];
			//write to global defect uy
			globalSolutionUy[vBatch.getDof(i)] += damp*localDefect[nVelocityDofs + i];
			#ifdef __3D__
			//write to global defect uz
			globalSolutionUz[vBatch.getDof(i)] += damp*localDefect[2*nVelocityDofs + i];
			#endif
		}
		for(size_t i=0; i < nPressureDofs ;i++){
			#ifdef __3D__
			//write to global defect p
			globalSolutionP[pBatch.getDof(i)] += damp*localDefect[3*nVelocityDofs + i];
			#else
			//write to global defect p
			globalSolutionP[pBatch.getDof(i)] += damp*localDefect[2*nVelocityDofs + i];
			#endif
		}
	}

	/*! Set up a the local matrix belonging to system systemIndex for NSTYPE 4.
	 *  After construction the matrix is simply pushed back to the localMatrices_ vector.
	 *  So make sure this is the place where you want it to go.
	 *
	 * @param[in] systemIndex
	 */
	void VankaSmoother_NSE4::setUpLocalMatrix(size_t systemIndex){
		//Catch the case that the back of localMatrices_ is not the place the matrix should go.
		if(localMatrices_.size() != systemIndex){
			ErrMsg("Local system index and size of localMatrices_ don't match!")
		}

		//Two references for convenience
		const DofBatch& veloBatch = velocityBatches_.at(systemIndex);
		const DofBatch& presBatch = pressureBatches_.at(systemIndex);
		size_t nLocalVeloDofs = veloBatch.getSize();
		size_t nLocalPresDofs = presBatch.getSize();

		//Construct the local matrix.
		#ifdef __3D__
		FortranStyleMatrix localMatrix(3*nLocalVeloDofs + nLocalPresDofs,
				3*nLocalVeloDofs + nLocalPresDofs);
		#else
		FortranStyleMatrix localMatrix(2*nLocalVeloDofs + nLocalPresDofs,
				2*nLocalVeloDofs + nLocalPresDofs);
		#endif

		//Initialise pointers for all matrices to reduce number of function calls.
		// A blocks
		TMatrix* A11 = multiGridLevel_.getA11(); //A11
		double* entriesA11 = A11->GetEntries();
		TMatrix* A12 = multiGridLevel_.getA12(); //A12
		double* entriesA12 = A12->GetEntries();
		TMatrix* A21 = multiGridLevel_.getA21(); //A21
		double* entriesA21 = A21->GetEntries();
		TMatrix* A22 = multiGridLevel_.getA22(); //A22
		double* entriesA22 = A22->GetEntries();
		#ifdef __3D__
		TMatrix* A13 = multiGridLevel_.getA13(); //A13
		double* entriesA13 = A13->GetEntries();
		TMatrix* A23 = multiGridLevel_.getA23(); //A23
		double* entriesA23 = A23->GetEntries();
		TMatrix* A31 = multiGridLevel_.getA31(); //A31
		double* entriesA31 = A31->GetEntries();
		TMatrix* A32 = multiGridLevel_.getA32(); //A32
		double* entriesA32 = A32->GetEntries();
		TMatrix* A33 = multiGridLevel_.getA33(); //A33
		double* entriesA33 = A33->GetEntries();
		#endif
		int* rowPtrA = A11->GetRowPtr(); //structure is the same
		int* KColA = A11->GetKCol();

		//BT blocks
		TMatrix* B1T = (TMatrix*)multiGridLevel_.getB1T(); //B1T
		double* entriesB1T = B1T->GetEntries();
		TMatrix* B2T = (TMatrix*)multiGridLevel_.getB2T(); //B2T
		double* entriesB2T = B2T->GetEntries();
		#ifdef __3D__
		TMatrix* B3T = (TMatrix*)multiGridLevel_.getB3T(); //B3T
		double* entriesB3T = B3T->GetEntries();
		#endif
		int* rowPtrBT = B1T->GetRowPtr();
		int* KColBT = B1T->GetKCol(); //structure is the same

		//B blocks
		TMatrix* B1 = (TMatrix*)multiGridLevel_.getB1(); //B1
		double* entriesB1 = B1->GetEntries();
		TMatrix* B2 = (TMatrix*)multiGridLevel_.getB2(); //B2
		double* entriesB2 = B2->GetEntries();
		#ifdef __3D__
		TMatrix* B3 = (TMatrix*)multiGridLevel_.getB3(); //B3
		double* entriesB3 = B3->GetEntries();
		#endif
		int* rowPtrB = B1->GetRowPtr();
		int* KColB = B1->GetKCol(); //structure is the same

		//loop for all A blocks
		size_t ansatzIndex = 0; //the index in the ansatz space dof batch
		for (auto blockGlobalRow : veloBatch){ //loop through ansatz batch

			//define pointers which can be used to iterate over the corresponding row of the global matrix
			int* kColSegmentBegin = &KColA[rowPtrA[blockGlobalRow]];
			int* kColSegmentEnd = &KColA[rowPtrA[blockGlobalRow+1]];
			int* kColSegmentIterator = kColSegmentBegin;

			double* entriesA11It = &entriesA11[rowPtrA[blockGlobalRow]];
			double* entriesA12It = &entriesA12[rowPtrA[blockGlobalRow]];
			double* entriesA21It = &entriesA21[rowPtrA[blockGlobalRow]];
			double* entriesA22It = &entriesA22[rowPtrA[blockGlobalRow]];
			#ifdef __3D__
			double* entriesA13It = &entriesA13[rowPtrA[blockGlobalRow]];
			double* entriesA23It = &entriesA23[rowPtrA[blockGlobalRow]];
			double* entriesA31It = &entriesA31[rowPtrA[blockGlobalRow]];
			double* entriesA32It = &entriesA32[rowPtrA[blockGlobalRow]];
			double* entriesA33It = &entriesA33[rowPtrA[blockGlobalRow]];
			#endif

			//loop through test batch (columns)
			size_t testIndex = 0; //the index in the test space dof batch

			for(auto blockGlobalCol : veloBatch){ //loop through test batch

				while  (*kColSegmentIterator < blockGlobalCol && kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
					++kColSegmentIterator;
					++entriesA11It;
					++entriesA12It;
					++entriesA21It;
					++entriesA22It;
					#ifdef __3D__
					++entriesA13It;
					++entriesA23It;
					++entriesA31It;
					++entriesA32It;
					++entriesA33It;
					#endif
				}
				if (kColSegmentIterator == kColSegmentEnd){ //end of segment is reached
					break; //break for loop and work in next row
				} else if (kColSegmentIterator > kColSegmentEnd) { //column is not in sparsity pattern
					// do nothing!
				} else if (*kColSegmentIterator == blockGlobalCol){ // column is in sparsity pattern
					//set all new entries from A blocks in the local matrix
					size_t shift = nLocalVeloDofs;
					localMatrix.setEntry(ansatzIndex, testIndex, *entriesA11It);
					localMatrix.setEntry(ansatzIndex, testIndex + shift, *entriesA12It);
					localMatrix.setEntry(ansatzIndex+ shift, testIndex, *entriesA21It);
					localMatrix.setEntry(ansatzIndex + shift, testIndex + shift, *entriesA22It);
					#ifdef __3D__
					size_t doubleshift = 2*nLocalVeloDofs;
					localMatrix.setEntry(ansatzIndex, testIndex + doubleshift, *entriesA13It);
					localMatrix.setEntry(ansatzIndex + shift, testIndex + doubleshift, *entriesA23It);
					localMatrix.setEntry(ansatzIndex + doubleshift, testIndex, *entriesA31It);
					localMatrix.setEntry(ansatzIndex + doubleshift, testIndex + shift, *entriesA32It);
					localMatrix.setEntry(ansatzIndex + doubleshift, testIndex + doubleshift, *entriesA33It);
					#endif
				}
				++testIndex;

			} //end loop through test batch
			++ansatzIndex;

		}	// end loop through ansatz batch

		//loop for BT blocks
		ansatzIndex = 0; //the index in the ansatz space dof batch
		for (auto blockGlobalRow : veloBatch){ //loop through ansatz batch

			//define pointers which can be used to iterate over the corresponding row of the global matrix
			int* kColSegmentBegin = &KColBT[rowPtrBT[blockGlobalRow]];
			int* kColSegmentEnd = &KColBT[rowPtrBT[blockGlobalRow+1]];
			int* kColSegmentIterator = kColSegmentBegin;

			double* entriesB1TIt = &entriesB1T[rowPtrBT[blockGlobalRow]];
			double* entriesB2TIt = &entriesB2T[rowPtrBT[blockGlobalRow]];
			#ifdef __3D__
			double* entriesB3TIt = &entriesB3T[rowPtrBT[blockGlobalRow]];
			#endif
			//loop through test batch (columns)
			size_t testIndex = 0; //the index in the test space dof batch

			for(auto blockGlobalCol : presBatch){ //loop through test batch

				while  (*kColSegmentIterator < blockGlobalCol && kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
					++kColSegmentIterator;
					++entriesB1TIt;
					++entriesB2TIt;
					#ifdef __3D__
					++entriesB3TIt;
					#endif
				}
				if (kColSegmentIterator == kColSegmentEnd){ //end of segment is reached
					break; //break for loop and work in next row
				} else if (kColSegmentIterator > kColSegmentEnd) { //column is not in sparsity pattern
					// do nothing!
				} else if (*kColSegmentIterator == blockGlobalCol){ // column is in sparsity pattern
					//set new entries in the local matrix
					size_t shift = nLocalVeloDofs;
					size_t doubleshift = 2*nLocalVeloDofs;
					#ifdef __3D__
					size_t tripleshift = 3*nLocalVeloDofs;
					localMatrix.setEntry(ansatzIndex , testIndex + tripleshift, *entriesB1TIt);
					localMatrix.setEntry(ansatzIndex + shift, testIndex + tripleshift, *entriesB2TIt);
					localMatrix.setEntry(ansatzIndex + doubleshift , testIndex + tripleshift, *entriesB3TIt);
					#else
					localMatrix.setEntry(ansatzIndex, testIndex+doubleshift , *entriesB1TIt);
					localMatrix.setEntry(ansatzIndex+shift, testIndex+doubleshift, *entriesB2TIt);
					#endif
				}
				++testIndex;

			} //end loop through test batch
			++ansatzIndex;

		}	// end loop through ansatz batch

		//loop for B blocks
		ansatzIndex = 0; //the index in the ansatz space dof batch
		for (auto blockGlobalRow : presBatch){ //loop through ansatz batch

			//define pointers which can be used to iterate over the corresponding row of the global matrix
			int* kColSegmentBegin = &KColB[rowPtrB[blockGlobalRow]];
			int* kColSegmentEnd = &KColB[rowPtrB[blockGlobalRow+1]];
			int* kColSegmentIterator = kColSegmentBegin;

			double* entriesB1It = &entriesB1[rowPtrB[blockGlobalRow]];
			double* entriesB2It = &entriesB2[rowPtrB[blockGlobalRow]];
			#ifdef __3D__
			double* entriesB3It = &entriesB3[rowPtrB[blockGlobalRow]];
			#endif

			//loop through test batch (columns)
			size_t testIndex = 0; //the index in the test space dof batch

			for(auto blockGlobalCol : veloBatch){ //loop through test batch

				while  (*kColSegmentIterator < blockGlobalCol && kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
					++kColSegmentIterator;
					++entriesB1It;
					++entriesB2It;
					#ifdef __3D__
					++entriesB3It;
					#endif
				}
				if (kColSegmentIterator == kColSegmentEnd){ //end of segment is reached
					break; //break for loop and work in next row
				} else if (kColSegmentIterator > kColSegmentEnd) { //column is not in sparsity pattern
					// do nothing!
				} else if (*kColSegmentIterator == blockGlobalCol){ // column is in sparsity pattern
					//set new entries in the local matrix
					size_t shift = nLocalVeloDofs;
					size_t doubleshift = 2*nLocalVeloDofs;
					#ifdef __3D__
					size_t tripleshift = 3*nLocalVeloDofs;
					localMatrix.setEntry(ansatzIndex+tripleshift, testIndex, *entriesB1It);
					localMatrix.setEntry(ansatzIndex+tripleshift, testIndex + shift, *entriesB2It);
					localMatrix.setEntry(ansatzIndex+tripleshift, testIndex + doubleshift, *entriesB3It);
					#else
					localMatrix.setEntry(ansatzIndex+doubleshift, testIndex, *entriesB1It);
					localMatrix.setEntry(ansatzIndex+doubleshift, testIndex + shift, *entriesB2It);
					#endif
				}
				++testIndex;

			} //end loop through test batch
			++ansatzIndex;

		}	// end loop through ansatz batch

		//Do some kind of a Dirichlet row correction -
		//this is necessary if it is not done in the assembling of the global matrices
		// set a 0 in those rows which come from dirichlet velo dofs
		//loop over all velo dofs
		for (size_t index = 0; index < nLocalVeloDofs ; ++index){
			if(veloBatch.getDof(index) >= multiGridLevel_.GetUSpace()->GetN_ActiveDegrees()){//we deal with a Dirichlet velo dof
				int rowB1T = index;
				int rowB2T = nLocalVeloDofs + index;
				#ifdef __3D__
				int rowB3T = 2*nLocalVeloDofs + index;
				#endif
				//loop through the entire row

				for (size_t column = 0; column < localMatrix.getNColumns(); ++column){
					//correct the row in B1T
					localMatrix.setEntry(rowB1T, column,0);
					//correct the row in B2T
					localMatrix.setEntry(rowB2T, column,0);
					#ifdef __3D__
					localMatrix.setEntry(rowB3T, column,0);
					#endif
				}//end loop through row
				//set the diagonal element 1
				localMatrix.setEntry(rowB1T, rowB1T,1);
				localMatrix.setEntry(rowB2T, rowB2T,1);
				#ifdef __3D__
				localMatrix.setEntry(rowB3T, rowB3T,1);
				#endif
			}
		} //end loop over velo dofs

		// Do the LU factorization.
		localMatrix.decomposeLU();

		//Store the local system by push back to localMatrices_.
		localMatrices_.push_back(localMatrix);
	}

	/*!
	 * @brief Sets up the right hand side for a local system.
	 *
	 * @param[in] systemIndex Index of the local system.
	 * @param[in] globalSolution The current global solution.
	 * @param[in] globalRhs The current global right hand side.
	 * @param[out] localDefect Writes the output (the defect in the rows of the local system) there.
	 */
	void VankaSmoother_NSE4::setUpLocalRhs(size_t systemIndex, const double* const globalSolution,
			const double* const globalRhs,	double* localDefect) const {
		//Two references to the relevant dof batchs for convenience
		const DofBatch& veloBatch = velocityBatches_.at(systemIndex);
		const DofBatch& presBatch = pressureBatches_.at(systemIndex);
		size_t nLocalVeloDofs = veloBatch.getSize();
		size_t nLocalPresDofs = presBatch.getSize();

		//Pluck apart the solutions and rhs array
		const double* const globalSolutionUx = globalSolution;
		const double* const globalSolutionUy = &globalSolution[multiGridLevel_.GetN_UDOF()];
		#ifdef __3D__
		const double* const globalSolutionUz = &globalSolution[2*multiGridLevel_.GetN_UDOF()];
		const double* const globalSolutionP = &globalSolution[3*multiGridLevel_.GetN_UDOF()];
		#else
		const double* const globalSolutionP = &globalSolution[2*multiGridLevel_.GetN_UDOF()];
		#endif

		const double* const globalRhsUx = globalRhs;
		const double* const globalRhsUy = &globalRhs[multiGridLevel_.GetN_UDOF()];
		#ifdef __3D__
		const double* const globalRhsUz = &globalRhs[2*multiGridLevel_.GetN_UDOF()];
		const double* const globalRhsP = &globalRhs[3*multiGridLevel_.GetN_UDOF()];
		#else
		const double* const globalRhsP = &globalRhs[2*multiGridLevel_.GetN_UDOF()];
		#endif

		// Initialise the vector entries for the right hand side.
		for(size_t i=0; i < nLocalVeloDofs; ++i){
			localDefect[i]= globalRhsUx[veloBatch.getDof(i)];
			localDefect[nLocalVeloDofs + i] = globalRhsUy[veloBatch.getDof(i)];
			#ifdef __3D__
			localDefect[2*nLocalVeloDofs + i] = globalRhsUz[veloBatch.getDof(i)];
			#endif
		}
		for(size_t i=0; i < nLocalPresDofs; ++i){
			#ifdef __3D__
			localDefect[3*nLocalVeloDofs + i] = globalRhsP[presBatch.getDof(i)];
			#else
			localDefect[2*nLocalVeloDofs + i] = globalRhsP[presBatch.getDof(i)];
			#endif
		}

		//Initialise pointers for all matrices to reduce number of function calls.
		// A blocks
		TMatrix* A11 = multiGridLevel_.getA11(); //A11
		double* entriesA11 = A11->GetEntries();
		TMatrix* A12 = multiGridLevel_.getA12(); //A12
		double* entriesA12 = A12->GetEntries();
		TMatrix* A21 = multiGridLevel_.getA21(); //A21
		double* entriesA21 = A21->GetEntries();
		TMatrix* A22 = multiGridLevel_.getA22(); //A22
		double* entriesA22 = A22->GetEntries();
		#ifdef __3D__
		TMatrix* A13 = multiGridLevel_.getA13(); //A13
		double* entriesA13 = A13->GetEntries();
		TMatrix* A23 = multiGridLevel_.getA23(); //A23
		double* entriesA23 = A23->GetEntries();
		TMatrix* A31 = multiGridLevel_.getA31(); //A31
		double* entriesA31 = A31->GetEntries();
		TMatrix* A32 = multiGridLevel_.getA32(); //A32
		double* entriesA32 = A32->GetEntries();
		TMatrix* A33 = multiGridLevel_.getA33(); //A33
		double* entriesA33 = A33->GetEntries();
		#endif
		int* rowPtrA = A11->GetRowPtr(); //structure is the same for all
		int* KColA = A11->GetKCol();

		//BT blocks.
		TMatrix* B1T = (TMatrix*)multiGridLevel_.getB1T(); //B1T
		double* entriesB1T = B1T->GetEntries();
		TMatrix* B2T = (TMatrix*)multiGridLevel_.getB2T(); //B2T
		double* entriesB2T = B2T->GetEntries();
		#ifdef __3D__
		TMatrix* B3T = (TMatrix*)multiGridLevel_.getB3T(); //B3T
		double* entriesB3T = B3T->GetEntries();
		#endif
		int* rowPtrBT = B1T->GetRowPtr(); //structure is the same for all
		int* KColBT = B1T->GetKCol();

		//B blocks
		TMatrix* B1 = (TMatrix*)multiGridLevel_.getB1(); //B1
		double* entriesB1 = B1->GetEntries();
		TMatrix* B2 = (TMatrix*)multiGridLevel_.getB2(); //B2
		double* entriesB2 = B2->GetEntries();
		#ifdef __3D__
		TMatrix* B3 = (TMatrix*)multiGridLevel_.getB3(); //B3
		double* entriesB3 = B3->GetEntries();
		#endif
		int* rowPtrB = B1->GetRowPtr(); //structure is the same for all
		int* KColB = B1->GetKCol();

		//loop for all A blocks
		size_t ansatzIndex = 0; //the index in the ansatz space dof batch
		for (auto blockGlobalRow : veloBatch){ //loop through ansatz batch

			//define pointers which can be used to iterate over the corresponding row of the global matrix
			int* kColSegmentBegin = &KColA[rowPtrA[blockGlobalRow]];
			int* kColSegmentEnd = &KColA[rowPtrA[blockGlobalRow+1]];
			int* kColSegmentIterator = kColSegmentBegin;

			double* entriesA11It = &entriesA11[rowPtrA[blockGlobalRow]];
			double* entriesA12It = &entriesA12[rowPtrA[blockGlobalRow]];
			double* entriesA21It = &entriesA21[rowPtrA[blockGlobalRow]];
			double* entriesA22It = &entriesA22[rowPtrA[blockGlobalRow]];
			#ifdef __3D__
			double* entriesA13It = &entriesA13[rowPtrA[blockGlobalRow]];
			double* entriesA23It = &entriesA23[rowPtrA[blockGlobalRow]];
			double* entriesA31It = &entriesA31[rowPtrA[blockGlobalRow]];
			double* entriesA32It = &entriesA32[rowPtrA[blockGlobalRow]];
			double* entriesA33It = &entriesA33[rowPtrA[blockGlobalRow]];
			#endif

			while  (kColSegmentIterator != kColSegmentEnd) {
				//update right hand side (localDefect)
				#ifdef __3D__
				localDefect[ansatzIndex]-= (*entriesA11It*globalSolutionUx[*kColSegmentIterator]
											 +*entriesA12It*globalSolutionUy[*kColSegmentIterator]
											 +*entriesA13It*globalSolutionUz[*kColSegmentIterator]);
				localDefect[nLocalVeloDofs + ansatzIndex]-= (*entriesA21It*globalSolutionUx[*kColSegmentIterator]
															  +*entriesA22It*globalSolutionUy[*kColSegmentIterator]
															  +*entriesA23It*globalSolutionUz[*kColSegmentIterator]);
				localDefect[2*nLocalVeloDofs + ansatzIndex]-= (*entriesA31It*globalSolutionUx[*kColSegmentIterator]
																+*entriesA32It*globalSolutionUy[*kColSegmentIterator]
																+*entriesA33It*globalSolutionUz[*kColSegmentIterator]);
				#else
				localDefect[ansatzIndex]-= (*entriesA11It*globalSolutionUx[*kColSegmentIterator]
											 +*entriesA12It*globalSolutionUy[*kColSegmentIterator]);
				localDefect[nLocalVeloDofs + ansatzIndex]-= (*entriesA21It*globalSolutionUx[*kColSegmentIterator]
															  +*entriesA22It*globalSolutionUy[*kColSegmentIterator]);
				#endif
				//count up the pointers in the segment arrays
				++kColSegmentIterator;
				++entriesA11It;
				++entriesA12It;
				++entriesA21It;
				++entriesA22It;
				#ifdef __3D__
				++entriesA13It;
				++entriesA23It;
				++entriesA31It;
				++entriesA32It;
				++entriesA33It;
				#endif
			}
			++ansatzIndex;
		}	// end loop through ansatz batch

		//Loop for BT blocks
		ansatzIndex = 0; //the index in the ansatz space dof batch
		for (auto blockGlobalRow : veloBatch){ //loop through ansatz batch

			//define pointers which can be used to iterate over the corresponding row of the global matrix
			int* kColSegmentBegin = &KColBT[rowPtrBT[blockGlobalRow]];
			int* kColSegmentEnd = &KColBT[rowPtrBT[blockGlobalRow+1]];
			int* kColSegmentIterator = kColSegmentBegin;

			double* entriesB1TIt = &entriesB1T[rowPtrBT[blockGlobalRow]];
			double* entriesB2TIt = &entriesB2T[rowPtrBT[blockGlobalRow]];
			#ifdef __3D__
			double* entriesB3TIt = &entriesB3T[rowPtrBT[blockGlobalRow]];
			#endif

			while  (kColSegmentIterator != kColSegmentEnd) {
				//update right hand side (localDefect)
				localDefect[ansatzIndex]-= (*entriesB1TIt*globalSolutionP[*kColSegmentIterator]);
				localDefect[nLocalVeloDofs + ansatzIndex]-= (*entriesB2TIt*globalSolutionP[*kColSegmentIterator]);
				#ifdef __3D__
				localDefect[2*nLocalVeloDofs + ansatzIndex]-= (*entriesB3TIt*globalSolutionP[*kColSegmentIterator]);
				#endif
				//count up the pointers in the segment arrays
				++kColSegmentIterator;
				++entriesB1TIt;
				++entriesB2TIt;
				#ifdef __3D__
				++entriesB3TIt;
				#endif
			}
		++ansatzIndex;
		}	// end loop through ansatz batch


		//loop for B blocks
		ansatzIndex = 0; //the index in the ansatz space dof batch
		for (auto blockGlobalRow : presBatch){ //loop through ansatz batch

			//define pointers which can be used to iterate over the corresponding row of the global matrix
			int* kColSegmentBegin = &KColB[rowPtrB[blockGlobalRow]];
			int* kColSegmentEnd = &KColB[rowPtrB[blockGlobalRow+1]];
			int* kColSegmentIterator = kColSegmentBegin;

			double* entriesB1It = &entriesB1[rowPtrB[blockGlobalRow]];
			double* entriesB2It = &entriesB2[rowPtrB[blockGlobalRow]];
			#ifdef __3D__
			double* entriesB3It = &entriesB3[rowPtrB[blockGlobalRow]];
			#endif

			//loop through test batch (columns)
			size_t testIndex = 0; //the index in the test space dof batch

			for(auto blockGlobalCol : veloBatch){ //loop through test batch

				while  (*kColSegmentIterator < blockGlobalCol && kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
					//update right hand side (localDefect), part for B1, B2 (and B3)
					#ifdef __3D__
					localDefect[3*nLocalVeloDofs + ansatzIndex] -= (*entriesB1It*globalSolutionUx[*kColSegmentIterator]
																	 +*entriesB2It*globalSolutionUy[*kColSegmentIterator]
																	 +*entriesB3It*globalSolutionUz[*kColSegmentIterator]);

					#else
					localDefect[2*nLocalVeloDofs + ansatzIndex] -= (*entriesB1It*globalSolutionUx[*kColSegmentIterator]
																+ *entriesB2It*globalSolutionUy[*kColSegmentIterator]);
					#endif

					++kColSegmentIterator;
					++entriesB1It;
					++entriesB2It;
					#ifdef __3D__
					++entriesB3It;
					#endif
				}
				++testIndex;
			} //end loop through test batch
			while  (kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
				#ifdef __3D__
				localDefect[3*nLocalVeloDofs + ansatzIndex] -= (*entriesB1It*globalSolutionUx[*kColSegmentIterator]
																+*entriesB2It*globalSolutionUy[*kColSegmentIterator]
														        +*entriesB3It*globalSolutionUz[*kColSegmentIterator]);
				#else
				localDefect[2*nLocalVeloDofs + ansatzIndex] -= (*entriesB1It*globalSolutionUx[*kColSegmentIterator]
																+ *entriesB2It*globalSolutionUy[*kColSegmentIterator]);
				#endif
				++kColSegmentIterator;
				++entriesB1It;
				++entriesB2It;
				#ifdef __3D__
				++entriesB3It;
				#endif
			}
			++ansatzIndex;
		}	// end loop through ansatz batch


//		//Correction in Dirichlet rows. Commented out, for not doing this seems to run even better.
//		//loop over all velo dofs
//		for (size_t index = 0; index < nLocalVeloDofs ; ++index){
//			if(veloBatch.getDof(index) >= multiGridLevel_.GetUSpace()->GetN_ActiveDegrees()){//we deal with a Dirichlet velo dof
//				int rowB1T = index;
//				int rowB2T = nLocalVeloDofs + index;
//				//Dirichlet correction in right hand side vector (localDefect)
//				localDefect[rowB1T]=0;
//				localDefect[rowB2T]=0;
//				}
//
//		} //end loop over velo dofs
	}

	/*!
	 * @brief Sets up the matrix and right hand side for a local system.
	 * The matrix is concatenated to localMatrices_ by  push_back.
	 *
	 * @param[in] systemIndex Index of the local system.
	 * @param[in] globalSolution The current global solution.
	 * @param[in] globalRhs The current global right hand side.
	 * @param[out] localDefect Writes the output (the defect in the rows of the local system) there.
	 */
	void VankaSmoother_NSE4::setUpLocalMatrixAndRhs(size_t systemIndex,
			const double* const globalSolution, const double* const globalRhs,
			double* localDefect) {

		//Two references to the relevant dof batchs for convenience
		const DofBatch& veloBatch = velocityBatches_.at(systemIndex);
		const DofBatch& presBatch = pressureBatches_.at(systemIndex);
		size_t nLocalVeloDofs = veloBatch.getSize();
		size_t nLocalPresDofs = presBatch.getSize();

		//Pluck apart the solutions and rhs array
		const double* const globalSolutionUx = globalSolution;
		const double* const globalSolutionUy  = &globalSolution[multiGridLevel_.GetN_UDOF()];
		#ifdef __3D__
		const double* const globalSolutionUz = &globalSolution[2*multiGridLevel_.GetN_UDOF()];
		const double* const globalSolutionP = &globalSolution[3*multiGridLevel_.GetN_UDOF()];
		#else
		const double* const globalSolutionP = &globalSolution[2*multiGridLevel_.GetN_UDOF()];
		#endif

		const double* const globalRhsUx = globalRhs;
		const double* const globalRhsUy = &globalRhs[multiGridLevel_.GetN_UDOF()];
		#ifdef __3D__
		const double* const globalRhsUz = &globalRhs[2*multiGridLevel_.GetN_UDOF()];
		const double* const globalRhsP = &globalRhs[3*multiGridLevel_.GetN_UDOF()];
		#else
		const double* const globalRhsP = &globalRhs[2*multiGridLevel_.GetN_UDOF()];
		#endif

		//Construct the local matrix...
		#ifdef __3D__
		FortranStyleMatrix localMatrix(3*nLocalVeloDofs + nLocalPresDofs,
				3*nLocalVeloDofs + nLocalPresDofs);
		#else
		FortranStyleMatrix localMatrix(2*nLocalVeloDofs + nLocalPresDofs,
				2*nLocalVeloDofs + nLocalPresDofs);
		#endif

		// Initialise the vector entries for the right hand side.
		for(size_t i=0; i < nLocalVeloDofs; ++i){
			localDefect[i]= globalRhsUx[veloBatch.getDof(i)];
			localDefect[nLocalVeloDofs + i] = globalRhsUy[veloBatch.getDof(i)];
			#ifdef __3D__
			localDefect[2*nLocalVeloDofs + i] = globalRhsUz[veloBatch.getDof(i)];
			#endif
		}
		for(size_t i=0; i < nLocalPresDofs; ++i){
			#ifdef __3D__
			localDefect[3*nLocalVeloDofs + i] = globalRhsP[presBatch.getDof(i)];
			#else
			localDefect[2*nLocalVeloDofs + i] = globalRhsP[presBatch.getDof(i)];
			#endif
		}

		//Initialise pointers for all matrices to reduce number of function calls.
		// A blocks
		TMatrix* A11 = multiGridLevel_.getA11(); //A11
		double* entriesA11 = A11->GetEntries();
		TMatrix* A12 = multiGridLevel_.getA12(); //A12
		double* entriesA12 = A12->GetEntries();
		TMatrix* A21 = multiGridLevel_.getA21(); //A21
		double* entriesA21 = A21->GetEntries();
		TMatrix* A22 = multiGridLevel_.getA22(); //A22
		double* entriesA22 = A22->GetEntries();
		#ifdef __3D__
		TMatrix* A13 = multiGridLevel_.getA13(); //A13
		double* entriesA13 = A13->GetEntries();
		TMatrix* A23 = multiGridLevel_.getA23(); //A23
		double* entriesA23 = A23->GetEntries();
		TMatrix* A31 = multiGridLevel_.getA31(); //A31
		double* entriesA31 = A31->GetEntries();
		TMatrix* A32 = multiGridLevel_.getA32(); //A32
		double* entriesA32 = A32->GetEntries();
		TMatrix* A33 = multiGridLevel_.getA33(); //A33
		double* entriesA33 = A33->GetEntries();
		#endif
		int* rowPtrA = A11->GetRowPtr(); //structure is the same
		int* KColA = A11->GetKCol();

		//BT blocks
		TMatrix* B1T = (TMatrix*)multiGridLevel_.getB1T(); //B1T
		double* entriesB1T = B1T->GetEntries();
		TMatrix* B2T = (TMatrix*)multiGridLevel_.getB2T(); //B2T
		double* entriesB2T = B2T->GetEntries();
		#ifdef __3D__
		TMatrix* B3T = (TMatrix*)multiGridLevel_.getB3T(); //B3T
		double* entriesB3T = B3T->GetEntries();
		#endif
		int* rowPtrBT = B1T->GetRowPtr();
		int* KColBT = B1T->GetKCol(); //structure is the same

		//B blocks
		TMatrix* B1 = (TMatrix*)multiGridLevel_.getB1(); //B1
		double* entriesB1 = B1->GetEntries();
		TMatrix* B2 = (TMatrix*)multiGridLevel_.getB2(); //B2
		double* entriesB2 = B2->GetEntries();
		#ifdef __3D__
		TMatrix* B3 = (TMatrix*)multiGridLevel_.getB3(); //B3
		double* entriesB3 = B3->GetEntries();
		#endif
		int* rowPtrB = B1->GetRowPtr();
		int* KColB = B1->GetKCol(); //structure is the same

		//loop for all A blocks
		size_t ansatzIndex = 0; //the index in the ansatz space dof batch
		for (auto blockGlobalRow : veloBatch){ //loop through ansatz batch

			//define pointers which can be used to iterate over the corresponding row of the global matrix
			int* kColSegmentBegin = &KColA[rowPtrA[blockGlobalRow]];
			int* kColSegmentEnd = &KColA[rowPtrA[blockGlobalRow+1]];
			int* kColSegmentIterator = kColSegmentBegin;

			double* entriesA11It = &entriesA11[rowPtrA[blockGlobalRow]];
			double* entriesA12It = &entriesA12[rowPtrA[blockGlobalRow]];
			double* entriesA21It = &entriesA21[rowPtrA[blockGlobalRow]];
			double* entriesA22It = &entriesA22[rowPtrA[blockGlobalRow]];
			#ifdef __3D__
			double* entriesA13It = &entriesA13[rowPtrA[blockGlobalRow]];
			double* entriesA23It = &entriesA23[rowPtrA[blockGlobalRow]];
			double* entriesA31It = &entriesA31[rowPtrA[blockGlobalRow]];
			double* entriesA32It = &entriesA32[rowPtrA[blockGlobalRow]];
			double* entriesA33It = &entriesA33[rowPtrA[blockGlobalRow]];
			#endif

			//loop through test batch (columns)
			size_t testIndex = 0; //the index in the test space dof batch
			for(auto blockGlobalCol : veloBatch){ //loop through test batch
				while  (*kColSegmentIterator < blockGlobalCol && kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
				#ifdef __3D__
				localDefect[ansatzIndex]-= (*entriesA11It*globalSolutionUx[*kColSegmentIterator]
											 +*entriesA12It*globalSolutionUy[*kColSegmentIterator]
											 +*entriesA13It*globalSolutionUz[*kColSegmentIterator]);
				localDefect[nLocalVeloDofs + ansatzIndex]-= (*entriesA21It*globalSolutionUx[*kColSegmentIterator]
															  +*entriesA22It*globalSolutionUy[*kColSegmentIterator]
															  +*entriesA23It*globalSolutionUz[*kColSegmentIterator]);
				localDefect[2*nLocalVeloDofs + ansatzIndex]-= (*entriesA31It*globalSolutionUx[*kColSegmentIterator]
																+*entriesA32It*globalSolutionUy[*kColSegmentIterator]
																+*entriesA33It*globalSolutionUz[*kColSegmentIterator]);
				#else
				localDefect[ansatzIndex]-= (*entriesA11It*globalSolutionUx[*kColSegmentIterator]
											 +*entriesA12It*globalSolutionUy[*kColSegmentIterator]);
				localDefect[nLocalVeloDofs + ansatzIndex]-= (*entriesA21It*globalSolutionUx[*kColSegmentIterator]
															  +*entriesA22It*globalSolutionUy[*kColSegmentIterator]);
				#endif
					//count up the pointers in the segment arrays
					++kColSegmentIterator;
					++entriesA11It;
					++entriesA12It;
					++entriesA21It;
					++entriesA22It;
					#ifdef __3D__
					++entriesA13It;
					++entriesA23It;
					++entriesA31It;
					++entriesA32It;
					++entriesA33It;
					#endif
				}
				if (kColSegmentIterator == kColSegmentEnd){ //end of segment is reached
					break; //break for loop and work in next row
				} else if (*kColSegmentIterator == blockGlobalCol){ //local column is in sparsity pattern
					//set all four new entries in the local matrix
					size_t shift = nLocalVeloDofs;
					localMatrix.setEntry(ansatzIndex, testIndex, *entriesA11It);
					localMatrix.setEntry(ansatzIndex, testIndex + shift, *entriesA12It);
					localMatrix.setEntry(ansatzIndex + shift, testIndex , *entriesA21It);
					localMatrix.setEntry(ansatzIndex + shift, testIndex + shift, *entriesA22It);
					#ifdef __3D__
					size_t doubleshift = 2*nLocalVeloDofs;
					localMatrix.setEntry(ansatzIndex, testIndex + doubleshift, *entriesA13It);
					localMatrix.setEntry(ansatzIndex + shift, testIndex + doubleshift, *entriesA23It);
					localMatrix.setEntry(ansatzIndex + doubleshift, testIndex, *entriesA31It);
					localMatrix.setEntry(ansatzIndex + doubleshift, testIndex + shift, *entriesA32It);
					localMatrix.setEntry(ansatzIndex + doubleshift, testIndex + doubleshift, *entriesA33It);
					#endif
				}
				++testIndex;
			} //end loop through test batch
			while  (kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
				//update right hand side (localDefect)
				#ifdef __3D__
				localDefect[ansatzIndex]-= (*entriesA11It*globalSolutionUx[*kColSegmentIterator]
											 +*entriesA12It*globalSolutionUy[*kColSegmentIterator]
											 +*entriesA13It*globalSolutionUz[*kColSegmentIterator]);
				localDefect[nLocalVeloDofs + ansatzIndex]-= (*entriesA21It*globalSolutionUx[*kColSegmentIterator]
															  +*entriesA22It*globalSolutionUy[*kColSegmentIterator]
															  +*entriesA23It*globalSolutionUz[*kColSegmentIterator]);
				localDefect[2*nLocalVeloDofs + ansatzIndex]-= (*entriesA31It*globalSolutionUx[*kColSegmentIterator]
																+*entriesA32It*globalSolutionUy[*kColSegmentIterator]
																+*entriesA33It*globalSolutionUz[*kColSegmentIterator]);
				#else
				localDefect[ansatzIndex]-= (*entriesA11It*globalSolutionUx[*kColSegmentIterator]
											 +*entriesA12It*globalSolutionUy[*kColSegmentIterator]);
				localDefect[nLocalVeloDofs + ansatzIndex]-= (*entriesA21It*globalSolutionUx[*kColSegmentIterator]
															  +*entriesA22It*globalSolutionUy[*kColSegmentIterator]);
				#endif
				//count up the pointers in the segment arrays
				++kColSegmentIterator;
				++entriesA11It;
				++entriesA12It;
				++entriesA21It;
				++entriesA22It;
				#ifdef __3D__
				++entriesA13It;
				++entriesA23It;
				++entriesA31It;
				++entriesA32It;
				++entriesA33It;
				#endif
			}

			++ansatzIndex;

		}	// end loop through ansatz batch

		//loop for BT blocks
		ansatzIndex = 0; //the index in the ansatz space dof batch
		for (auto blockGlobalRow : veloBatch){ //loop through ansatz batch

			//define pointers which can be used to iterate over the corresponding row of the global matrix
			int* kColSegmentBegin = &KColBT[rowPtrBT[blockGlobalRow]];
			int* kColSegmentEnd = &KColBT[rowPtrBT[blockGlobalRow+1]];
			int* kColSegmentIterator = kColSegmentBegin;

			double* entriesB1TIt = &entriesB1T[rowPtrBT[blockGlobalRow]];
			double* entriesB2TIt = &entriesB2T[rowPtrBT[blockGlobalRow]];
			#ifdef __3D__
			double* entriesB3TIt = &entriesB3T[rowPtrBT[blockGlobalRow]];
			#endif

			//loop through test batch (columns)
			size_t testIndex = 0; //the index in the test space dof batch

			for(auto blockGlobalCol : presBatch){ //loop through test batch

				while  (*kColSegmentIterator < blockGlobalCol && kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
					//update right hand side (localDefect)
					localDefect[ansatzIndex]-= (*entriesB1TIt*globalSolutionP[*kColSegmentIterator]);
					localDefect[nLocalVeloDofs + ansatzIndex]-= (*entriesB2TIt*globalSolutionP[*kColSegmentIterator]);
					#ifdef __3D__
					localDefect[2*nLocalVeloDofs + ansatzIndex]-= (*entriesB3TIt*globalSolutionP[*kColSegmentIterator]);
					#endif
					//count up the pointers in the segment arrays
					++kColSegmentIterator;
					++entriesB1TIt;
					++entriesB2TIt;
					#ifdef __3D__
					++entriesB3TIt;
					#endif
				}
				if (kColSegmentIterator == kColSegmentEnd){ //end of segment is reached
					break; //break for loop and work in next row
				} else if (kColSegmentIterator > kColSegmentEnd) { //column is not in sparsity pattern
					// do nothing!
				} else if (*kColSegmentIterator == blockGlobalCol){ // column is in sparsity pattern
					//set new entries in the local matrix
					size_t shift = nLocalVeloDofs;
					size_t doubleshift = 2*nLocalVeloDofs;
					#ifdef __3D__
					size_t tripleshift = 3*nLocalVeloDofs;
					localMatrix.setEntry(ansatzIndex , testIndex + tripleshift, *entriesB1TIt);
					localMatrix.setEntry(ansatzIndex + shift, testIndex + tripleshift, *entriesB2TIt);
					localMatrix.setEntry(ansatzIndex + doubleshift , testIndex + tripleshift, *entriesB3TIt);
					#else
					localMatrix.setEntry(ansatzIndex, testIndex+doubleshift , *entriesB1TIt);
					localMatrix.setEntry(ansatzIndex+shift, testIndex+doubleshift, *entriesB2TIt);
					#endif
				}
				++testIndex;

			} //end loop through test batch
			while  (kColSegmentIterator != kColSegmentEnd) { //"nach-iterieren" for defect claculation
				//update right hand side (localDefect)
				localDefect[ansatzIndex]-= (*entriesB1TIt*globalSolutionP[*kColSegmentIterator]);
				localDefect[nLocalVeloDofs + ansatzIndex]-= (*entriesB2TIt*globalSolutionP[*kColSegmentIterator]);
				#ifdef __3D__
				localDefect[2*nLocalVeloDofs + ansatzIndex]-= (*entriesB3TIt*globalSolutionP[*kColSegmentIterator]);
				#endif
				//count up the pointers in the segment arrays
				++kColSegmentIterator;
				++entriesB1TIt;
				++entriesB2TIt;
				#ifdef __3D__
				++entriesB3TIt;
				#endif
			}
			++ansatzIndex;

		}	// end loop through ansatz batch

		//loop for B blocks
		ansatzIndex = 0; //the index in the ansatz space dof batch
		for (auto blockGlobalRow : presBatch){ //loop through ansatz batch

			//define pointers which can be used to iterate over the corresponding row of the global matrix
			int* kColSegmentBegin = &KColB[rowPtrB[blockGlobalRow]];
			int* kColSegmentEnd = &KColB[rowPtrB[blockGlobalRow+1]];
			int* kColSegmentIterator = kColSegmentBegin;

			double* entriesB1It = &entriesB1[rowPtrB[blockGlobalRow]];
			double* entriesB2It = &entriesB2[rowPtrB[blockGlobalRow]];
			#ifdef __3D__
			double* entriesB3It = &entriesB3[rowPtrB[blockGlobalRow]];
			#endif

			//loop through test batch (columns)
			size_t testIndex = 0; //the index in the test space dof batch

			for(auto blockGlobalCol : veloBatch){ //loop through test batch

				while  (*kColSegmentIterator < blockGlobalCol && kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
					//update right hand side (localDefect), part for B1, B2 (and B3)
					#ifdef __3D__
					localDefect[3*nLocalVeloDofs + ansatzIndex] -= (*entriesB1It*globalSolutionUx[*kColSegmentIterator]
																	 +*entriesB2It*globalSolutionUy[*kColSegmentIterator]
																	 +*entriesB3It*globalSolutionUz[*kColSegmentIterator]);

					#else
					localDefect[2*nLocalVeloDofs + ansatzIndex] -= (*entriesB1It*globalSolutionUx[*kColSegmentIterator]
																+ *entriesB2It*globalSolutionUy[*kColSegmentIterator]);
					#endif

					++kColSegmentIterator;
					++entriesB1It;
					++entriesB2It;
					#ifdef __3D__
					++entriesB3It;
					#endif
				}
				if (*kColSegmentIterator == blockGlobalCol){ // column is in sparsity pattern
					//set new entries in the local matrix
					size_t shift = nLocalVeloDofs;
					size_t doubleshift = 2*nLocalVeloDofs;
					#ifdef __3D__
					size_t tripleshift = 3*nLocalVeloDofs;
					localMatrix.setEntry(ansatzIndex+tripleshift, testIndex, *entriesB1It);
					localMatrix.setEntry(ansatzIndex+tripleshift, testIndex + shift, *entriesB2It);
					localMatrix.setEntry(ansatzIndex+tripleshift, testIndex + doubleshift, *entriesB3It);
					#else
					localMatrix.setEntry(ansatzIndex+doubleshift, testIndex, *entriesB1It);
					localMatrix.setEntry(ansatzIndex+doubleshift, testIndex + shift, *entriesB2It);
					#endif
				}

				++testIndex;
			} //end loop through test batch
			while  (kColSegmentIterator != kColSegmentEnd) { //assume ordering of KCol!
				#ifdef __3D__
				localDefect[3*nLocalVeloDofs + ansatzIndex] -= (*entriesB1It*globalSolutionUx[*kColSegmentIterator]
																+*entriesB2It*globalSolutionUy[*kColSegmentIterator]
																+*entriesB3It*globalSolutionUz[*kColSegmentIterator]);
				#else
				localDefect[2*nLocalVeloDofs + ansatzIndex] -= (*entriesB1It*globalSolutionUx[*kColSegmentIterator]
																+ *entriesB2It*globalSolutionUy[*kColSegmentIterator]);
				#endif
				++kColSegmentIterator;
				++entriesB1It;
				++entriesB2It;
				#ifdef __3D__
				++entriesB3It;
				#endif
			}

			++ansatzIndex;
		}	// end loop through ansatz batch


		//Do some kind of a Dirichlet value correction -
		//this is necessary if it is not done in the assembling of the global matrices
		// set a 0 in those rows of the local analogue of the B1T and B2T blocks which come from dirichlet velo dofs
		//loop over all velo dofs
		for (size_t index = 0; index < nLocalVeloDofs ; ++index){
			if(veloBatch.getDof(index) >= multiGridLevel_.GetUSpace()->GetN_ActiveDegrees()){//we deal with a Dirichlet velo dof
				int rowB1T = index;
				int rowB2T = nLocalVeloDofs + index;
				#ifdef __3D__
				int rowB3T = 2*nLocalVeloDofs + index;
				#endif
				//loop through the entire row
				for (size_t column = 0; column < 2*nLocalVeloDofs +nLocalPresDofs; ++column){
					//correct the row in B1T
					localMatrix.setEntry(rowB1T, column,0);
					//correct the row in B2T
					localMatrix.setEntry(rowB2T, column,0);
					#ifdef __3D__
					localMatrix.setEntry(rowB3T, column,0);
					#endif
				}//end loop through row
				//set the diagonal element 1
				localMatrix.setEntry(rowB1T, rowB1T,1);
				localMatrix.setEntry(rowB2T, rowB2T,1);
				#ifdef __3D__
				localMatrix.setEntry(rowB3T, rowB3T,1);
				#endif
//				//Dirichlet correction in right hand side vector. Note: better convergence without this!
//				localDefect[rowB1T]=0;
//				localDefect[rowB2T]=0;
				}

		} //end loop over velo dofs

		//Do the LU decomposition
		localMatrix.decomposeLU();
		//Store the local system
		localMatrices_.push_back(localMatrix);
	}

/*!
 *
 * NodalVanka.h
 *
 * @date Mar 31, 2015
 * @author Clemens Bartsch
 *
 * Declaration of an abstract base class for Vanka smoothers.
 * These are used in the multigrid solving of large linear systems.
 */

#ifndef VANKASMOOTHER_H_
#define VANKASMOOTHER_H_

#include <MooNMD_Io.h>

#include <DofBatch.h>
#include <FortranStyleMatrix.h>

//Forward declaration
class TMatrix;
class TFESpace;

//! @brief Enum class of Vanka types that can be used.
//! The type "NONE" means a slim object, which is constructed but not needed.
enum class VankaType{CELL, NODAL, CELLBATCH, NONE};

/*!
 * @brief Abstract base class for the Nodal Vanka local smoothing procedure.
 *
 * The class declares and defines all the data fields and all the functions which can be used by all
 * implementations of Vanka smoothers, without regarding which matrix type the smoother will
 * be applied to. Thus, in the implementation of derived classes for every NSE type, particular care must be taken
 * for declaring and defining those methods which depend on the specific NSE type.
 */
class VankaSmoother{

public:
	//An std::vector of "DofBatch"es is renamed to batchVec.
	typedef std::vector<DofBatch> batchVec;

public:
	//! Constructor.
	VankaSmoother(size_t nBatches, size_t level);

	/*! Get number of pressure batches.
	 * @return The number of stored pressure Batches.
	 */
	size_t getNPressureBatches() const {
		return pressureBatches_.size();
	}

	/*! Get number of pressure batches.
	 * @return The number of stored pressure Batches.
	 */
	size_t getNVelocityBatches() const {
		return velocityBatches_.size();
	}

	//! Initialization method.
	virtual void initialize() = 0;

	//! Solve all local systems.
	virtual void solveLocalSystems( double* currentSolution,
			const double* const currentRHS) = 0;

protected:
	/***************
	 * Protected methods to set up the Vanka Object.
	 ***************/
	//!Add a pressure batch.
	void addPressureBatch(const DofBatch& batch);

	//!Add a velocity batch.
	void addVelocityBatch(const DofBatch& batch);

	//! Put together batches of pressure dofs.
	void assortPressureBatches(const TFESpace& pressureSpace);

	//! Put together batches of velocity dofs.
	void assortVelocityBatches(const TMatrix& pressureVelocityMatrix,
			const TFESpace& velocitySpace);

	//! Method to solve one local system.
	virtual void solveOneLocalSystem(size_t systemIndex,
						double* globalSolution, double* localDefect) const = 0;

	//! Set up local matrix and rhs for system nr systemIndex at once.
	virtual void setUpLocalMatrixAndRhs(size_t systemIndex,
			const double* const globalSolution, const double* const globalRhs,
			double* localDefect) = 0;

	//! Assemble a local matrix and store it.
	virtual void setUpLocalMatrix(size_t systemIndex) = 0;

	//! Calculate the defect of a system, but only in currently needed rows.
	virtual void setUpLocalRhs(size_t systemIndex,
			const double* const globalSolution, const double* const globalRhs,
			double* localDefect) const = 0;

protected:
	/***************
	 * Data members.
	 ***************/
	//! A vector containing all pressure dof batches.
	batchVec pressureBatches_;
	//! A vector containing all velocity dof batches.
	batchVec velocityBatches_;

	//! Whether local Matrices get stored or not.
	bool storesMatrices_;
	//! Which Vanka type is used.
	VankaType type_;

	//! List of local matrices in Fortran style.
	std::vector<FortranStyleMatrix> localMatrices_;
};


#endif /* VANKASMOOTHER_H_ */

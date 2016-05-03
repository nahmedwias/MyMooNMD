/*
 * VankaSmoother_NSE4.h
 *
 * @date 2015/06/03
 * @author Clemens Bartsch
 *
 * Declaration of class VankaSmoother_NSE4.
 * Vanka smoother object used in multigrid with NSTYPE 4 matrices.
 *
 */

#ifndef VANKASMOOTHER_NSE4_H_
#define VANKASMOOTHER_NSE4_H_

#include <VankaSmoother.h>

//Forward declaration
class TNSE_MGLevel4;

/*!
 * @brief Class for the Nodal Vanka local smoothing procedure for MooNMD matrices of NSTYPE 4.
 *
 * All data members and all methods which need no specific
 * knowledge of the matrix type are implemented in the base class.
 */
class VankaSmoother_NSE4 : public VankaSmoother {

public:
	//! Constructor.
	VankaSmoother_NSE4(const TNSE_MGLevel4& multiGridLevel,
			size_t nBatches, size_t levelNo);

	//! Initialization method.
	virtual void initialize() override;

	//! Solve all local systems.
	virtual void solveLocalSystems( double* currentSolution,
			const double* const currentRHS) override;

protected:
	//! Method to solve one local system.
	virtual void solveOneLocalSystem(size_t systemIndex,
						double* globalSolution, double* localDefect) const override;

	//! Set up local matrix and rhs for system nr systemIndex at once.
	virtual void setUpLocalMatrixAndRhs(size_t systemIndex,
			const double* const globalSolution, const double* const globalRhs,
			double* localDefect) override;

	//! Assemble a local matrix and store it.
	void setUpLocalMatrix(size_t systemIndex) override;

	//! Calculate the defect of a system, but only in currently needed rows.
	void setUpLocalRhs(size_t systemIndex, const double* const globalSolution,
			const double* const globalRhs,	double* localDefect) const override;

	//! A constant reference to the MGLevel this smoother is working on.
	const TNSE_MGLevel4& multiGridLevel_;

};

#endif /* VANKASMOOTHER_NSE4_H_ */

/*****************************************************************************
 *  @name CoupledReaction.h
 *	@brief Hold the reaction part of a single CDR equation in a coupled CDR system,
 *	so the part where the actual coupling happens.
 *
 *	The class is intended to know the reaction term and its derivative and
 *	be able to evaluate both with given input data, plus be able to assemble
 *	matrices and vectors which hold entries coming from either the derivative or
 *	the reaction function itself.
 *
 *	The classe does not know what it is used for. So clients have to be careful
 *	which assembling function they call, depending on the solving strategy they use.
 *
 *  @date May 8, 2015
 *  @author Clemens Bartsch
 *****************************************************************************/

#ifndef COUPLEDREACTION_H_
#define COUPLEDREACTION_H_

#include <Constants.h>
#include <CDR_2D_System.h>
#include <stdexcept>
#include <SystemMatScalar2D.h>

//Forward declarations.
class TFESpace2D;

class CoupledReaction {
public:

public:

	/*! Constructor, to be used for linearized_decoupled solving strategy.. */
	CoupledReaction(AssembleFctParam2D* rhsAssemblingFct, ParamFct* paramFunction,
			size_t nCoupled, TFESpace2D* const rhsFESpace);

	/*! Standard destructor. */
	~CoupledReaction();

	/*! @brief Prohibit copy construction and copy assignment. */
	CoupledReaction(const CoupledReaction &obj) = delete;
	CoupledReaction& operator=( const CoupledReaction& obj ) = delete;

	/*! Assembles the right hand side for the linearized_decoupled strategy.
	 * @param latestSolutions The TFEFunctions pointer Array to be handed to the aux Object.
	 * @param coupledTerm The Abbildungsvorschrift of the coupling.*/
	void assembleLinearDecoupled(TFEFunction2D** latestSolutions);

	/*! Returns a pointer to the right hand side vector.
	 * @return A pointer to the right hand side. */
	double* getRightHandSide() const{
		return rightHandSide_;
	}

//	/*!
//	 * This nested class contains a trinity of functions that are used in the process of assembling
//	 * the rhs for the "linearized_decoupled" strategy.
//	 * Because these have to be adjusted to each other, we decided to group them in a common namespace.
//	 */
//	struct Linearized_Decoupled{
//		//! Handed to the aux object, responsible for the output order.
//		static ParamFct* ParameterFunction;
//		//! Handed to the discrete form, responsible for arranging the coefficients.
//		static CoeffFct2D* CoefficientFunction;
//		//! Handed to the discrete form, contains the actual assembling direction.
//		static AssembleFctParam2D* AssemblingFunction;
//	};

private:

	/*! @brief The number of variables in the system this coupled term belongs to.
	 * Example: If the system consists of 3 equations for c1, c2, c3
	 * and this is a term F(c2,c3), the number is 3.*/
	size_t nCoupled_;

	AssembleFctParam2D* rhsAssemblingFct_;

	ParamFct* paramFunction_;

	/*! @brief The part of the right hand side vector which is due to coupling.*/
	double* rightHandSide_;

//	/*!
//	 * @brief A list storing Matrices which are due to coupling.
//	 *
//	 * If this list is filled and what with depends on the solving strategy.
//	 */
//	std::vector<TSystemMatScalar2D*> matrices_;

	//! The  FE space the right hand side and the matrix entries refer to.
	TFESpace2D* feSpace_;

	//! Boundary Condition Function and Values ensuring 0 bdry conditions.
	static void BoundCondition(int BdComp, double t, BoundCond &cond)
	{
		cond = DIRICHLET;
	}
	static void BoundValue(int BdComp, double Param, double &value)
	{
		value = 0.0;
	}

};



#endif /* COUPLEDREACTION_H_ */

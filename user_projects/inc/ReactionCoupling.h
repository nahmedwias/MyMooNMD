/*****************************************************************************
 *  @name ReactionCoupling.h
 *	@brief Hold the reaction part of a single CDR equation in a coupled CDR system,
 *	the part where the actual coupling happens.
 *
 *	The class is intended to know the reaction term and its derivative and
 *	be able to evaluate both with given input data, plus be able to assemble
 *	matrices and vectors which hold entries coming from either the derivative or
 *	the reaction function itself.
 *
 *	So far the class is used for the linearized_decoupled solving strategy only.
 *	Implementing other solving strategies (which is a to do) will require modifications
 *	of the interface, too.
 *
 *  @date May 8, 2015; Oct 29, 2015
 *  @author Clemens Bartsch
 *****************************************************************************/

#ifndef ReactionCoupling_H_
#define ReactionCoupling_H_

#include <Constants.h>
#include <CDR_2D_System.h>
#include <stdexcept>
#include <BlockMatrixCD2D.h>

//Forward declaration.
class TFESpace2D;

class ReactionCoupling {
public:

public:

	/** @brief Constructor. So far the interface only supports linearized_decoupled strategy.
	 *	@param[in] strategy The solving strategy, so far only linear_decoupled is supported.
	 *	@param[in] rhsAssemblingFct An assembling function for the right hand side.
	 *	@param[in] paramFunction A parameter function ("in-out function").
	 *	@param[in] nCoupled
	 *	@param[in] rhsFESpace The FE space used for the right hand side.
	 */
	ReactionCoupling(CDR_2D_System::SolvingStrategy strategy,
			AssembleFctParam2D* rhsAssemblingFct, ParamFct* paramFunction,
			size_t nCoupled, const TFESpace2D&  rhsFESpace);

	/** @brief Assembles the right hand side for the linearized_decoupled strategy.
	 * @param latestSolutions The TFEFunctions pointer Array to be handed to the aux Object.
	 * @param coupledTerm The Abbildungsvorschrift of the coupling.
	 * */
	void assembleLinearDecoupled(TFEFunction2D** latestSolutions);


	//Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
	ReactionCoupling(const ReactionCoupling&) = default;

    //! Default move constructor.
	ReactionCoupling(ReactionCoupling&&) = default;

    //! Default copy assignment operator. Performs deep copy.
	ReactionCoupling& operator=(const ReactionCoupling&) = default;

    //! Default move assignment operator
	ReactionCoupling& operator=(ReactionCoupling&&) = default;

    //! Default destructor.
    ~ReactionCoupling() = default;


	//Getter methods.

	const TFESpace2D& getFeSpace() const {
		return feSpace_;
	}

	size_t getCoupled() const {
		return nCoupled_;
	}

	ParamFct* getParamFunction() const {
		return paramFunction_;
	}

	AssembleFctParam2D* getRhsAssemblingFct() const {
		return rhsAssemblingFct_;
	}

	const BlockVector& getRightHandSide() const {
		return rightHandSide_;
	}

protected:

	/*! @brief The number of variables in the system this coupled term belongs to.
	 * Example: If the system consists of 3 equations for c1, c2, c3
	 * and this is a term F(c2,c3), the number is 3.*/
	size_t nCoupled_;

	/*! A function pointer to the assembling function for the right hand side,
	 *  hard coded in the example file. */
	AssembleFctParam2D* rhsAssemblingFct_;

	/*! A function pointer to the param function (only for the right hand side,?)
	 *  hard coded in the example file. */
	ParamFct* paramFunction_;

	/** The  FE space the right hand side and the matrix entries refer to.
	 * 	@note The object is supposed to be managed by the CD class this belongs to,
	 * 	so just hold a const reference.
	 */
	const TFESpace2D& feSpace_;

	//! The part of the right hand side vector which is due to coupling.
	BlockVector rightHandSide_;

private:
	//! @brief Boundary Condition Function and Values ensuring 0 dirichlet bdry conditions.
	static void DirichletBoundCondition(int BdComp, double t, BoundCond &cond)
	{
		cond = DIRICHLET;
	}
	//! @brief Boundary Condition Function and Values ensuring 0 dirichlet bdry conditions.
	static void ZeroBoundValue(int BdComp, double Param, double &value)
	{
		value = 0.0;
	}

};



#endif /* ReactionCoupling_H_ */

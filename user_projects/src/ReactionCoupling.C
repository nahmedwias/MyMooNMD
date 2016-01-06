/*
 * ReactionCoupling.C
 *
 * @brief Implements class ReactionCoupling declared in ReactionCoupling.h.
 *
 * @date: May 15, 2015
 * @author: Clemens Bartsch
 */
#include <cstring> //For prototype "memset"

#include <ReactionCoupling.h>

#include <FESpace2D.h>
#include <Enumerations.h>
#include <DiscreteForm2D.h>
#include <AuxParam2D.h>
#include <Assemble2D.h>

//Constructor.
ReactionCoupling::ReactionCoupling(CoupledCDR_2D::SolvingStrategy strategy, AssembleFctParam2D* rhsAssemblingFct, ParamFct* paramFunction,
		size_t nCoupled, const TFESpace2D& rhsFESpace) :
		nCoupled_(nCoupled), rhsAssemblingFct_(rhsAssemblingFct),
		paramFunction_(paramFunction), feSpace_(rhsFESpace), rightHandSide_(feSpace_.GetN_DegreesOfFreedom())

		{
			switch (strategy) {
			case CoupledCDR_2D::SolvingStrategy::linear_decoupled:{
			  break;
			  }
			default:
				ErrMsg("Unknown or unimplemented solving strategy!");
			  break;
			}


		}

// Assembling routine for the "linearized_decoupled" solution strategy.
void ReactionCoupling::assembleLinearDecoupled(TFEFunction2D** latestSolutions){

	// ****** Start constructing the LocalAssembling2D object ******
	int myN_Terms = 1; //only 1 term to assemble
	std::vector<int> myFESpaceNumber({0}); //for this term space "0" is used
	std::vector<MultiIndex2D> myDerivatives({D00}); //No derivatives of the ansatz functions used.
	std::vector<int> myRowSpace; //No matrices assembled here.
	std::vector<int> myColumnSpace; //No matrices assembled here.
	std::vector<int> myRhsSpace({0}); //just one right hand side is assembled

	CoeffFct2D* myCoeffs = nullptr; // The unused coefficient function
	AssembleFctParam2D* myAssembleParam = rhsAssemblingFct_; // The static assembling function

	ManipulateFct2D* myManipulate = nullptr; //nobody uses the manipulate function
	int myN_Matrices = 0; //this is of almost no importance anyway
	int myN_Rhs = 1; //this is of almost no importance anyway

	int myN_ParamFct = 1; //use only one parameter function
	std::vector<ParamFct*> myParameterFct({paramFunction_}); //The static parameter function.
	std::vector<int> myBeginParameter({0}); //The 1 parameter function begins working at 0.
	int myN_Parameters = nCoupled_; //the number of parameters equals the number of functions in the coupling
	TFEFunction2D **myFEFunctions2D = latestSolutions; //The FE functions to be evaluated are the latest solutions
	int myN_FEValues = nCoupled_; //here all parameters are FE values, and their number equals nCoupled_
	std::vector<int> myFEValue_FctIndex(myN_FEValues);
	for (size_t i =0; i < nCoupled_;++i) myFEValue_FctIndex[i]=i; //each FE_Value[i] comes from FE_Function[i]
	std::vector<MultiIndex2D> myFEValue_MultiIndex(myN_FEValues);
	for (size_t i =0; i < nCoupled_;++i) myFEValue_MultiIndex[i]=D00;  //which is to say, from its "underived" version

	//Construct the LocalAssembling2D object from the above values.
	LocalAssembling2D localAssembler(
			myN_Terms, myDerivatives, myFESpaceNumber, myRowSpace, myColumnSpace, myRhsSpace,
			myCoeffs, myAssembleParam, myManipulate,
			myN_Matrices, myN_Rhs,
			myN_ParamFct, myParameterFct, myBeginParameter, myN_Parameters,
			myFEFunctions2D, myN_FEValues, myFEValue_FctIndex, myFEValue_MultiIndex);

	// ****** End construction LocalAssembling2D object. ******

	// Reset the rhs vector to zero, Assemble2D(...) adds the new values.
	rightHandSide_.reset();

	// The feSpace_ has to be wrapped up as a const TFESpace2D**
	const TFESpace2D* feSpaceWrap[1] = {&feSpace_};
	//The rhs has to be wrapped up as a double**
	double* rhsWrap[1] = {rightHandSide_.get_entries()};
	//Wrap up the zero dirichlet bdry condition and values used for the coupled part.
	BoundCondFunct2D* bdryCondPtrArray[1] = {&DirichletBoundCondition};
	BoundValueFunct2D* bdryDataPtrArray[1] = {&ZeroBoundValue};

	//Kick off the assembling routine.
	Assemble2D(
			1, //Only one fe space involved.
			feSpaceWrap, //The 1 fe space used is stored in rhs vector..
			0, NULL, 0, NULL, //No assembling of matrices (neither rectangular nor square)
			1, // Only one part of the rhs gets assembled here.
			rhsWrap, //This is where to write the assembled values.
			feSpaceWrap, // There is only one fe space to use as rhs space. NOTE: What exactly is the difference to the 2nd argument?
			bdryCondPtrArray, // The Dirichlet bdry conditions
			bdryDataPtrArray, // The Dirichlet bdry values (0)
			localAssembler); //The local assembling object.
}

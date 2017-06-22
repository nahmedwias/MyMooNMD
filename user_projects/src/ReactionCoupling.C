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
ReactionCoupling::ReactionCoupling(
		CoupledCDR_2D::SolvingStrategy strategy,
		AssembleFctParam2D* rhsAssemblingFct,
		ParamFct* paramFunction,
		size_t nCoupled,
		const TFESpace2D& rhsFESpace,
		bool axisymmetric) :
		nCoupled_(nCoupled), rhsAssemblingFct_(rhsAssemblingFct),
		paramFunction_(paramFunction), feSpace_(rhsFESpace),
		rightHandSide_(feSpace_.GetN_DegreesOfFreedom()),
		axisymmetric_(axisymmetric)
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
void ReactionCoupling::assembleLinearDecoupled(
    TFEFunction2D** latestSolutions){

  size_t n_incoming_functions = nCoupled_;
  TFEFunction2D* incoming_functions[n_incoming_functions];

  for(size_t i = 0; i < n_incoming_functions; ++i)
      incoming_functions[i] = latestSolutions[i];


	// ****** Start constructing the LocalAssembling2D object ******
	int myN_Terms = 1; //only 1 evaluated fe function is used ("c", 0th derivative)
	std::vector<int> myFESpaceNumber({0}); //for this evaluated fe function the space "0" is used
	std::vector<MultiIndex2D> myDerivatives({D00}); //of the evaluated fe function 0th derivative is needed.
	std::vector<int> myRowSpace; //empty - no matrices assembled here.
	std::vector<int> myColumnSpace; //empty - no matrices assembled here.
	std::vector<int> myRhsSpace({0}); //just one right hand side is assembled

	CoeffFct2D* myCoeffs = nullptr; // We need no coefficient function here, it's all in the rhsAssemblingFct_
	AssembleFctParam2D* myAssembleParam = rhsAssemblingFct_; // The assembling function

	ManipulateFct2D* myManipulate = nullptr; //nobody uses the manipulate function
	int myN_Matrices = 0;   //no matrices assembled
	int myN_Rhs = 1; 		//1 right hand side vector gets assembled

	int myN_ParamFct = 1; //use only one parameter function
	std::vector<ParamFct*> myParameterFct({paramFunction_}); //The static parameter function.
	std::vector<int> myBeginParameter({0}); //The 1 parameter function begins working at 0.
	int myN_Parameters = n_incoming_functions; //the number of parameters
	if(axisymmetric_)
		myN_Parameters += 1; // in axisymmetric case there is one extra parameter: the radius 'r' (= y value)
	TFEFunction2D **myFEFunctions2D = incoming_functions; //The FE functions to be evaluated are the latest solutions
	int myN_FEValues = n_incoming_functions; //n_incoming_functions of the parameters stem from the evaluation of fe functions
	std::vector<int> myFEValue_FctIndex(myN_FEValues);
	for (int i =0; i < myN_FEValues;++i) myFEValue_FctIndex[i]=i; //for each 'i', FE_Value[i] comes from FE_Function[i]
	std::vector<MultiIndex2D> myFEValue_MultiIndex(myN_FEValues);
	for (int i =0; i < myN_FEValues;++i) myFEValue_MultiIndex[i]=D00;  //which is to say, from its "underived" version

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
	//TODO Shouldn't we rather take the bdry of the feSpace_?
	BoundCondFunct2D* bdryCondPtrArray[1] = {&DirichletBoundCondition};
	BoundValueFunct2D* bdryDataPtrArray[1] = {&ZeroBoundValue};

	//Kick off the assembling routine.
	Assemble2D(
			1, //Only one fe space involved.
			feSpaceWrap, //The 1 fe space used is stored in rhs vector..
			0, NULL, 0, NULL, //No assembling of matrices (neither rectangular nor square)
			1, // Only one rhs vector gets assembled here.
			rhsWrap, //This is where to write the assembled values.
			feSpaceWrap, // There is only one fe space to use as rhs space
			bdryCondPtrArray, // The bdry conditions
			bdryDataPtrArray, // The bdry values (0)
			localAssembler); //The local assembling object.
}

/*
 * CoupledReaction.C
 *
 * @brief Implements class CoupledReaction declared in CoupledReaction.h.
 *
 * @date: May 15, 2015
 * @author: Clemens Bartsch
 */
#include <cstring> //For prototype "memset"

#include <CoupledReaction.h>
#include <CDR_2D_System.h>

#include <FESpace2D.h>
#include <Enumerations.h>
#include <DiscreteForm2D.h>
#include <AuxParam2D.h>
#include <Assemble2D.h>

//Constructor.
/*! Constructor, to be used for linearized_decoupled solving strategy.. */
CoupledReaction::CoupledReaction(AssembleFctParam2D* rhsAssemblingFct, ParamFct* paramFunction,
		size_t nCoupled, TFESpace2D* const rhsFESpace) :
		feSpace_(rhsFESpace), rhsAssemblingFct_(rhsAssemblingFct), paramFunction_(paramFunction),
		nCoupled_(nCoupled)

		{
			// Construct the right hand side vector and fill with zeroes.
			rightHandSide_ = new double[feSpace_->GetN_DegreesOfFreedom()]();

//			//When constructing the matrices - how to deal with boundary conditions here??
//			//Initialize the matrices as system matrices constructed with FESpace.
//			matrices_.resize(nCoupled_);
//			for(size_t i =0; i<nCoupled_ ; ++i){
//				TSystemMatScalar2D* matrix = new TSystemMatScalar2D(feSpace_);
//				matrices_.push_back(matrix);
//			}
		}

//Destructor.
CoupledReaction::~CoupledReaction(){
	delete[] rightHandSide_;
	//for (auto odm : matrices_) delete odm;
}

// Assembling routine for the "linearized_decoupled" solution strategy.
void CoupledReaction::assembleLinearDecoupled(TFEFunction2D** latestSolutions){

	// Construct the LocalAssembling2D object which wraps up the former TDiscreteForm2D and TAuxParam2D
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

	// End construction LocalAssembling2D object.

	// Reset the rhs vector to zero, Assemble2D(...) adds the new values.
	memset(rightHandSide_, 0, feSpace_->GetN_DegreesOfFreedom()*SizeOfDouble);

	// The feSpace_ has to be wrapped up as a TFESpace2D**
	TFESpace2D* feSpaceWrap[1] = {feSpace_};
	//The rhs has to be wrapped up as a double**
	double* rhsWrap[1] = {rightHandSide_};
	//Wrap up bdry condition and values.
	BoundCondFunct2D* bdryCondPtrArray[1] = {&BoundCondition};
	BoundValueFunct2D* bdryDataPtrArray[1] = {&BoundValue};

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

//void CoupledReaction::Linearized_Decoupled::ParameterFunction(double* in, double* out){
//	// Just skip the first two entries - this is where TAuxParam2D->GetParameters() places the x and y value,
//	// and pass on as many values as nCoupled_.
//	for (size_t i = 0; i<CDRDatabase::CURRENT_ASSEMBLING_N_UNKNOWNS ; ++i){
//		out[i] = in [i+2];
//	}
//}
//
//void CoupledReaction::Linearized_Decoupled::CoefficientFunction(
//		int n_points, double *X, double *Y,
//		double **parameters, double **coeffs){
//	//There is no coefficient handling, they come with the couplingTerm_.
//	return;
//}
//
//void CoupledReaction::Linearized_Decoupled::AssemblingFunction(
//		double Mult, double *coeff, double *param,
//		double hK, double **OrigValues, int *N_BaseFuncts,
//		double ***LocMatrices, double **LocRhs)
//{
//	// For convenience: rename the place to write to (Rhs) and the place to read from (Orig,
//	// which suposedly contains values of test functions)
//	double* Rhs = LocRhs[0];
//	double* Orig = OrigValues[0];
//
//	// Calculate the value of the coupling term at the quad point.
//	// This relies on the correct interaction with the ParameterFunction and the
//	// input order of fe functions to the aux object.
//	// Make sure to set the global parameter to the right value before calling.
//	double coupledTerm = CDRDatabase::CURRENT_ASSEMBLING_COUPLING_FUNCTION(param);
//				// coupled term evaluated at the quad point this
//				// AssembleFctParam2D is called upon
//
//	//Loop over all local base functions.
//	for(int i=0;i<N_BaseFuncts[0];i++)
//	{
//		Rhs[i] -= Mult*coupledTerm*Orig[i]; //(Mult contains quad weigth, and maybe even the determinant due to trafo)
//	}
//}
//
////Definition of the static data members of CDRDatabase, so the linker knows them.
//DoubleFunctScalar* CDRDatabase::CURRENT_ASSEMBLING_COUPLING_FUNCTION = nullptr;
//size_t CDRDatabase::CURRENT_ASSEMBLING_N_UNKNOWNS=0;

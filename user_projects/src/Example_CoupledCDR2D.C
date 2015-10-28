#include <Example_CoupledCDR2D.h>

#include <FEDatabase2D.h>
#include <Database.h>
#include <MainUtilities.h>


/* examples */

namespace constant
{
#include "Coupled_CDR_2D/Constant_functions.h"
}

Example_CoupledCDR2D::Example_CoupledCDR2D() : Example2D()
{
	switch( TDatabase::ParamDB->EXAMPLE )
	{
	case 0: //Constant function example.
		using namespace constant;

		/** exact_solution */
		exact_solution.push_back( ExactC1 );
		exact_solution.push_back( ExactC2 );

		/** boundary condition */
		boundary_conditions.push_back( BoundCondition );

		/** boundary values */
		boundary_data.push_back( BoundValueC1 );
		boundary_data.push_back( BoundValueC2 );

		/** coefficient functions */
		bilinCoeffs_.push_back(BilinearCoeffsC1);
		bilinCoeffs_.push_back(BilinearCoeffsC2);

		/** assembling functions */
		rhsAssemblingFunctions_.push_back(AssemblingFunctionC1);
		rhsAssemblingFunctions_.push_back(AssemblingFunctionC2);

		/** parameter functions */
		parameterFunctions_.push_back(ParameterFunction);
		parameterFunctions_.push_back(ParameterFunction);

		/** number of equations */
		nEquations_ = 2;

		ExampleFile();
		break;
	default:
		ErrMsg("Unknown name of the convection-diffusion-reaction (CDR) " <<
				"example!");
		exit(0);
	}
}

/** ************************************************************************ */
Example_CD2D* Example_CoupledCDR2D::getDecoupledExample(size_t n) const
{
	if(n >= nEquations_){
		throw std::runtime_error("There is no example with this index. ");
	}

	if(TDatabase::ParamDB->SC_VERBOSE > 1)
		OutPut("extracting coupled CDR with constant functions\n");

	std::vector <DoubleFunct2D*> exact_coupled;
	std::vector <BoundCondFunct2D*> bc_coupled;
	std::vector <BoundValueFunct2D*> bd_coupled;

	//when using the parts, ignore the exact solution - it's the exact solution of the coupled
	//system, not of any single equation.
	exact_coupled.push_back(exact_solution.at(n)); // n = 0 : "concentration"
	// n = 1 : "temperature"

	bc_coupled.push_back(boundary_conditions.at(0)); // boundary is the same for both

	bd_coupled.push_back(boundary_data.at(n)); // n = 0 : "concentration"
	// n = 1 : "temperature"
	//Construct the example and return a pointer to it.
	return new Example_CD2D(exact_coupled, bc_coupled, bd_coupled,
			bilinCoeffs_.at(n));
}

//! Return function pointer to coefficient function of equation nr. equationIndex.
CoeffFct2D* Example_CoupledCDR2D::getCoeffFct(size_t equationIndex) const{
	return bilinCoeffs_.at(equationIndex);
}

//! Return function pointer to assembling function of equation nr. equationIndex.
AssembleFctParam2D* Example_CoupledCDR2D::getAssemblingFct(size_t equationIndex) const{
	return rhsAssemblingFunctions_.at(equationIndex);
}

//! Return function pointer to parameter function of equation nr. equationIndex.
ParamFct* Example_CoupledCDR2D::getParamFct(size_t equationIndex) const{
	return parameterFunctions_.at(equationIndex);
}

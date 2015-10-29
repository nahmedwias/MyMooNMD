#include <Example_CoupledCDR2D.h>

#include <FEDatabase2D.h>
#include <Database.h>
#include <MainUtilities.h>


// Namespaces for the hard coded examples.
namespace constant
{
#include "Constant_functions.h"
}

// Implementation of member methods
Example_CoupledCDR2D::Example_CoupledCDR2D() : Example2D()
{
	switch( TDatabase::ParamDB->EXAMPLE )
	{
	case 0: //Constant function example.
		using namespace constant;

		/** exact_solution */
		exact_solution.push_back( ExactC1 );
		exact_solution.push_back( ExactC2 );

		/** boundary condition - here it's twice the same*/
		boundary_conditions.push_back( BoundCondition );
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
		ErrMsg("Unknown index of the convection-diffusion-reaction (CDR) " <<
				"example!");
		exit(0);
	}

	// For all examples: fill list of decoupled examples.
	generateDecoupledExamples();

}

const Example_CD2D& Example_CoupledCDR2D::getDecoupledExample(size_t n) const
{
	// return a const ref to the example n, using range checked vector method "at"
	return decoupledExamples_.at(n);
}

CoeffFct2D* Example_CoupledCDR2D::getCoeffFct(size_t equationIndex) const{
	return bilinCoeffs_.at(equationIndex);
}

AssembleFctParam2D* Example_CoupledCDR2D::getAssemblingFct(size_t equationIndex) const{
	return rhsAssemblingFunctions_.at(equationIndex);
}

ParamFct* Example_CoupledCDR2D::getParamFct(size_t equationIndex) const{
	return parameterFunctions_.at(equationIndex);
}

void Example_CoupledCDR2D::generateDecoupledExamples() {

  decoupledExamples_.reserve(nEquations_); //reserve space to avoid reallocations

  for(size_t n = 0; n<nEquations_;++n){

    std::vector <DoubleFunct2D*> exact_coupled;
	std::vector <BoundCondFunct2D*> bc_coupled;
	std::vector <BoundValueFunct2D*> bd_coupled;

	//when using the parts, ignore the exact solution - it's the exact solution of the coupled
	//system, not of any single equation.
	exact_coupled.push_back(exact_solution.at(n));

	// Boundary conditions and boundary data.
	bc_coupled.push_back(boundary_conditions.at(n));
	bd_coupled.push_back(boundary_data.at(n));

	//Construct the example and push back a copy of it.
	decoupledExamples_.push_back(Example_CD2D(exact_coupled, bc_coupled, bd_coupled,
			bilinCoeffs_.at(n)));
	}
}










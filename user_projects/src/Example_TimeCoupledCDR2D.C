/**
 * Implementation of class Example_TimeCoupledCDR2D.
 *
 * @date 2016/06/06
 * @author Clemens Bartsch
 */

#include <Example_TimeCoupledCDR2D.h>
#include <MainUtilities.h>
#include <Database.h>

#include <cmath>


// Namespaces for the hard coded examples.
namespace test_time
{
#include "Test_Time.h"
}

namespace asa_crystallizer
{
#include "ASA_crystallizer.h"
}

/// Zero function R^2 -> R, f must have 4 double allocated.
/// On output: f[0] function value, f[1] to f[3]: x to z derivatives.
void zero_function_2D(double x, double y, double* f)
{
  f[0] = 0;
  f[1] = 0;
  f[2] = 0;
  f[3] = 0;
}

// Implementation of member methods
Example_TimeCoupledCDR2D::Example_TimeCoupledCDR2D(const ParameterDatabase & db)
: Example_NonStationary2D(db)
{
  int example_code = example_database["example"];
  switch( example_code )
  {
    case 0: //cobbled together test example
    {
      using namespace test_time;

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

      /** parameter function - twice the same*/
      parameterFunctions_.push_back(ParameterFunction);
      parameterFunctions_.push_back(ParameterFunction);

      /** initial conditions */
      initialCondition.push_back(InitialConditionC1);
      initialCondition.push_back(InitialConditionC2);

      /** number of equations */
      nEquations_ = 2;

      ExampleFile();
      break;
    }
    case 1: // Example 1 for PhD Thesis Clemens: Eder et al. (2010)
    {
      using namespace asa_crystallizer;

      /** Unknwon exact solutions, put to zero. */
      exact_solution.push_back( zero_function_2D );
      exact_solution.push_back( zero_function_2D );

      /** Boundary conditions.*/
      boundary_conditions.push_back( BoundCond_T );
      boundary_conditions.push_back( BoundCond_C_ASA );

      /** Boundary values */
      boundary_data.push_back( BoundValue_T );
      boundary_data.push_back( BoundValue_C_ASA );

      /** coefficient functions */
      bilinCoeffs_.push_back(Coefficients_T);
      bilinCoeffs_.push_back(Coefficients_C_ASA);

      /** Assembling functions */
      rhsAssemblingFunctions_.push_back(RhsAssemblingFunction_T);
      rhsAssemblingFunctions_.push_back(RhsAssemblingFunction_C_ASA);

      /** parameter function - twice the same*/
      parameterFunctions_.push_back(ParameterFunction);
      parameterFunctions_.push_back(ParameterFunction);

      /** initial conditions */
      initialCondition.push_back(InitialCondition_T);
      initialCondition.push_back(InitialCondition_C_ASA);

      /** number of equations */
      nEquations_ = 2;

      ExampleFile();
      break;
    }
    default:
      ErrMsg("Unknown index of the convection-diffusion-reaction (CDR) " <<
             "example!");
      exit(0);

  }

  // For all examples: fill list of decoupled examples.
  generateDecoupledExamples();

}

const Example_TimeCD2D& Example_TimeCoupledCDR2D::getDecoupledExample(size_t n) const
{
  // return a const ref to the example n, using range checked vector method "at"
  return decoupledExamples_.at(n);
}

CoeffFct2D* Example_TimeCoupledCDR2D::getCoeffFct(size_t equationIndex) const{
  return bilinCoeffs_.at(equationIndex);
}

AssembleFctParam2D* Example_TimeCoupledCDR2D::getAssemblingFct(size_t equationIndex) const{
  return rhsAssemblingFunctions_.at(equationIndex);
}

ParamFct* Example_TimeCoupledCDR2D::getParamFct(size_t equationIndex) const{
  return parameterFunctions_.at(equationIndex);
}

CoeffFct2D* Example_TimeCoupledCDR2D::get_coeffs() const
{
  ErrThrow("Do not use get_coeffs for Example_TimeCoupledCDR2D. "
      "That class stores more than one coefficient function, "
      "use getCoeffFct(size_t equationIndex) instead.");
  return problem_coefficients; //is nullptr, just get rid of a compiler warning
}

void Example_TimeCoupledCDR2D::generateDecoupledExamples() {

  decoupledExamples_.reserve(nEquations_); //reserve space to avoid reallocations

  for(size_t n = 0; n<nEquations_;++n){

    std::vector <DoubleFunct2D*> exact_coupled;
    std::vector <BoundCondFunct2D*> bc_coupled;
    std::vector <BoundValueFunct2D*> bd_coupled;
    std::vector <DoubleFunct2D*> init_cond;

    // Exact solution, is just set to zero
    exact_coupled.push_back( zero_function_2D );

    // Boundary conditions and boundary data.
    bc_coupled.push_back(boundary_conditions.at(n));
    bd_coupled.push_back(boundary_data.at(n));

    //Initial conditions.
    init_cond.push_back(initialCondition.at(n));

    //Construct the example and push back a copy of it.
    decoupledExamples_.push_back(
        Example_TimeCD2D(exact_coupled, bc_coupled, bd_coupled,
                     bilinCoeffs_.at(n), true, true, init_cond ));
    // TODO Think whether "true, true" is always correct here,
    // or not better the example should be asked!
  }
}












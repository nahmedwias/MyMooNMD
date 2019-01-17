#include "Example_CoupledNS_Stress.h"
#include "Database.h"
#include "MainUtilities.h"

namespace simp_ns_stress
{
 #include "simple_example_stressNS.h"
}

Example_CoupledNS_Stress::Example_CoupledNS_Stress(const ParameterDatabase& param_db)
:Example2D(param_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case 0:
    /** exact_solution */
    exact_solution.push_back( simp_ns_stress::Stress_XX );
    exact_solution.push_back( simp_ns_stress::Stress_XY );
    exact_solution.push_back( simp_ns_stress::Stress_YY );
    exact_solution.push_back( simp_ns_stress::ExactU1 );
    exact_solution.push_back( simp_ns_stress::ExactU2 );
    exact_solution.push_back( simp_ns_stress::ExactP );
    
    /** boundary condition */
    boundary_conditions.push_back( BoundConditionNoBoundCondition );
    boundary_conditions.push_back( BoundConditionNoBoundCondition );
    boundary_conditions.push_back( BoundConditionNoBoundCondition );
    boundary_conditions.push_back( simp_ns_stress::BoundConditionNS );
    boundary_conditions.push_back( simp_ns_stress::BoundConditionNS );    
    boundary_conditions.push_back( BoundConditionNoBoundCondition );
    
    /** boundary values */
    boundary_data.push_back( BoundaryValueHomogenous );
    boundary_data.push_back( BoundaryValueHomogenous );
    boundary_data.push_back( BoundaryValueHomogenous );
    boundary_data.push_back( simp_ns_stress::U1BoundValue );
    boundary_data.push_back( simp_ns_stress::U2BoundValue );
    boundary_data.push_back( BoundaryValueHomogenous );
    
    /** coefficients */
    problem_coefficients = simp_ns_stress::LinCoeffs;
    
    // Set dimensionless viscosity
    double nu = this->example_database["reynolds_number"];    
    simp_ns_stress::DIMENSIONLESS_VISCOSITY = 1./nu;
    
    // boundary conditions
    simp_ns_stress::ExampleFile();
    break;
  }
}

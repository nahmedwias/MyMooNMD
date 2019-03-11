#include "Example_NSPR_NSE2D.h"
#include <Database.h>
#include <FEDatabase2D.h>
#include <SquareMatrix2D.h>
#include <string.h>
#include <MainUtilities.h>

namespace pot_flow_1
{
#include "potential_flow_ex1.h"
}

Example_NSPR_NSE2D::Example_NSPR_NSE2D(const ParameterDatabase& param_db) : Example2D(param_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case 0:
    {
       /** exact_solution */
      exact_solution.push_back( pot_flow_1::ExactU1 );
      exact_solution.push_back( pot_flow_1::ExactU2 );
      exact_solution.push_back( pot_flow_1::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( pot_flow_1::BoundCondition );
      boundary_conditions.push_back( pot_flow_1::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( pot_flow_1::U1BoundValue );
      boundary_data.push_back( pot_flow_1::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = pot_flow_1::LinCoeffs;
      
      pot_flow_1::ExampleFile();
    }
    break;
  }
}

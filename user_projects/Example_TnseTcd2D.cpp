#include <Example_TnseTcd2D.h>
#include <Database.h>
#include <MainUtilities.h>
#include <FEDatabase2D.h>
#include "TimeNavierStokes.h"
#include <Example_TimeCD2D.h>
#include "AuxParam2D.h" //used in MixingLayerSlipSmallSquares.h
#include <string>
#include <FEFunction2D.h>

namespace lsp_tnse_tcd    
{
 #include "linear_space_time_mphase.h"
}

Example_TnseTcd2D::Example_TnseTcd2D(
  const ParameterDatabase& user_input_parameter_db)
 : Example_NonStationary2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case 1:
      /** exact_solution */
      exact_solution.push_back( lsp_tnse_tcd::ExactU1 );
      exact_solution.push_back( lsp_tnse_tcd::ExactU2 );
      exact_solution.push_back( lsp_tnse_tcd::ExactP );
      exact_solution.push_back( lsp_tnse_tcd::ExactT );
      
      /** boundary condition */
      boundary_conditions.push_back( lsp_tnse_tcd::BoundCondition );
      boundary_conditions.push_back( lsp_tnse_tcd::BoundCondition );
      boundary_conditions.push_back( lsp_tnse_tcd::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( lsp_tnse_tcd::U1BoundValue );
      boundary_data.push_back( lsp_tnse_tcd::U2BoundValue );
      boundary_data.push_back( lsp_tnse_tcd::TBoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = lsp_tnse_tcd::LinCoeffs;
      
      initialCondition.push_back(lsp_tnse_tcd::InitialU1);
      initialCondition.push_back(lsp_tnse_tcd::InitialU2);
      initialCondition.push_back(lsp_tnse_tcd::InitialT);
      
      lsp_tnse_tcd::ExampleFile();
      
     this->timeDependentRhs = lsp_tnse_tcd::rhs_depends_on_time;
     this->timeDependentCoeffs=lsp_tnse_tcd::coefficients_depend_on_time;
     break;
      
         default:
      ErrThrow("Unknown name of the transient Navier-Stokes & convection-diffusion example!", 
               example_code);
      
  }
}

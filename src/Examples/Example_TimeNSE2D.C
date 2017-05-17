#include <Example_TimeNSE2D.h>
#include <Database.h>
#include <MainUtilities.h>
#include <FEDatabase2D.h>
#include <Time_NSE2D.h>

#include <string>

namespace bsp1
{
 #include "TNSE_2D/Bsp1.h"
}

namespace lin_space_time
{
#include "TNSE_2D/linear_space_time.h"
}

namespace sincosexp
{
#include "TNSE_2D/SinCosExp.h"
}

namespace flow_around_cylinder_steady_inflow
{
#include "flow_around_cylinder_steady_inflow.h"
}

Example_TimeNSE2D::Example_TimeNSE2D(
  const ParameterDatabase& user_input_parameter_db)
 : Example_NonStationary2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case 0:
      /** exact_solution */
      exact_solution.push_back( bsp1::ExactU1 );
      exact_solution.push_back( bsp1::ExactU2 );
      exact_solution.push_back( bsp1::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( bsp1::BoundCondition );
      boundary_conditions.push_back( bsp1::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( bsp1::U1BoundValue );
      boundary_data.push_back( bsp1::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = bsp1::LinCoeffs;
      
      /** initial condition */
      initialCondition.push_back(bsp1::InitialU1);
      initialCondition.push_back(bsp1::InitialU2);
      bsp1::ExampleFile();
      break;
    case 1:
      /** exact_solution */
      exact_solution.push_back( lin_space_time::ExactU1 );
      exact_solution.push_back( lin_space_time::ExactU2 );
      exact_solution.push_back( lin_space_time::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( lin_space_time::U1BoundValue );
      boundary_data.push_back( lin_space_time::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = lin_space_time::LinCoeffs;
      
      initialCondition.push_back(lin_space_time::InitialU1);
      initialCondition.push_back(lin_space_time::InitialU2);
      
      lin_space_time::ExampleFile();
      break;
      
    case 2: // SinCosExp
      /** exact_solution */
      exact_solution.push_back(sincosexp::ExactU1 );
      exact_solution.push_back(sincosexp::ExactU2 );
      exact_solution.push_back(sincosexp::ExactP );

      /** boundary condition */
      boundary_conditions.push_back(sincosexp::BoundCondition );
      boundary_conditions.push_back(sincosexp::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back(sincosexp::U1BoundValue );
      boundary_data.push_back(sincosexp::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients =sincosexp::LinCoeffs;

      initialCondition.push_back(sincosexp::InitialU1);
      initialCondition.push_back(sincosexp::InitialU2);

     sincosexp::ExampleFile();
    
    case 3:
      /** exact_solution */
      exact_solution.push_back( flow_around_cylinder_steady_inflow::ExactU1 );
      exact_solution.push_back( flow_around_cylinder_steady_inflow::ExactU2 );
      exact_solution.push_back( flow_around_cylinder_steady_inflow::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( flow_around_cylinder_steady_inflow::BoundCondition );
      boundary_conditions.push_back( flow_around_cylinder_steady_inflow::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( flow_around_cylinder_steady_inflow::U1BoundValue );
      boundary_data.push_back( flow_around_cylinder_steady_inflow::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = flow_around_cylinder_steady_inflow::LinCoeffs;
      
      initialCondition.push_back(flow_around_cylinder_steady_inflow::InitialU1);
      initialCondition.push_back(flow_around_cylinder_steady_inflow::InitialU2);
      
      // Set dimensionless viscosity
      flow_around_cylinder_steady_inflow::DIMENSIONLESS_VISCOSITY = get_nu();

      /**post processing - drag and lift calculation and output */
      post_processing_stat = flow_around_cylinder_steady_inflow::compute_drag_lift_pdiff;

      flow_around_cylinder_steady_inflow::ExampleFile();
      break;
    default:
      ErrThrow("Unknown time-dependent Example_TimeNSE2D example!");
  }
}

void Example_TimeNSE2D::do_post_processing(Time_NSE2D& tnse2d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(tnse2d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_TimeNSE2D","No post processing done for the current example.");
  }
}

double Example_TimeNSE2D::get_nu() const
{
  double inverse_reynolds = this->example_database["reynolds_number"];
  inverse_reynolds = 1/inverse_reynolds;
  return inverse_reynolds;
}

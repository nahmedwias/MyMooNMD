#include <Example_TimeNSE2D.h>
#include <Database.h>
#include <MainUtilities.h>
#include <Time_NSE2D.h>
#include <FEDatabase2D.h>

#include <string>

namespace bsp1
{
 #include "TNSE_2D/Bsp1.h"
}
namespace lin_space_time
{
#include "TNSE_2D/linear_space_time.h"
}

namespace flow_around_cylinder
{
#include "TNSE_2D/flow_around_cylinder.h"
}

namespace mixing_layer_slip
{
#include "TNSE_2D/MixingLayerSlip.h"
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
      initialCOndtion.push_back(bsp1::InitialU1);
      initialCOndtion.push_back(bsp1::InitialU2);
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
      
      initialCOndtion.push_back(lin_space_time::InitialU1);
      initialCOndtion.push_back(lin_space_time::InitialU2);
      
      lin_space_time::ExampleFile();
      lin_space_time::DIMENSIONLESS_VISCOSITY = this->get_nu();
      break;
    case 2:
      exact_solution.push_back( flow_around_cylinder::ExactU1 );
      exact_solution.push_back( flow_around_cylinder::ExactU2 );
      exact_solution.push_back( flow_around_cylinder::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( flow_around_cylinder::BoundCondition );
      boundary_conditions.push_back( flow_around_cylinder::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( flow_around_cylinder::U1BoundValue );
      boundary_data.push_back( flow_around_cylinder::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = flow_around_cylinder::LinCoeffs;
      
      initialCOndtion.push_back(flow_around_cylinder::InitialU1);
      initialCOndtion.push_back(flow_around_cylinder::InitialU2);
      
      flow_around_cylinder::ExampleFile();
      
      flow_around_cylinder::DIMENSIONLESS_VISCOSITY = this->get_nu();
      
      /**post processing - drag and lift calculation and output */
      post_processing_stat = flow_around_cylinder::compute_drag_and_lift;
      break;
    case 3:
      exact_solution.push_back( mixing_layer_slip::ExactU1 );
      exact_solution.push_back( mixing_layer_slip::ExactU2 );
      exact_solution.push_back( mixing_layer_slip::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( mixing_layer_slip::BoundCondition );
      boundary_conditions.push_back( mixing_layer_slip::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( mixing_layer_slip::U1BoundValue );
      boundary_data.push_back( mixing_layer_slip::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = mixing_layer_slip::LinCoeffs;
      
      initialCOndtion.push_back(mixing_layer_slip::InitialU1);
      initialCOndtion.push_back(mixing_layer_slip::InitialU2);
      
      mixing_layer_slip::ExampleFile();
      
      mixing_layer_slip::DIMENSIONLESS_VISCOSITY = this->get_nu();
      
      /**post processing - drag and lift calculation and output */
      post_processing_stat = mixing_layer_slip::EvaluateSolution;
      break;
    default:
      ErrThrow("Unknown Time dependent Example_TimeNSE2D example!");
  }
}

void Example_TimeNSE2D::do_post_processing(const Time_NSE2D& tnse2d) const
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

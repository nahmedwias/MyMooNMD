#include <Example_TimeNSE2D.h>
#include <Database.h>
#include <MainUtilities.h>
#include <FEDatabase2D.h>
#include <Time_NSE2D.h>
#include "AuxParam2D.h" // used in MixingLayerSlipSmallSquares.h
#include <string>

namespace bsp1              // case 0
{
 #include "TNSE_2D/Bsp1.h"
}

namespace lin_space_time   // case 1
{
#include "TNSE_2D/linear_space_time.h"
}

namespace sincosexp        // case 2
{
#include "TNSE_2D/SinCosExp.h"
}

namespace flow_around_cylinder_steady_inflow     // case 3
{
#include "flow_around_cylinder_steady_inflow.h"
}

namespace backward_facing_step  // case 4
{
#include "TNSE_2D/backward_facing_step.h"
}

namespace driven_cavity         // case 5
{
#include "TNSE_2D/DrivenCavity.h"
}

namespace mixing_layer_us       // case 6
{
#include "TNSE_2D/MixingLayerSlipSmallSquares.h"
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

      // Set dimensionless viscosity
      sincosexp::DIMENSIONLESS_VISCOSITY = get_nu();

      /** coefficients */
      problem_coefficients =sincosexp::LinCoeffs;

      initialCondition.push_back(sincosexp::InitialU1);
      initialCondition.push_back(sincosexp::InitialU2);

      sincosexp::ExampleFile();
      break;

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
    case 4:
      exact_solution.push_back( backward_facing_step::ExactU1 );
      exact_solution.push_back( backward_facing_step::ExactU2 );
      exact_solution.push_back( backward_facing_step::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( backward_facing_step::BoundCondition );
      boundary_conditions.push_back( backward_facing_step::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( backward_facing_step::U1BoundValue );
      boundary_data.push_back( backward_facing_step::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = backward_facing_step::LinCoeffs;
      
      initialCondition.push_back(backward_facing_step::InitialU1);
      initialCondition.push_back(backward_facing_step::InitialU2);
      initialCondition.push_back(backward_facing_step::InitialP);
      
      backward_facing_step::ExampleFile();
      backward_facing_step::DIMENSIONLESS_VISCOSITY = this->get_nu();
      break;
    case 5:
      exact_solution.push_back( driven_cavity::ExactU1 );
      exact_solution.push_back( driven_cavity::ExactU2 );
      exact_solution.push_back( driven_cavity::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( driven_cavity::BoundCondition );
      boundary_conditions.push_back( driven_cavity::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( driven_cavity::U1BoundValue );
      boundary_data.push_back( driven_cavity::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = driven_cavity::LinCoeffs;

      initialCondition.push_back(driven_cavity::InitialU1);
      initialCondition.push_back(driven_cavity::InitialU2);

      driven_cavity::DIMENSIONLESS_VISCOSITY = this->get_nu();
      driven_cavity::ExampleFile();
      break;
    case 6:
      exact_solution.push_back( mixing_layer_us::ExactU1 );
      exact_solution.push_back( mixing_layer_us::ExactU2 );
      exact_solution.push_back( mixing_layer_us::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( mixing_layer_us::BoundCondition );
      boundary_conditions.push_back( mixing_layer_us::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( mixing_layer_us::U1BoundValue );
      boundary_data.push_back( mixing_layer_us::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = mixing_layer_us::LinCoeffs;
      
      initialCondition.push_back(mixing_layer_us::InitialU1);
      initialCondition.push_back(mixing_layer_us::InitialU2);
      
      mixing_layer_us::ExampleFile();
      
      mixing_layer_us::DIMENSIONLESS_VISCOSITY = this->get_nu();
      
      /**post processing - drag and lift calculation and output */
      post_processing_stat = mixing_layer_us::EvaluateSolution;
      break;
    default:
      ErrThrow("Unknown time-dependent Example_TimeNSE2D example!");
  }
}

Example_TimeNSE2D::Example_TimeNSE2D(std::vector <DoubleFunct2D*> exact,
                   std::vector <BoundCondFunct2D*> bc,
                   std::vector <BoundValueFunct2D*> bd,
                   CoeffFct2D coeffs,
                   bool timedependentrhs, bool timedependentcoeffs,
                   std::vector <DoubleFunct2D*> init_cond)
  : Example_NonStationary2D(exact, bc, bd, coeffs, timedependentrhs,
                            timedependentcoeffs, init_cond)
  {

  };

void Example_TimeNSE2D::do_post_processing(Time_NSE2D& tnse2d, double& val) const
{
  if(post_processing_stat)
  {
    post_processing_stat(tnse2d, val);
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

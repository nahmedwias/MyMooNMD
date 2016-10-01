#include <Example_NSE2D.h>
#include <NSE2D.h>

#include <Database.h>
#include <BoundEdge.h>
#include <Example_NSE2D.h>
#include <FEDatabase2D.h>
#include <SquareMatrix2D.h>
#include <string.h>
#include <MainUtilities.h>

/* examples */

namespace poiseuille
{
  #include "NSE_2D/Poiseuille.h"  
}
namespace driven_cavity
{
  #include "NSE_2D/DrivenCavity.h"
}
namespace sine_cosine
{
  #include "NSE_2D/SinCos.h"
}
namespace flow_around_cylinder
{
  #include "NSE_2D/flow_around_cylinder.h"
}
// tests for the pressure robust methods
//=========================================
namespace zerosolution
{
#include "NSE_2D/StokesZeroSol.h"
}

namespace quad_pres
{
  #include "NSE_2D/quadratic_pressure.h"
}

namespace bsp1_pr1
{
#include "NSE_2D/Bsp_PR1.h"
}

//=========================================
namespace td_quad_pres
{
#include "TNSE_2D/stokes_quadratic_pressure.h"
}

namespace time_flow_around_cylinder
{
#include "TNSE_2D/flow_over_cylinder.h"
}

namespace ns_test_code1
{
#include "TNSE_2D/Navier_Stokes_Test_code1.h"
}

namespace potential_flow_ex1
{
#include "TNSE_2D/potential_flow_td.h"
}

namespace potential_flow_ex2
{
#include "TNSE_2D/potential_flow_td_ex2.h"
}
namespace cosine_sin
{
#include "TNSE_2D/cosine_sine.h"
}
//=========================================

Example_NSE2D::Example_NSE2D(const ParameterDatabase& user_input_parameter_db) 
 : Example2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch( example_code )
  {
    case 0:
      /** exact_solution */
      exact_solution.push_back( poiseuille::ExactU1 );
      exact_solution.push_back( poiseuille::ExactU2 );
      exact_solution.push_back( poiseuille::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( poiseuille::BoundCondition );
      boundary_conditions.push_back( poiseuille::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( poiseuille::U1BoundValue );
      boundary_data.push_back( poiseuille::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = poiseuille::LinCoeffs;
      
      poiseuille::ExampleFile();
      break;
    case 1:
      /** exact_solution */
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
      
      driven_cavity::ExampleFile();
      break;
    case 2:
      /** exact_solution */
      exact_solution.push_back( sine_cosine::ExactU1 );
      exact_solution.push_back( sine_cosine::ExactU2 );
      exact_solution.push_back( sine_cosine::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( sine_cosine::BoundCondition );
      boundary_conditions.push_back( sine_cosine::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sine_cosine::U1BoundValue );
      boundary_data.push_back( sine_cosine::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = sine_cosine::LinCoeffs;
      
      sine_cosine::ExampleFile();
      break;
    case 3:
      /** exact_solution */
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
      
      // Set dimensionless viscosity
      flow_around_cylinder::DIMENSIONLESS_VISCOSITY = get_nu();
      // post_processing = flow_around_cylinder::compute_drag_and_lift;
      
    case 40:
      /** exact_solution */
      exact_solution.push_back( zerosolution::ExactU1 );
      exact_solution.push_back( zerosolution::ExactU2 );
      exact_solution.push_back( zerosolution::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( zerosolution::BoundCondition );
      boundary_conditions.push_back( zerosolution::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( zerosolution::U1BoundValue );
      boundary_data.push_back( zerosolution::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = zerosolution::LinCoeffs;
      
      zerosolution::ExampleFile();
      break;
    case 41:
      /** exact_solution */
      exact_solution.push_back( quad_pres::ExactU1 );
      exact_solution.push_back( quad_pres::ExactU2 );
      exact_solution.push_back( quad_pres::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( quad_pres::BoundCondition );
      boundary_conditions.push_back( quad_pres::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( quad_pres::U1BoundValue );
      boundary_data.push_back( quad_pres::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = quad_pres::LinCoeffs;
      
      quad_pres::ExampleFile();
      break;
      
       case 42:
      /** exact_solution */
      exact_solution.push_back( bsp1_pr1::ExactU1 );
      exact_solution.push_back( bsp1_pr1::ExactU2 );
      exact_solution.push_back( bsp1_pr1::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( bsp1_pr1::BoundCondition );
      boundary_conditions.push_back( bsp1_pr1::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( bsp1_pr1::U1BoundValue );
      boundary_data.push_back( bsp1_pr1::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = bsp1_pr1::LinCoeffs;
      
      bsp1_pr1::ExampleFile();
      break;
      
//    case 101:
//      /** exact_solution */
//      exact_solution.push_back( bsp1::ExactU1 );
//      exact_solution.push_back( bsp1::ExactU2 );
//      exact_solution.push_back( bsp1::ExactP );
      
//      /** boundary condition */
//      boundary_conditions.push_back( bsp1::BoundCondition );
//      boundary_conditions.push_back( bsp1::BoundCondition );
//      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
//      /** boundary values */
//      boundary_data.push_back( bsp1::U1BoundValue );
//      boundary_data.push_back( bsp1::U2BoundValue );
//      boundary_data.push_back( BoundaryValueHomogenous );
      
//      /** coefficients */
//      problem_coefficients = bsp1::LinCoeffs;
      
//      /** initial condition */
//      initial_conditions.push_back(bsp1::InitialU1);
//      initial_conditions.push_back(bsp1::InitialU2);
//      bsp1::ExampleFile();
//      break;
//    case 102:
//      /** exact_solution */
//      exact_solution.push_back( lin_space_time::ExactU1 );
//      exact_solution.push_back( lin_space_time::ExactU2 );
//      exact_solution.push_back( lin_space_time::ExactP );
      
//      /** boundary condition */
//      boundary_conditions.push_back( lin_space_time::BoundCondition );
//      boundary_conditions.push_back( lin_space_time::BoundCondition );
//      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
//      /** boundary values */
//      boundary_data.push_back( lin_space_time::U1BoundValue );
//      boundary_data.push_back( lin_space_time::U2BoundValue );
//      boundary_data.push_back( BoundaryValueHomogenous );
      
//      /** coefficients */
//      problem_coefficients = lin_space_time::LinCoeffs;
      
//      initial_conditions.push_back(lin_space_time::InitialU1);
//      initial_conditions.push_back(lin_space_time::InitialU2);
      
//      lin_space_time::ExampleFile();
//      break;
//    case 103:
//      /** exact_solution */
//      exact_solution.push_back( td_quad_pres::ExactU1 );
//      exact_solution.push_back( td_quad_pres::ExactU2 );
//      exact_solution.push_back( td_quad_pres::ExactP );
      
//      /** boundary condition */
//      boundary_conditions.push_back( td_quad_pres::BoundCondition );
//      boundary_conditions.push_back( td_quad_pres::BoundCondition );
//      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
//      /** boundary values */
//      boundary_data.push_back( td_quad_pres::U1BoundValue );
//      boundary_data.push_back( td_quad_pres::U2BoundValue );
//      boundary_data.push_back( BoundaryValueHomogenous );
      
//      /** coefficients */
//      problem_coefficients = td_quad_pres::LinCoeffs;
      
//      initial_conditions.push_back(td_quad_pres::InitialU1);
//      initial_conditions.push_back(td_quad_pres::InitialU2);
      
//      td_quad_pres::ExampleFile();
//      break;
//    case 104:
//      /** exact_solution */
//      exact_solution.push_back( time_flow_around_cylinder::ExactU1 );
//      exact_solution.push_back( time_flow_around_cylinder::ExactU2 );
//      exact_solution.push_back( time_flow_around_cylinder::ExactP );
      
//      /** boundary condition */
//      boundary_conditions.push_back( time_flow_around_cylinder::BoundCondition );
//      boundary_conditions.push_back( time_flow_around_cylinder::BoundCondition );
//      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
//      /** boundary values */
//      boundary_data.push_back( time_flow_around_cylinder::U1BoundValue );
//      boundary_data.push_back( time_flow_around_cylinder::U2BoundValue );
//      boundary_data.push_back( BoundaryValueHomogenous );
      
//      /** initial conditions, in case of a non-stationary problem */
//      initial_conditions.push_back(time_flow_around_cylinder::InitialU1);
//      initial_conditions.push_back(time_flow_around_cylinder::InitialU2);
//      initial_conditions.push_back(time_flow_around_cylinder::InitialP);
      
//      /** post processing function to  compute drag and lift */
//      post_processing_time = time_flow_around_cylinder::compute_drag_and_lift;

//      /** coefficients */
//      problem_coefficients = time_flow_around_cylinder::LinCoeffs;
      
//      time_flow_around_cylinder::ExampleFile();
//      break;
//    case 105:
//      /** exact_solution */
//      exact_solution.push_back( ns_test_code1::ExactU1 );
//      exact_solution.push_back( ns_test_code1::ExactU2 );
//      exact_solution.push_back( ns_test_code1::ExactP );
      
//      /** boundary condition */
//      boundary_conditions.push_back( ns_test_code1::BoundCondition );
//      boundary_conditions.push_back( ns_test_code1::BoundCondition );
//      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
//      /** boundary values */
//      boundary_data.push_back( ns_test_code1::U1BoundValue );
//      boundary_data.push_back( ns_test_code1::U2BoundValue );
//      boundary_data.push_back( BoundaryValueHomogenous );
      
//      /** initial conditions, in case of a non-stationary problem */
//      initial_conditions.push_back(ns_test_code1::InitialU1);
//      initial_conditions.push_back(ns_test_code1::InitialU2);
//      initial_conditions.push_back(ns_test_code1::InitialP);
//      /** coefficients */
//      problem_coefficients = ns_test_code1::LinCoeffs;
      
//      ns_test_code1::ExampleFile();
//      break;
      
//    case 106:
//      exact_solution.push_back( potential_flow_ex1::ExactU1 );
//      exact_solution.push_back( potential_flow_ex1::ExactU2 );
//      exact_solution.push_back( potential_flow_ex1::ExactP );
      
//      /** boundary condition */
//      boundary_conditions.push_back( potential_flow_ex1::BoundCondition );
//      boundary_conditions.push_back( potential_flow_ex1::BoundCondition );
//      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
//      /** boundary values */
//      boundary_data.push_back( potential_flow_ex1::U1BoundValue );
//      boundary_data.push_back( potential_flow_ex1::U2BoundValue );
//      boundary_data.push_back( BoundaryValueHomogenous );
      
//      /** initial conditions, in case of a non-stationary problem */
//      initial_conditions.push_back(potential_flow_ex1::InitialU1);
//      initial_conditions.push_back(potential_flow_ex1::InitialU2);
//      initial_conditions.push_back(potential_flow_ex1::InitialP);
//      /** coefficients */
//      problem_coefficients = potential_flow_ex1::LinCoeffs;
      
//      potential_flow_ex1::ExampleFile();
//      break;
//      case 107:
//      /** exact_solution */
//      exact_solution.push_back( potential_flow_ex2::ExactU1 );
//      exact_solution.push_back( potential_flow_ex2::ExactU2 );
//      exact_solution.push_back( potential_flow_ex2::ExactP );
      
//      /** boundary condition */
//      boundary_conditions.push_back( potential_flow_ex2::BoundCondition );
//      boundary_conditions.push_back( potential_flow_ex2::BoundCondition );
//      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
//      /** boundary values */
//      boundary_data.push_back( potential_flow_ex2::U1BoundValue );
//      boundary_data.push_back( potential_flow_ex2::U2BoundValue );
//      boundary_data.push_back( BoundaryValueHomogenous );
      
//      /** initial conditions, in case of a non-stationary problem */
//      initial_conditions.push_back(potential_flow_ex2::InitialU1);
//      initial_conditions.push_back(potential_flow_ex2::InitialU2);
//      initial_conditions.push_back(potential_flow_ex2::InitialP);
//      /** coefficients */
//      problem_coefficients = potential_flow_ex2::LinCoeffs;
      
//      potential_flow_ex2::ExampleFile();
//      break;
//      case 108:
//      /** exact_solution */
//      exact_solution.push_back( cosine_sin::ExactU1 );
//      exact_solution.push_back( cosine_sin::ExactU2 );
//      exact_solution.push_back( cosine_sin::ExactP );
      
//      /** boundary condition */
//      boundary_conditions.push_back( cosine_sin::BoundCondition );
//      boundary_conditions.push_back( cosine_sin::BoundCondition );
//      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
//      /** boundary values */
//      boundary_data.push_back( cosine_sin::U1BoundValue );
//      boundary_data.push_back( cosine_sin::U2BoundValue );
//      boundary_data.push_back( BoundaryValueHomogenous );
      
//      /** initial conditions, in case of a non-stationary problem */
//      initial_conditions.push_back(cosine_sin::InitialU1);
//      initial_conditions.push_back(cosine_sin::InitialU2);
//      initial_conditions.push_back(cosine_sin::InitialP);
//      /** coefficients */
//      problem_coefficients = cosine_sin::LinCoeffs;
      
//      cosine_sin::ExampleFile();
//      break;
    default:
      ErrThrow("Unknown Navier-Stokes example!");
  }
}

void Example_NSE2D::do_post_processing(NSE2D& nse2d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(nse2d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_NSE2D","No post processing done for the current example.");
  }
}

double Example_NSE2D::get_nu() const
{
  double inverse_reynolds = this->example_database["reynolds_number"];
  inverse_reynolds = 1/inverse_reynolds;
  return inverse_reynolds;
}


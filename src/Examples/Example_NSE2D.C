#include <Example_NSE2D.h>

#include <Database.h>
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
namespace quad_pres
{
  #include "NSE_2D/quadratic_pressure.h"
}

//=========================================
// tests for the pressure robust methods
//=========================================
namespace zerosolution
{
#include "NSE_2D/StokesZeroSol.h"
}
//=========================================
// time dependent case 
namespace bsp1
{
 #include "TNSE_2D/Bsp1.h"
}
namespace lin_space_time
{
#include "TNSE_2D/linear_space_time.h"
}
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
//=========================================

Example_NSE2D::Example_NSE2D() : Example2D()
{
  switch( TDatabase::ParamDB->EXAMPLE ) 
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
      
      flow_around_cylinder::ExampleFile();
      break;
      
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
      
    case 101:
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
      initial_conditions.push_back(bsp1::InitialU1);
      initial_conditions.push_back(bsp1::InitialU2);
      bsp1::ExampleFile();
      break;
    case 102:
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
      
      initial_conditions.push_back(lin_space_time::InitialU1);
      initial_conditions.push_back(lin_space_time::InitialU2);
      
      lin_space_time::ExampleFile();
      break;
    case 103:
      /** exact_solution */
      exact_solution.push_back( td_quad_pres::ExactU1 );
      exact_solution.push_back( td_quad_pres::ExactU2 );
      exact_solution.push_back( td_quad_pres::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( td_quad_pres::BoundCondition );
      boundary_conditions.push_back( td_quad_pres::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( td_quad_pres::U1BoundValue );
      boundary_data.push_back( td_quad_pres::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = td_quad_pres::LinCoeffs;
      
      initial_conditions.push_back(td_quad_pres::InitialU1);
      initial_conditions.push_back(td_quad_pres::InitialU2);
      
      td_quad_pres::ExampleFile();
      break;    
    case 104:
      /** exact_solution */
      exact_solution.push_back( time_flow_around_cylinder::ExactU1 );
      exact_solution.push_back( time_flow_around_cylinder::ExactU2 );
      exact_solution.push_back( time_flow_around_cylinder::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( time_flow_around_cylinder::BoundCondition );
      boundary_conditions.push_back( time_flow_around_cylinder::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( time_flow_around_cylinder::U1BoundValue );
      boundary_data.push_back( time_flow_around_cylinder::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** initial conditions, in case of a non-stationary problem */
      initial_conditions.push_back(time_flow_around_cylinder::InitialU1);
      initial_conditions.push_back(time_flow_around_cylinder::InitialU2);
      initial_conditions.push_back(time_flow_around_cylinder::InitialP);
      
      /** post processing function to  compute drag and lift */
      post_processing_time = time_flow_around_cylinder::compute_drag_and_lift;

      /** coefficients */
      problem_coefficients = time_flow_around_cylinder::LinCoeffs;
      
      time_flow_around_cylinder::ExampleFile();
      break;
    case 105:
      /** exact_solution */
      exact_solution.push_back( ns_test_code1::ExactU1 );
      exact_solution.push_back( ns_test_code1::ExactU2 );
      exact_solution.push_back( ns_test_code1::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( ns_test_code1::BoundCondition );
      boundary_conditions.push_back( ns_test_code1::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( ns_test_code1::U1BoundValue );
      boundary_data.push_back( ns_test_code1::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** initial conditions, in case of a non-stationary problem */
      initial_conditions.push_back(ns_test_code1::InitialU1);
      initial_conditions.push_back(ns_test_code1::InitialU2);
      initial_conditions.push_back(ns_test_code1::InitialP);
      /** coefficients */
      problem_coefficients = ns_test_code1::LinCoeffs;
      
      ns_test_code1::ExampleFile();
      break;
    default:
      ErrThrow("Unknown Navier-Stokes example!");
  }
}

      

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

namespace potential_flow_ex3 {
#include "NSE_2D/potential_flow_ex3.h"
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
      post_processing_stat = flow_around_cylinder::compute_drag_lift_pdiff;

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
      
      zerosolution::DIMENSIONLESS_VISCOSITY = this->get_nu();
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
      
      quad_pres::DIMENSIONLESS_VISCOSITY = this->get_nu();
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
      bsp1_pr1::DIMENSIONLESS_VISCOSITY = this->get_nu();
      bsp1_pr1::ExampleFile();
      break;
  case 43:
      /** exact_solution */
      exact_solution.push_back( potential_flow_ex3::ExactU1 );
      exact_solution.push_back( potential_flow_ex3::ExactU2 );
      exact_solution.push_back( potential_flow_ex3::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( potential_flow_ex3::BoundCondition );
      boundary_conditions.push_back( potential_flow_ex3::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( potential_flow_ex3::U1BoundValue );
      boundary_data.push_back( potential_flow_ex3::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      potential_flow_ex3::DIMENSIONLESS_VISCOSITY = this->get_nu();

      /** coefficients */
      problem_coefficients = potential_flow_ex3::LinCoeffs;

      potential_flow_ex3::ExampleFile();
      break;
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


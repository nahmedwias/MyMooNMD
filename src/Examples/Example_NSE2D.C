#include <Example_NSE2D.h>
#include "NavierStokes.h"
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
namespace backward_facing_step
{
  #include "NSE_2D/backward_facing_step.h"
}
namespace exampleD3
{
  #include "NSE_2D/polynomial_solution.h"
}
namespace brinkman_poiseuille
{
  #include "NSE_2D/Brinkman_Poiseuille.h"
}
namespace brinkman_sincos_darcyflow
{
  #include "NSE_2D/Brinkman_SinCos_DarcyFlow.h"
}
namespace brinkman_discacciatiflow
{
#include "NSE_2D/Brinkman_DiscacciatiFlow.h"
}

//============================================================================

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
      
      // Set dimensionless viscosity
      driven_cavity::DIMENSIONLESS_VISCOSITY = get_nu();

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

      /**post processing - drag and lift calculation and output */
      post_processing_stat = flow_around_cylinder::compute_drag_lift_pdiff;

      flow_around_cylinder::ExampleFile();
      break;
    case 4:
      /** exact_solution */
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
      
      // Set dimensionless viscosity
      backward_facing_step::DIMENSIONLESS_VISCOSITY = get_nu();
      backward_facing_step::ExampleFile();
      break;
    case 5:
      /** exact_solution */
      exact_solution.push_back( exampleD3::ExactU1 );
      exact_solution.push_back( exampleD3::ExactU2 );
      exact_solution.push_back( exampleD3::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( exampleD3::BoundCondition );
      boundary_conditions.push_back( exampleD3::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( exampleD3::U1BoundValue );
      boundary_data.push_back( exampleD3::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = exampleD3::LinCoeffs;
      
      // Set dimensionless viscosity
      exampleD3::DIMENSIONLESS_VISCOSITY = get_nu();
      exampleD3::ExampleFile();
      break;
    case 6:
      /** exact_solution */
      exact_solution.push_back( brinkman_poiseuille::ExactU1 );
      exact_solution.push_back( brinkman_poiseuille::ExactU2 );
      exact_solution.push_back( brinkman_poiseuille::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( brinkman_poiseuille::BoundCondition );
      boundary_conditions.push_back( brinkman_poiseuille::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( brinkman_poiseuille::U1BoundValue );
      boundary_data.push_back( brinkman_poiseuille::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = brinkman_poiseuille::LinCoeffs;
      
      // Set dimensionless viscosity
      brinkman_poiseuille::effective_viscosity = get_nu();
      brinkman_poiseuille::sigma = get_inverse_permeability();
      brinkman_poiseuille::neumann_id = get_neumann_id();
      brinkman_poiseuille::nitsche_id = get_nitsche_id();
      
      brinkman_poiseuille::ExampleFile();
      break;
case 7:
      /** exact_solution */
      exact_solution.push_back( brinkman_sincos_darcyflow::ExactU1 );
      exact_solution.push_back( brinkman_sincos_darcyflow::ExactU2 );
      exact_solution.push_back( brinkman_sincos_darcyflow::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( brinkman_sincos_darcyflow::BoundCondition );
      boundary_conditions.push_back( brinkman_sincos_darcyflow::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( brinkman_sincos_darcyflow::U1BoundValue );
      boundary_data.push_back( brinkman_sincos_darcyflow::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = brinkman_sincos_darcyflow::LinCoeffs;
      
      // Set dimensionless viscosity
      brinkman_sincos_darcyflow::effective_viscosity = get_nu();
      brinkman_sincos_darcyflow::sigma = get_inverse_permeability();
      brinkman_sincos_darcyflow::ExampleFile();
      break;
case 8:
      /** exact_solution */
      exact_solution.push_back( brinkman_discacciatiflow::ExactU1 );
      exact_solution.push_back( brinkman_discacciatiflow::ExactU2 );
      exact_solution.push_back( brinkman_discacciatiflow::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( brinkman_discacciatiflow::BoundCondition );
      boundary_conditions.push_back( brinkman_discacciatiflow::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( brinkman_discacciatiflow::U1BoundValue );
      boundary_data.push_back( brinkman_discacciatiflow::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = brinkman_discacciatiflow::LinCoeffs;
      
      // Set dimensionless viscosity
      brinkman_discacciatiflow::effective_viscosity = get_nu();
      brinkman_discacciatiflow::sigma = get_inverse_permeability();
      brinkman_discacciatiflow::ExampleFile();
      break;

    default:
      ErrThrow("Unknown Navier-Stokes example!");
  }
}

Example_NSE2D::Example_NSE2D(std::vector<DoubleFunct2D *> exact, 
                             std::vector<BoundCondFunct2D *> bc,
                             std::vector<BoundValueFunct2D *> bd, 
                             CoeffFct2D coeffs, double nu) 
 : Example2D(exact, bc, bd, coeffs)
{
  this->example_database["reynolds_number"] = 1./nu;
}


void Example_NSE2D::do_post_processing(NavierStokes<2 >& nse2d) const
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

double Example_NSE2D::get_inverse_permeability() const
{
  return this->example_database["inverse_permeability"];
}

std::vector<size_t> Example_NSE2D::get_neumann_id() const
{
  return this->example_database["neumann_id"];
}


std::vector<size_t> Example_NSE2D::get_nitsche_id() const
{
  return this->example_database["nitsche_id"];
}

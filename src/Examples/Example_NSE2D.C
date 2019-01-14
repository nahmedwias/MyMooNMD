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
namespace brinkman_radial_flow_with_hole
{
#include "NSE_2D/Brinkman_Radial_Flow_with_hole.h"
}
namespace brinkman_circle_with_immersed_hole
{
#include "NSE_2D/Brinkman_circle_with_immersed_hole.h"
}
namespace brinkman_two_wells
{
#include "NSE_2D/Brinkman_two_wells.h"
}
namespace brinkman_expo_poiseuille
{
#include "NSE_2D/Brinkman_Exponential_Poiseuille.h"
}
namespace brinkman_darcy_varying_permeability
{
#include "NSE_2D/Brinkman_Darcy_varying_permeability.h"
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

      brinkman_sincos_darcyflow::neumann_id = get_neumann_id();
      brinkman_sincos_darcyflow::nitsche_id = get_nitsche_id();
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

      // boundary conditions
      brinkman_discacciatiflow::neumann_id = get_neumann_id();
      brinkman_discacciatiflow::nitsche_id = get_nitsche_id();
      brinkman_discacciatiflow::ExampleFile();
      break;

  case 9:
    /** exact_solution */
    exact_solution.push_back( brinkman_radial_flow_with_hole::ExactU1 );
    exact_solution.push_back( brinkman_radial_flow_with_hole::ExactU2 );
    exact_solution.push_back( brinkman_radial_flow_with_hole::ExactP );
    
    /** boundary condition */
    boundary_conditions.push_back( brinkman_radial_flow_with_hole::BoundCondition );
    boundary_conditions.push_back( brinkman_radial_flow_with_hole::BoundCondition );
    boundary_conditions.push_back( BoundConditionNoBoundCondition );
    
    /** boundary values */
    boundary_data.push_back( brinkman_radial_flow_with_hole::U1BoundValue );
    boundary_data.push_back( brinkman_radial_flow_with_hole::U2BoundValue );
    boundary_data.push_back( BoundaryValueHomogenous );
    
    /** coefficients */
    problem_coefficients = brinkman_radial_flow_with_hole::LinCoeffs;
    
    // Set dimensionless viscosity
    brinkman_radial_flow_with_hole::effective_viscosity = get_nu();
    brinkman_radial_flow_with_hole::sigma = get_inverse_permeability();
    
    // boundary conditions
    brinkman_radial_flow_with_hole::neumann_id = get_neumann_id();
    brinkman_radial_flow_with_hole::nitsche_id = get_nitsche_id();
    brinkman_radial_flow_with_hole::ExampleFile();
    break;
    
  case 10:
    /** exact_solution */
    exact_solution.push_back( brinkman_circle_with_immersed_hole::ExactU1 );
    exact_solution.push_back( brinkman_circle_with_immersed_hole::ExactU2 );
    exact_solution.push_back( brinkman_circle_with_immersed_hole::ExactP );
    
    /** boundary condition */
    boundary_conditions.push_back( brinkman_circle_with_immersed_hole::BoundCondition );
    boundary_conditions.push_back( brinkman_circle_with_immersed_hole::BoundCondition );
    boundary_conditions.push_back( BoundConditionNoBoundCondition );
    
    /** boundary values */
    boundary_data.push_back( brinkman_circle_with_immersed_hole::U1BoundValue );
    boundary_data.push_back( brinkman_circle_with_immersed_hole::U2BoundValue );
    boundary_data.push_back( BoundaryValueHomogenous );
    
    /** coefficients */
    problem_coefficients = brinkman_circle_with_immersed_hole::LinCoeffs;
    
    // Set dimensionless viscosity
    brinkman_circle_with_immersed_hole::effective_viscosity = get_nu();
    brinkman_circle_with_immersed_hole::sigma = get_inverse_permeability();
    
    // boundary conditions
    brinkman_circle_with_immersed_hole::neumann_id = get_neumann_id();
    brinkman_circle_with_immersed_hole::nitsche_id = get_nitsche_id();
    brinkman_circle_with_immersed_hole::ExampleFile();
    break;

  case 11:
    /** exact_solution */
    exact_solution.push_back( brinkman_two_wells::ExactU1 );
    exact_solution.push_back( brinkman_two_wells::ExactU2 );
    exact_solution.push_back( brinkman_two_wells::ExactP );
    
    /** boundary condition */
    boundary_conditions.push_back( brinkman_two_wells::BoundCondition );
    boundary_conditions.push_back( brinkman_two_wells::BoundCondition );
    boundary_conditions.push_back( BoundConditionNoBoundCondition );
    
    /** boundary values */
    boundary_data.push_back( brinkman_two_wells::U1BoundValue );
    boundary_data.push_back( brinkman_two_wells::U2BoundValue );
    boundary_data.push_back( BoundaryValueHomogenous );
    
    /** coefficients */
    problem_coefficients = brinkman_two_wells::LinCoeffs;
    
    // Set dimensionless viscosity
    brinkman_two_wells::effective_viscosity = get_nu();
    brinkman_two_wells::sigma = get_inverse_permeability();
    
    // boundary conditions
    brinkman_two_wells::neumann_id = get_neumann_id();
    brinkman_two_wells::nitsche_id = get_nitsche_id();
    brinkman_two_wells::ExampleFile();
    break;
    
  case 12:
      /** exact_solution */
      exact_solution.push_back( brinkman_expo_poiseuille::ExactU1 );
      exact_solution.push_back( brinkman_expo_poiseuille::ExactU2 );
      exact_solution.push_back( brinkman_expo_poiseuille::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( brinkman_expo_poiseuille::BoundCondition );
      boundary_conditions.push_back( brinkman_expo_poiseuille::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( brinkman_expo_poiseuille::U1BoundValue );
      boundary_data.push_back( brinkman_expo_poiseuille::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = brinkman_expo_poiseuille::LinCoeffs;

      // Set dimensionless viscosity
      brinkman_expo_poiseuille::effective_viscosity = get_nu();
      brinkman_expo_poiseuille::sigma = get_inverse_permeability();

      // boundary conditions
      brinkman_expo_poiseuille::neumann_id = get_neumann_id();
      brinkman_expo_poiseuille::nitsche_id = get_nitsche_id();
      brinkman_expo_poiseuille::ExampleFile();
      break;

  case 13:
      /** exact_solution */
      exact_solution.push_back( brinkman_darcy_varying_permeability::ExactU1 );
      exact_solution.push_back( brinkman_darcy_varying_permeability::ExactU2 );
      exact_solution.push_back( brinkman_darcy_varying_permeability::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( brinkman_darcy_varying_permeability::BoundCondition );
      boundary_conditions.push_back( brinkman_darcy_varying_permeability::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( brinkman_darcy_varying_permeability::U1BoundValue );
      boundary_data.push_back( brinkman_darcy_varying_permeability::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = brinkman_darcy_varying_permeability::LinCoeffs;

      // Set dimensionless viscosity
      brinkman_darcy_varying_permeability::effective_viscosity = get_nu();
      brinkman_darcy_varying_permeability::sigma = get_inverse_permeability();

      // boundary conditions
      brinkman_darcy_varying_permeability::neumann_id = get_neumann_id();
      brinkman_darcy_varying_permeability::nitsche_id = get_nitsche_id();
      brinkman_darcy_varying_permeability::ExampleFile();
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
  int example_code = this->example_database["example"];
  double nu=-1.;
  switch( example_code )
  {
  case 0:
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
    {
      nu = this->example_database["reynolds_number"];
      nu = 1./nu;
      break;
    }
  case 6:
  case 7:
  case 8:
  case 9:
  case 10:
  case 11:
  case 12:
  case 13:
    {
      nu = this->example_database["effective_viscosity"];
      break;
    }
  default:
    ErrThrow("Unknown Navier-Stokes example!");
  }
  return nu;
}


double Example_NSE2D::get_inverse_permeability() const
{
  return this->example_database["inverse_permeability"];
}

std::vector<size_t> Example_NSE2D::get_neumann_id() const
{
  std::vector<size_t> empty_vector;
  int n_neumann_bd = this->example_database["n_neumann_bd"];
  if (n_neumann_bd){
    return this->example_database["neumann_id"];
  } else {
    return empty_vector;
  }
}


std::vector<size_t> Example_NSE2D::get_nitsche_id() const
{
  std::vector<size_t> empty_vector;
  int n_nitsche_bd = this->example_database["n_nitsche_bd"];
  if (n_nitsche_bd){
    return this->example_database["nitsche_id"];
  } else {
    return empty_vector;
  }
}

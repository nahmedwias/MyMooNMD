#include <Example_NSE2D.h>
#include <NSE2D.h>

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




/***** BELOW THIS LINE, EXAMPLES ARE SPECIFIC TO USER PROJCET *****/
namespace example_layout
{
#include "0_ExampleLayout.h"
}

namespace test_example1_poiseuille
{
#include "1_TestExample1.h"
}

namespace test_example2_sincos
{
#include "2_TestExample2.h"
}

namespace test_example3_drivencavity
{
#include "3_TestExample3.h"
}

namespace test_example4_velocityforconvection
{
#include "4_TestExample4.h"
}

namespace variable_viscosity
{
#include "../../user_projects/include/Examples/NSE2D/variable_viscosity.h"
}
// ********* END OF USER PROJECT CODE





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





      /**************************************************************/
      /* BELOW THIS LINE, THE EXAMPLES ARE SPECIFIC TO USER PROJECT */
    case 10:
      /** exact_solution */
      exact_solution.push_back( example_layout::ExactU1 );
      exact_solution.push_back( example_layout::ExactU2 );
      exact_solution.push_back( example_layout::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( example_layout::BoundCondition );
      boundary_conditions.push_back( example_layout::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( example_layout::U1BoundValue );
      boundary_data.push_back( example_layout::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = example_layout::LinCoeffs;

      // Set dimensionless viscosity
      example_layout::REYNOLDS_number = get_nu();
      example_layout::USER_parameter1 = this->example_database["user_parameter1"];
      example_layout::USER_parameter2 = this->example_database["user_parameter2"];

      example_layout::ExampleFile();
      break;
    case 11:  // TestExample 1 = Poiseuille
      /** exact_solution */
      exact_solution.push_back( test_example1_poiseuille::ExactU1 );
      exact_solution.push_back( test_example1_poiseuille::ExactU2 );
      exact_solution.push_back( test_example1_poiseuille::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( test_example1_poiseuille::BoundCondition );
      boundary_conditions.push_back( test_example1_poiseuille::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example1_poiseuille::U1BoundValue );
      boundary_data.push_back( test_example1_poiseuille::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = test_example1_poiseuille::LinCoeffs;

      // Set dimensionless viscosity
      test_example1_poiseuille::REYNOLDS_number = get_nu();
      test_example1_poiseuille::USER_parameter1 = this->example_database["user_parameter1"];
      test_example1_poiseuille::USER_parameter2 = this->example_database["user_parameter2"];

      test_example1_poiseuille::ExampleFile();
      break;
    case 12:  // TestExample 2 = SinCos
      /** exact_solution */
      exact_solution.push_back( test_example2_sincos::ExactU1 );
      exact_solution.push_back( test_example2_sincos::ExactU2 );
      exact_solution.push_back( test_example2_sincos::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( test_example2_sincos::BoundCondition );
      boundary_conditions.push_back( test_example2_sincos::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example2_sincos::U1BoundValue );
      boundary_data.push_back( test_example2_sincos::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = test_example2_sincos::LinCoeffs;

      // Set dimensionless viscosity
      test_example2_sincos::REYNOLDS_number = get_nu();
      test_example2_sincos::USER_parameter1 = this->example_database["user_parameter1"];
      test_example2_sincos::USER_parameter2 = this->example_database["user_parameter2"];

      test_example2_sincos::ExampleFile();
      break;
    case 13:  // TestExample 3 = DrivenCavity
      /** exact_solution */
      exact_solution.push_back( test_example3_drivencavity::ExactU1 );
      exact_solution.push_back( test_example3_drivencavity::ExactU2 );
      exact_solution.push_back( test_example3_drivencavity::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( test_example3_drivencavity::BoundCondition );
      boundary_conditions.push_back( test_example3_drivencavity::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example3_drivencavity::U1BoundValue );
      boundary_data.push_back( test_example3_drivencavity::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = test_example3_drivencavity::LinCoeffs;

      // Set dimensionless viscosity
      test_example3_drivencavity::REYNOLDS_number = get_nu();
      test_example3_drivencavity::USER_parameter1 = this->example_database["user_parameter1"];
      test_example3_drivencavity::USER_parameter2 = this->example_database["user_parameter2"];

      test_example3_drivencavity::ExampleFile();
      break;
    case 14:  // TestExample 4 = Generates the velocity field for the convection example "Two interior Layers"
      /** exact_solution */
      exact_solution.push_back( test_example4_velocityforconvection::ExactU1 );
      exact_solution.push_back( test_example4_velocityforconvection::ExactU2 );
      exact_solution.push_back( test_example4_velocityforconvection::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( test_example4_velocityforconvection::BoundCondition );
      boundary_conditions.push_back( test_example4_velocityforconvection::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example4_velocityforconvection::U1BoundValue );
      boundary_data.push_back( test_example4_velocityforconvection::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = test_example4_velocityforconvection::LinCoeffs;

      // Set dimensionless viscosity
      test_example4_velocityforconvection::REYNOLDS_number = get_nu();
      test_example4_velocityforconvection::USER_parameter1 = this->example_database["user_parameter1"];
      test_example4_velocityforconvection::USER_parameter2 = this->example_database["user_parameter2"];

      test_example4_velocityforconvection::ExampleFile();
      break;
    case 15:  // variable viscosity = Generates variable viscosity problem
      /** exact_solution */
      exact_solution.push_back( variable_viscosity::ExactU1 );
      exact_solution.push_back( variable_viscosity::ExactU2 );
      exact_solution.push_back( variable_viscosity::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( variable_viscosity::BoundCondition );
      boundary_conditions.push_back( variable_viscosity::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( variable_viscosity::U1BoundValue );
      boundary_data.push_back( variable_viscosity::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = variable_viscosity::LinCoeffs;

      // Set dimensionless viscosity
//      variable_viscosity::REYNOLDS_number = get_nu();
//      variable_viscosity::USER_parameter1 = this->example_database["user_parameter1"];
//      variable_viscosity::USER_parameter2 = this->example_database["user_parameter2"];

      variable_viscosity::ExampleFile();
      break;

      /**************************************************************/
      /* ABOVE THIS LINE, THE EXAMPLES ARE SPECIFIC TO USER PROJECT */



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


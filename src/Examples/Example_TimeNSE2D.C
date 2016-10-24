#include <Example_TimeNSE2D.h>
#include <Database.h>
#include <MainUtilities.h>

#include <string>

namespace bsp1
{
 #include "TNSE_2D/Bsp1.h"
}
namespace lin_space_time
{
#include "TNSE_2D/linear_space_time.h"
}


/***** BELOW THIS LINE, EXAMPLES ARE SPECIFIC TO USER PROJECT *****/
namespace test_example1_sincos
{
#include "1_TestExampleTNSE.h"
}
namespace test_example2_forconvectiontest1
{
#include "2_TestExampleForTCDTestExample1.h"
}
namespace test_example3_forconvectiontest2
{
#include "3_TestExampleForTCDTestExample2.h"
}
namespace test_4 // for dam break, tcd example 15
{
#include "4_TestTNSE2D.h"
}
namespace test_5 // for semicircle, tcd example 16
{
#include "5_TestTNSE2D.h"
}
// ********* END OF USER PROJECT CODE



//=========================================
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
      break;







      /***** BELOW THIS LINE, EXAMPLES ARE SPECIFIC TO USER PROJECT *****/
    case 11:                // TestExample1 = Bsp1 = sincos
      /** exact_solution */
      exact_solution.push_back( test_example1_sincos::ExactU1 );
      exact_solution.push_back( test_example1_sincos::ExactU2 );
      exact_solution.push_back( test_example1_sincos::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( test_example1_sincos::BoundCondition );
      boundary_conditions.push_back( test_example1_sincos::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example1_sincos::U1BoundValue );
      boundary_data.push_back( test_example1_sincos::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = test_example1_sincos::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(test_example1_sincos::InitialU1);
      initialCOndtion.push_back(test_example1_sincos::InitialU2);
      test_example1_sincos::REYNOLDS_number = get_nu();
      test_example1_sincos::USER_parameter1 = this->example_database["user_parameter1"];
      test_example1_sincos::USER_parameter2 = this->example_database["user_parameter2"];



      test_example1_sincos::ExampleFile();
      break;

    case 12:                // TestExample2 = Generate the velocity field for Convection Test Example 1 (1,-1)
      /** exact_solution */
      exact_solution.push_back( test_example2_forconvectiontest1::ExactU1 );
      exact_solution.push_back( test_example2_forconvectiontest1::ExactU2 );
      exact_solution.push_back( test_example2_forconvectiontest1::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( test_example2_forconvectiontest1::BoundCondition );
      boundary_conditions.push_back( test_example2_forconvectiontest1::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example2_forconvectiontest1::U1BoundValue );
      boundary_data.push_back( test_example2_forconvectiontest1::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = test_example2_forconvectiontest1::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(test_example2_forconvectiontest1::InitialU1);
      initialCOndtion.push_back(test_example2_forconvectiontest1::InitialU2);
      test_example2_forconvectiontest1::REYNOLDS_number = get_nu();
      test_example2_forconvectiontest1::USER_parameter1 = this->example_database["user_parameter1"];
      test_example2_forconvectiontest1::USER_parameter2 = this->example_database["user_parameter2"];



      test_example2_forconvectiontest1::ExampleFile();
      break;

    case 13:                // TestExample3 = Generate the velocity field for Convection Test Example 2 (1,2)
      /** exact_solution */
      exact_solution.push_back( test_example3_forconvectiontest2::ExactU1 );
      exact_solution.push_back( test_example3_forconvectiontest2::ExactU2 );
      exact_solution.push_back( test_example3_forconvectiontest2::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( test_example3_forconvectiontest2::BoundCondition );
      boundary_conditions.push_back( test_example3_forconvectiontest2::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example3_forconvectiontest2::U1BoundValue );
      boundary_data.push_back( test_example3_forconvectiontest2::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = test_example3_forconvectiontest2::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(test_example3_forconvectiontest2::InitialU1);
      initialCOndtion.push_back(test_example3_forconvectiontest2::InitialU2);
      test_example3_forconvectiontest2::REYNOLDS_number = get_nu();
      test_example3_forconvectiontest2::USER_parameter1 = this->example_database["user_parameter1"];
      test_example3_forconvectiontest2::USER_parameter2 = this->example_database["user_parameter2"];

      test_example3_forconvectiontest2::ExampleFile();
      break;


    case 14:                // Test4 = Navier Stokes for the test 5 of TCD2D (sort of dam break)
      /** exact_solution */
      exact_solution.push_back( test_4::ExactU1 );
      exact_solution.push_back( test_4::ExactU2 );
      exact_solution.push_back( test_4::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( test_4::BoundCondition );
      boundary_conditions.push_back( test_4::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( test_4::U1BoundValue );
      boundary_data.push_back( test_4::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = test_4::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(test_4::InitialU1);
      initialCOndtion.push_back(test_4::InitialU2);
      test_4::REYNOLDS_number = get_nu();
      test_4::USER_parameter1 = this->example_database["user_parameter1"];
      test_4::USER_parameter2 = this->example_database["user_parameter2"];

      test_4::ExampleFile();
      break;

    case 15:                // Test5 = Navier Stokes for the test 6 of TCD2D (semi-circle)
      /** exact_solution */
      exact_solution.push_back( test_5::ExactU1 );
      exact_solution.push_back( test_5::ExactU2 );
      exact_solution.push_back( test_5::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( test_5::BoundCondition );
      boundary_conditions.push_back( test_5::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( test_5::U1BoundValue );
      boundary_data.push_back( test_5::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = test_5::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(test_5::InitialU1);
      initialCOndtion.push_back(test_5::InitialU2);
      test_5::REYNOLDS_number = get_nu();
      test_5::USER_parameter1 = this->example_database["user_parameter1"];
      test_5::USER_parameter2 = this->example_database["user_parameter2"];

      test_5::ExampleFile();
      break;


    default:
      ErrThrow("Unknown Time dependent Example_TimeNSE2D example!");
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

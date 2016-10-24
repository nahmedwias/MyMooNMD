#include <Example_TimeCD2D.h>
#include <Database.h>
#include <FEFunction2D.h>

#include<string>

namespace linear_space_time
{
#include "TCD_2D/linear_space_time.h"
}
namespace exp_sin_cos
{
#include "TCD_2D/exp.h"
}
namespace sin_sin_sin
{
#include "TCD_2D/Sin3.h"
}

namespace sin_cos
{
#include "TCD_2D/SinCos1.h"
}

namespace rotating_bodies_1
{
#include "TCD_2D/Rotating_Bodies.h"
}


/***** BELOW THIS LINE, EXAMPLES ARE SPECIFIC TO USER PROJECT *****/
namespace test_example1_expsincos
{
#include "1_TestExampleTCD2D.h"
}
namespace test_example2_sinsinsin
{
#include "2_TestExampleTCD2D.h"
}
namespace test_example3_sincos
{
#include "3_TestExampleTCD2D.h"
}
namespace test_example4_linear
{
#include "4_TestExampleTCD2D.h"
}
namespace test_2phase
{
#include "5_TestTCD2D.h"
}
namespace test_semicircle
{
#include "6_TestTCD2D.h"
}

// ********* END OF USER PROJECT CODE



Example_TimeCD2D::Example_TimeCD2D(
  const ParameterDatabase& user_input_parameter_db)
 : Example_NonStationary2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case -1:
    {
      /**Exact solution"**/
      exact_solution.push_back(linear_space_time::Exact);
      /** boundary condition */
      boundary_conditions.push_back( linear_space_time::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( linear_space_time::BoundValue );
      
      /** coefficients */
      problem_coefficients = linear_space_time::BilinearCoeffs;
      
      /** Initial condition*/
     initialCOndtion.push_back(linear_space_time::InitialCondition);
     linear_space_time::ExampleFile();

     this->timeDependentRhs = linear_space_time::rhs_depends_on_time;
     this->timeDependentCoeffs=linear_space_time::coefficients_depend_on_time;
    }
    break;
    case 0:
      /**Exact solution"**/
      exact_solution.push_back(exp_sin_cos::Exact);
      /** boundary condition */
      boundary_conditions.push_back( exp_sin_cos::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( exp_sin_cos::BoundValue );
      
      /** coefficients */
      problem_coefficients = exp_sin_cos::BilinearCoeffs;
      
      /** Initial condition*/
     initialCOndtion.push_back(exp_sin_cos::InitialCondition);
     exp_sin_cos::ExampleFile();
     
     this->timeDependentRhs = exp_sin_cos::rhs_depends_on_time;
     this->timeDependentCoeffs=exp_sin_cos::coefficients_depend_on_time;
     break;
    case 1:
      /**Exact solution"**/
       exact_solution.push_back(sin_sin_sin::Exact);
      /** boundary condition */
      boundary_conditions.push_back( sin_sin_sin::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sin_sin_sin::BoundValue );
      
      /** coefficients */
      problem_coefficients = sin_sin_sin::BilinearCoeffs;
      
      /** Initial condition*/
      initialCOndtion.push_back(sin_sin_sin::InitialCondition);
      sin_sin_sin::ExampleFile();
      
      this->timeDependentRhs = sin_sin_sin::rhs_depends_on_time;
     this->timeDependentCoeffs=sin_sin_sin::coefficients_depend_on_time;
      break;
    case 2:
      /**Exact solution"**/
      exact_solution.push_back(sin_cos::Exact);

      /** boundary condition */
      boundary_conditions.push_back( sin_cos::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sin_cos::BoundValue );
      
      /** coefficients */
      problem_coefficients = sin_cos::BilinearCoeffs;
      
      /** Initial condition*/
      initialCOndtion.push_back(sin_cos::InitialCondition);
      sin_cos::ExampleFile();      
      break;
    case 3:
      /**Exact solution"**/
      exact_solution.push_back(rotating_bodies_1::Exact);

      /** boundary condition */
      boundary_conditions.push_back( rotating_bodies_1::BoundCondition );

      /** boundary values */
      boundary_data.push_back( rotating_bodies_1::BoundValue );

      /** coefficients */
      problem_coefficients = rotating_bodies_1::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(rotating_bodies_1::InitialCondition);

      // Print some example specific information.
      rotating_bodies_1::ExampleFile();
      break;






      /***** BELOW THIS LINE, EXAMPLES ARE SPECIFIC TO USER PROJECT *****/
    case 11:                // 1_TestExample Exp-Sin-Cos
      /**Exact solution"**/
      exact_solution.push_back(test_example1_expsincos::Exact);
      /** boundary condition */
      boundary_conditions.push_back( test_example1_expsincos::BoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example1_expsincos::BoundValue );

      /** coefficients */
      problem_coefficients = test_example1_expsincos::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(test_example1_expsincos::InitialCondition);
      test_example1_expsincos::ExampleFile();

      this->timeDependentRhs = test_example1_expsincos::rhs_depends_on_time;
      this->timeDependentCoeffs=test_example1_expsincos::coefficients_depend_on_time;
      break;

    case 12:                // 2_TestExample Sin-Sin-Sin
      /**Exact solution"**/
      exact_solution.push_back(test_example2_sinsinsin::Exact);
      /** boundary condition */
      boundary_conditions.push_back( test_example2_sinsinsin::BoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example2_sinsinsin::BoundValue );

      /** coefficients */
      problem_coefficients = test_example2_sinsinsin::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(test_example2_sinsinsin::InitialCondition);
      test_example2_sinsinsin::ExampleFile();

      this->timeDependentRhs = test_example2_sinsinsin::rhs_depends_on_time;
      this->timeDependentCoeffs=test_example2_sinsinsin::coefficients_depend_on_time;
      break;

    case 13:                // 3_TestExample Sin-Cos
      /**Exact solution"**/
      exact_solution.push_back(test_example3_sincos::Exact);
      /** boundary condition */
      boundary_conditions.push_back( test_example3_sincos::BoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example3_sincos::BoundValue );

      /** coefficients */
      problem_coefficients = test_example3_sincos::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(test_example3_sincos::InitialCondition);
      test_example3_sincos::ExampleFile();

      this->timeDependentRhs = test_example3_sincos::rhs_depends_on_time;
      this->timeDependentCoeffs=test_example3_sincos::coefficients_depend_on_time;
      break;

    case 14:                // 4_TestExample Linear Space Time
      /**Exact solution"**/
      exact_solution.push_back(test_example4_linear::Exact);
      /** boundary condition */
      boundary_conditions.push_back( test_example4_linear::BoundCondition );

      /** boundary values */
      boundary_data.push_back( test_example4_linear::BoundValue );

      /** coefficients */
      problem_coefficients = test_example4_linear::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(test_example4_linear::InitialCondition);
      test_example4_linear::ExampleFile();

      this->timeDependentRhs = test_example4_linear::rhs_depends_on_time;
      this->timeDependentCoeffs=test_example4_linear::coefficients_depend_on_time;
      break;

    case 15:                // 5_Test 2 phase in a square box
      /**Exact solution"**/
      exact_solution.push_back(test_2phase::Exact);
      /** boundary condition */
      boundary_conditions.push_back( test_2phase::BoundCondition );

      /** boundary values */
      boundary_data.push_back( test_2phase::BoundValue );

      /** coefficients */
      problem_coefficients = test_2phase::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(test_2phase::InitialCondition);
      test_2phase::ExampleFile();

      this->timeDependentRhs = test_2phase::rhs_depends_on_time;
      this->timeDependentCoeffs=test_2phase::coefficients_depend_on_time;
      break;

    case 16:                // 6_Test rotating semi-circle
      /**Exact solution"**/
      exact_solution.push_back(test_semicircle::Exact);
      /** boundary condition */
      boundary_conditions.push_back( test_semicircle::BoundCondition );

      /** boundary values */
      boundary_data.push_back( test_semicircle::BoundValue );

      /** coefficients */
      problem_coefficients = test_semicircle::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(test_semicircle::InitialCondition);
      test_semicircle::ExampleFile();

      this->timeDependentRhs = test_semicircle::rhs_depends_on_time;
      this->timeDependentCoeffs=test_semicircle::coefficients_depend_on_time;
      break;




    default:
      ErrThrow("Unknown name of the transient-convection-diffusion (Time_CD2D) example!", 
               example_code);
  }
}

void Example_TimeCD2D::do_post_processing(Time_CD2D& tcd2d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(tcd2d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_TimeCD2D","No post processing done for the current example.");
  }
}

double Example_TimeCD2D::get_nu() const
{
  double diffusion_coefficient = this->example_database["diffusion_coefficient"];
  return diffusion_coefficient;
}


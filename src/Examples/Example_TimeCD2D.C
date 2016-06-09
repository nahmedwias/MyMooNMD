#include <Example_TimeCD2D.h>
#include <Database.h>
#include <FEFunction2D.h>

#include<string>

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

Example_TimeCD2D::Example_TimeCD2D(int example_code)
{
  switch(example_code)
  {
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
      initialCondition.push_back(exp_sin_cos::InitialCondition);
     exp_sin_cos::ExampleFile();
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
      initialCondition.push_back(sin_sin_sin::InitialCondition);
      sin_sin_sin::ExampleFile();
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
      initialCondition.push_back(sin_cos::InitialCondition);
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
      initialCondition.push_back(rotating_bodies_1::InitialCondition);

      // Print some example specific information.
      rotating_bodies_1::ExampleFile();
      break;
    default:
      ErrThrow("Unknown name of the transient-convection-diffusion (Time_CD2D) example!", 
               example_code);
  }
}

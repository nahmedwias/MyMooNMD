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

namespace Geothermal_Energy_TCD2D
{
  #include "TCD_2D/Geothermal_Energy_TCD2D.h"
}

namespace Tube2D
{
  #include "TCD_2D/Tube2D.h"
}


//=======================================================

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
      initialCondition.push_back(linear_space_time::InitialCondition);
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
      initialCondition.push_back(exp_sin_cos::InitialCondition);
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
      initialCondition.push_back(sin_sin_sin::InitialCondition);
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
 case 4:
      /**Exact solution"**/
      exact_solution.push_back(Geothermal_Energy_TCD2D::Exact);

      /** boundary condition */
      boundary_conditions.push_back( Geothermal_Energy_TCD2D::BoundCondition );

      /** boundary values */
      boundary_data.push_back( Geothermal_Energy_TCD2D::BoundValue );

      /** coefficients */
      problem_coefficients = Geothermal_Energy_TCD2D::BilinearCoeffs;

      /** Initial condition*/
      initialCondition.push_back(Geothermal_Energy_TCD2D::InitialCondition);
 
      // Print some example specific information.
      Geothermal_Energy_TCD2D::ExampleFile();
      break;
      
  case 5:
      /**Exact solution"**/
      exact_solution.push_back(Tube2D::Exact);

      /** boundary condition */
      boundary_conditions.push_back( Tube2D::BoundCondition );

      /** boundary values */
      boundary_data.push_back( Tube2D::BoundValue );

      /** coefficients */
      problem_coefficients = Tube2D::BilinearCoeffs;

      /** Initial condition*/
      initialCondition.push_back(Tube2D::InitialCondition);

      Tube2D::diffusion_coefficient = get_nu();
      
      // Print some example specific information.
      Tube2D::ExampleFile();
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


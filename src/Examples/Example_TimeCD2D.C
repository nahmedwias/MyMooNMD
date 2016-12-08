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
namespace example10_sincos_tcd2d            // example 10
{
#include "../../user_projects/include/Examples/Time_CD2D/10_SinCos_TCD2D.h"
}
namespace example20_coupling_expsincos      //example 20
{
#include "../../user_projects/include/Examples/Time_CD2D/20_CouplingNSE_CD_ExpSinCos.h"
}
namespace example21_coupling_sinsinsin      //example 21
{
#include "../../user_projects/include/Examples/Time_CD2D/21_CouplingNSE_CD_SinSinSin.h"
}
namespace example22_coupling_semicircle       // example 2
{
#include "../../user_projects/include/Examples/Time_CD2D/22_SemiCircleNSE_CD.h"
}
namespace example30_linear_variablevisco    // example 30
{
#include "../../user_projects/include/Examples/Time_CD2D/30_Variable_viscosity_linear.h"
}
namespace example31_expon_variablevisco     // example 31
{
#include "../../user_projects/include/Examples/Time_CD2D/31_Variable_viscosity_exponential.h"
}
namespace example32_variablevisco_beltrami  // example 32
{
#include "../../user_projects/include/Examples/Time_CD2D/32_Variable_viscosity_Beltrami.h"
}
namespace example40_dambreak_cd_nse         // example 40
{
#include "../../user_projects/include/Examples/Time_CD2D/40_DamBreakCD_NSE.h"
}
namespace example42_rayleightaylor_cd_nse   // example 42
{
#include "../../user_projects/include/Examples/Time_CD2D/42_RayleighTaylorCD_NSE.h"
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
    case 10:                // Example 10 Standard TCD2D example Sin-Cos, 0 coupling
      /**Exact solution"**/
      exact_solution.push_back(example10_sincos_tcd2d::Exact);
      /** boundary condition */
      boundary_conditions.push_back( example10_sincos_tcd2d::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example10_sincos_tcd2d::BoundValue );

      /** coefficients */
      problem_coefficients = example10_sincos_tcd2d::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(example10_sincos_tcd2d::InitialCondition);
      example10_sincos_tcd2d::ExampleFile();

      this->timeDependentRhs = example10_sincos_tcd2d::rhs_depends_on_time;
      this->timeDependentCoeffs=example10_sincos_tcd2d::coefficients_depend_on_time;
      break;

    case 20:                // 20_CouplingNSE_CD Exp-Sin-Cos
      /**Exact solution"**/
      exact_solution.push_back(example20_coupling_expsincos::Exact);
      /** boundary condition */
      boundary_conditions.push_back( example20_coupling_expsincos::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example20_coupling_expsincos::BoundValue );

      /** coefficients */
      problem_coefficients = example20_coupling_expsincos::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(example20_coupling_expsincos::InitialCondition);
      example20_coupling_expsincos::ExampleFile();

      this->timeDependentRhs = example20_coupling_expsincos::rhs_depends_on_time;
      this->timeDependentCoeffs=example20_coupling_expsincos::coefficients_depend_on_time;
      break;

    case 21:                // 21_CouplingNSE_CD Sin-Sin-Sin
      /**Exact solution"**/
      exact_solution.push_back(example21_coupling_sinsinsin::Exact);
      /** boundary condition */
      boundary_conditions.push_back( example21_coupling_sinsinsin::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example21_coupling_sinsinsin::BoundValue );

      /** coefficients */
      problem_coefficients = example21_coupling_sinsinsin::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(example21_coupling_sinsinsin::InitialCondition);
      example21_coupling_sinsinsin::ExampleFile();

      this->timeDependentRhs = example21_coupling_sinsinsin::rhs_depends_on_time;
      this->timeDependentCoeffs=example21_coupling_sinsinsin::coefficients_depend_on_time;
      break;

    case 22:                // Example 22: 1 way coupling: rotating semi-circle
      /**Exact solution"**/
      exact_solution.push_back(example22_coupling_semicircle::Exact);
      /** boundary condition */
      boundary_conditions.push_back( example22_coupling_semicircle::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example22_coupling_semicircle::BoundValue );

      /** coefficients */
      problem_coefficients = example22_coupling_semicircle::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(example22_coupling_semicircle::InitialCondition);
      example22_coupling_semicircle::ExampleFile();

      this->timeDependentRhs = example22_coupling_semicircle::rhs_depends_on_time;
      this->timeDependentCoeffs=example22_coupling_semicircle::coefficients_depend_on_time;
      break;

    case 30:                // Example 30 CouplingCD_NSE_ Poiseuille Variable Viscosity Linear
      /**Exact solution"**/
      exact_solution.push_back(example30_linear_variablevisco::Exact);
      /** boundary condition */
      boundary_conditions.push_back( example30_linear_variablevisco::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example30_linear_variablevisco::BoundValue );

      /** coefficients */
      problem_coefficients = example30_linear_variablevisco::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(example30_linear_variablevisco::InitialCondition);
      example30_linear_variablevisco::ExampleFile();

      this->timeDependentRhs = example30_linear_variablevisco::rhs_depends_on_time;
      this->timeDependentCoeffs=example30_linear_variablevisco::coefficients_depend_on_time;
      break;

    case 31:                // Example 31: Coupling CD>NSE Variable Viscosity Exponential
      /**Exact solution"**/
      exact_solution.push_back(example31_expon_variablevisco::Exact);
      /** boundary condition */
      boundary_conditions.push_back( example31_expon_variablevisco::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example31_expon_variablevisco::BoundValue );

      /** coefficients */
      problem_coefficients = example31_expon_variablevisco::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(example31_expon_variablevisco::InitialCondition);
      example31_expon_variablevisco::ExampleFile();

      this->timeDependentRhs = example31_expon_variablevisco::rhs_depends_on_time;
      this->timeDependentCoeffs=example31_expon_variablevisco::coefficients_depend_on_time;
      break;

    case 32:                // Example 32: Coupling CD>NSE Variable Viscosity Beltrami
      /**Exact solution"**/
      exact_solution.push_back(example32_variablevisco_beltrami::Exact);
      /** boundary condition */
      boundary_conditions.push_back( example32_variablevisco_beltrami::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example32_variablevisco_beltrami::BoundValue );

      /** coefficients */
      problem_coefficients = example32_variablevisco_beltrami::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(example32_variablevisco_beltrami::InitialCondition);
      example32_variablevisco_beltrami::ExampleFile();

      this->timeDependentRhs = example32_variablevisco_beltrami::rhs_depends_on_time;
      this->timeDependentCoeffs=example32_variablevisco_beltrami::coefficients_depend_on_time;
      break;

    case 40:                // Example 40: 2 way-coupling for Dam Break (2 phases in a square box)
      /**Exact solution"**/
      exact_solution.push_back(example40_dambreak_cd_nse::Exact);
      /** boundary condition */
      boundary_conditions.push_back( example40_dambreak_cd_nse::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example40_dambreak_cd_nse::BoundValue );

      /** coefficients */
      problem_coefficients = example40_dambreak_cd_nse::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(example40_dambreak_cd_nse::InitialCondition);
      example40_dambreak_cd_nse::ExampleFile();

      this->timeDependentRhs = example40_dambreak_cd_nse::rhs_depends_on_time;
      this->timeDependentCoeffs=example40_dambreak_cd_nse::coefficients_depend_on_time;
      break;


    case 42:                // Example 42: 2-WAY-COUPLING for Rayleigh-Taylor Instability
      /**Exact solution"**/
      exact_solution.push_back(example42_rayleightaylor_cd_nse::Exact);
      /** boundary condition */
      boundary_conditions.push_back( example42_rayleightaylor_cd_nse::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example42_rayleightaylor_cd_nse::BoundValue );

      /** coefficients */
      problem_coefficients = example42_rayleightaylor_cd_nse::BilinearCoeffs;

      /** Initial condition*/
      initialCOndtion.push_back(example42_rayleightaylor_cd_nse::InitialCondition);
      example42_rayleightaylor_cd_nse::ExampleFile();

      this->timeDependentRhs = example42_rayleightaylor_cd_nse::rhs_depends_on_time;
      this->timeDependentCoeffs=example42_rayleightaylor_cd_nse::coefficients_depend_on_time;
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


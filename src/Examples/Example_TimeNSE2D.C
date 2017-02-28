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
namespace example10_sincos_tnse2d   //
{
#include "../../user_projects/include/Examples/Time_NSE2D/10_SinCos_TNSE2D.h"
}
namespace example20_coupling_nse_cd // test coupling NSE>CD = Generate the velocity field (1,-1)
{
#include "../../user_projects/include/Examples/Time_NSE2D/20_CouplingNSE_CD.h"
}
namespace example21_coupling_nse_cd // test coupling NSE>CD = Generate the velocity field (1,2)
{
#include "../../user_projects/include/Examples/Time_NSE2D/21_CouplingNSE_CD.h"
}
namespace example22_semicircle_nse_cd // 1 way coupling for rotating semi circle
{
#include "../../user_projects/include/Examples/Time_NSE2D/22_SemiCircleNSE_CD.h"
}
namespace example30_poiseuille_variablevisco // Coupling CD>NSE Poiseuille with variable viscosity
{
#include "../../user_projects/include/Examples/Time_NSE2D/30_CouplingCD_NSE_Poiseuille_variableviscosity.h"
}
namespace example31_poiseuille_variablevisco // Coupling CD>NSE Poiseuille with variable viscosity exponential
{
#include "../../user_projects/include/Examples/Time_NSE2D/31_CouplingCD_NSE_Exponential_variableviscosity.h"
}
namespace example32_SinCosExp // Coupling CD>NSE with SinCosExp
{
#include "../../user_projects/include/Examples/Time_NSE2D/32_CouplingCD_NSE_SinCosExp.h"
}
namespace example40_dambreak_nse_cd // 2 way coupling for dam break
{
#include "../../user_projects/include/Examples/Time_NSE2D/40_DamBreakNSE_CD.h"
}
namespace example41_rayleightaylor_nse_cd   // Rayleigh-Taylor instability (Fraigneau 2001)
{
#include "../../user_projects/include/Examples/Time_NSE2D/41_RayleighTaylorNSE_CD.h"
}
namespace example42_rayleightaylor2_nse_cd   // Rayleigh-Taylor instability (Pochet 2013)
{
#include "../../user_projects/include/Examples/Time_NSE2D/42_RayleighTaylor2NSE_CD.h"
}
namespace example50_gasstirring_nse_cd   // 2 way coupling for Gas Stirring test
{
#include "../../user_projects/include/Examples/Time_NSE2D/50_GasStirringTestNSE_CD.h"
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
    case 10:                // Bsp1 = sincos  = SAME AS IN TNSE2D
      /** exact_solution */
      exact_solution.push_back( example10_sincos_tnse2d::ExactU1 );
      exact_solution.push_back( example10_sincos_tnse2d::ExactU2 );
      exact_solution.push_back( example10_sincos_tnse2d::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( example10_sincos_tnse2d::BoundCondition );
      boundary_conditions.push_back( example10_sincos_tnse2d::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( example10_sincos_tnse2d::U1BoundValue );
      boundary_data.push_back( example10_sincos_tnse2d::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = example10_sincos_tnse2d::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(example10_sincos_tnse2d::InitialU1);
      initialCOndtion.push_back(example10_sincos_tnse2d::InitialU2);
      example10_sincos_tnse2d::REYNOLDS_number = get_nu();
      example10_sincos_tnse2d::USER_parameter1 = this->example_database["user_parameter1"];
      example10_sincos_tnse2d::USER_parameter2 = this->example_database["user_parameter2"];

      example10_sincos_tnse2d::ExampleFile();
      break;

    case 20:                // Test coupling NSE>CD = Generate the velocity field (1,-1) for Convection example 20
      /** exact_solution */
      exact_solution.push_back( example20_coupling_nse_cd::ExactU1 );
      exact_solution.push_back( example20_coupling_nse_cd::ExactU2 );
      exact_solution.push_back( example20_coupling_nse_cd::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( example20_coupling_nse_cd::BoundCondition );
      boundary_conditions.push_back( example20_coupling_nse_cd::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( example20_coupling_nse_cd::U1BoundValue );
      boundary_data.push_back( example20_coupling_nse_cd::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = example20_coupling_nse_cd::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(example20_coupling_nse_cd::InitialU1);
      initialCOndtion.push_back(example20_coupling_nse_cd::InitialU2);
      example20_coupling_nse_cd::REYNOLDS_number = get_nu();
      example20_coupling_nse_cd::USER_parameter1 = this->example_database["user_parameter1"];
      example20_coupling_nse_cd::USER_parameter2 = this->example_database["user_parameter2"];

      example20_coupling_nse_cd::ExampleFile();
      break;

    case 21:                // CouplingNSE_CD = Generate the velocity field (1,2) for Convection Example 21
      /** exact_solution */
      exact_solution.push_back( example21_coupling_nse_cd::ExactU1 );
      exact_solution.push_back( example21_coupling_nse_cd::ExactU2 );
      exact_solution.push_back( example21_coupling_nse_cd::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( example21_coupling_nse_cd::BoundCondition );
      boundary_conditions.push_back( example21_coupling_nse_cd::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( example21_coupling_nse_cd::U1BoundValue );
      boundary_data.push_back( example21_coupling_nse_cd::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = example21_coupling_nse_cd::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(example21_coupling_nse_cd::InitialU1);
      initialCOndtion.push_back(example21_coupling_nse_cd::InitialU2);
      example21_coupling_nse_cd::REYNOLDS_number = get_nu();
      example21_coupling_nse_cd::USER_parameter1 = this->example_database["user_parameter1"];
      example21_coupling_nse_cd::USER_parameter2 = this->example_database["user_parameter2"];

      example21_coupling_nse_cd::ExampleFile();
      break;

    case 22:                // Example 22= 1 way coupling for rotating semi-circle
      /** exact_solution */
      exact_solution.push_back( example22_semicircle_nse_cd::ExactU1 );
      exact_solution.push_back( example22_semicircle_nse_cd::ExactU2 );
      exact_solution.push_back( example22_semicircle_nse_cd::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( example22_semicircle_nse_cd::BoundCondition );
      boundary_conditions.push_back( example22_semicircle_nse_cd::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( example22_semicircle_nse_cd::U1BoundValue );
      boundary_data.push_back( example22_semicircle_nse_cd::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = example22_semicircle_nse_cd::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(example22_semicircle_nse_cd::InitialU1);
      initialCOndtion.push_back(example22_semicircle_nse_cd::InitialU2);
      example22_semicircle_nse_cd::REYNOLDS_number = get_nu();
      example22_semicircle_nse_cd::USER_parameter1 = this->example_database["user_parameter1"];
      example22_semicircle_nse_cd::USER_parameter2 = this->example_database["user_parameter2"];

      example22_semicircle_nse_cd::ExampleFile();
      break;

    case 30:                // Coupling CD_NSE Poiseuille with variable viscosity
      /** exact_solution */
      exact_solution.push_back( example30_poiseuille_variablevisco::ExactU1 );
      exact_solution.push_back( example30_poiseuille_variablevisco::ExactU2 );
      exact_solution.push_back( example30_poiseuille_variablevisco::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( example30_poiseuille_variablevisco::BoundCondition );
      boundary_conditions.push_back( example30_poiseuille_variablevisco::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( example30_poiseuille_variablevisco::U1BoundValue );
      boundary_data.push_back( example30_poiseuille_variablevisco::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = example30_poiseuille_variablevisco::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(example30_poiseuille_variablevisco::InitialU1);
      initialCOndtion.push_back(example30_poiseuille_variablevisco::InitialU2);
      example30_poiseuille_variablevisco::REYNOLDS_number = get_nu();
      example30_poiseuille_variablevisco::USER_parameter1 = this->example_database["user_parameter1"];
      example30_poiseuille_variablevisco::USER_parameter2 = this->example_database["user_parameter2"];

      example30_poiseuille_variablevisco::ExampleFile();
      break;


    case 31:                // Example31 Coupling CD>NSE = Poiseuille with variable viscosity exponential
      /** exact_solution */
      exact_solution.push_back( example31_poiseuille_variablevisco::ExactU1 );
      exact_solution.push_back( example31_poiseuille_variablevisco::ExactU2 );
      exact_solution.push_back( example31_poiseuille_variablevisco::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( example31_poiseuille_variablevisco::BoundCondition );
      boundary_conditions.push_back( example31_poiseuille_variablevisco::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( example31_poiseuille_variablevisco::U1BoundValue );
      boundary_data.push_back( example31_poiseuille_variablevisco::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = example31_poiseuille_variablevisco::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(example31_poiseuille_variablevisco::InitialU1);
      initialCOndtion.push_back(example31_poiseuille_variablevisco::InitialU2);
      example31_poiseuille_variablevisco::REYNOLDS_number = get_nu();
      example31_poiseuille_variablevisco::USER_parameter1 = this->example_database["user_parameter1"];
      example31_poiseuille_variablevisco::USER_parameter2 = this->example_database["user_parameter2"];

      example31_poiseuille_variablevisco::ExampleFile();
      break;

    case 32:                // Example32 Coupling CD>NSE = SinCosExp
      /** exact_solution */
      exact_solution.push_back( example32_SinCosExp::ExactU1 );
      exact_solution.push_back( example32_SinCosExp::ExactU2 );
      exact_solution.push_back( example32_SinCosExp::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( example32_SinCosExp::BoundCondition );
      boundary_conditions.push_back( example32_SinCosExp::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( example32_SinCosExp::U1BoundValue );
      boundary_data.push_back( example32_SinCosExp::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = example32_SinCosExp::LinCoeffs;

      /** initial condition */
      initialCOndtion.push_back(example32_SinCosExp::InitialU1);
      initialCOndtion.push_back(example32_SinCosExp::InitialU2);
//      example32_SinCosExp::REYNOLDS_number = get_nu();
//      example32_SinCosExp::USER_parameter1 = this->example_database["user_parameter1"];
//      example32_SinCosExp::USER_parameter2 = this->example_database["user_parameter2"];

      example32_SinCosExp::ExampleFile();
      break;


    case 40:                // Example40 = 2-way-coupling :sort of dam break
         /** exact_solution */
         exact_solution.push_back( example40_dambreak_nse_cd::ExactU1 );
         exact_solution.push_back( example40_dambreak_nse_cd::ExactU2 );
         exact_solution.push_back( example40_dambreak_nse_cd::ExactP );

         /** boundary condition */
         boundary_conditions.push_back( example40_dambreak_nse_cd::BoundCondition );
         boundary_conditions.push_back( example40_dambreak_nse_cd::BoundCondition );
         boundary_conditions.push_back( BoundConditionNoBoundCondition );

         /** boundary values */
         boundary_data.push_back( example40_dambreak_nse_cd::U1BoundValue );
         boundary_data.push_back( example40_dambreak_nse_cd::U2BoundValue );
         boundary_data.push_back( BoundaryValueHomogenous );

         /** coefficients */
         problem_coefficients = example40_dambreak_nse_cd::LinCoeffs;

         /** initial condition */
         initialCOndtion.push_back(example40_dambreak_nse_cd::InitialU1);
         initialCOndtion.push_back(example40_dambreak_nse_cd::InitialU2);
         example40_dambreak_nse_cd::REYNOLDS_number = get_nu();
         example40_dambreak_nse_cd::USER_parameter1 = this->example_database["user_parameter1"];
         example40_dambreak_nse_cd::USER_parameter2 = this->example_database["user_parameter2"];

         example40_dambreak_nse_cd::ExampleFile();
         break;

       case 41:                // Example 41 = 2-WAY-COUPLING FOR RAYLEIGH TAYLOR INSTABILITY
         /** exact_solution */
         exact_solution.push_back( example41_rayleightaylor_nse_cd::ExactU1 );
         exact_solution.push_back( example41_rayleightaylor_nse_cd::ExactU2 );
         exact_solution.push_back( example41_rayleightaylor_nse_cd::ExactP );

         /** boundary condition */
         boundary_conditions.push_back( example41_rayleightaylor_nse_cd::BoundCondition );
         boundary_conditions.push_back( example41_rayleightaylor_nse_cd::BoundCondition );
         boundary_conditions.push_back( BoundConditionNoBoundCondition );

         /** boundary values */
         boundary_data.push_back( example41_rayleightaylor_nse_cd::U1BoundValue );
         boundary_data.push_back( example41_rayleightaylor_nse_cd::U2BoundValue );
         boundary_data.push_back( BoundaryValueHomogenous );

         /** coefficients */
         problem_coefficients = example41_rayleightaylor_nse_cd::LinCoeffs;

         /** initial condition */
         initialCOndtion.push_back(example41_rayleightaylor_nse_cd::InitialU1);
         initialCOndtion.push_back(example41_rayleightaylor_nse_cd::InitialU2);
         example41_rayleightaylor_nse_cd::REYNOLDS_number = get_nu();
         example41_rayleightaylor_nse_cd::USER_parameter1 = this->example_database["user_parameter1"];
         example41_rayleightaylor_nse_cd::USER_parameter2 = this->example_database["user_parameter2"];

         example41_rayleightaylor_nse_cd::ExampleFile();
         break;

       case 42:                // Example 42 = 2-WAY-COUPLING FOR RAYLEIGH TAYLOR INSTABILITY
                               // From Pochet et al 2013
         /** exact_solution */
         exact_solution.push_back( example42_rayleightaylor2_nse_cd::ExactU1 );
         exact_solution.push_back( example42_rayleightaylor2_nse_cd::ExactU2 );
         exact_solution.push_back( example42_rayleightaylor2_nse_cd::ExactP );

         /** boundary condition */
         boundary_conditions.push_back( example42_rayleightaylor2_nse_cd::BoundCondition );
         boundary_conditions.push_back( example42_rayleightaylor2_nse_cd::BoundCondition );
         boundary_conditions.push_back( BoundConditionNoBoundCondition );

         /** boundary values */
         boundary_data.push_back( example42_rayleightaylor2_nse_cd::U1BoundValue );
         boundary_data.push_back( example42_rayleightaylor2_nse_cd::U2BoundValue );
         boundary_data.push_back( BoundaryValueHomogenous );

         /** coefficients */
         problem_coefficients = example42_rayleightaylor2_nse_cd::LinCoeffs;

         /** initial condition */
         initialCOndtion.push_back(example42_rayleightaylor2_nse_cd::InitialU1);
         initialCOndtion.push_back(example42_rayleightaylor2_nse_cd::InitialU2);
         example42_rayleightaylor2_nse_cd::REYNOLDS_number = get_nu();
         example42_rayleightaylor2_nse_cd::USER_parameter1 = this->example_database["user_parameter1"];
         example42_rayleightaylor2_nse_cd::USER_parameter2 = this->example_database["user_parameter2"];

         example42_rayleightaylor2_nse_cd::ExampleFile();
         break;

       case 50:                // Example 50 = 2-WAY-COUPLING FOR Gas Stirring Test
         /** exact_solution */
         exact_solution.push_back( example50_gasstirring_nse_cd::ExactU1 );
         exact_solution.push_back( example50_gasstirring_nse_cd::ExactU2 );
         exact_solution.push_back( example50_gasstirring_nse_cd::ExactP );

         /** boundary condition */
         boundary_conditions.push_back( example50_gasstirring_nse_cd::BoundCondition );
         boundary_conditions.push_back( example50_gasstirring_nse_cd::BoundCondition );
         boundary_conditions.push_back( BoundConditionNoBoundCondition );

         /** boundary values */
         boundary_data.push_back( example50_gasstirring_nse_cd::U1BoundValue );
         boundary_data.push_back( example50_gasstirring_nse_cd::U2BoundValue );
         boundary_data.push_back( BoundaryValueHomogenous );

         /** coefficients */
         problem_coefficients = example50_gasstirring_nse_cd::LinCoeffs;

         /** initial condition */
         initialCOndtion.push_back(example50_gasstirring_nse_cd::InitialU1);
         initialCOndtion.push_back(example50_gasstirring_nse_cd::InitialU2);
         example50_gasstirring_nse_cd::REYNOLDS_number = get_nu();
         example50_gasstirring_nse_cd::USER_parameter1 = this->example_database["user_parameter1"];
         example50_gasstirring_nse_cd::USER_parameter2 = this->example_database["user_parameter2"];

         example50_gasstirring_nse_cd::ExampleFile();
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

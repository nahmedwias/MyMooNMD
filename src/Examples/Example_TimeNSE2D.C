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

namespace  cosine_sin {
#include "TNSE_2D/cosine_sine.h"
}

namespace potential_flow_ex1
{
#include "TNSE_2D/potential_flow_td.h"
}

namespace potential_flow_ex2
{
#include "TNSE_2D/potential_flow_td_ex2.h"
}

namespace potential_flow_ex3 {
#include "TNSE_2D/potential_flow_ex3.h"
}

namespace potential_flow_ex4 {
#include "TNSE_2D/potential_flow_ex4.h"
}

namespace potential_flow_ex5 {
#include "TNSE_2D/potential_flow_td_ex5.h"
}

namespace potential_flow_ex6 {
#include "TNSE_2D/potential_flow_td_ex6.h"
}

namespace cosine_sin_poly
{
#include "TNSE_2D/cosine_sine_poly.h"
}

namespace cori_betaplane
{
#include "TNSE_2D/coriolis_betaplane.h"
}

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

      bsp1::DIMENSIONLESS_VISCOSITY=this->get_nu();
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
      
      lin_space_time::DIMENSIONLESS_VISCOSITY=this->get_nu();
      lin_space_time::ExampleFile();
      break;
    case 2:
      /** exact_solution */
      exact_solution.push_back( cosine_sin::ExactU1 );
      exact_solution.push_back( cosine_sin::ExactU2 );
      exact_solution.push_back( cosine_sin::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( cosine_sin::BoundCondition );
      boundary_conditions.push_back( cosine_sin::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( cosine_sin::U1BoundValue );
      boundary_data.push_back( cosine_sin::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** initial conditions, in case of a non-stationary problem */
      initialCOndtion.push_back(cosine_sin::InitialU1);
      initialCOndtion.push_back(cosine_sin::InitialU2);
      initialCOndtion.push_back(cosine_sin::InitialP);
      /** coefficients */
      problem_coefficients = cosine_sin::LinCoeffs;

      cosine_sin::DIMENSIONLESS_VISCOSITY=this->get_nu();

      cosine_sin::ExampleFile();
      break;
  case 3:
      exact_solution.push_back( potential_flow_ex1::ExactU1 );
      exact_solution.push_back( potential_flow_ex1::ExactU2 );
      exact_solution.push_back( potential_flow_ex1::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( potential_flow_ex1::BoundCondition );
      boundary_conditions.push_back( potential_flow_ex1::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( potential_flow_ex1::U1BoundValue );
      boundary_data.push_back( potential_flow_ex1::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** initial conditions, in case of a non-stationary problem */
      initialCOndtion.push_back(potential_flow_ex1::InitialU1);
      initialCOndtion.push_back(potential_flow_ex1::InitialU2);
      initialCOndtion.push_back(potential_flow_ex1::InitialP);
      /** coefficients */
      problem_coefficients = potential_flow_ex1::LinCoeffs;

      potential_flow_ex1::DIMENSIONLESS_VISCOSITY=this->get_nu();
      potential_flow_ex1::ExampleFile();
      break;
  case 4:
      /** exact_solution */
      exact_solution.push_back( potential_flow_ex2::ExactU1 );
      exact_solution.push_back( potential_flow_ex2::ExactU2 );
      exact_solution.push_back( potential_flow_ex2::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( potential_flow_ex2::BoundCondition );
      boundary_conditions.push_back( potential_flow_ex2::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( potential_flow_ex2::U1BoundValue );
      boundary_data.push_back( potential_flow_ex2::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** initial conditions, in case of a non-stationary problem */
      initialCOndtion.push_back(potential_flow_ex2::InitialU1);
      initialCOndtion.push_back(potential_flow_ex2::InitialU2);
      initialCOndtion.push_back(potential_flow_ex2::InitialP);
      /** coefficients */
      problem_coefficients = potential_flow_ex2::LinCoeffs;

      potential_flow_ex2::DIMENSIONLESS_VISCOSITY=this->get_nu();
      potential_flow_ex2::ExampleFile();
      break;
  case 5:
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

      /** initial conditions, in case of a non-stationary problem */
      initialCOndtion.push_back(potential_flow_ex3::InitialU1);
      initialCOndtion.push_back(potential_flow_ex3::InitialU2);
      initialCOndtion.push_back(potential_flow_ex3::InitialP);
      /** coefficients */
      problem_coefficients = potential_flow_ex3::LinCoeffs;

      potential_flow_ex3::DIMENSIONLESS_VISCOSITY=this->get_nu();
      potential_flow_ex3::ExampleFile();
      break;
    case 6:
    {
      using namespace potential_flow_ex4;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** initial conditions, in case of a non-stationary problem */
      initialCOndtion.push_back(InitialU1);
      initialCOndtion.push_back(InitialU2);
      initialCOndtion.push_back(InitialP);
      /** coefficients */
      problem_coefficients = potential_flow_ex4::LinCoeffs;

      potential_flow_ex4::DIMENSIONLESS_VISCOSITY=this->get_nu();
      ExampleFile();
      break;
    }
    case 7:
    {
      using namespace potential_flow_ex5;
      /** exact_solution */
      exact_solution.push_back(ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** initial conditions, in case of a non-stationary problem */
      initialCOndtion.push_back(InitialU1);
      initialCOndtion.push_back(InitialU2);
      initialCOndtion.push_back(InitialP);
      /** coefficients */
      problem_coefficients = potential_flow_ex5::LinCoeffs;

      potential_flow_ex5::DIMENSIONLESS_VISCOSITY=this->get_nu();
      ExampleFile();
      break;
    }
    case 8:
    {
      using namespace potential_flow_ex6;
      /** exact_solution */
      exact_solution.push_back(ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** initial conditions, in case of a non-stationary problem */
      initialCOndtion.push_back(InitialU1);
      initialCOndtion.push_back(InitialU2);
      initialCOndtion.push_back(InitialP);
      /** coefficients */
      problem_coefficients = potential_flow_ex6::LinCoeffs;

      potential_flow_ex6::DIMENSIONLESS_VISCOSITY=this->get_nu();
      ExampleFile();
      break;
    }
    case 9:
    {
      using namespace cosine_sin_poly;
      exact_solution.push_back(ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** initial conditions, in case of a non-stationary problem */
      initialCOndtion.push_back(InitialU1);
      initialCOndtion.push_back(InitialU2);
      initialCOndtion.push_back(InitialP);
      /** coefficients */
      problem_coefficients = cosine_sin_poly::LinCoeffs;

      cosine_sin_poly::DIMENSIONLESS_VISCOSITY=this->get_nu();
      ExampleFile();
      break;
    }
    case 10:
    {
      using namespace cori_betaplane;
      exact_solution.push_back(ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** initial conditions, in case of a non-stationary problem */
      initialCOndtion.push_back(InitialU1);
      initialCOndtion.push_back(InitialU2);
      initialCOndtion.push_back(InitialP);
      /** coefficients */
      problem_coefficients = cori_betaplane::LinCoeffs;

      cori_betaplane::DIMENSIONLESS_VISCOSITY=this->get_nu();
      ExampleFile();
      break;
    }
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

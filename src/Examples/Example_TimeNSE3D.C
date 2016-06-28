#include <Example_TimeNSE3D.h>
#include <Time_NSE3D.h>
#include <FEDatabase3D.h>
#include <Database.h>
#include <MainUtilities.h>

namespace lin_space_time
{
  #include "TNSE_3D/linear_space_time.h"  // 0
}
namespace AnsatzLinConst
{
  #include "TNSE_3D/AnsatzLinConst.h"     // 1
}
namespace Bsp0
{
  #include "TNSE_3D/Bsp0.h"   // 2
}
namespace Bsp1
{
  #include "TNSE_3D/Bsp1.h"   // 3
}
namespace Bsp2
{
 #include "TNSE_3D/Bsp2.h"    // 4
}
namespace Bsp3
{
  #include "TNSE_3D/Bsp3.h"   // 5
}
namespace tunnel_tdp1
{
#include "TNSE_3D/Tunnel1.h"  // 6
}
namespace flow_around_cylinder_instationary
{
#include "TNSE_3D/FlowAroundCylinder_instat.h"   // 7
}

//=========================================================
Example_TimeNSE3D::Example_TimeNSE3D(int example_code,
                                     const ParameterDatabase& user_input_parameter_db)
{
  this->example_database.merge(user_input_parameter_db,false);
  switch(example_code)
  {
    case 0:
    {
      using namespace lin_space_time;
      /** exact_solution */
      exact_solution.push_back( lin_space_time::ExactU1 );
      exact_solution.push_back( lin_space_time::ExactU2 );
      exact_solution.push_back( lin_space_time::ExactU3 );
      exact_solution.push_back( lin_space_time::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( lin_space_time::U1BoundValue );
      boundary_data.push_back( lin_space_time::U2BoundValue );
      boundary_data.push_back( lin_space_time::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = lin_space_time::LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( lin_space_time::InitialU1 );
      initialCondtion.push_back( lin_space_time::InitialU2 );
      initialCondtion.push_back( lin_space_time::InitialU3 );

      ExampleFile();
      break;
    }
    case 1:
    {
      using namespace AnsatzLinConst;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    case 2:
    {
      using namespace Bsp0;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    case 3:
    {
      using namespace Bsp1;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }    
    case 4:
    {
      using namespace Bsp2;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    case 5:
    {
      using namespace Bsp3;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    case 6:
    {
      using namespace tunnel_tdp1;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    case 7:
    {
      using namespace flow_around_cylinder_instationary;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      /**post processing - drag and lift calculation and output */
      post_processing_stat = compute_drag_lift_pdiff;

      nu = flow_around_cylinder_instationary::get_nu();
      ExampleFile();
      break;
    }
    default:
      ErrThrow("Unknown Example_TimeNSE3D example!");
  }
}

void Example_TimeNSE3D::do_post_processing(Time_NSE3D& tnse3d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(tnse3d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_TimeNSE3D","No post processing done for the current example.");
  }
}

double Example_TimeNSE3D::get_nu() const
{
  if(nu==-1)
    ErrThrow("Kinematic viscosity is not set in this example!");
  
  return nu;
}


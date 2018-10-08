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

namespace flow_around_cylinder_instationary
{
#include "TNSE_3D/FlowAroundCylinder_instat.h"   // 6
}

namespace cylinder
{
#include "TNSE_3D/Cylinder.h"   // 7
}

//=========================================================
Example_TimeNSE3D::Example_TimeNSE3D(
  const ParameterDatabase& user_input_parameter_db)
 : Example_NonStationary3D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case 0:
    {
      using namespace lin_space_time;
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

      /** some variables to change values in the example */
      lin_space_time::DIMENSIONLESS_VISCOSITY = this->get_nu();

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

      /** some variables to change values in the example */
      AnsatzLinConst::DIMENSIONLESS_VISCOSITY = this->get_nu();

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

      /** some variables to change values in the example */
      Bsp0::DIMENSIONLESS_VISCOSITY = this->get_nu();

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

      /** some variables to change values in the example */
      Bsp1::DIMENSIONLESS_VISCOSITY = this->get_nu();

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

      /** some variables to change values in the example */
      Bsp2::DIMENSIONLESS_VISCOSITY = this->get_nu();

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

      /** some variables to change values in the example */
      Bsp3::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 6:
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
      post_processing_stat = flow_around_cylinder_instationary::compute_drag_lift_pdiff;

      /** some variables to change values in the example */
      flow_around_cylinder_instationary::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 7:
    {
      using namespace cylinder;
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

 
      /** some variables to change values in the example */
      cylinder::DIMENSIONLESS_VISCOSITY = this->get_nu();

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
  double inverse_reynolds = this->example_database["reynolds_number"];
  inverse_reynolds = 1/inverse_reynolds;
  return inverse_reynolds;
}


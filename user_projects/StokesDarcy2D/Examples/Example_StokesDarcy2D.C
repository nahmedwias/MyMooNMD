#include <Example_StokesDarcy2D.h>

#include <FEDatabase2D.h>
#include <Database.h>
#include <MainUtilities.h>


/* examples */

namespace polynomial_ut0
{
  #include "Examples/Polynomial_ut0.h"
}
namespace sin_cos_polynomial_BJS
{
  #include "Examples/SinCosPolynomial_BJS.h"
}
namespace riverbed
{
  #include "Examples/riverbed.h"
}
namespace  ana_sol_Stokes_Darcy3
{
  #include "Examples/ana_sol_Stokes_Darcy.h"
}
namespace porousObstacle
{
#include "Examples/porousObstacle.h"
}

Example_StokesDarcy2D::Example_StokesDarcy2D() : Example2D()
{
  switch( TDatabase::ParamDB->EXAMPLE ) 
  {
    case 0:
      /** exact_solution */
      exact_solution.push_back( sin_cos_polynomial_BJS::ExactU1_NSE );
      exact_solution.push_back( sin_cos_polynomial_BJS::ExactU2_NSE );
      exact_solution.push_back( sin_cos_polynomial_BJS::ExactP_NSE );
      // only mixed
      exact_solution.push_back( sin_cos_polynomial_BJS::ExactU1_Darcy ); 
      // only mixed
      exact_solution.push_back( sin_cos_polynomial_BJS::ExactU2_Darcy ); 
      exact_solution.push_back( sin_cos_polynomial_BJS::ExactP_Darcy );
  
      /** boundary condition */
      boundary_conditions.push_back( sin_cos_polynomial_BJS::BoundCondition_NSE );
      boundary_conditions.push_back( sin_cos_polynomial_BJS::BoundCondition_NSE );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      boundary_conditions.push_back( sin_cos_polynomial_BJS::BC_Velocity_Darcy);//mixed
      boundary_conditions.push_back( BoundConditionNoBoundCondition ); // mixed
      //primal
      boundary_conditions.push_back( sin_cos_polynomial_BJS::BC_Darcy );
  
      /** boundary values */
      boundary_data.push_back( sin_cos_polynomial_BJS::U1BoundValue_NSE );
      boundary_data.push_back( sin_cos_polynomial_BJS::U2BoundValue_NSE );
      boundary_data.push_back( BoundaryValueHomogenous );
      //mixed
      boundary_data.push_back(sin_cos_polynomial_BJS::BoundValueVelocity_Darcy);
      boundary_data.push_back( BoundaryValueHomogenous ); // mixed
      // primal
      boundary_data.push_back( sin_cos_polynomial_BJS::BoundValue_Darcy ); 
      
      /** coefficients */
      problem_coefficients = sin_cos_polynomial_BJS::LinCoeffs_NSE;
      mixed_darcy_coeffs   = sin_cos_polynomial_BJS::LinCoeffs_DarcyMixed;
      primal_darcy_coeffs  = sin_cos_polynomial_BJS::LinCoeffs_Darcy;
      
      sin_cos_polynomial_BJS::ExampleFile();
      break;
    case 1:
      /** exact_solution */
      exact_solution.push_back( polynomial_ut0::ExactU1_NSE );
      exact_solution.push_back( polynomial_ut0::ExactU2_NSE );
      exact_solution.push_back( polynomial_ut0::ExactP_NSE );
      exact_solution.push_back( polynomial_ut0::ExactU1_Darcy ); // only mixed
      exact_solution.push_back( polynomial_ut0::ExactU2_Darcy ); // only mixed
      exact_solution.push_back( polynomial_ut0::ExactP_Darcy );
  
      /** boundary condition */
      boundary_conditions.push_back( polynomial_ut0::BoundCondition_NSE );
      boundary_conditions.push_back( polynomial_ut0::BoundCondition_NSE );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      boundary_conditions.push_back( polynomial_ut0::BC_Velocity_Darcy);//mixed
      boundary_conditions.push_back( BoundConditionNoBoundCondition ); // mixed
      boundary_conditions.push_back( polynomial_ut0::BC_Darcy ); // primal
  
      /** boundary values */
      boundary_data.push_back( polynomial_ut0::U1BoundValue_NSE );
      boundary_data.push_back( polynomial_ut0::U2BoundValue_NSE );
      boundary_data.push_back( BoundaryValueHomogenous );
      boundary_data.push_back(polynomial_ut0::BoundValueVelocity_Darcy);//mixed
      boundary_data.push_back( BoundaryValueHomogenous ); // mixed
      boundary_data.push_back( polynomial_ut0::BoundValue_Darcy ); // primal
      
      /** coefficients */
      problem_coefficients = polynomial_ut0::LinCoeffs_NSE;
      mixed_darcy_coeffs   = polynomial_ut0::LinCoeffs_DarcyMixed;
      primal_darcy_coeffs  = polynomial_ut0::LinCoeffs_Darcy;
      
      polynomial_ut0::ExampleFile();
      break;
    case 2:
      /** exact_solution */
      exact_solution.push_back( riverbed::ExactU1_NSE );
      exact_solution.push_back( riverbed::ExactU2_NSE );
      exact_solution.push_back( riverbed::ExactP_NSE );
      exact_solution.push_back( riverbed::ExactU1_Darcy ); // only mixed
      exact_solution.push_back( riverbed::ExactU2_Darcy ); // only mixed
      exact_solution.push_back( riverbed::ExactP_Darcy );
  
      /** boundary condition */
      boundary_conditions.push_back( riverbed::BoundCondition_NSE );
      boundary_conditions.push_back( riverbed::BoundCondition_NSE );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      boundary_conditions.push_back( riverbed::BC_Velocity_Darcy);//mixed
      boundary_conditions.push_back( BoundConditionNoBoundCondition ); // mixed
      boundary_conditions.push_back( riverbed::BC_Darcy ); // primal
  
      /** boundary values */
      boundary_data.push_back( riverbed::U1BoundValue_NSE );
      boundary_data.push_back( riverbed::U2BoundValue_NSE );
      boundary_data.push_back( BoundaryValueHomogenous );
      boundary_data.push_back(riverbed::BoundValueVelocity_Darcy);//mixed
      boundary_data.push_back( BoundaryValueHomogenous ); // mixed
      boundary_data.push_back( riverbed::BoundValue_Darcy ); // primal
      
      /** coefficients */
      problem_coefficients = riverbed::LinCoeffs_NSE;
      mixed_darcy_coeffs   = riverbed::LinCoeffs_DarcyMixed;
      primal_darcy_coeffs  = riverbed::LinCoeffs_Darcy;
      
      riverbed::ExampleFile();
      break;
    case 3:
      /** exact_solution */
      exact_solution.push_back( ana_sol_Stokes_Darcy3::ExactU1_NSE );
      exact_solution.push_back( ana_sol_Stokes_Darcy3::ExactU2_NSE );
      exact_solution.push_back( ana_sol_Stokes_Darcy3::ExactP_NSE );
      exact_solution.push_back( ana_sol_Stokes_Darcy3::ExactU1_Darcy ); // mixed
      exact_solution.push_back( ana_sol_Stokes_Darcy3::ExactU2_Darcy ); // mixed
      exact_solution.push_back( ana_sol_Stokes_Darcy3::ExactP_Darcy );
  
      /** boundary condition */
      boundary_conditions.push_back( ana_sol_Stokes_Darcy3::BoundCondition_NSE);
      boundary_conditions.push_back( ana_sol_Stokes_Darcy3::BoundCondition_NSE);
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      boundary_conditions.push_back( ana_sol_Stokes_Darcy3::BC_Velocity_Darcy);//mixed
      boundary_conditions.push_back( BoundConditionNoBoundCondition ); // mixed
      boundary_conditions.push_back( ana_sol_Stokes_Darcy3::BC_Darcy );// primal
  
      /** boundary values */
      boundary_data.push_back( ana_sol_Stokes_Darcy3::U1BoundValue_NSE );
      boundary_data.push_back( ana_sol_Stokes_Darcy3::U2BoundValue_NSE );
      boundary_data.push_back( BoundaryValueHomogenous );
      //mixed
      boundary_data.push_back(ana_sol_Stokes_Darcy3::BoundValueVelocity_Darcy);
      boundary_data.push_back( BoundaryValueHomogenous ); // mixed
      // primal
      boundary_data.push_back( ana_sol_Stokes_Darcy3::BoundValue_Darcy ); 
      
      /** coefficients */
      problem_coefficients = ana_sol_Stokes_Darcy3::LinCoeffs_NSE;
      mixed_darcy_coeffs   = ana_sol_Stokes_Darcy3::LinCoeffs_DarcyMixed;
      primal_darcy_coeffs  = ana_sol_Stokes_Darcy3::LinCoeffs_Darcy;
      
      ana_sol_Stokes_Darcy3::ExampleFile();
      break;
    case 4:
      /** exact_solution */
      exact_solution.push_back( porousObstacle::ExactU1_NSE );
      exact_solution.push_back( porousObstacle::ExactU2_NSE );
      exact_solution.push_back( porousObstacle::ExactP_NSE );
      exact_solution.push_back( porousObstacle::ExactU1_Darcy ); // mixed
      exact_solution.push_back( porousObstacle::ExactU2_Darcy ); // mixed
      exact_solution.push_back( porousObstacle::ExactP_Darcy );
  
      /** boundary condition */
      boundary_conditions.push_back( porousObstacle::BoundCondition_NSE);
      boundary_conditions.push_back( porousObstacle::BoundCondition_NSE);
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      boundary_conditions.push_back( porousObstacle::BC_Velocity_Darcy);//mixed
      boundary_conditions.push_back( BoundConditionNoBoundCondition ); // mixed
      boundary_conditions.push_back( porousObstacle::BC_Darcy );// primal
  
      /** boundary values */
      boundary_data.push_back( porousObstacle::U1BoundValue_NSE );
      boundary_data.push_back( porousObstacle::U2BoundValue_NSE );
      boundary_data.push_back( BoundaryValueHomogenous );
      //mixed
      boundary_data.push_back(porousObstacle::BoundValueVelocity_Darcy);
      boundary_data.push_back( BoundaryValueHomogenous ); // mixed
      // primal
      boundary_data.push_back( porousObstacle::BoundValue_Darcy ); 
      
      /** coefficients */
      problem_coefficients = porousObstacle::LinCoeffs_NSE;
      mixed_darcy_coeffs   = porousObstacle::LinCoeffs_DarcyMixed;
      primal_darcy_coeffs  = porousObstacle::LinCoeffs_Darcy;
      
      porousObstacle::ExampleFile();
      break;
    default:
      ErrMsg("Unknown name of the Stokes--Darcy' example!");
      exit(0);
  }
}

/** ************************************************************************ */
std::shared_ptr<Example_NSE2D> Example_StokesDarcy2D::get_stokes_example() const
{
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
    OutPut("extracting Stokes example from Stokes--Darcy example\n");
  
  std::vector <DoubleFunct2D*> exact_stokes;
  std::vector <BoundCondFunct2D*> bc_stokes;
  std::vector <BoundValueFunct2D*> bd_stokes;
  
  exact_stokes.push_back(exact_solution[0]); // velocity first component
  exact_stokes.push_back(exact_solution[1]); // velocity second component
  exact_stokes.push_back(exact_solution[2]); // pressure
  
  bc_stokes.push_back(boundary_conditions[0]); // velocity first component
  bc_stokes.push_back(boundary_conditions[1]); // velocity second component
  bc_stokes.push_back(boundary_conditions[2]); // pressure
  
  bd_stokes.push_back(boundary_data[0]); // velocity first component
  bd_stokes.push_back(boundary_data[1]); // velocity second component
  bd_stokes.push_back(boundary_data[2]); // pressure
  
  return std::make_shared<Example_NSE2D>(exact_stokes, bc_stokes, bd_stokes,
                                         problem_coefficients);
}

/** ************************************************************************ */
std::shared_ptr<Example_Darcy2D> 
 Example_StokesDarcy2D::get_mixed_darcy_example() const
{
  if(TDatabase::ParamDB->SC_VERBOSE > 1) 
    OutPut("extracting mixed Darcy example from Stokes--Darcy example\n");
  
  std::vector <DoubleFunct2D*> exact_darcy;
  std::vector <BoundCondFunct2D*> bc_darcy;
  std::vector <BoundValueFunct2D*> bd_darcy;
  
  exact_darcy.push_back(exact_solution[3]); // velocity first component
  exact_darcy.push_back(exact_solution[4]); // velocity second component
  exact_darcy.push_back(exact_solution[5]); // pressure
  
  bc_darcy.push_back(boundary_conditions[3]); // boundary condition for velocity space
  bc_darcy.push_back(boundary_conditions[4]); // bounadry condition for pressure space (dummy)
  
  bd_darcy.push_back(boundary_data[3]); // flux (Dirichlet) or pressure (Neumann)
  bd_darcy.push_back(boundary_data[4]); // boundary data for pressure (dummy)
  
  return std::make_shared<Example_Darcy2D>(exact_darcy, bc_darcy, bd_darcy,
                                           mixed_darcy_coeffs);  
}

/** ************************************************************************ */
std::shared_ptr<Example_CD2D> Example_StokesDarcy2D::get_primal_darcy_example() 
 const
{
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
    OutPut("extracting primal Darcy example from Stokes--Darcy example\n");
  
  std::vector <DoubleFunct2D*> exact_darcy;
  std::vector <BoundCondFunct2D*> bc_darcy;
  std::vector <BoundValueFunct2D*> bd_darcy;
  
  exact_darcy.push_back(exact_solution[5]);
  bc_darcy.push_back(boundary_conditions[5]);
  bd_darcy.push_back(boundary_data[5]);
  
  return std::make_shared<Example_CD2D>(exact_darcy, bc_darcy, bd_darcy,
                                        primal_darcy_coeffs);  
}

      
#include <Example_CD3D.h>

#include <Database.h>
#include <FEFunction3D.h>
#include <FEDatabase3D.h>
#include <SquareMatrix3D.h>
#include <string.h>

/* examples */

namespace sine_laplace_3D
{
  #include "CD_3D/Laplace.h"
}

namespace hemker_3d
{
#include "CD_3D/Hemker_3D.h"
}

namespace test_p0_zero
{
  #include <test_p0.h>
}

namespace test_p1
{
  #include "CD_3D/test_p1.h"
}

namespace test_p2
{
  #include "CD_3D/test_p2.h"
}

//=========================================================================

Example_CD3D::Example_CD3D(const ParameterDatabase& user_input_parameter_db)
 : Example3D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch( example_code )
  {
    //steady-state problems
    case 0:
    {
      using namespace sine_laplace_3D;
      /** exact_solution */
      exact_solution.push_back( Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( BoundValue );
      
      /** coefficients */
      problem_coefficients = BilinearCoeffs;
      
      sine_laplace_3D::PECLET_NUMBER = this->get_nu();
      ExampleFile();
    }
      break;
      
    case 1:
    {
      using namespace hemker_3d;
      /** exact_solution */
      exact_solution.push_back( hemker_3d::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( hemker_3d::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( hemker_3d::BoundValue );
      
      /** coefficients */
      problem_coefficients = hemker_3d::BilinearCoeffs;
      
      hemker_3d::PECLET_NUMBER = this->get_nu();
      
      ExampleFile();
      break;
    }

    //negative integers are reserved for pure test examples
    case -1:
    {//constant zero solution example
      using namespace test_p0_zero;
      exact_solution.push_back( Exact );
      boundary_conditions.push_back( BoundCondition );
      boundary_data.push_back( BoundValue );
      problem_coefficients = BilinearCoeffs;
      ExampleFile();
      break;
    }
    case -2:
    {//linear solution example
      using namespace test_p1;
      exact_solution.push_back( Exact );
      boundary_conditions.push_back( BoundCondition );
      boundary_data.push_back( BoundValue );
      problem_coefficients = BilinearCoeffs;
      ExampleFile();
      break;
    }
    case -3:
    {//quadratic solution example
      using namespace test_p2;
      exact_solution.push_back( Exact );
      boundary_conditions.push_back( BoundCondition );
      boundary_data.push_back( BoundValue );
      problem_coefficients = BilinearCoeffs;
      ExampleFile();
      break;
    }
    default:
      ErrThrow("Unknown name of the convection-diffusion-reaction example CD3D or Time_CD3D!");
  }
}

void Example_CD3D::do_post_processing(ConvectionDiffusion<3>& cd3d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(cd3d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_CD3D","No post processing done for the current example.");
  }
}

double Example_CD3D::get_nu() const
{
  double diffusion_coefficient = this->example_database["diffusion_coefficient"];
  return diffusion_coefficient;
}


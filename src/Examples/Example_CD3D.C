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

Example_CD3D::Example_CD3D(int example_code,
                           const ParameterDatabase& user_input_parameter_db) : Example3D()
{
  this->example_database.merge(user_input_parameter_db,false);

  switch( example_code )
  {
    //steady-state problems
    case 0:
      /** exact_solution */
      exact_solution.push_back( sine_laplace_3D::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( sine_laplace_3D::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sine_laplace_3D::BoundValue );
      
      /** coefficients */
      problem_coefficients = sine_laplace_3D::BilinearCoeffs;
      
      sine_laplace_3D::ExampleFile();
      break;

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

void Example_CD3D::do_post_processing(CD3D& cd3d) const
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
  double inverse_reynolds = this->example_database["reynolds_number"];
  inverse_reynolds = 1/inverse_reynolds;
  return inverse_reynolds;
}


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

Example_CD3D::Example_CD3D() : Example3D()
{
  switch( TDatabase::ParamDB->EXAMPLE ) 
  {
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
    default:
      ErrMsg("Unknown name of the convection-diffusion (CD3D) example!");
      throw("Unknown name of the convection-diffusion (CD3D) example!");
  }
}

      

#include <Example_TimeCD3D.h>
#include <Database.h>
#include <string.h>

namespace linear_space_time  //-2
{
  #include "TCD_3D/linear_space_time.h"
}

namespace quad_space_time  //-1
{
#include "TCD_3D/quadratic_space_time.h"
}

namespace concentration  //0
{
#include "TCD_3D/concentrationOfSpecies_3d.h"
}


/********* BELOW THIS LINE, USER-SPECIFIC EXAMPLES *******/
namespace quadratic   // 10
{
#include "Time_CD3D/10_Quadratic_TCD3D.h"
}



Example_TimeCD3D::Example_TimeCD3D(
  const ParameterDatabase& user_input_parameter_db)
 : Example_NonStationary3D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case -2:
    {
      // linear space and time solution example 
      using namespace linear_space_time;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initialCondtion.push_back(InitialCondition);
      ExampleFile();
    }
      break;
    case -1:
    {
      using namespace quad_space_time;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initialCondtion.push_back(InitialCondition);
      ExampleFile();
    }
      break;
    case 0:
      using namespace concentration;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initialCondtion.push_back(InitialCondition);
      
      concentration::PECLET_NUMBER = this->get_nu();
      ExampleFile();
      break;

    /********* BELOW THIS LINE, USER-SPECIFIC EXAMPLES *******/
    case 10:        // quadratic space time, like example -1
    {
      using namespace quadratic;
      exact_solution.push_back(quadratic::Exact);
      boundary_conditions.push_back(quadratic::BoundCondition);
      boundary_data.push_back(quadratic::BoundValue);
      problem_coefficients = quadratic::BilinearCoeffs;
      initialCondtion.push_back(quadratic::InitialCondition);
      quadratic::ExampleFile();
      break;
    }
    default:
      ErrThrow("Unknown Example_TimeCD3D example!");
  }
}

void Example_TimeCD3D::do_post_processing(Time_CD3D& tcd3d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(tcd3d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_TimeCD3D","No post processing done for the current example.");
  }
}

double Example_TimeCD3D::get_nu() const
{
  double diffusion_coefficient= this->example_database["diffusion_coefficient"];
  return diffusion_coefficient;
}

#include <Example_TimeCD3D.h>
#include <Database.h>
#include <string.h>

namespace linear_space_time
{
  #include "TCD_3D/linear_space_time.h"
}

namespace quad_space_time
{
#include "TCD_3D/quadratic_space_time.h"
}

namespace concentration
{
#include "TCD_3D/concentrationOfSpecies_3d.h"
}


////CB PROJECT include the wiedmeyer fluidized bed crystallizer example here
#include <WiedmeyerBatchCrystallizer.h>
////END PROJECT

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

    //CB PROJECT
    case 71:
    {// c: mass balance in batch crystallizer
      using namespace wiedmeyer_batch_crystallizer;
      ConcentrationProperties::initialise_inlet_values(
          user_input_parameter_db["initial_concentration"],
          user_input_parameter_db["fluid_cycle_time"],
          TDatabase::TimeDB->TIMESTEPLENGTH);
      exact_solution.push_back(Exact_cALUM);
      boundary_conditions.push_back(BoundCondition_cALUM);
      boundary_data.push_back(BoundValue_cALUM);
      problem_coefficients = BilinearCoeffs_cALUM;
      initialCondtion.push_back(InitialCondition_cALUM);
      break;
    }
    case 72:
    {// T: energy balance in batch crystallizer
      using namespace wiedmeyer_batch_crystallizer;
      TemperatureConditions::set_T_start(user_input_parameter_db["T_start"]);
      TemperatureConditions::set_T_end(user_input_parameter_db["T_end"]);
      TemperatureConditions::set_t_start( TDatabase::TimeDB->STARTTIME );
      TemperatureConditions::set_t_end(TDatabase::TimeDB->ENDTIME);

      exact_solution.push_back(Exact_T);
      boundary_conditions.push_back(BoundCondition_T);
      boundary_data.push_back(BoundValue_T);
      problem_coefficients = BilinearCoeffs_T;
      initialCondtion.push_back(InitialCondition_T);
      break;
    }
    //END PROJECT

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

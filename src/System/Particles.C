#include "Particles.h"
#include "Database.h"

void Particles::update_particle(const std::vector< double >& fluid_velocity, 
                                double fluid_density)
{
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double reynold_pr;
  double hdynamic_force_coeff;
}

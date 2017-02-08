#include "Particle_Tracer.h"

Particle_Tracer::Particle_Tracer(std::shared_ptr< std::vector< Particles > >& p, 
 std::shared_ptr< std::vector< TFEFunction2D > >& funtct, const TDomain& dom)
: particles(p), u(funtct), domain(dom)
{
}

void Particle_Tracer::check_boundary_collision()
{
  
}

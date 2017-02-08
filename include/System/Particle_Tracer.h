/** ***************************************************************************
 * 
 * @name Particle_Tracer
 * @brief 
 * 
 * @author Naveed Ahmed, Felix Anker
 * @date  2017/02/08
** ***************************************************************************/


#ifndef PARTICLE_TRACER_H
#define PARTICLE_TRACER_H

#include <Particles.h>
#include <FEFunction2D.h>
#include <Domain.h>
#include <memory>

class Particle_Tracer
{
protected:
  /// 
  std::shared_ptr<std::vector<Particles>> particles;
  /// 
  std::shared_ptr<std::vector<TFEFunction2D>> u;  
  ///
  const TDomain domain;
public:
  // constructor
  Particle_Tracer(std::shared_ptr<std::vector<Particles>> &p, 
       std::shared_ptr<std::vector<TFEFunction2D>> &funtct, const TDomain &dom);
  
  /// 
  void check_boundary_collision();
};

#endif // PARTICLE_TRACER_H

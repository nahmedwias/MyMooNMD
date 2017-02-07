/** ***************************************************************************
 * 
 * @name Particles
 * @brief 
 * 
 * @author Naveed Ahmed, Felix Anker
 * @date  2017/02/06
/** ***************************************************************************/

#ifndef PARTICLES_H
#define PARTICLES_H

#include <Point.h>

class Particles
{
protected:
  /** @brief positions of particles at each time point*/
  std::vector<Point> positions;
  /** @brief diameters of particles at each time point*/
  std::vector<double> diameters;
  /** @brief density of the particle*/
  double rho;
  /** @brief velocity of the particles*/
  std::vector<double> velocity;
  /** @brief average particle slip velocity */
  double slip_velo;
  
public:
  /// default constructor
  Particles() = default;
  
  void update_particle(const std::vector<double> &fluid_velocity, 
                       double fluid_density);
  
  std::vector<double> get_velocity() const
  {
    return velocity;
  }
};

#endif // PARTICLES_H
#ifndef MEANVELOCITY_H
#define MEANVELOCITY_H

#include <vector>
#include<stddef.h>

#ifdef __2D__

#include <TimeNavierStokes.h>

#endif

namespace MeanVelocity
{
  static std::vector<double> xlayers, ylayers;
  static size_t n_xlayers, n_ylayers;
  static std::vector<double>xDofs, yDofs;
  static size_t nDofs;
  
  static std::vector< double > temporal_mean_u1; 
  static std::vector< double > temporal_mean_u2; 
  static std::vector< double > R11; 
  static std::vector< double > R22;
#ifdef __2D__  
  void fill_arrays(const TimeNavierStokes<2>& tnse2d );
  void compute_mean_velocity(const TimeNavierStokes<2>& tnse2d);
  
  void compute_mean_velocity_on_points(const TimeNavierStokes<2>& tnse2d,
                          const std::vector<double>& vec, const std::vector<TBaseCell *> cells);
  
  void  compute_velocity_on_points(const TimeNavierStokes<2>& tnse2d,
                          const std::vector<double>& vec_x_y, const std::vector<TBaseCell *> cells,
                          std::string basename);
#endif
};

#endif // MEANVELOCITY_H

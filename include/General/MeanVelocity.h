#ifndef MEANVELOCITY_H
#define MEANVELOCITY_H

#include <vector>

#ifdef __2D__

#include <Time_NSE2D_Merged.h>

class Time_NSE2D_Merged;
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
  void fill_arrays(const Time_NSE2D_Merged& tnse2d );
  void  compute_mean_velocity(const Time_NSE2D_Merged& tnse2d);
  void  compute_mean_velocity_on_points(const Time_NSE2D_Merged& tnse2d, const std::vector<double>& vec,
                                        const std::vector<TBaseCell *> cells);
  void  compute_velocity_on_points(const Time_NSE2D_Merged& tnse2d, const std::vector<double>& vec,
                                        const std::vector<TBaseCell *> cells);
#endif
};

#endif // MEANVELOCITY_H

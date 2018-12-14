#ifndef _PrePost_Cylinder_Square_
#define _PrePost_Cylinder_Square_

#include <ParameterDatabase.h>
#include <Collection.h>
#include <TimeNavierStokes.h>

class TCollection;
namespace Cylinder_Square
{
  /// needed data
  static double dimensionless_viscosity;
  /// set parameters used in the corresponding example
  void setParameters(ParameterDatabase& db_);
  
  /// computation of drag, lift coefficient
  void get_Drag_Lift(TFEFunction3D *u1, TFEFunction3D *u2,
             TFEFunction3D *u3, TFEFunction3D *p,
             TFEFunction3D *u1old, TFEFunction3D *u2old,
             double &cdrag, double &clift, const double dt);
  // compute pressure difference
  double get_p_diff(const std::array< double, int(3) >& point_A, const std::array< double, int(3) >& point_B, const TFEFunction3D& p);

  void compute_drag_lift_pdiff(TimeNavierStokes<3>& tnse3d);
  
  /// variable used to compute the ..???
  static size_t n_center_velo;
  static std::vector<double> center_velo;
  // this variable is set in the preparation and used in
  // computing the centerline velocities
  static int counter_av;
  void PrepareCenterlineVelocities(TCollection* coll);
  // compute velocities at centerline
  void CenterlineVelocities(TimeNavierStokes<3>& tnse3d);
  
  static size_t n_cyl_velo;
  static std::vector<double> cyl_velo;
  static int counter_av_single;
  void PrepareVelocityAtCylinder(TCollection* coll);
  // compute velocity at cyliner
  void VelocityAtCylinder(TimeNavierStokes<3>& tnse3d);
  
  static size_t n_pres_nodes;
  static std::vector<double> press_cyl;
  void PreparePressureAtCylinder(const TCollection* coll);
  // compute pressure at cylinder 
  static int counter_av_pres;
  void PressureAtCylinder(TimeNavierStokes<3>& tnse3d);
  
  /// set no pentration values 
  void SetNoPenetrationValues(TimeNavierStokes<3>& tnse3d);
  
  /// compute fraction velocities
  static int count_fric_vel;
  static std::vector<double> velo_friction;
  void ComputeFrictionVelocities(const TimeNavierStokes<3>& tnse3d);
}

#endif // _PrePost_Cylinder_Square_

/**
 * @file This file holds routines used to set the coordinates, periodic joints etc.
 * 
 * And is used within the project "channel_rey180" and can/will be extended
 * to the wind_tunnel model
 * 
 * followed from the twisted pipe project
 */

#ifndef _ChannelFlowRoutines_
#define _ChannelFlowRoutines_

#include <vector>
#include <ParameterDatabase.h>
#include <Time_NSE3D.h>

class TCollection;

namespace ChannelFlowRoutines
{
  /// set parameters used for the corresponding example
  void setParameters(ParameterDatabase& db);
  /// set the z coordinates for a channel flow problem
  /// with reynolds no 180
  /// periodic b.c. : x = -1 and x = 1 and z = 0 and z = 2
  void setZCoordinates(TCollection* Coll, int level);

  /// check the coordinates
  void checkZCoordinates(TCollection* Coll, int level);

  /// set the refinement descriptor
  void setRefineDesc(TCollection* coll);

  /// set the periodic joints
  void setPeriodicFaceJoints(TCollection* Coll);

  /// degrees of freedoms for x, y and z coordinates
  static std::vector<double> xDofs;
  static std::vector<double> yDofs;
  static std::vector<double> zDofs;
  /// array with coordinates of layers in z-direction
  static std::vector<double> zLayers;
  /// number of dof-layers in z-direction
  static size_t nZLayers;
  // number of basis functions
  static size_t nBasisFunction;

  /// routine for computing the coordinates of the d.o.f. for a finite
  void GetCoordinatesOfDof(const Time_NSE3D& tnse3d);

  /// routine for computing mean velocity
  /// currently implemented for only Channel Flow with
  /// Reynolds number=180 and Turbulent Model used is
  /// the Smagorinsky
  void computeMeanVelocity(Time_NSE3D& tnse3d);
  /// compute the average velocity
  /// @Output: list of velocity gradients for all three components
  /// of velocity with respect to x, y and z

  void computeAverageVelocity(
          std::array< std::vector< double >, int(9) > velo, 
          int ndofs, Time_NSE3D& tnse3d);
  /// no of summations in the computations of mean velocity
  static std::vector<int> sum_layer_dofs;
  /// 
  void summation(size_t length);
  /// computations of spatial mean velocity at current time
 std::vector<double> spatialMean(std::vector< double > in);
 /// mean velocity profile (arithmatic mean)
 std::vector<double> meanVelocity(std::vector< double > vecin);
 /// compute the product of two vectors
 std::vector<double> product(std::vector<double>vec1, 
                             std::vector<double>vec2);
 /// A_ij: contribution from the eddy-viscosity model
 /// maximum six vectors to be returned
 /// eddy_xx, eddy_yy, eddy_zz, eddy_xy, eddy_xz, eddy_yz
 ///TODO: not considered yet 
 /// @input: 
 void eddy_viscosity(std::array< std::vector< double>, int(6)> eddy, 
                     Time_NSE3D& tnse3d);
 /// u_tau = frictionvelocity
 double getFrictionVelocity(std::vector<double> vec);
 /// root mean square velocities
 std::deque<std::vector< double>> getrms(std::deque<std::vector<double>> reynoldStress, 
                            std::deque<std::vector<double>> meanvelo);
 /// out put for plotting
 void saveData(std::deque<std::vector<double>> m, std::deque<std::vector<double>> rms, 
             std::deque<std::vector<double>> R);
 
}

#endif

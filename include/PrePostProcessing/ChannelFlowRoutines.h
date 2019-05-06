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
#include <deque>
#include <array>
#include <FESpace3D.h>
#include <TimeNavierStokes.h>

class TCollection;
template <int d> class TimeNavierStokes; //forward declaration
namespace ChannelTau180
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
  
  // inverse of the Reynolds number
  static double reynolds_number;

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
  
  /// number of layers
  static int n_layers;

  ///@brief computing the coordinates of the d.o.f. 
  /// for a finite element: currently only Q2-element
  /// @param[in]
  void GetCoordinatesOfDof(const TimeNavierStokes<3>& tnse3d);
  
  /// @brief mean velocities computed at each time step
  static std::deque<std::vector<double>> MeanVelocity;
  /// @brief Reynold Stress
  static std::deque<std::vector<double>> ReynoldsStress;
  /// @brief derivative of the mean velocity
  /// only the first component of velocity
  static std::vector<double> DerivmeanVelo;
  /// @brief set the memory 
  void set_up_memory();
  /// routine for computing mean velocity
  /// currently implemented for only Channel Flow with
  /// Reynolds number=180 and Turbulent Model used is
  /// the Smagorinsky
  void computeMeanVelocity(const TimeNavierStokes<3>& tnse3d); // replace the names when completed
  /// @brief this computes the temporal mean 
  void temporalMean(std::vector< double > spatial_mean, std::vector< double >& temporal_mean);
  ///TODO: not used and therefore not completed yet 
  /// need to modify also when started implementation
  /// this function is used in the computations of "eddy_viscosity"
  /// update accordingly:
  /// @brief compute the average velocity
  /// @param[out] list of velocity gradients for all three components
  /// of velocity with respect to x, y and z
  /// @param[in] TimeNavierStokes<3> object
  void computeAverageVelocity( std::array< std::vector< double >, int(9) > velo, 
                               const TimeNavierStokes<3>& tnse3d);
  /// compute the summation of layers
  void count_dofs_per_layer(std::vector< int >& summ, std::shared_ptr<const TFESpace3D> fesp);
  /// computations of spatial mean velocity at current time
 std::vector<double> compute_sum_of_velocity(std::vector< double > in, std::shared_ptr<const TFESpace3D> space);
 /// A_ij: contribution from the eddy-viscosity model
 /// maximum six vectors to be returned
 /// eddy_xx, eddy_yy, eddy_zz, eddy_xy, eddy_xz, eddy_yz
 ///TODO: not considered yet 
 /// @param[out] eddy eddy viscosity
 /// @param[in]  TimeNavierStokes<3> class object 
 void eddy_viscosity(std::array< std::vector< double >, int(6) > eddy, 
                     const TimeNavierStokes<3>& tnse3d);
 /// @brief computes the friction viscosity u_tau
 /// @param[in] vec mean velocity vector for the first component only
 /// @param[in out], derivative of the mean velcoity
 double getFrictionVelocity(std::vector< double > vec, std::vector< double >& meanDeriv);
 ///@brief compute the root mean square velocities
 ///@param[in] reynoldStress all components of Reynold Stress 
 /// (spatial,Temporal  averaged)
 ///@param[in] meanVelo all components: mean velcoity (spatial, Temporal averaged)
 /// @return root mean square velcoity all components 
 std::deque<std::vector< double>> getrms(std::deque<std::vector<double>> reynoldStress, 
                            std::deque<std::vector<double>> meanvelo);
 ///@brief print out the interested quantities
 void print_quantity_of_interest(std::deque<std::vector<double>> m, std::deque<std::vector<double>> rms, 
             std::deque<std::vector<double>> R);
 
}

#endif

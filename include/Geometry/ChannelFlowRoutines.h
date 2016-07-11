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

class TCollection;

namespace ChannelFlowRoutines
{
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
};

#endif
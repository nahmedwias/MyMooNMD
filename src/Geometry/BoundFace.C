// =======================================================================
// @(#)BoundFace.C        1.2 10/18/99
// 
// Class:       TBoundFace
// Purpose:     face on a boundary component
//
// Author:      Volker Behns  28.08.97
//
// =======================================================================

#ifndef __3D__
#define __3D__
#endif

#include <BoundFace.h>

// Constructors
TBoundFace::TBoundFace(TBoundComp3D *bdcomp, double *param1, double *param2)
 : TJoint()
{
  int i;
  ID = BoundaryFace;

  BoundComp = bdcomp;
  for(i=0;i<4;i++)
  {
    Param1[i] = param1[i];
    Param2[i] = param2[i];
  }
}

TBoundFace::TBoundFace(TBoundComp3D *bdcomp) : TJoint()
{
  ID = BoundaryFace;

  BoundComp = bdcomp;
  Param1[0] = 0.0;
  Param2[0] = 0.0;
  Param1[1] = 1.0;
  Param2[1] = 0.0;
  Param1[2] = 1.0;
  Param2[2] = 1.0;
  Param1[3] = 0.0;
  Param2[3] = 1.0;
}

// Methods
int TBoundFace::CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp)
{
  Tmp.Filled = false;
  return 0;
}

// create a new instance of this class
TJoint *TBoundFace::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TBoundFace(BoundComp);
}

TJoint *TBoundFace::NewInst()
{
  return new TBoundFace(BoundComp);
}

void TBoundFace::SetParameters(double *param1, double *param2)
{
  int i;

  for(i=0;i<4;i++)
  {
    Param1[i] = param1[i];
    Param2[i] = param2[i];
  }
}

void TBoundFace::GetParameters(double *param1, double *param2)
{
  int i;

  for(i=0;i<4;i++)
  {
    param1[i] = Param1[i];
    param2[i] = Param2[i];
  }
}

//NEW LB 11.10.18
/**
 @brief set the local index of this joint in given neighbor
 @todo write a more general function that does not take index as input
 */
void TBoundFace::set_index_in_neighbour(TBaseCell *neigh, int index)
{
    if(neigh == Neighb0)
        IndexInNeighbor[0] = index;
    else
    {
        ErrMsg("ERROR in TBoundFace::set_index_in_neighbour");
        exit(-4711);
    }
}


/** get the index of this joint in given neighbor */
int TBoundFace::get_index_in_neighbour(const TBaseCell*const neigh) const
{
    if(neigh == Neighb0)
        return IndexInNeighbor[0];
    else
    {
        ErrMsg("ERROR in TBoundFace::get_index_in_neighbour");
        exit(-4711);
    }
}

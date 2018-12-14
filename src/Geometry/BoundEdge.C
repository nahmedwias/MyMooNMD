// =======================================================================
// @(#)BoundEdge.C        1.1 10/30/98
// 
// Class:       TBoundEdge
// Purpose:     edge on a boundary component
//
// Author:      Volker Behns  02.08.97
//
// =======================================================================

#include <BoundEdge.h>

// Constructors
TBoundEdge::TBoundEdge(TBoundComp2D *bdcomp, double t_0, double t_1)
{
  ID = BoundaryEdge;

  BoundComp = bdcomp;
  T_0 = t_0;
  T_1 = t_1;
}

// Methods
int TBoundEdge::CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp)
{
  Tmp.Filled = false;
  return 0;
}

#ifdef __2D__
/** update parameters according to the new vertex positions */
void TBoundEdge::UpdateParameters(const TVertex *Begin, const TVertex *End)
{
  double x1, y1, x2, y2;
  double t1, t2;

#ifdef __2D__
  Begin->GetCoords(x1, y1);
  End->GetCoords(x2, y2);
#else
  double z1, z2;
  Begin->GetCoords(x1, y1, z1);
  End->GetCoords(x2, y2, z2);
#endif

  BoundComp->GetTofXY(x1, y1, t1);
  BoundComp->GetTofXY(x2, y2, t2);

  T_0 = t1;
  T_1 = t2;
}
#endif

// create a new instance of this class
TJoint *TBoundEdge::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TBoundEdge(BoundComp, T_0 + newT_0*(T_1 - T_0),
                        T_0 + newT_1*(T_1 - T_0));
}

TJoint *TBoundEdge::NewInst()
{
  return new TBoundEdge(BoundComp, T_0, T_1);
}

/**
 @brief set the local index of this joint in given neighbor 
 @todo write a more general function that does not take index as input
 */
void TBoundEdge::set_index_in_neighbour(TBaseCell *neigh, int index)
{
    if(neigh == Neighb0)
        IndexInNeighbor[0] = index;
    else
    {
        ErrMsg("ERROR, TInnerInterfaceJoint::SetIndexInNeighbor !!!!!!!!");
        exit(-4711);
    }
}


/** get the index of this joint in given neighbor */
int TBoundEdge::get_index_in_neighbour(const TBaseCell*const neigh) const
{
    if(neigh == Neighb0)
        return IndexInNeighbor[0];
    else
    {
        ErrMsg("ERROR, TBoundEdge::GetIndexInNeighbor !!!!!!!!");
        exit(-4711);
    }
}


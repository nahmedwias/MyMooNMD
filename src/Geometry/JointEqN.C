// =======================================================================
// @(#)JointEqN.C        1.8 11/15/99
// 
// Class:       TJointEqN
// Purpose:     connects two cells
//
// Author:      Volker Behns  23.07.97
//
// =======================================================================

#include <BaseCell.h>
#include <Constants.h>
#include <Database.h>
#include <JointEqN.h>
#include <RefDesc.h>

// Constructors
TJointEqN::TJointEqN(TBaseCell *neighb0) : TJoint()
{
  ID = JointEqN;

  Neighb0 = neighb0;
}

TJointEqN::TJointEqN(TBaseCell *neighb0, TBaseCell *neighb1) : TJoint()
{
  ID = JointEqN;

  Neighb0 = neighb0;
  Neighb1 = neighb1;
}

// Methods

int TJointEqN::CheckMatchingRef(TBaseCell *Me, int J_i, struct StoreGeom &Tmp)
{
  TBaseCell *Neighb, *Parent = NULL;
  int i, J_j, N_, MaxLen1, MaxLen2, auxi, aux;
  const int *TmpValues1, *TmpValues2;
  TRefDesc *NeighbRefDesc, *ParRefDesc;
  TJoint *ParJoint;
#ifdef __2D__
  Refinements NeibEdgeRef, MyEdgeRef;
  const int *TmpoEnE, *TmpEC, *TmpECI, *TmpoEnV, *TmpVC, *TmpVCI;
#else
  Refinements NeibFaceRef, MyFaceRef;
  const int *TmpLen1, *TmpLen2, *TmpLen3, *TmpLen4;
  const int *TmpoFnF, *TmpFC, *TmpFCI, *TmpoFnV, *TmpVC, *TmpVCI;
  const int *MapVerts, *MapFaces;
  int MaxLen3;
#endif

  Tmp.Filled = FALSE;

  if (Neighb0 == Me)
    Neighb = Neighb1;
  else
    Neighb = Neighb0;

  if (!Neighb)
  {
    Parent = Me;
    while ((Parent = Parent->GetParent()))
    {
      N_ = Parent->GetN_Children();
      for (i=0;i<N_;i++)
        if (Parent->GetChild(i) == Me) break;

      ParRefDesc = Parent->GetRefDesc();
#ifdef __2D__
      ParRefDesc->GetChildEdge(TmpValues1, MaxLen1);
      ParRefDesc->GetNewEdgeOldEdge(TmpValues2);
#else
      ParRefDesc->GetChildFace(TmpValues1, MaxLen1);
      ParRefDesc->GetNewFaceOldFace(TmpValues2);
#endif

      J_i = TmpValues2[TmpValues1[MaxLen1*i + J_i]];

      // we can go deeper only if nothing happens (NoRef)
#ifdef __2D__
      if (Parent->GetEdgeRef(J_i) != NoRef) return 0;
#else
      if (Parent->GetFaceRef(J_i) != NoRef) return 0;
#endif

      if ((Neighb = Parent->GetJoint(J_i)->GetNeighbour(Parent)))
        break;
    }

    if (!Parent) return 0;
    ParJoint = Parent->GetJoint(J_i);
  }

  if (Neighb->ExistChildren())
  {

#ifdef __2D__
    N_ = Neighb->GetN_Edges();
#else
    N_ = Neighb->GetN_Faces();
#endif
    if (Parent)
    {
      for (J_j=0;J_j<N_;J_j++)
        if (Neighb->GetJoint(J_j) == ParJoint) break;
    }
    else
      for (J_j=0;J_j<N_;J_j++)
        if (Neighb->GetJoint(J_j) == this) break;

#ifdef __2D__
    NeibEdgeRef = Neighb->GetEdgeRef(J_j);
    MyEdgeRef = Me->GetEdgeRef(J_i);

  #ifdef __MORTAR__
    // fix problem after mortar refinement layer
    if (MyEdgeRef != NoRef)
    {
      while (NeibEdgeRef == NoRef && Neighb->ExistChildren())
      {
        Neighb->GetRefDesc()->GetOldEdgeNewEdge(TmpoEnE, TmpValues1, MaxLen1);
        Neighb->GetRefDesc()->GetEdgeChild(TmpEC, TmpValues1, MaxLen2);
        Neighb->GetRefDesc()->GetEdgeChildIndex(TmpECI, TmpValues1, MaxLen2);

        J_j = TmpoEnE[MaxLen1 * J_j];
        Neighb = Neighb->GetChild(TmpEC[MaxLen2 * J_j]);
        J_j = TmpECI[MaxLen2 * J_j];
        NeibEdgeRef = Neighb->GetEdgeRef(J_j);
      }

      if (NeibEdgeRef == NoRef) return 0;
    }
  #endif // __MORTAR__

  
      if (NeibEdgeRef != MyEdgeRef)
       cout << "test JointEqn" << endl;
  
    if (NeibEdgeRef != MyEdgeRef)
      return -2;
   
    Tmp.Filled = TRUE;

    NeighbRefDesc = Neighb->GetRefDesc();

    NeighbRefDesc->GetOldEdgeNewEdge(TmpoEnE, TmpValues1, MaxLen1);
    N_ = TmpValues1[J_j] + 1;

    NeighbRefDesc->GetEdgeChild(TmpEC, TmpValues1, MaxLen2);
    NeighbRefDesc->GetEdgeChildIndex(TmpECI, TmpValues1, MaxLen2);

    for (i=0;i<N_;i++)
    {
      auxi = TmpoEnE[J_j * MaxLen1 + i] * MaxLen2;
      Tmp.Joints[i] = Neighb->GetChild(TmpEC[auxi])->
                              GetJoint(TmpECI[auxi]);
    }

    NeighbRefDesc->GetOldEdgeNewVertex(TmpoEnV, TmpValues1, MaxLen1);
    NeighbRefDesc->GetVertexChild(TmpVC, TmpValues1, MaxLen2);
    NeighbRefDesc->GetVertexChildIndex(TmpVCI, TmpValues1, MaxLen2);

    auxi = J_j * MaxLen1 + 1;
    for (i=1;i<N_;i++)
    {
      aux = TmpoEnV[auxi++] * MaxLen2;
      Tmp.Vertices[i] = Neighb->GetChild(TmpVC[aux])->
                          GetVertex(TmpVCI[aux]);
    }
#else
    NeibFaceRef = Neighb->GetFaceRef(J_j);
    MyFaceRef = Me->GetFaceRef(J_i);

    if (NeibFaceRef != MyFaceRef)
      return -2;

    Tmp.Filled = TRUE;

    NeighbRefDesc = Neighb->GetRefDesc();

    NeighbRefDesc->GetOldFaceNewFace(TmpoFnF, TmpLen1, MaxLen1);
    N_ = TmpLen1[J_j];

    NeighbRefDesc->GetFaceChild(TmpFC, TmpLen2, MaxLen2);
    NeighbRefDesc->GetFaceChildIndex(TmpFCI, TmpLen2, MaxLen2);

    NeighbRefDesc->GetOldFaceNewVertex(TmpoFnV, TmpLen3, MaxLen3);

    Neighb->GetJoint(J_j)->GetMapperRef(MapVerts, MapFaces);

    auxi = J_j * MaxLen1;
    for (i=0;i<N_;i++)
    {
      aux = TmpoFnF[auxi++] * MaxLen2;
      Tmp.Joints[MapFaces[i]] = Neighb->
            GetChild(TmpFC[aux])->GetJoint(TmpFCI[aux]);
    }

    NeighbRefDesc->GetVertexChild(TmpVC, TmpValues1, MaxLen2);
    NeighbRefDesc->GetVertexChildIndex(TmpVCI, TmpValues1, MaxLen2);

    N_ = TmpLen3[J_j];
    auxi = J_j * MaxLen3;
    for (i=0;i<N_;i++)
    {
      aux = TmpoFnV[auxi++] * MaxLen2;
      Tmp.Vertices[MapVerts[i]] = Neighb->GetChild(TmpVC[aux])->
                          GetVertex(TmpVCI[aux]);
    }
#endif // __2D__
  }

  return 0;
}

#ifdef __MORTAR__

int TJointEqN::CheckMatchingRef(TBaseCell *Me, int J_i, StoreGeomMortar &Tmp)
{
  TBaseCell *Neighb;
  int i, J_j, N_, MaxLen1, MaxLen2, auxi, aux;
  const int *TmpValue;
  const int *TmpoEnE, *TmpEC,  *TmpECI, *TmpoEnV, *TmpVC, *TmpVCI;
  TRefDesc *NeighbRefDesc;
  Refinements NeibEdgeRef, MyEdgeRef;

  Tmp.Filled = FALSE;

  if (Neighb0 == Me)
    Neighb = Neighb1;
  else
    Neighb = Neighb0;

  if (!Neighb) return 0;

  if (Neighb->ExistChildren())
  {
    N_ = Neighb->GetN_Edges();
    for (J_j=0;J_j<N_;J_j++)
      if (Neighb->GetJoint(J_j) == this) break;

    NeibEdgeRef = Neighb->GetEdgeRef(J_j);
    MyEdgeRef = Me->GetEdgeRef(J_i);

    if (NeibEdgeRef != MyEdgeRef)
      return -2;

    Tmp.Filled = TRUE;

    NeighbRefDesc = Neighb->GetRefDesc();

    NeighbRefDesc->GetOldEdgeNewEdge(TmpoEnE, TmpValue, MaxLen1);
    N_ = TmpValue[J_j] + 1;

    NeighbRefDesc->GetEdgeChild(TmpEC, TmpValue, MaxLen2);
    NeighbRefDesc->GetEdgeChildIndex(TmpECI, TmpValue, MaxLen2);

    for (i=0;i<N_;i++)
    {
      auxi = TmpoEnE[J_j * MaxLen1 + i] * MaxLen2;
      Tmp.Joints[i] = Neighb->GetChild(TmpEC[auxi])->
                              GetJoint(TmpECI[auxi]);
    }

    NeighbRefDesc->GetOldEdgeNewVertex(TmpoEnV, TmpValue, MaxLen1);
    NeighbRefDesc->GetVertexChild(TmpVC, TmpValue, MaxLen2);
    NeighbRefDesc->GetVertexChildIndex(TmpVCI, TmpValue, MaxLen2);

    auxi = J_j * MaxLen1 + 1;
    for (i=1;i<N_;i++)
    {
      aux = TmpoEnV[auxi++];
      Tmp.Vertices[i] = Neighb->GetChild(TmpVC[aux])->
                          GetVertex(TmpVCI[aux]);
    }
  }

  return 0;
}
#endif


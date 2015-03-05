// =======================================================================
// @(#)Joint.C        1.7 11/18/99
// 
// Class:       TJoint
// Purpose:     superclass for all joints
//
// Author:      Volker Behns  23.09.98
//
// =======================================================================

#include <BaseCell.h>
#include <Database.h>
#include <Joint.h>
#include <stdlib.h>

// Constructors
TJoint::TJoint()
{
  ID = Joint;

  Neighb0 = NULL;
  Neighb1 = NULL;

  NeibSubDomainLocalJointNo = -1;
#ifdef __3D__
  MapType = -1;
#endif
}

// Methods
int TJoint::SetNeighbour(TBaseCell *Neighb)
{
  switch (this->GetType())
  {
    case JointEqN:
    case MortarBaseJoint:
    case InterfaceJoint:
    case IsoInterfaceJoint:
#ifdef __3D__
    case InterfaceJoint3D:
    case IsoInterfaceJoint3D:
#endif
    case PeriodicJoint:
      if (Neighb0)
        Neighb1 = Neighb;
      else
        Neighb0 = Neighb;

      return 0;

    default:
      return -1;
  }
}

TBaseCell *TJoint::GetNeighbour(TBaseCell *Me)
{
  if (Neighb0 == Me)
    return Neighb1;
  else
    return Neighb0;
}

int TJoint::SetNeighbour(int i, TBaseCell *Neighb)
{
  switch (this->GetType())
  {
    case JointEqN:
    case MortarBaseJoint:
    case InterfaceJoint:
    case IsoInterfaceJoint:
#ifdef __3D__
    case InterfaceJoint3D:
    case IsoInterfaceJoint3D:
#endif
    case PeriodicJoint:
      if (i)
        Neighb1 = Neighb;
      else
        Neighb0 = Neighb;

      return 0;

    default:
      return -1;
  }
}

TBaseCell *TJoint::GetNeighbour(int i)
{
  if (i)
    return Neighb1;
  else
    return Neighb0;
}

void TJoint::Delete(TBaseCell *Me)
{
  if (Neighb0 == Me)
    Neighb0 = NULL;
  else
    Neighb1 = NULL;
}

#ifdef __3D__
void TJoint::SetMapType()
{
  int N_, LocJoint0, LocJoint1, MaxLen, aux;
  const int *TmpFV, *TmpLen;
  TVertex *Vert;

  if (Neighb0 && Neighb1)
  {
    N_ = Neighb0->GetN_Faces();
    for (LocJoint0=0;LocJoint0<N_;LocJoint0++)
      if (Neighb0->GetJoint(LocJoint0) == this) break;

    N_ = Neighb1->GetN_Faces();
    for (LocJoint1=0;LocJoint1<N_;LocJoint1++)
      if (Neighb1->GetJoint(LocJoint1) == this) break;

    Neighb0->GetRefDesc()->GetShapeDesc()->
             GetFaceVertex(TmpFV, TmpLen, MaxLen);

    Vert = Neighb0->GetVertex(TmpFV[LocJoint0 * MaxLen]);

    Neighb1->GetRefDesc()->GetShapeDesc()->
             GetFaceVertex(TmpFV, TmpLen, MaxLen);

    N_ = TmpLen[LocJoint1];
    aux = LocJoint1 * MaxLen;
    for (MapType=0;MapType<N_;MapType++)
      if (Neighb1->GetVertex(TmpFV[aux + MapType]) == Vert) break;

    if (MapType == N_)
    {
      /*
      int i;
      N_ = Neighb0->GetN_Vertices();
      for (i=0;i<N_;i++)
        cout << " test 0:" << i << Neighb0->GetVertex(i) << "  " << (int)
                Neighb0->GetVertex(i) << endl;

      N_ = Neighb0->GetN_Vertices();
      for (i=0;i<N_;i++)
        cout << " test 1:" << i << Neighb1->GetVertex(i) << "  " << (int)
                Neighb1->GetVertex(i) << endl;
      */

      cerr << "Error in SetMapType: could not find vertex" << endl;
      exit (-1);
      return ;
    }
  }
}

void TJoint::GetMapperRef(const int *&MapVerts, const int *&MapFaces)
{
  if (MapType != -1)
    switch (Neighb0->GetType())
    {
      case Tetrahedron: TDatabase::MapperDB[MapTriReg0 + MapType]->
                          GetMapperRef(MapVerts, MapFaces);
                        break;

      case Brick:
      case Hexahedron: TDatabase::MapperDB[MapQuadReg0 + MapType]->
                         GetMapperRef(MapVerts, MapFaces);
                       break;
    }
  else
  {
    cerr << "Error in GetMapperRef: wrong MapType " << MapType << endl;
    exit (-1);
  }
}

void TJoint::GetMapperOrig(const int *&MapVerts, const int *&MapEdges)
{
  if (MapType != -1)
    switch (Neighb0->GetType())
    {
      case Tetrahedron: TDatabase::MapperDB[MapTriReg0 + MapType]->
                          GetMapperOrig(MapVerts, MapEdges);
                        break;

      case Brick:
      case Hexahedron: TDatabase::MapperDB[MapQuadReg0 + MapType]->
                         GetMapperOrig(MapVerts, MapEdges);
                       break;
    }
  else
  {
    cerr << "Error in GetMapperOrig: wrong MapType " << MapType << endl;
    exit (-1);
  }
}
#endif

// Destructor
TJoint::~TJoint()
{
  if(Neighb0)
  { Neighb0 = NULL;}
  
  if(Neighb1)
  { Neighb1 = NULL;}
  
}




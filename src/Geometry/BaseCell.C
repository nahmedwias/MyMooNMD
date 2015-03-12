
#ifdef _MPI
#  include "mpi.h"
#endif

#include <BaseCell.h>
#include <Joint.h>
#include <stdlib.h>

#include <BoundComp3D.h>
#include <BoundFace.h>
#include <Edge.h>


// Constructor
TBaseCell::TBaseCell(TRefDesc *refdesc)
{
 int i, N_;

  RefDesc = refdesc;

  N_ = RefDesc->GetShapeDesc()->GetN_Joints();  
  Joints = new TJoint*[N_];

  for (i=0;i<N_;i++)
   Joints[i] = NULL;

  ClipBoard = 0;
  Phase_ID = 0;
  CellIndex = -1;
  region = 0;
  LayerCell = 0;  
  
#ifdef __3D__   
  N_ = RefDesc->GetN_OrigEdges();
  Edges = new TEdge*[N_];
  for (i=0;i<N_;i++)
   Edges[i] = NULL;
#endif

#ifdef  _MPI
  ClipBoard_Par  = -1;
  SubDomainNumber = 0;
  GlobalCellNo = -1;
  SubDomainLocalCellNo = -1;

  OwnCell=FALSE;
  HaloCell=FALSE;
  SubDomainInterfaceCell=FALSE;
  DependentCell=FALSE;
  N_NeibProcesses = 0;
  NeibProcessesIds = NULL;
#endif
}

// Destructor
TBaseCell::~TBaseCell()
{
  int i, N_;
  TJoint *CurrJoint;

#ifdef __2D__
  N_ = RefDesc->GetN_OrigEdges();
#else
  N_ = RefDesc->GetN_OrigFaces();
#endif
  for (i=0;i<N_;i++)
  {
    CurrJoint = Joints[i];
    switch (CurrJoint->GetType())
    {
      case JointEqN:
      case MortarBaseJoint:
      case InterfaceJoint:
      case IsoInterfaceJoint:
#ifdef __3D__
      case InterfaceJoint3D:
      case IsoInterfaceJoint3D:
#endif
        if (CurrJoint->GetNeighbour(this))
	{ CurrJoint->Delete(this);}
        else
	{ delete CurrJoint;}
        break;

      default:
        delete CurrJoint;
        break;
    }
  }

  delete Joints;
}

// added 25.04.2010 for fixing refinement problem
#ifdef __3D__
void TBaseCell::CorrectBoundaryVertices(TVertex **NewVertices, TJoint **NewJoints)
{
  int N_NewFaces = RefDesc->GetN_Faces();
  const int *FaceVertex;
  int MaxN_VpF;
  
  TBoundComp3D *BoundComp;
  TVertex *BoundVertices[3], *Vertex;
  double x, y, z, t, s;

  for (int i=0;i<N_NewFaces;++i)
  {
    if ( NewJoints[i]->GetType() == BoundaryFace )
    {
      BoundComp = ((TBoundFace*) NewJoints[i])->GetBoundComp();
      
      if ( BoundComp->GetType() == Plane ) continue;
      
      RefDesc->GetFaceVertex(FaceVertex, MaxN_VpF);
      
      for (int j=0;j<MaxN_VpF;++j)
      {
	Vertex = NewVertices[FaceVertex[MaxN_VpF*i+j]];
	
	Vertex->GetCoords(x ,y, z);
	
	BoundComp->GetTSofXYZ(x, y, z, t, s);
	BoundComp->GetXYZofTS(t, s, x, y, z);
	
	Vertex->SetCoords(x, y, z);
      }     
    }
    else continue;
  }
  
}
#endif


// Methods
#ifdef  _MPI
void TBaseCell::SetNeibProcessesIds(int *Neiblist)
{
 int i;

 if(N_NeibProcesses==0)
  {
   printf(" Set the N_NeibProcesses for this cell first !!!! \n");
   MPI_Finalize();
   exit(0);
  }

 if(NeibProcessesIds) delete [] NeibProcessesIds;

 NeibProcessesIds = new int[N_NeibProcesses];

 for(i=0;i<N_NeibProcesses;i++)
  NeibProcessesIds[i] = Neiblist[i];

}
#endif

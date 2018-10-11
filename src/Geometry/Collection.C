// =======================================================================
// @(#)Collection.C        1.2 08/12/99
//
// Class:       TCollection
// Purpose:     store cells in an array
//              used by cell iterators
//
// Author:      Gunar Matthies  14.10.97
//
// History:     14.10.97 Starting implementation
// =======================================================================

#include <algorithm>
#include <Collection.h>
#include <BaseCell.h>
#include <string.h>
#include <JointCollection.h>   
#include <IsoBoundEdge.h>
#include <BoundComp.h>
#include <array>
#include <BoundEdge.h>
#include <BoundFace.h>

#ifdef _MPI
#include <mpi.h>
#endif

/** constructor */
TCollection::TCollection(int n_cells, TBaseCell **cells)
{
  N_Cells = n_cells;
  Cells = cells;

  SortedCells = nullptr;
  Index = nullptr;

  #ifdef  _MPI
  N_OwnCells = 0;
  GlobalIndex = new int[N_Cells];

  for(int i=0; i<N_Cells; i++)
   GlobalIndex[i] = Cells[i]->GetGlobalCellNo();
  #endif
}

void TCollection::GenerateSortedArrays()
{
  if(!SortedCells)
  {
    SortedCells = new TBaseCell*[N_Cells];
    Index = new int[N_Cells];

    memcpy(SortedCells, Cells, N_Cells*sizeof(TBaseCell*));
    std::sort(SortedCells, SortedCells+N_Cells);

    for(int i=0;i<N_Cells;i++)
      Index[GetSortedIndex(Cells[i])] = i;
  }
}

/** destructor: delete arrays */
TCollection::~TCollection()
{
  if(Index) delete [] Index;
  if(SortedCells) delete [] SortedCells;

  if(Cells) delete [] Cells;
}

/** get maximal and minimal diameter */
int TCollection::GetHminHmax(double *hmin, double *hmax)
{
  int i;
  double h_min = 1e10 , h_max= 0, h;
  TBaseCell *cell;
    
  for (i=0;i<N_Cells;i++)
  {
    cell = GetCell(i);
    h = cell->GetDiameter();
    if (h<h_min)
      h_min = h;
    if (h>h_max)
      h_max = h;
  }
  *hmin = h_min;
  *hmax = h_max;
  return(0);
}

/** return Index of cell in SortedCells-array */
int TCollection::GetSortedIndex(TBaseCell *cell)
{
  int Left = 0;
  int Mid;
  int Right = N_Cells - 1;
  int ret = -1;
  TBaseCell *b;

  while (Left <= Right )
  {
    Mid = Left + ((Right - Left) / 2);

    b = SortedCells[Mid];

    if (b == cell)
    {
      ret = Mid;
      break;
    }

    if(b > cell)
      Right = Mid - 1;
    else
      Left = Mid + 1;
  }

  return ret;
}

/** return Index of cell in SortedCells-array */
int TCollection::GetIndex(TBaseCell *cell)
{
  int ret, gsi;

  GenerateSortedArrays();

  gsi = GetSortedIndex(cell);
  ret = (gsi==-1)?(N_Cells+1):Index[gsi];

  return ret;
}

 /** return the Index of the vertex in the sorted array */
int TCollection::GetIndex(TVertex **Array, int Length, TVertex *Element)
{
  int l=0, r=Length, m=(r+l)/2;
  TVertex *Mid;

  Mid=Array[m];
  while(Mid!=Element)
  {
    if(Mid>Element)
    {
      l=m;
    }
    else
    {
      r=m;
    }
    m=(r+l)/2;
    Mid=Array[m];
  }

  return m;
}


/** mark the vertices that are on the boundary */
int TCollection::MarkBoundaryVertices()
{
  int i, j, N_, comp;
  double x0, y0, x2, y2, t0, t1, eps=1e-8;
  const int *TmpEdVer;  
  TBaseCell *cell;
  TVertex *vertex0, *vertex1;
  TRefDesc *refdesc;
  TBoundEdge *boundedge;
  TJoint *joint;

  // initialization 
  // loop over all mesh cells
  for (i=0;i<N_Cells;i++)
  {
    cell = GetCell(i);
    N_ = cell->GetN_Vertices();
    // loop over the vertices
    for (j=0;j<N_;j++)
    {
     vertex1 = cell->GetVertex(j);
     vertex1->SetClipBoard(-1);
    }
  }

  // set ClipBoard for vertices on the boundary
  for (i=0;i<N_Cells;i++)
  {
    cell = GetCell(i);
    // get refinement descriptor
    refdesc=cell->GetRefDesc();                   
    // get information to compute vertices from edges
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
    // number of edges
    N_ = cell->GetN_Edges();
    for (j=0;j<N_;j++)
    {
      joint = cell->GetJoint(j);
	// if boundary edge
	if (joint->GetType() == BoundaryEdge ||
	    joint->GetType() == IsoBoundEdge)
	{
	  boundedge=(TBoundEdge *)joint;
	  auto BoundComp=boundedge->GetBoundComp(); // get boundary component
	  comp=BoundComp->GetID();                  // boundary id
	  boundedge->GetParameters(t0, t1);         // parameter interval
	  vertex0 = cell->GetVertex(TmpEdVer[2*j]);
	  x0 = vertex0->GetX();
	  y0 = vertex0->GetY();
	  vertex1 = cell->GetVertex(TmpEdVer[2*j+1]);
	  boundedge->GetXYofT(t0,x2,y2);
	  if ((fabs(x0-x2)<eps)&&(fabs(y0-y2)<eps))
	    {
	       vertex0->SetClipBoard(comp+t0);
	       vertex1->SetClipBoard(comp+t1);
	    }
	  else
	    {
	       vertex0->SetClipBoard(comp+t1);
	       vertex1->SetClipBoard(comp+t0);
	    }
	  //OutPut(t0 << " " << vertex0->GetClipBoard() << " " << vertex1->GetClipBoard() << endl);
	}
    }
  }
  
  return(0);
  
}

// methods for TJointCollection, 03.11.09  (Sashi)
static void Sort(TJoint **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, *rr, len;
  TJoint *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete [] rr;

}

static void Sort(TVertex **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, *rr, len;
  TVertex *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete [] rr;

}


/** return Index of joints in Cells-array */
TJointCollection *TCollection::GetJointCollection()
{
 int i, j, N_Joints, N, N_RootJoints;
 TBaseCell *Me;
 TJoint **joints, **RootJoints, *Last, *Current;
 TJointCollection *JointColl;


 #ifdef __3D__ 
 joints = new TJoint*[6*N_Cells];
 #else 
 joints = new TJoint*[4*N_Cells];
 #endif

 N=0;
 for(i=0; i<N_Cells; i++)
  {
   Me = Cells[i];
   Me->SetCellIndex(i); // needed for DG matrices assembling
   N_Joints = Me->GetN_Joints();
   
   for(j=0; j<N_Joints; j++) 
    {
     joints[N] = Me->GetJoint(j);
     N++;
    } //for(j=0; j<N_Joints; j++) 
  } //for(i=0; i<N_Cells; i++)


  N--;
  // sort the Vertices array
  Sort(joints, N);
  N++;

  Last=nullptr;
  N_RootJoints=0;
  for(i=0;i<N;i++)
   {
    Current=joints[i];
    if(Current!=Last)
    {
      N_RootJoints++;
      Last=Current;
    }
   }

  RootJoints =  new TJoint*[N_RootJoints];
  Last=nullptr;
  N_RootJoints=0;
  for(i=0;i<N;i++)
   {
    Current=joints[i];
    if(Current!=Last)
    {
      RootJoints[N_RootJoints] = Current;
      Last = Current;
      N_RootJoints++;
    }
   }

  JointColl = new TJointCollection(N_RootJoints, RootJoints);

 return JointColl;
}


// for operator-split nodal point collection, 14.07.2010 (Sashi)
void TCollection::GenerateCellVertNeibs()
{
 int i, j, k, m, N, N_VertInCell, N_LocVertices;
 int Max_N_VertInCell, N_RootVertices, *NumberVertex, *VertexNumbers;
 TVertex *Current, *Last, **Vertices;

   Max_N_VertInCell = 0;
   for(i=0;i<N_Cells;i++)
    if(Max_N_VertInCell<Cells[i]->GetN_Vertices())
      Max_N_VertInCell=Cells[i]->GetN_Vertices();

   Vertices=new TVertex*[Max_N_VertInCell*N_Cells];

   N=0;
   for(i=0;i<N_Cells;i++)
    {
     N_VertInCell = Cells[i]->GetN_Vertices();

     for(j=0;j<N_VertInCell;j++)
      {
       Vertices[N]=Cells[i]->GetVertex(j);
       N++;
      } // j
    } // i
  N_LocVertices = N;
  N--;
  NumberVertex =new int[N_LocVertices];
  VertexNumbers= new int[N_LocVertices];
  // sort the Vertices array based on vertices pointer values
  Sort(Vertices, N);

  Last=nullptr;
  N_RootVertices=-1;
  for(i=0;i<N_LocVertices;i++)
   {
    if((Current=Vertices[i])!=Last)
    {
      N_RootVertices++;
      Last=Current;
    }
    NumberVertex[i]=N_RootVertices;
   }
  N_RootVertices++;

  m=0;
  for(i=0;i<N_Cells;i++)
   {
    k=Cells[i]->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=Cells[i]->GetVertex(j);
      N=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[N];
      m++;
    } // endfor j
  } //endfor i



}


/// @brief create lists with vertex coordinates and element ids
int TCollection::createElementLists()
{
    
  int nVertexPerFace;
  int nBoundaryFaces;

  #ifdef __2D__
  nVertexPerFace = 2;
#else
  if ( Cells[0]->GetType() == Tetrahedron) {
    nVertexPerFace = 3; // vertex per face
  } else {
    nVertexPerFace = 4;
  }
#endif

  // if arrays have been created before
  // free and recreate (to handle multigrid levels)
  if (ElementNodes.size()) {
    ElementNodes.clear();
  }
  if (NodesCoords.size()) {
    NodesCoords.resize(0);
  }
  if (BdFacesNodes.size()) {
    BdFacesNodes.resize(0);
  }

  // create a list with all local vertices
  std::vector<TVertex*> localVertices;
  localVertices.resize(0);
  for(int i=0;i<N_Cells;i++) {
    int k = Cells[i]->GetN_Vertices();
    for(int j=0; j<k; j++) {
      localVertices.push_back(Cells[i]->GetVertex(j));
    }
  }
  NLocVertices = localVertices.size();
  std::sort(localVertices.begin(),localVertices.end());
  // remove duplicate
  auto it = std::unique(localVertices.begin(), localVertices.end());
  localVertices.resize(std::distance(localVertices.begin(), it));
  unsigned int nPoints = localVertices.size();

  // fill the array with nodes coordinates
#ifdef __2D__
  NodesCoords.resize(2*nPoints);
  NodesReferences.resize(nPoints,0.);
  int N_=0;
  for(unsigned int i=0;i<localVertices.size();i++)
  {
    localVertices[i]->GetCoords(NodesCoords[N_],NodesCoords[N_+1]);
    N_ += 2;
  }
#else
  NodesCoords.resize(3*nPoints);
  NodesReferences.resize(nPoints,0.);
  int N_=0;
  for(unsigned int i=0;i<localVertices.size();i++)
  {
    localVertices[i]->GetCoords(NodesCoords[N_],NodesCoords[N_+1],NodesCoords[N_+2]);
    N_ += 3;
  }
#endif
  
  /*
     @attention in the .mesh file, numbering of the vertices within an elements
     starts from 1, i.e. first (see 'VERTEX OFFEST' below)
  */
  // elements array
  ElementNodes.resize(N_Cells);
  ElementReferences.resize(N_Cells);
  
  for(int i=0; i<N_Cells; i++)
  {
    ElementReferences[i] = Cells[i]->GetReference_ID();
    ElementNodes[i].resize(Cells[i]->GetN_Vertices());

    for (int j=0; j<Cells[i]->GetN_Vertices();j++)
    {
      TVertex *current = Cells[i]->GetVertex(j);
      for (unsigned int s=0; s<localVertices.size(); s++)
      {
        if(current == localVertices[s])
        {
          ElementNodes[i][j] = s; // VERTEX OFFSET
          break;
        }
      }  
    }
  }
 

#ifdef __2D__
  nBoundaryFaces = 0;
  int nInterfaceFaces = 0;
  for(int i=0;i<N_Cells;i++){
    for (int j=0; j<Cells[i]->GetN_Edges(); j++) {
      
      TJoint *joint = Cells[i]->GetJoint(j);
      if(!(joint->InnerJoint())) {
	nBoundaryFaces++;	
      } else {
	///@todo check for internal/interface joints
	/*
	// it is a inner joint but it could be an InterfaceJoint
	// get neighbor and properties

	TBaseCell *cell_neighbor = joint->GetNeighbour(Cells[i]);
	int neighbor_reference = cell_neighbor->GetReference_ID();
	int neighbor_index = getIndexInCollection(cell_neighbor); 
	if (neighbor_index == -1) { // neighbor does not belong to this collection
	  nBoundaryFaces++;
	  cout << " interface joint" << endl;
	} else { // neighbor has a different reference
	  int cell_reference = Cells[i]->GetReference_ID();
	  if (cell_reference != neighbor_reference) {
	    if (neighbor_index>i) { // to do this check only once
	      cout << " internal/interface joint" << endl;
	      nBoundaryFaces++;
	      nInternalBdFaces++;
	    }
	  }
	}
	*/
      } 
    } //for (int j=0;j<N_Joints;j++) {

  } //for(int i=0;i<nElements;i++){

  nBoundaryFaces = nBoundaryFaces+nInterfaceFaces;
  BdFacesReferences.resize(nBoundaryFaces);
  BdFacesNodes.resize(nVertexPerFace*nBoundaryFaces);
  
#else
  
  BdFacesReferences.clear();
  for(int i=0;i<N_Cells;i++)
  {
    for (int j=0;j<Cells[i]->GetN_Faces();j++)
    {
      TJoint *joint = Cells[i]->GetJoint(j);
      if(!(joint->InnerJoint()))
      {
	///@todo assign a meaningful reference to the boundary face
	//int local_reference = 1;//Cells[i]->GetJointReference(j);
	TBoundFace* boundface = (TBoundFace *)joint;
	int local_reference = boundface->GetBoundComp()->get_physical_id();
	BdFacesReferences.push_back(local_reference);
      }
    }
  }
  nBoundaryFaces = BdFacesReferences.size();
  BdFacesNodes.resize(nVertexPerFace*nBoundaryFaces);
  
#endif

  int count_boundary_elements = 0;

  for(int i=0;i<N_Cells;i++)
  {
    
#ifdef __2D__

    int nJoints = Cells[i]->GetN_Edges();
    for (int j=0;j<nJoints;j++)
    {
      bool foundVertex1= false;
      bool foundVertex2 = false;
      TJoint *joint = Cells[i]->GetJoint(j);
      if(!(joint->InnerJoint()))
      {
	// 2D: joint n. j is a straight line between vertex j and vertex j+1
	TVertex *v1 = Cells[i]->GetVertex(j);
	TVertex *v2 = Cells[i]->GetVertex((j+1)%nJoints);
	for (unsigned int k=0;k<localVertices.size(); k++)
	{
	  if (v1==localVertices[k])
	  {
	    BdFacesNodes[nVertexPerFace*count_boundary_elements]=k+1;
	    foundVertex1 = true;
	  }
	  if (v2==localVertices[k])
	  {
	    BdFacesNodes[nVertexPerFace*count_boundary_elements+1]=k+1;
	    foundVertex2 = true;
	  }
	  if (foundVertex1&&foundVertex2)
	    break;
	}
	count_boundary_elements++;
	
      } else { //if(!(joint->InnerJoint())) {
	
	///@todo inner/interface joints still to be finished
	// it is a inner joint but it could be an InterfaceJoint
	// get neighbor
	/*TBaseCell *cell_neighbor = joint->GetNeighbour(Cells[i]);
	  int cell_reference = Cells[i]->GetReference_ID();
	  int neighbor_reference = cell_neighbor->GetReference_ID();
	  if (cell_reference < neighbor_reference) {
	  double xE1 =  Cells[i] ->GetVertex(j)->GetX();
	  double yE1 =  Cells[i] ->GetVertex(j)->GetY();
	  double xE2 =  Cells[i] ->GetVertex((j+1)%nJoints)->GetX();
	  double yE2 =  Cells[i] ->GetVertex((j+1)%nJoints)->GetY();
	    
	  for (int k=0;k<nPoints;k++) {
	  double xV = NodesCoords[dim*k];
	  double yV = NodesCoords[dim*k+1];
	  
	  double d1 = (xV-xE1)*(xV-xE1)+(yV-yE1)*(yV-yE1);
	  if (d1<1e-10) {
	  BdFacesNodes[nVertexPerFace*ibd]=k+1;
	  }
	  }
	  for (int k=0;k<nPoints;k++) {
	  double xV = NodesCoords[dim*k];
	  double yV = NodesCoords[dim*k+1];
	  
	  double d1 = (xV-xE2)*(xV-xE2)+(yV-yE2)*(yV-yE2);
	  if (d1<1e-10) {
	  BdFacesNodes[nVertexPerFace*ibd+1]=k+1;
	  }
	  }
	  ibd++;
	  }
	*/
      }
      
    } //for (int j=0;j<N_Joints;j++) {
    
#else
    // 3D case
    
    // ParMooN numbering: face j-th <-> local vertices
    int FaceVertexTetra[][3] = { {0, 1, 2} ,{0, 3, 1}, {2, 1, 3}, {0, 2, 3} };
    int FaceVertexHexa[][4] = { {0, 1, 2, 3} ,{0, 1, 5, 4}, {1, 2, 6, 5},
				{3, 2, 6, 7}, {0, 3, 7, 4}, {4, 5, 6, 7} };
    
    TBaseCell *cell = Cells[i]; 
    for (int j=0;j< cell->GetN_Faces();j++)
    {
      TJoint *joint = cell->GetJoint(j);
      if(!(joint->InnerJoint()))
      {
	
	// loop over the vertices of each face and find
	// the corresponding index in the global vertices list
	for (int kvertex = 0; kvertex<nVertexPerFace; kvertex++)
	{
	  int localVertexIndex;
	  // local index of j-th vertex (cell numeration)
	  if (nVertexPerFace==3)
	    localVertexIndex = FaceVertexTetra[j][kvertex];
	  else
	    localVertexIndex = FaceVertexHexa[j][kvertex];

	  TVertex *v_on_face = cell->GetVertex(localVertexIndex);
	  for (unsigned int k=0; k<localVertices.size(); k++)
	  {
	    if (v_on_face==localVertices[k])
	    {
	      BdFacesNodes[nVertexPerFace*count_boundary_elements + kvertex] =
		k; //+ VERTEX_OFFSET
	    }
	  } 
	}
	
	count_boundary_elements++;
      } 
	  
    }//for (int j=0;j<n_faces
    
#endif

  } // loop over cells
  return 0;
}

//############################################################
//access the data which is generated on createElementList
//############################################################

///@brief Get number of vertices
unsigned int TCollection::GetN_Vertices()
{
  if(this->NodesReferences.size()==0)
    this->createElementLists();
  return this->NodesReferences.size();
}

///@brief Get number of boundary faces
unsigned int TCollection::GetN_BdFaces()
{
    if(this->NodesReferences.size()==0)
      this->createElementLists();
    return this->BdFacesReferences.size();  
}

///@brief direct access to the vector NodesCoords
double TCollection::GetCoord(unsigned int vert)
{
    if(this->NodesReferences.size()==0)
      this->createElementLists();
    return this->NodesCoords.at(vert);
}

///@brief direct access to the vector BdFacesNodes
unsigned int TCollection::GetBdFacesNode(unsigned int node)
{
  if(this->NodesReferences.size()==0)
    this->createElementLists();
  return this->BdFacesNodes[node];
}

///@brief direct access to global number of the jth vertex of the ith cell
unsigned int TCollection::GetGlobalVerNo(unsigned int cell, unsigned int locvert)
{
  if(this->NodesReferences.size()==0)
    this->createElementLists();
  return this->ElementNodes[cell][locvert];
}

///@brief Get the sum of the numbers of local vertices over all cells
unsigned int TCollection::GetNLocVertices()
{
  if(this->NodesReferences.size()==0)
    this->createElementLists();
  return this->NLocVertices;
}


//#############################################################

/** @brief return the index of cell in SortedCells-array 
    1.- now during construction the array GlobalIndex[] is filled
    2.- cell index is compared with the ones in GlobalIndex[]
    3.- return -1 if cell does not belong to collection
    @warning The previous function GetIndex() does not seem to work properly in general
*/
int TCollection::getIndexInCollection(TBaseCell *cell)
{
  for (int i=0; i<N_Cells; i++) {
    if (GlobalIndex[i]==cell->GetCellIndex()) {
      return i;
    }
  }
  return -1;
}

/**
   Notes: 
   (1) we always write a .mesh file with three dimensions. This
   is done for visualization purposes (with medit).
   (2) In two dimensions, mixed meshes (triangles+quadrilaterals) are supported
   (3) In 3D, we only write boundary faces and volume elements
*/
int TCollection::writeMesh(const char *meshFileName)
{
  int dim = 3;
#ifdef __2D__
  dim = 2;
#endif
  
  //index change for .geo output
  unsigned int VERTEX_OFFSET=1;

  int nVertexPerFace;
  if (dim==2) {
    nVertexPerFace = 2;
  } else {
    if ( Cells[0]->GetType() == Tetrahedron) {
      nVertexPerFace = 3; // vertex per face
    } else {
      nVertexPerFace = 4;
    }
  }
  int nVertexPerElement = Cells[0]->GetN_Vertices();

  createElementLists();

  std::ofstream MESHfile; 
  // header of .mesh file
  MESHfile.open(meshFileName);
  MESHfile << "MeshVersionFormatted 1" << endl;
  MESHfile << endl;
  MESHfile << "Dimension 3" << endl; //note: dim always 3 (for visualization)
  MESHfile << endl;
  MESHfile << "Vertices" << endl;
  unsigned int nPoints = NodesReferences.size();
  MESHfile << nPoints << endl;
  for (unsigned int i=0; i<nPoints; i++) {
    for (int j=0; j<dim; j++) {
     MESHfile << NodesCoords[i*dim+j] << "  ";
    }
    if (dim==2) {
      MESHfile << "0.0000  " ;
    }
    MESHfile << NodesReferences[i] << endl;
  }
  MESHfile << endl;

  // faces (egdes in 2D, triangles/quads in 3D)
  unsigned int nBoundaryFaces = BdFacesReferences.size();

  // write elements: edges + surface elements in 2D, boundary faces + volume el. in 3D
  if (dim==2) {
    MESHfile << "Edges" << endl;
    MESHfile << nBoundaryFaces << endl;
    for (unsigned int i=0; i<nBoundaryFaces; i++) {
      MESHfile << BdFacesNodes[nVertexPerFace*i] << "  " 
	       << BdFacesNodes[nVertexPerFace*i+1] << " "  
	       << BdFacesReferences[i] << endl;
    }
    MESHfile << endl;
    
    // write cells:
    // In order to support mixed (tria+quad) meshes
    // we need to store trias and quads in two tmp
    // separate lists, then write them on the file
    std::vector<int> nodesTria;
    std::vector<int> nodesQuad;
    
    // Note: we store nodes (3 or 4) + reference
    for (int i=0; i<N_Cells; i++) {
      if (Cells[i]->GetN_Vertices() == 3) {

	for (int j=0;j<3;j++) 
	  nodesTria.push_back(ElementNodes[i][j]+VERTEX_OFFSET);
	nodesTria.push_back(ElementReferences[i]);
	
      } else if (Cells[i]->GetN_Vertices() == 4) {
	
	for (int j=0;j<4;j++) 
	  nodesQuad.push_back(ElementNodes[i][j]+VERTEX_OFFSET);
	nodesQuad.push_back(ElementReferences[i]);
      }
    }
    unsigned int nTrias = nodesTria.size()/4;
    unsigned int nQuads = nodesQuad.size()/5;
    cout << " I found " << nTrias << " and " << nQuads << " quads " << endl;
    // now we write the list of triangles
    if (nodesTria.size()) {
      MESHfile << "Triangles" << endl;
      MESHfile << nTrias << endl;
      for (unsigned int i=0; i<nTrias; i++) {
	for (int j=0;j<4;j++) {
	  MESHfile << nodesTria[i*4+j] << "  ";
	}
	MESHfile << endl;
      } 
    }

    if (nodesQuad.size()) {
      MESHfile << "Quadrilaterals" << endl;
      MESHfile << nQuads << endl;
      for (unsigned int i=0; i<nQuads; i++) {
	for (int j=0;j<5;j++) {
	  MESHfile << nodesQuad[i*5+j] << "  ";
	}
	MESHfile << endl;
      } 
    }
    

    MESHfile << endl;
    MESHfile << "End" << endl;
    
  } else { // dim=3

    // surface elements: mixed meshes not allowed in 3D
    if (nVertexPerFace == 3) {
      MESHfile << "Triangles" << endl;
    } else {
      MESHfile << "Quadrilaterals" << endl;
    }
    MESHfile << nBoundaryFaces << endl;
    for (unsigned int i=0; i<nBoundaryFaces; i++) {
      for (int j=0;j<nVertexPerFace;j++) {
	MESHfile << BdFacesNodes[nVertexPerFace*i+j] << "  ";
      }
      MESHfile << BdFacesReferences[i] << endl;
    }// for (int i=0; i<nBdFaces; i++) {
    MESHfile << endl;

    // write cells
    if (nVertexPerElement == 4) {
      MESHfile << "Tetrahedra" << endl;
    } else {
      MESHfile << "Hexahedra" << endl;
    }
    MESHfile << N_Cells << endl;
    for (int i=0; i<N_Cells ; i++) {
      for (int j=0;j<nVertexPerElement;j++) {
	MESHfile << ElementNodes[i][j]+VERTEX_OFFSET << "  ";
      }
      MESHfile << ElementReferences[i] << endl;
    } // for (int i=0; i<nElements; i++) {
    MESHfile << endl;
    MESHfile << "End" << endl;
  
  } // if dim==2

  MESHfile.close();
  cout << "TCollection::writeMesh mesh written on " << meshFileName << endl;

  
  return 0;
  
}

void TCollection::get_edge_list_on_component(int id,std::vector<TBoundEdge*> &edges)
{
  edges.clear();
  for(int i=0;i<this->N_Cells; i++)
    {
      TBaseCell *cell = this->Cells[i];
      for(int j=0;  j < cell->GetN_Joints(); j++)
        {
    TJoint *joint= cell->GetJoint(j);
    if (joint->GetType()==BoundaryEdge)
            {
        TBoundEdge *boundedge = (TBoundEdge *)joint;
        const TBoundComp *BoundComp = boundedge->GetBoundComp();
     if (BoundComp->GetID() == id)
                {
    ///@todo set the boundedge properties in the function MakeGrid
    boundedge->SetNeighbour(cell);
    boundedge->set_index_in_neighbour(cell,j);
    edges.push_back(boundedge);                 
                }
            }
        }
    }
}

void TCollection::get_boundary_edge_list(std::vector<TBoundEdge*> &edges)
{
  edges.clear();
  for(int i = 0; i < this->N_Cells; i++)
  {
    TBaseCell *cell = this->Cells[i];
    for(int j = 0; j < cell->GetN_Joints(); j++)
    {
      TJoint *joint = cell->GetJoint(j);
      if (joint->GetType() == BoundaryEdge)
      {
        TBoundEdge *boundedge = (TBoundEdge *)joint;
          ///@todo set the boundedge properties in the function MakeGrid
          boundedge->SetNeighbour(cell);
          boundedge->set_index_in_neighbour(cell,j);
          edges.push_back(boundedge);                 
      }
    }
  }
}

//New LB 11.10.18
#ifdef __3D__
void TCollection::get_face_list_on_component(int boundary_component_id, std::vector<TBoundFace*> &faces)
{
  faces.clear();
  for (int i = 0; i < this->N_Cells; i++)
  {
    TBaseCell *cell = this->Cells[i];

    for(size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        if (boundface->GetBoundComp()->get_physical_id() == boundary_component_id)
        {
          ///@todo set the boundedge properties in the function MakeGrid
          //boundface->SetNeighbour(cell); //todo check if ok
          //boundface->set_index_in_neighbour(cell, joint_id); //todo check if ok
          faces.push_back(boundface);
        }
      }
    }
  }
}
#endif

#ifdef _MPI
int TCollection::find_process_of_point(double x, double y, double z) const
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_rank, size;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &my_rank);

  bool found_on_proc = false;
  for (int i=0;i<N_Cells;i++)
  {
    TBaseCell* cell = GetCell(i);
    if(!cell->IsHaloCell())
    {
      if(cell->PointInCell(x,y,z))
      {
        found_on_proc = true;
        break;
      }
    }
  }

  // each process sends either its rank (when found) or "size",
  // the sent values get reduced to their minimum
  int sendbuf[1] = {size};
  if (found_on_proc)
    sendbuf[0] = my_rank;
  int recvbuf[1] = {0};

  MPI_Allreduce(sendbuf, recvbuf, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  int containing_proc = recvbuf[0];

  if(containing_proc == size)
  {
    ErrThrow("Point (",x,",",y,",",z,") not found on any process!");
  }

  return containing_proc;


}



#endif

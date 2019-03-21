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

  #ifdef  _MPI
  N_OwnCells = 0;
  #endif
}

/** destructor: delete arrays */
TCollection::~TCollection()
{
  if(Cells) delete [] Cells;
}

/** get maximal and minimal diameter */
int TCollection::GetHminHmax(double *hmin, double *hmax) const
{
  double h_min = 1e10;
  double h_max = 0;
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = GetCell(i);
    double h = cell->GetDiameter();
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
int TCollection::get_cell_index(const TBaseCell *cell) const
{
  if(cell_to_index_map.empty()) // only done once
  {
    for(int i = 0; i < N_Cells; ++i)
    {
      auto it = cell_to_index_map.insert({Cells[i], i});
      if(it.second == false)
      {
        ErrThrow("cell ", Cells[i],
                 " appears more than once in this collection. Indices ",
                 it.first->second, " and ", i);
      }
    }
  }
  return cell_to_index_map.at(cell);
}

/// @brief create lists with vertex coordinates and element ids
void TCollection::createElementLists() const
{
    
  int nVertexPerFace;
  int nBoundaryFaces;

#ifdef __2D__
  nVertexPerFace = 2;
#else  
  nVertexPerFace = Cells[0]->GetType() == Tetrahedron ? 3 : 4;
#endif

  // if arrays have been created before free and recreate (to handle multigrid
  // levels). Replace by new (empty) object.
  element_lists = ElementLists();

  // create a list with all local vertices
  std::vector<TVertex*> localVertices;
  localVertices.resize(0);
  for(int i=0;i<N_Cells;i++) {
    int k = Cells[i]->GetN_Vertices();
    for(int j=0; j<k; j++) {
      localVertices.push_back(Cells[i]->GetVertex(j));
    }
  }
  element_lists.NLocVertices = localVertices.size();
  std::sort(localVertices.begin(),localVertices.end());
  // remove duplicate
  auto it = std::unique(localVertices.begin(), localVertices.end());
  localVertices.resize(std::distance(localVertices.begin(), it));
  unsigned int nPoints = localVertices.size();

  // fill the array with nodes coordinates
#ifdef __2D__
  element_lists.NodesCoords.resize(2*nPoints);
  element_lists.NodesReferences.resize(nPoints,0.);
  int N_=0;
  for(unsigned int i=0;i<localVertices.size();i++)
  {
    localVertices[i]->GetCoords(element_lists.NodesCoords[N_],
                                element_lists.NodesCoords[N_+1]);
    N_ += 2;
  }
#else
  element_lists.NodesCoords.resize(3*nPoints);
  element_lists.NodesReferences.resize(nPoints,0.);
  int N_=0;
  for(unsigned int i=0;i<localVertices.size();i++)
  {
    localVertices[i]->GetCoords(element_lists.NodesCoords[N_],
                                element_lists.NodesCoords[N_+1],
                                element_lists.NodesCoords[N_+2]);
    N_ += 3;
  }
#endif
  
  /*
     @attention in the .mesh file, numbering of the vertices within an elements
     starts from 1, i.e. first (see 'VERTEX OFFEST' below)
  */
  // elements array
  element_lists.ElementNodes.resize(N_Cells);
  element_lists.ElementReferences.resize(N_Cells);
  
  for(int i=0; i<N_Cells; i++)
  {
    element_lists.ElementReferences[i] = Cells[i]->GetReference_ID();
    element_lists.ElementNodes[i].resize(Cells[i]->GetN_Vertices());

    for (int j=0; j<Cells[i]->GetN_Vertices();j++)
    {
      TVertex *current = Cells[i]->GetVertex(j);
      for (unsigned int s=0; s<localVertices.size(); s++)
      {
        if(current == localVertices[s])
        {
          element_lists.ElementNodes[i][j] = s; // VERTEX OFFSET
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
      
      auto joint = Cells[i]->GetJoint(j);
      if(!(joint->InnerJoint()))
      {
        nBoundaryFaces++;	
      }
      else
      {
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
  element_lists.BdFacesReferences.resize(nBoundaryFaces);
  element_lists.BdFacesNodes.resize(nVertexPerFace*nBoundaryFaces);
  
#else
  
  //BdFacesReferences.clear();
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
        element_lists.BdFacesReferences.push_back(local_reference);
      }
    }
  }
  nBoundaryFaces = element_lists.BdFacesReferences.size();
  element_lists.BdFacesNodes.resize(nVertexPerFace*nBoundaryFaces);
  
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
	    element_lists.BdFacesNodes[nVertexPerFace*count_boundary_elements]=k+1;
	    foundVertex1 = true;
	  }
	  if (v2==localVertices[k])
	  {
	    element_lists.BdFacesNodes[nVertexPerFace*count_boundary_elements+1]=k+1;
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
	      element_lists.BdFacesNodes[nVertexPerFace*count_boundary_elements + kvertex] =
		k; //+ VERTEX_OFFSET
	    }
	  } 
	}
	
	count_boundary_elements++;
      } 
	  
    }//for (int j=0;j<n_faces
    
#endif

  } // loop over cells
}

//############################################################
//access the data which is generated on createElementList
//############################################################

///@brief Get number of vertices
unsigned int TCollection::GetN_Vertices() const
{
  if(this->element_lists.empty())
    this->createElementLists();
  return this->element_lists.NodesReferences.size();
}

///@brief Get number of boundary faces
unsigned int TCollection::GetN_BdFaces() const
{
    if(this->element_lists.empty())
      this->createElementLists();
    return this->element_lists.BdFacesReferences.size();  
}

///@brief direct access to the vector NodesCoords
double TCollection::GetCoord(unsigned int vert) const
{
    if(this->element_lists.empty())
      this->createElementLists();
    return this->element_lists.NodesCoords.at(vert);
}

///@brief direct access to the vector BdFacesNodes
unsigned int TCollection::GetBdFacesNode(unsigned int node) const
{
  if(this->element_lists.empty())
    this->createElementLists();
  return this->element_lists.BdFacesNodes[node];
}

///@brief direct access to global number of the jth vertex of the ith cell
unsigned int TCollection::GetGlobalVerNo(unsigned int cell, unsigned int locvert) const
{
  if(this->element_lists.empty())
    this->createElementLists();
  return this->element_lists.ElementNodes[cell][locvert];
}

///@brief Get the sum of the numbers of local vertices over all cells
unsigned int TCollection::GetNLocVertices() const
{
  if(this->element_lists.empty())
    this->createElementLists();
  return this->element_lists.NLocVertices;
}


//#############################################################

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
  unsigned int nPoints = element_lists.NodesReferences.size();
  MESHfile << nPoints << endl;
  for (unsigned int i=0; i<nPoints; i++) {
    for (int j=0; j<dim; j++) {
     MESHfile << element_lists.NodesCoords[i*dim+j] << "  ";
    }
    if (dim==2) {
      MESHfile << "0.0000  " ;
    }
    MESHfile << element_lists.NodesReferences[i] << endl;
  }
  MESHfile << endl;

  // faces (egdes in 2D, triangles/quads in 3D)
  unsigned int nBoundaryFaces = element_lists.BdFacesReferences.size();

  // write elements: edges + surface elements in 2D, boundary faces + volume el. in 3D
  if (dim==2) {
    MESHfile << "Edges" << endl;
    MESHfile << nBoundaryFaces << endl;
    for (unsigned int i=0; i<nBoundaryFaces; i++) {
      MESHfile << element_lists.BdFacesNodes[nVertexPerFace*i] << "  " 
	       << element_lists.BdFacesNodes[nVertexPerFace*i+1] << " "  
	       << element_lists.BdFacesReferences[i] << endl;
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
	  nodesTria.push_back(element_lists.ElementNodes[i][j]+VERTEX_OFFSET);
	nodesTria.push_back(element_lists.ElementReferences[i]);
	
      } else if (Cells[i]->GetN_Vertices() == 4) {
	
	for (int j=0;j<4;j++) 
	  nodesQuad.push_back(element_lists.ElementNodes[i][j]+VERTEX_OFFSET);
	nodesQuad.push_back(element_lists.ElementReferences[i]);
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
	MESHfile << element_lists.BdFacesNodes[nVertexPerFace*i+j] << "  ";
      }
      MESHfile << element_lists.BdFacesReferences[i] << endl;
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
	MESHfile << element_lists.ElementNodes[i][j]+VERTEX_OFFSET << "  ";
      }
      MESHfile << element_lists.ElementReferences[i] << endl;
    } // for (int i=0; i<nElements; i++) {
    MESHfile << endl;
    MESHfile << "End" << endl;
  
  } // if dim==2

  MESHfile.close();
  cout << "TCollection::writeMesh mesh written on " << meshFileName << endl;

  
  return 0;
  
}

void TCollection::get_edge_list_on_component(
  int id, std::vector<TBoundEdge*> &edges) const
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

void TCollection::get_boundary_edge_list(std::vector<TBoundEdge*> &edges) const
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
void TCollection::get_face_list_on_component(int boundary_component_id,
                                             std::vector<TBoundFace*> &faces)
const
{
  faces.clear();
  for (int i = 0; i < this->N_Cells; i++)
  {
    TBaseCell *cell = this->Cells[i];
//cout << "Jooooooo" << endl;
    for(size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);
     // cout << "joint_id: " << joint_id << endl;
      if ((joint->GetType() == BoundaryFace) || (joint->GetType() == IsoBoundFace))
      {
        //cout << "joint_id on boundary: " << joint_id << endl;
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        if (boundface->GetBoundComp()->get_physical_id() == boundary_component_id)
        {
          ///@todo set the boundedge properties in the function MakeGrid
          boundface->SetNeighbour(cell); //todo check if ok

          boundface->set_index_in_neighbour(cell, joint_id); //todo check if ok
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

#include <fstream>

#include "Mesh.h"
#include "MooNMD_Io.h"
#include <assert.h>

// default initialization of the mesh (dimension = 0, no elements)
Mesh::Mesh() {
  dimension = -1;
  vertex.resize(0);
  edge.resize(0);
  triangle.resize(0);
  quad.resize(0);
  tetra.resize(0);
  hexa.resize(0);
  meshTrifaceHash.resize(0);
  faceToTetra.resize(0);
  n_boundary_faces = 0;
  boundaryFacesMarker.resize(0);
  hasBothTriaAndQuads = false;
}

// initialize from a file
// note: dimension is set to 2, but it will be changed (if necessary) reading the file
Mesh::Mesh(std::string filename) {
  dimension = -1;
  vertex.resize(0);
  edge.resize(0);
  triangle.resize(0);
  quad.resize(0);
  tetra.resize(0);
  hexa.resize(0);
  meshTrifaceHash.resize(0);
  faceToTetra.resize(0);
  n_boundary_faces = 0;
  boundaryFacesMarker.resize(0);
  hasBothTriaAndQuads = false;
  
  readFromFile(filename);
}

Mesh::Mesh(std::string filename,std::string filenameBoundary) {
  dimension = -1;
  vertex.resize(0);
  edge.resize(0);
  triangle.resize(0);
  quad.resize(0);
  tetra.resize(0);
  hexa.resize(0);
  hasBothTriaAndQuads = false;
  meshTrifaceHash.resize(0);
  faceToTetra.resize(0);
  n_boundary_faces = 0;
  boundaryFacesMarker.resize(0);
  
  readFromFile(filename);
  setBoundary(filenameBoundary);
}


// read the mesh data from a file (.mesh)
void Mesh::readFromFile(std::string filename)
{
  std::ifstream ifile;
  ifile.open(filename.c_str());
  if (!ifile)
  {
    Output::print(" *** Error(Mesh::readFromFile) I could not open ",filename);
    exit(-1);
  };

  std::string line;
  // read header and dimension
  do {
    // read the whole line up to [enter], if the line is empty or starting with #
    getline(ifile,line,'\n');
    stripSpace(line); 
  } while (line.find("Dimension")==std::string::npos);

  line.erase(0,9);
  if (line.length()>0)
    sscanf(line.c_str(),"%d",&dimension);
  else
    ifile >> dimension;
  
  Output::print<4>("Read dimension: ",dimension);
  if (dimension==-1)
  {
    Output::print<1>(" ** ERROR: I could not read the mesh dimension. Check the .mesh file **");
    exit(1);
  }
  // read nodes
  unsigned int numberOfNodes;
  do {
    getline(ifile,line,'\n');
    stripSpace(line);     
  } while (line.find("Vertices")==std::string::npos);
  ifile >> numberOfNodes;
  Output::print<4>("Read n of nodes: ",numberOfNodes);

  
  vertex.resize(numberOfNodes);  
  for (unsigned int i=0; i<numberOfNodes; i++) {
    if (dimension==2) {
      ifile >> vertex[i].x >> vertex[i].y >> vertex[i].reference;
    } else {
      ifile >> vertex[i].x >> vertex[i].y >> vertex[i].z >> vertex[i].reference;
    }
  }
  int current_line_number = ifile.tellg();


  // read edges
  unsigned int numberOfEdges;
  do {
    getline(ifile,line,'\n');
    stripSpace(line); // borra espacios en blanco
  } 
  while (line.find("Edges")==std::string::npos && !ifile.eof());

  if (line.find("Edges")!=std::string::npos) {
    ifile >> numberOfEdges;
   
    if (numberOfEdges) {
      edge.resize(numberOfEdges);
      for (unsigned int i=0; i<numberOfEdges; i++) {
	ifile >> edge[i].nodes[0] >> edge[i].nodes[1] >> edge[i].reference;
      }
    }
  }

  ifile.clear();
  ifile.seekg(current_line_number);
  
  // Triangles
  unsigned int numberOfTriangles = 0;
  do {
    getline(ifile,line,'\n');
    stripSpace(line); 
    if (line.find("#")!=std::string::npos)
      line="####COMMENT###";
  } while (line.find("Triangles")==std::string::npos && !ifile.eof());
  
  if (line.find("Triangles")!=std::string::npos) {
    ifile >> numberOfTriangles;
    
    if (numberOfTriangles) {
      triangle.resize(numberOfTriangles);
      for (unsigned int i=0;i<numberOfTriangles;i++) {
	for (unsigned int k=0;k<3; k++) 
	  ifile >> triangle[i].nodes[k];
	ifile >> triangle[i].reference;
      }
    }
    
  }
  ifile.clear();
  ifile.seekg(current_line_number);

  // Quads
  unsigned int numberOfQuads = 0;
  do {
    getline(ifile,line,'\n');
    stripSpace(line); 
    if (line.find("#")!=std::string::npos)
      line="####COMMENT###";
  } while (line.find("Quadrilaterals")==std::string::npos && !ifile.eof());
  
  if (line.find("Quadrilaterals")!=std::string::npos) {
    ifile >> numberOfQuads;
   
    if (numberOfQuads) {
      quad.resize(numberOfQuads);
      for (unsigned int i=0;i<numberOfQuads;i++) {
	for (unsigned int k=0;k<4; k++) 
	  ifile >> quad[i].nodes[k];
	ifile >> quad[i].reference;
      }
    }
    
  }
  ifile.clear();
  ifile.seekg(current_line_number);

  if (quad.size() && triangle.size() )
  {
    hasBothTriaAndQuads = true;
  }
  
  // tetra
  unsigned int numberOfTetra = 0;
  do {
    getline(ifile,line,'\n');
    stripSpace(line); 
    if (line.find("#")!=std::string::npos)
      line="####COMMENT###";
  } while (line.find("Tetrahedra")==std::string::npos && !ifile.eof());
  
  if (line.find("Tetrahedra")!=std::string::npos) {
    ifile >> numberOfTetra;
    
    if (numberOfTetra) {
      tetra.resize(numberOfTetra);
      for (unsigned int i=0;i<numberOfTetra;i++) {
	for (unsigned int k=0;k<4; k++) 
	  ifile >> tetra[i].nodes[k];
	ifile >> tetra[i].reference;
      }
    }
    
  }
  ifile.clear();
  ifile.seekg(current_line_number);



  // hexa
  unsigned int numberOfHexa = 0;
  do {
    getline(ifile,line,'\n');
    stripSpace(line); 
    if (line.find("#")!=std::string::npos)
      line="####COMMENT###";
  } while (line.find("Hexahedra")==std::string::npos && !ifile.eof());
  
  if (line.find("Hexahedra")!=std::string::npos) {
    ifile >> numberOfHexa;
    
    if (numberOfHexa) {
      hexa.resize(numberOfHexa);
      for (unsigned int i=0;i<numberOfHexa;i++) {
	for (unsigned int k=0;k<8; k++) 
	  ifile >> hexa[i].nodes[k];
	ifile >> hexa[i].reference;
      }
    }
    
  }
  ifile.clear();
  ifile.seekg(current_line_number);
  ifile.close();
  
  if(numberOfHexa != 0)
  {
    // why is this not implemented???
    ErrThrow("Missing implementation for reading mesh files with Hexahedra");
  }
  
  if(numberOfTetra + numberOfHexa == 0 && dimension == 3)
  {
    // this should be 2D and is reset here after a check.
    // All vertices should have the same z-component
    // we don't support 2D meshes embedded in 3D (yet)
    for (unsigned int i=1; i<numberOfNodes; i++) {
      if(vertex[i].z != vertex[0].z)
        ErrThrow("Cannot read a 3D mesh with no tetrahedra and no hexahedra.");
    }
    dimension = 2;
  }
}


/** 
    This function creates a map between
    face(a,b,c) -> hash = a+b+c
    the map is contained in a vector<int*> meshTrifaceHash, where
    * meshTrifaceHash[hash][0] is the number of faces with a+b+c = hast
    * meshTrifaceHash[hash][1],...,meshTrifaceHash[hash][p] point to these faces
    
    

    */
void Mesh::hashTriFaces()
{
  int n_points = vertex.size();
  meshTrifaceHash.resize(3*n_points);
  for(auto &e : meshTrifaceHash)
  {
    e=nullptr;
  }
  std::vector<int> BucketCount(3*n_points);
  
  for(size_t i=0;i<triangle.size();++i)
  {
    int hash = triangle[i].nodes[0]+triangle[i].nodes[1]+triangle[i].nodes[2];
    (BucketCount[hash])++;
  }

  for(size_t i=0;i<triangle.size();++i)
  {
    int hash = triangle[i].nodes[0]+triangle[i].nodes[1]+triangle[i].nodes[2];

    if(meshTrifaceHash.at(hash) == nullptr)
    {
      meshTrifaceHash.at(hash) = new int[BucketCount.at(hash)+1];
      meshTrifaceHash.at(hash)[0] = 1;
      meshTrifaceHash.at(hash)[1] = i;
    }
    else
    {
      int pos = ++(meshTrifaceHash.at(hash)[0]);
      meshTrifaceHash.at(hash)[pos] = i;
    }
  }

  /*
  // printout the whole hash vector
  for(int i=0;i<triangle.size();++i)
  {
    int hash = triangle[i].nodes[0]+triangle[i].nodes[1]+triangle[i].nodes[2];
    cout << "hash = " << hash;
    if (meshTrifaceHash.at(hash)==nullptr) {
      cout << " no faces " << endl;
    } else {
      cout << " total: " << meshTrifaceHash.at(hash)[0] << " -> "; 
      for (int k=1; k<=meshTrifaceHash.at(hash)[0]; k++){
	cout << meshTrifaceHash.at(hash)[k] << " ";
      }
      cout << endl;
    }
  }
  */
}

// fill the list mapping each face to the neighbors tetra
// -1 indicates that a face is on the boundary
///@todo this function could be extended to hexa as well
void Mesh::createFaceToTetrahedraMap()
{

  // create the meshTrifaceHash vector (used in findTriFace)
  if (meshTrifaceHash.size()==0) {
    Output::print("Mesh::createFaceTetrahedraMap() is creating the hash vector");
    this->hashTriFaces();
  }

  // this vector describes the order in which the local vertices of a tetrahedra
  // (0,1,2,3) appear on the faces
  int FaceVertex[][3] = { {0, 1, 2} ,{0, 3, 1}, {2, 1, 3}, {0, 2, 3} };

  this->faceToTetra.resize(this->triangle.size());

  for(unsigned int i=0;i<this->triangle.size();++i)
  {
    // set all values to -1
    this->faceToTetra[i].resize(2);
    this->faceToTetra[i][0] = -1;
    this->faceToTetra[i][1] = -1;
  }
 
  for(unsigned int i=0;i<this->tetra.size();++i)
  {
    for(int face=0;face<4;++face)
    {
      // find a triface with the same vertices
      int triface = findTriFace(tetra[i].nodes[FaceVertex[face][0]],
				tetra[i].nodes[FaceVertex[face][1]],
				tetra[i].nodes[FaceVertex[face][2]]);
      assert (triface != -1);
      // if the element [triface,0] is still = -1, set it = to the current tetrahedra (i)
      if(this->faceToTetra[triface][0] == -1)
      {
	this->faceToTetra[triface][0] = i;
      }
      else
      {
	// if the face has been already associated to a tetrahedra,
	// e.g. faceToTetra[triface][0]=k, then set faceToTetra[triface][1] equal to i
	assert(this->faceToTetra[triface][1] == -1);
        this->faceToTetra[triface][1]= i;
	Output::print<3>("tetrahedra ",i," faces: ",
		      this->faceToTetra[triface][0]," ",this->faceToTetra[triface][1]);
      }
      
    }
  }

}


// find the index of the face with vertices a,b,c
// return -1 if face is not found
///@todo extend this to quadrilateral faces
int Mesh::findTriFace(int a, int b, int c)
{
  int hash = a+b+c;

  if (this->meshTrifaceHash.size()==0) {
    Output::print("Mesh::findTriFace() is creating the hash vector");
    this->hashTriFaces();
  }

  assert(meshTrifaceHash[hash] != nullptr);
  int count = meshTrifaceHash[hash][0];

  for(int i=1;i<=count;++i)
  {
    int triface = meshTrifaceHash[hash][i];
    int found = 0;
    for(int j=0;j<3;++j)
    {
      int vertex = triangle[triface].nodes[j];

      if(a == vertex || b == vertex || c == vertex)
        ++found;
    }
    if(found == 3)
      return triface;
  }
  Output::print("Mesh::findTriFace() ** WARNING ** no face with vertices ",
		a,",",b,",",c," found.");
  return -1;
}

void Mesh::computeNumberOfBoundaryFaces()
{
  int n_inner_faces = 0;
  this->n_boundary_faces = 0;

  if (faceToTetra.size()==0)
  {
    Output::print("Mesh::computeNumberOfBoundaryFaces() is creating the faceToTetra map");
    this->createFaceToTetrahedraMap();
  }

  for (unsigned int i=0; i<faceToTetra.size(); i++)
  {
    // increase the number of boundary faces if face[i] is on the boundary
    if ( (faceToTetra[i][0] < 0)  || (faceToTetra[i][1] < 0) ) 
      this->n_boundary_faces++;
    else
      n_inner_faces++;
  }

  Output::print("Mesh::computeNumberOfBoundaryFaces() - I found ",
		this->n_boundary_faces, " boundary faces and ",
		n_inner_faces, " inner faces");

  // check if faces are missing
  if ( triangle.size() != (2*tetra.size()+n_boundary_faces/2) )
  {
    Output::print("Mesh::computeNumberOfBoundaryFaces() ** WARNING: it looks like some face is missing, ",
		  " n.triangles = ", triangle.size(), ", expected = ",2*tetra.size()+n_boundary_faces/2, " **");
  }
  
}

/**
   This function checks if the mesh read from file contains all the faces
   (inner+boundary). If not, the missing faces are created.

   Note: this function create a local hash vector (to handle the faces
   read from the original file). This vector is deleted at the end, and
   the meshTrifaceHash has to be created later.
   
   This function should be called only once.
 */
void Mesh::createInnerFaces()
{
  Output::print("Mesh::createInnerFaces()");
  
  std::vector< std::vector<int> > localHash;
  int n_points = vertex.size();
  localHash.resize(3*n_points);
  for(unsigned int i=0; i<localHash.size(); i++)
    localHash[i].resize(0);
  
  for(size_t i=0;i<triangle.size();++i)
  {
    int hash = triangle[i].nodes[0]+triangle[i].nodes[1]+triangle[i].nodes[2];

    if( localHash[hash].size()==0 )
    {
      localHash[hash].resize(2);
      localHash[hash][0] = 1;
      localHash[hash][1] = i;
    }
    else
    {
      localHash[hash][0]++;
      localHash[hash].push_back(i);
    }
  }
  // now check if there are missing faces
  int n_created_faces = 0;
  // map vertex id -> face
  int FaceVertex[][3] = { {0, 1, 2} ,{0, 3, 1}, {2, 1, 3}, {0, 2, 3} };
  // loop over tetra
  for (unsigned int i=0; i<tetra.size(); i++)
  {
    
    // loop over faces
    for (unsigned int face=0; face<4; face++)
    {
      int a = tetra[i].nodes[FaceVertex[face][0]];
      int b = tetra[i].nodes[FaceVertex[face][1]];
      int c = tetra[i].nodes[FaceVertex[face][2]];
      int hash = a + b + c;
      // check if face j is in the list of triangles
      if(localHash[hash].size() == 0)
      {
	// face is missing
	n_created_faces++;
	meshTriangle new_triangle;
	new_triangle.nodes[0] = a;
	new_triangle.nodes[1] = b;
	new_triangle.nodes[2] = c;
	new_triangle.reference = 0;
	triangle.push_back(new_triangle);
	Output::print<4>(" Create face ",a," ",b," ",c);
	// update the hash vector (to avoid creating twice the same face)
	localHash[hash].resize(2);
	localHash[hash][0] = 1;
	localHash[hash][1] = triangle.size()-1;
      }
      else
      {
	int count = localHash[hash][0];

	int faces_found = 0;
	for(int i=1;i<=count;++i)
	{
	  int triface = localHash[hash][i];
	  int found = 0;
	  for(int j=0;j<3;++j)
	  {
	    int vertex = triangle[triface].nodes[j];
	    if(a == vertex || b == vertex || c == vertex)
	      ++found;
	  }
	  if(found == 3)
	    faces_found++;	    
	}
	if (faces_found==0)
	{
	  // face is missing
	  n_created_faces++;
	  Output::print<4>(" Create face ",a," ",b," ",c);
	  meshTriangle new_triangle;
	  new_triangle.nodes[0] = a;
	  new_triangle.nodes[1] = b;
	  new_triangle.nodes[2] = c;
	  new_triangle.reference = 0;
	  triangle.push_back(new_triangle);
	  // update the hash vector (to avoid creating twice the same face)
	  localHash[hash][0]++;
	  localHash[hash].push_back(triangle.size()-1);//
	}
      }
      
    } // loop over faces
  } // loop over tetra
  Output::print("Mesh::createInnerFaces() finished. Created ",
		n_created_faces, " triangles");
  
}
// read a PRM file into a Boundary class
void Mesh::setBoundary(std::string PRM)
{
  boundary.initFromFile(PRM);
}

// write the mesh onto a .mesh file
void Mesh::writeToMesh(std::string filename)
{
  std::ofstream ofile;
  ofile.open(filename.c_str(),std::ios::out);
  if (!ofile) {
    Output::print(" *** Error(Mesh::writeToFile) I could not open ",filename);
    exit(-1);
  }

  Output::print(" writing mesh on ", filename);
  
  // header (.mesh)
  ofile << "MeshVersionFormatted 1" << endl;
  ofile << "Dimension" << endl;
  ofile << dimension << endl;
  ofile << "Vertices" << endl;
  ofile << vertex.size() << endl;

  // write nodes
  for (unsigned int i=0; i<vertex.size(); i++){
    ofile << vertex[i].x << " " << vertex[i].y << " ";
    if (dimension>2) {
      ofile << vertex[i].z << " ";
    }
    ofile << vertex[i].reference << endl;
  }

  ///@brief write edges only in dimension 2
  if (dimension==2){
    ofile << endl;
    ofile << "Edges" << endl;
    ofile << edge.size() << endl;
    for (unsigned int i=0; i<edge.size(); i++) {
      ofile << edge[i].nodes[0] << " " <<  edge[i].nodes[1]
	     << " " << edge[i].reference << endl;
    }
  }

  
  // write triangles
  ofile << endl;
  ofile << "Triangles" << endl;
  ofile <<  triangle.size() << endl;
  for (unsigned int i=0; i<triangle.size(); i++) {
    for (unsigned int k=0; k<3; k++) 
      ofile << triangle[i].nodes[k] << " ";
    ofile << triangle[i].reference << endl;
  }

  // write quadrilaterals
  ofile << endl;
  ofile << "Quadrilaterals" << endl;
  ofile <<  quad.size() << endl;
  for (unsigned int i=0; i<quad.size(); i++) {
    for (unsigned int k=0; k<4; k++) 
      ofile << quad[i].nodes[k] << " ";
    ofile << quad[i].reference << endl;
  }


  if (dimension==3){
    // write tetra
    ofile << endl;
    ofile << "Tetrahedra" << endl;
    ofile << tetra.size() << endl;
    for (unsigned int i=0; i<tetra.size(); i++) {
      for (unsigned int k=0; k<4; k++) 
	ofile << tetra[i].nodes[k] << " ";
      ofile << tetra[i].reference << endl;
    }

    // write hexa
    ofile << endl;
    ofile << "Hexahedra" << endl;
    ofile << hexa.size() << endl;
    for (unsigned int i=0; i<hexa.size(); i++) {
      for (unsigned int k=0; k<8; k++) 
	ofile << hexa[i].nodes[k] << " ";
      ofile << hexa[i].reference << endl;
    }
    
  }


}

// write the mesh onto a .(x)GEO file
// note: in 2D it needs the boundary description. If the boundary
// has not been initialized, the function returns an error
void Mesh::writeToGEO(std::string geoFilename)
{
  bool writeXgeo = false;
  std::ofstream geofile;
  geofile.open(geoFilename.c_str(),std::ios::out);
  if (!geofile) {
    Output::print(" *** Error(Mesh::writeToGEO) I could not open ",geoFilename);
    exit(-1);
  }

  // check errors
  if ((dimension==2)&(boundary.parts.size()==0)) {
    Output::print(" *** Error(Mesh::writeToGEO) I need a boundary description (PRM file) in order to write a 2D GEO file");
    exit(-1);
  }

    
  // write output (x)GEO file
  Output::print(" writing mesh on ", geoFilename);

  // check if input file is an extended geo file (.xGEO)
  unsigned int nn=0;
  while (geoFilename[nn] != 0) {++nn;}
  
  if (geoFilename[nn-4]=='x') {
    Output::print("  Mesh::writeToGEO: writing mesh in xGEO format (with physical references) ");
    writeXgeo = true;
  }
  unsigned int numberOfElements=0, maxNVertexPerElem=0;
  if (triangle.size()) {
    numberOfElements = triangle.size();
    maxNVertexPerElem = 3;
  }
  if (quad.size()) {
    numberOfElements += quad.size();
    maxNVertexPerElem = 4;
  }
  
  if ( (quad.size()==0)&&(triangle.size()==0) ) {
    Output::print(" Mesh:Write2GEO Error: I cannot write a mesh without elements");
    exit(-1);
  }

  unsigned int nBoundParts = boundary.parts.size();

  geofile << "Grid 2D #generated by Mesh::writeToGEO " << endl;
  geofile << "n_elem  n_vert ignored  max_n_vert_per_elem  n_bound_parts" << endl;
  geofile <<  numberOfElements << " " <<  vertex.size() << " "
	  <<  " 0 " << maxNVertexPerElem  << " " << nBoundParts << " " << endl;
  geofile << "DCORVG" << endl;

  // vertices
  /**
     @note in 2D, The vertex belonging to the boundary are written
     according to their local parametrization. I.e., we write
     compID+t 0
     where compID is the ID of the boundary component of the vertex, and t is its 
     local parametrization.
     Inner vertices are written simply using their coordinates.
   */
  double localParam=-1;
  int partID;
  for (unsigned int i=0; i<vertex.size(); i++) {

    // check if the vertex is on the boundary
    partID = boundary.isOnComponent(vertex[i].x,vertex[i].y,localParam);
    if (partID>=0) {
      Output::print("Vertex ",vertex[i].x," ",vertex[i].y," is on part ", partID+1, "; t = ", localParam, " --> ", localParam, " 0.0");
      geofile << localParam << " 0.000000000 " << endl;
      // note: a boundary vertex takes the reference of its boundary component+1
      vertex[i].reference = partID+1;
    } else {
      geofile << vertex[i].x << " " << vertex[i].y << endl;
      // inner nodes are marked with reference 0
      vertex[i].reference = 0;
    }  
  }
  
  geofile << "KVERT" << endl;
  // elements
  for (unsigned int i=0; i<triangle.size(); i++) {
    for (unsigned int k=0; k<3; k++) {
      geofile << triangle[i].nodes[k] << " " ;
    }
    // note: for mixed meshes, we weite a 0 as triangles 4th elements
    if (maxNVertexPerElem==4) {
      geofile << 0 << " ";
    }
    if (writeXgeo) {
      geofile << triangle[i].reference;
      }
    geofile << endl;
  }

  for (unsigned int i=0; i<quad.size(); i++) {
    for (unsigned int k=0; k<4; k++) { 
      geofile << quad[i].nodes[k] << " " ;
    }
    if (writeXgeo) {
      geofile << quad[i].reference;
    }
    geofile << endl;
  }
  // for elements
  

  geofile << "KNPR" << endl;
  for (unsigned int i=0; i<vertex.size(); i++) {

    geofile << vertex[i].reference << " ";
    
    if ( ((i+1) % 10)==0 ) geofile << endl;

  }
  
  geofile << endl;
  geofile << "KMM";
  geofile << endl;
  ///@todo what does this last row mean?
  geofile << "1  2  3  4" << endl;
  geofile.close();

  
}


void Mesh::stripSpace(std::string &str) 
{
  for (size_t i=0;i<str.length();i++)
    {
      if (str[i]==' ') 
	{
	  str.erase(i,1);
	  i--;
	}
    }
}


void Mesh::info()
{
  Output::print(" ---- MESH INFORMATION ---- ");
  Output::print(" -- Dimension: ",dimension);
  Output::print(" --  N. Nodes: ",vertex.size());
  Output::print(" --  N. Edges: ",edge.size());
  Output::print(" --  N. Triangles: ",triangle.size());
  Output::print(" --  N. Quandrilaterals: ",quad.size());
  Output::print(" --  N. Tetrahedra: ",tetra.size());
  Output::print(" --  N. Hexahedra: ",hexa.size());
  
}

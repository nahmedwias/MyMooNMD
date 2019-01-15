/** ************************************************************************ 
*
* @class     Mesh
* @brief     stores mesh arrays and allows for conversion
*     
*            Store the mesh using several lists
*            (vertices, tria, quad, etc.), which are
*            implemented as separate struct (in this file)
*            Moreover, it contains a Boundary object that
*            can store information about the boundary (i.e. a PRM file)
*            Note: the Boundary must be initialized in order to
*            write the geometry in .GEO format.
* 
* @author    Alfonso Caiazzo 
* @date      21.03.16
 ************************************************************************  */

#include <vector>
#include <array>
#include <string>
#include "Boundary.h"

#ifndef __MESH__
#define __MESH__

struct meshNode
{
  double x,y,z;
  int reference;
};

///@brief a structure for edges 
struct meshEdge
{
  std::array<int,2> nodes;
  int reference;
};

///@brief a structure for triangles
struct meshTriangle
{
  std::array<int,3> nodes;
  int reference;
};

///@brief a structure for quadrilaterals
struct meshQuad
{
  std::array<int,4> nodes;
  int reference;
};

///@brief a structure for tetrahedra
struct meshTetrahedron
{
  std::array<int,4> nodes;
  int reference;
};

///@brief a structure for hexahedra
struct meshHexahedron
{
  std::array<int,8> nodes;
  int reference;
};


/**  @brief mesh details */
class Mesh 
{
 protected:
  
 public:

  /**
   *   @brief dimension
   *   @attention this is not necessarily the geometrical dimension
   *  (a mesh file could be written in dimension 3 even if
   * the domain lies on a plane)
  */
  unsigned int dimension;
  
  /** @brief nodes */    
  std::vector<meshNode> vertex;

  /** @brief edges */    
  std::vector<meshEdge> edge;

  /** @brief triangles */    
  std::vector<meshTriangle> triangle;

  /** @brief quadrilateral */    
  std::vector<meshQuad> quad;

  /** 
      @brief flag for the case of two dimensional mesh with both
      triangles and quadrilaterals
  */    
  bool hasBothTriaAndQuads;
  
  /** @brief tetrahedron */    
  std::vector<meshTetrahedron> tetra;

  /** @brief hexahedron */    
  std::vector<meshHexahedron> hexa;

  ///@brief boundary handler class
  Boundary boundary;

  /**
     @brief number of boundary faces
  */
  int n_boundary_faces;
  
  /**
     @brief meshTrifaceHash[k] contains the id of the faces (a,b,c) 
     such that k=a+b+c

     meshTrifaceHash is a vector<int*> of size = 3*total n. of points
     -- meshTrifaceHash[hash][0] containts the number of faces with a+b+c = hash
     -- meshTrifaceHast[hash][i] contains the id of the i-th triangles s.t. a+b+c = hash

  */
  std::vector<int*> meshTrifaceHash;

  /**
     @brief map face to neighboring tetra
    
     faceToTetra[i][0] and faceToTetra[i][1] contains the the indices
     of two elements sharing face i
     faceToTetra[i][0]=faceToTetra[i][1]=-1 if the face is on the boundary
  */
  std::vector< std::vector<int> > faceToTetra;

  ///@brief boundaryFacesMarker: 0 = inner face, > 0 = boundary face
  std::vector<int> boundaryFacesMarker;
  /**
   *   @brief Emtpy constructor initialize an empty mesh,
   *   a single file initialize the mesh structures (.mesh)
   *   while passing also a boundary file initializes also
   *   the boundary description.
   */
  Mesh();
  explicit Mesh(std::string f);
  Mesh(std::string filename,std::string filenameBoundary);
  // Destructor
  ~Mesh(){};

  /**
     @brief read mesh from a file
     @note supported formats: .mesh
  */
  void readFromFile(std::string filename);

  ///@brief write mesh to a file .mesh
  void writeToMesh(std::string filename);

  /**
     @brief write mesh to a file .xGEO (extended ParMooN format)
     @param prmfile is an input prm file describing the boundary
     @attention the Boundary object must have been initialized
     (from a  prm file consistent with the geometry)
     @warning it works only in 2D at the moment
  */
  void writeToGEO(std::string filename);

  /**
     @brief initialize the Boundary class reading a PRM file 
     @warning it works only in 2D at the moment
   */
  void setBoundary(std::string PRM);


  /** 
      @brief count the faces with the same hash 
      
      This function is used to speed up the search of a nodes' face
      Once the "hash" array has been created, given three vertices (a,b,c), 
      we can search the ID of the corresponding face by looking up only the 
      faces with the same hash (a+b+c).
  */
  void hashTriFaces();

  /**
     @brief create a list mapping each face to the neighboring tetra

   */
  void createFaceToTetrahedraMap();

  /**
     @brief compute the number of boundary faces of the mesh

     Note: this functions uses the vector faceToTetra.
     If this has not been filled yet, it will be created by the function
   */
  void computeNumberOfBoundaryFaces();

  /**
     @brief create inner faces is these are not written in the mesh file
   */
  void createInnerFaces();

  /**
     @brief find the index of the triangular face with indices a,b,c
     
     Given the vertex indices (a,b,c), this function uses the hash a+b+c
     to find the corresponing face.
     Note: the meshTrifaceHash vector must have been created before. If not,
     the funcion creates it.
     Return -1 (and a warning) if no face is found.
  */
  int findTriFace(int a, int b, int c);
  
  ///@brief remove empry spaces in a string (general utility)
  void stripSpace(std::string &str);

  ///@brief display some info on screen
  void info();

  /**
     @brief number of nodes in the mesh
   */
  unsigned int nPoints() {
    return vertex.size();
  };
  
  /**
     @brief number of (inner) elements in the mesh
   */
  unsigned int nElements() {
    if (hexa.size()+tetra.size())
      return hexa.size()+tetra.size();
    else
      return triangle.size()+quad.size();
  };

};

#endif

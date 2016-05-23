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
   *   @brief Emtpy constructor initialize an empty mesh,
   *   a single file initialize the mesh structures (.mesh)
   *   while passing also a boundary file initializes also
   *   the boundary description.
   */
  Mesh();
  Mesh(std::string f);
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

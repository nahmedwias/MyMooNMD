/** ************************************************************************ 
*
* @class     Mesh
* @brief     stores Mesh-structures (Pts, Elements)
* @author    Alfonso Caiazzo 
* @date      21.03.16
* @History 
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
  /** @brief dimension */
  unsigned int dimension;
  
  /** @brief nodes */    
  std::vector<meshNode> vertex;

  
  /** @brief edges */    
  std::vector<meshEdge> edge;

  /** @brief triangles */    
  std::vector<meshTriangle> triangle;

  /** @brief quadrilateral */    
  std::vector<meshQuad> quad;

  /** @brief tetrahedron */    
  std::vector<meshTetrahedron> tetra;

  /** @brief hexahedron */    
  std::vector<meshHexahedron> hexa;

  ///@brief boundary handler class
  Boundary boundary;
  
  // Constructors
  Mesh();
  Mesh(std::string f);
  
  // Destructor
  ~Mesh(){};

  /**
     @brief read mesh from a file
     @note supported formats: .mesh
  */
  void readFromFile(std::string filename);

  ///@brief write mesh to a file .mesh
  void writeToMesh(std::string filename);

  /**@brief write mesh to a file .xGEO (extended ParMooN format)
     @param filename is the output geofile
     @param prmfile is an input prm file describing the boundary
     @attention the prm file must be consistent with the geometry
     @warning it works only for 2D at the moment
  **/
  void writeToGEO(std::string filename);

  /**
     @brief initialize the boundary reading a PRM file
   */
  void setBoundary(std::string PRM);
  
  ///@brief remove empry spaces in a string
  void stripSpace(std::string &str);

  ///@brief display some info on screen
  void info();
};

#endif

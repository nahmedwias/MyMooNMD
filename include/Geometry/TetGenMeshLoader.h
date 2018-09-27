/** ************************************************************************ 
*
* @class TTetGenMeshLoader  
* @date  05.05.00
* @brief handle smesh and mesh input (surface) meshes and creates a 3D mesh using tetgen
* @author Naveed Ahmed, Alfonso Caiazzo

************************************************************************  */
#ifndef _TTetGenMeshLoader
#define _TTetGenMeshLoader
/// using 1.4 version of the TetGen
//#define __TETGEN_14X__ 

#include <tetgen.h>

#include <Constants.h>
#include <MooNMD_Io.h>
#include <ParameterDatabase.h>

#include <iostream>
#include <memory>
#include <vector>

class TBoundPart;
class TBoundComp3D;
class TVertex;
class TBaseCell;
class TJoint;
class TDomain;

class TTetGenMeshLoader
{
  
public:

  /// @brief Name of the mesh file  
  std::string meshFileName;
  
  /// @brief counts how often face (a,b,c) with same hash=a+b+c occurs 
  std::vector<int*> meshTrifaceHash;
  
  /// @brief ParMooN Parameter database
  ParameterDatabase db;
  
  // three different options can be used    
  bool plc; // piecewise linear complex
  bool reconstruct;
  bool insertpoints;
  
  /// @brief TetGen in object
  tetgenio meshTetGenIn;

  /// @brief TetGen out object
  tetgenio meshTetGenOut;
#ifdef __TETGEN_14X__
    tetgenio meshTetAddIn;
#endif

  
  /** 
   * This constructor will only merge the current database 
   * with the ParMooN database.
   * 
   * @param filename copies the file name to meshFileName
   * @param parmoon_db merged with the default database
   */
  TTetGenMeshLoader(std::string filename, ParameterDatabase& parmoon_db);
  
  /// default Destructor 
  ~TTetGenMeshLoader()=default;
  
  /**
   * read mesh file in tetgen format and Generate 
   * the tetrahedral mesh.
   * 
   * This function will read the mesh file ".smesh, or .." 
   * and call tetgen to generate the volume mesh
   */
  void GenerateMesh();

  /// create mesh with options
  void Tetgen();


  ///@todo These members are not used now. To be removed.
  /// @brief count the faces with the same hash and store it in meshTrifaceHash vector
  void hashTriFaces();
  
  /** 
   * create adjacent tetrahedra to the faces of trifacelist
   * 
   * @return number of number of boudary components "N_BoundComp"
   * @attention this functions does two things
   */
  int nBoundaryComponents;
  int CreateAdjacency();

  ///find the face of triangle
  int findTriFace(int a, int b, int c);


};

#endif // _TTetGenMeshLoader

//===================================================================================================#include <TetGenMeshLoader_Test.h>

// =======================================================================
// @(#)Domain.C        1.23 05/05/00
// 
// Class:       TTetGenMeshLoader
// Purpose:     creates domain with tetgen from a .smesh file
//
// Author:      Andreas Hahn 26.04.2010
//              Naveed Ahmed, 10.06.2016  
// History: 
//
// =======================================================================
#ifndef _TTetGenMeshLoader
#define _TTetGenMeshLoader
/// using 1.4 version of the TetGen
#define __TETGEN_14X__ 

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
protected:
  /// Name of the mesh file  
  std::string meshFileName;
  
  /// list of all boundary components
  std::vector<TBoundComp3D*> meshBoundComps;
  
  /// list of all vertices
  std::vector<TVertex*> meshVertices;
  
  /// list of all joints
  std::vector<TJoint*> meshJoints;
  
  /// counts how often face with 
  /// same hash occurs 
  std::vector<int*> meshTrifaceHash;
  
  /// ParMooN Parameter database
  ParameterDatabase db;
  
  ///three different options can be used    
  bool plc; // piecewise linear complex
  bool reconstruct;
  bool insertpoints;
  
  /// TetGen in and out objects
  tetgenio meshTetGenIn;
  tetgenio meshTetGenOut;
#ifdef __TETGEN_14X__
    tetgenio meshTetAddIn;
#endif

public:  
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
   * the ParMooN mesh.
   * 
   * This function will read the mesh file ".smesh, or .." 
   * and modiefy the domain to ressemble the mesh, 
   * the boundary, and build the mesh which will be 
   * used within ParMooN.
   * 
   * @param[out] domain domain object where the mesh is loaded
   */
    void Generate(TDomain& domain);
  
protected:
  /// create mesh with options
  void Tetgen();
  
  /// build the boundary
  void buildBoundary(TBoundPart**& BdParts, int& N_BoundParts, 
            int& N_BoundComps, int*& StartBdCompID, int*& Interfaces);

  /// this creates the "meshTrifaceHash"
  void hashTriFaces();

  /** 
   * create adjacent tetrahedra to the faces of trifacelist
   * 
   * @return number of number of boudary components "N_BoundComp" 
   */
  int CreateAdjacency();

  ///find the face of triangle
  int findTriFace(int a, int b, int c);

  ///this function build the mesh which 
  void buildParMooNMDMesh(TBaseCell**& CellTree, int& N_RootCells,
                           double& StartX, double& StartY, double& StartZ, 
                           double& BoundX, double& BoundY, double& BoundZ);
  
  /// set vertices and bounds 
  void setVertices(double &StartX, double &StartY, double &StartZ,
                    double &BoundX, double &BoundY, double &BoundZ);
  
  /// allocate roor cells
  void allocRootCells(TBaseCell**& CellTree, int& N_RootCells);
 
  /// set the joints, boundary face params and map type
  void distributeJoints(TBaseCell **CellTree);  
};

#endif // _TTetGenMeshLoader

//===================================================================================================#include <TetGenMeshLoader_Test.h>

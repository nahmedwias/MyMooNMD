#include <BoundPart.h>
#include <BoundComp3D.h>
#include <BdPlane.h>
#include <BdCylinder.h>
#include <Vertex.h>
#include <MacroCell.h>
#include <Database.h>
#include <JointEqN.h>
#include <BoundFace.h>
#include <MooNMD_Io.h>
#include <TetGenMeshLoader.h>

ParameterDatabase get_default_param_database()
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  // merge default database of TetGen
  db.merge(ParameterDatabase::default_tetgen_database(), true);
  db.set_name("TetGen mesh loader default database ");
  return db;
}

TTetGenMeshLoader::TTetGenMeshLoader(std::string filename, 
                             ParameterDatabase &parmoon_db)
 : meshFileName(filename), db(get_default_param_database()),
   plc(false), reconstruct(false), insertpoints(false)  
{
  this->db.merge(parmoon_db);
  this->GenerateMesh();
}


void TTetGenMeshLoader::GenerateMesh()
{
  std::string remove_extension;
  // remove the extension b/c load_poly need the file without extension
  size_t lastindex = meshFileName.find_last_of("."); 
  std::string rawname = meshFileName.substr(0, lastindex); 

  // checking first the file extension
  //if(meshFileName.compare(rawname))
  if (meshFileName.substr(meshFileName.find_last_of(".")+1) == "smesh")
  {
    // set true for smesh files 
    plc = true;
    // load the .smesh file 
    // convert std::string to char* b/c load_poly(char*)
#ifdef __TETGEN_14X__
    meshTetGenIn.load_poly((char*)rawname.c_str());
#else
    ErrThrow("missing tetgen library, unable to generate a mesh ");
#endif
    this->Tetgen();
  }
  else
  {
    ErrThrow("TTetGenMeshLoader::GenerateMesh() - only .smesh files are supported");
    exit(1);
  }  

}

void TTetGenMeshLoader::Tetgen()
{
  std::string options;
  std::ostringstream osq, osa;
    
  double tetgen_quality = db["tetgen_quality"];
  double tetgen_volume = db["tetgen_volume"];
  
  osq << tetgen_quality;
  osa << tetgen_volume;
  
  if(plc)
  {
    options += "p"; // Tetrahedralizes a piecewise linear comple.
  }
  
  if(reconstruct)
  {
    options += "r"; // reconstruct mesh
  }
  
  if(tetgen_quality >= 1.0)
  {
    options += "q";
    options += osq.str();
  }
  
  if(tetgen_volume > 0.0 /*&& plc*/)
  {
    options += "a";
    options += osa.str();
  }
  
  if(insertpoints)
  {
    options += "i";
  }
  options += "z"; // Numbers all output items starting from zero.
  options += "f"; // Outputs faces (including non-boundary faces).
  options += "g"; // write output in .mesh format, do not write .node etc.
  
  if(db["tetgen_steiner"].is(1) && plc)
  {
    options += "Y"; // Suppresses boundary facets/segments splitting.
  }
  options += "C"; // Checks the consistency of the final mesh.
  if(db["tetgen_quiet"].is(0))
    options += "Q"; // Quiet
  
  if(db["tetgen_merge_colplaner"].is(0))
    options += "M"; // Does not merge coplanar facets.
  if((size_t)db["verbosity"]> 1)
  {
    options += "V"; // Verbose: Detailed information, more terminal output.
    options += "n"; // Outputs tetrahedra neighbors.
  }
  
  options += "A"; // assigns region attributes

  
  tetgenbehavior TetBeh;
  // call tetgen with switch options
#ifdef __TETGEN_14X__
  
  if(TetBeh.parse_commandline((char*) options.c_str()))
  {
    if(insertpoints)
    {
      Output::print("adding points");
      tetrahedralize(&TetBeh, &meshTetGenIn, &meshTetGenOut, &meshTetAddIn);
    }
    else
    {
      Output::print(" calling: tetgen -",options);
      tetrahedralize(&TetBeh, &meshTetGenIn, &meshTetGenOut);
    }
  }
  else
  {
    ErrThrow("parse_commandline failed !");    
  }
#else
  ErrThrow("missing tetgen library, unable to generate a mesh");
  //tetrahedralize((char*) options.c_str(), &meshTetGenIn,&meshTetGenOut);
#endif

  Output::print(" -- TTetGenMeshLoader::TetGen() completed --");
  
  // max 1 attribute (reference) per element allowed
  if(meshTetGenOut.numberoftetrahedronattributes <= 1)
  {
    Output::print<2>("number of tetrahedron attributes: ",
		     meshTetGenOut.numberoftetrahedronattributes);
  } else {
    Output::print("number of tetrahedron attributes: ",
		  meshTetGenOut.numberoftetrahedronattributes);
    ErrThrow("multiple attributes (>1) is not yet implemented");
  }

  // check consistence of the tetrahedral mesh
  if(meshTetGenOut.numberofcorners != 4)
  {
    ErrThrow("Wrong number of corners !");
  }

}

/**
   meshTrifaceHash is a vector<int*> of size = total n. of points
   - each element is initialized to a nullptr
   - for each triangle (a,b,c) in the meshTetGenOtut, we take hash = a+b+c
   - first, we coumpute the amount of triangles with the same hash
   (stored in the vector<int>BucketCount)
   - then, 
   -- meshTrifaceHash[hash][0] containts the number of faces with a+b+c = hash
   -- meshTrifaceHast[hash][i] contains the id of the i-th triangles s.t. a+b+c = hash

   @todo This function is not used now (rewritten in Mesh.C)
 */
void TTetGenMeshLoader::hashTriFaces()
{
  int a, b, c, hash;
  int n_points = meshTetGenOut.numberofpoints;
  meshTrifaceHash.resize(3*n_points);
  for(auto &e : meshTrifaceHash)
  {
    e=nullptr;
  }
  std::vector<int> BucketCount(3*n_points);
  
  for(int i=0;i<meshTetGenOut.numberoftrifaces;++i)
  {
    a = meshTetGenOut.trifacelist[3*i  ];
    b = meshTetGenOut.trifacelist[3*i+1];
    c = meshTetGenOut.trifacelist[3*i+2];

    hash = a+b+c;
    (BucketCount[hash])++;
    //cout << " BucketCount [ " << hash << "]=" << BucketCount[hash] << endl;
  }

  for(int i=0;i<meshTetGenOut.numberoftrifaces;++i)
  {
    a = meshTetGenOut.trifacelist[3*i  ];
    b = meshTetGenOut.trifacelist[3*i+1];
    c = meshTetGenOut.trifacelist[3*i+2];

    hash = a+b+c;

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
}

/// @todo This function is not used now (rewritten in Mesh.C)
int TTetGenMeshLoader::CreateAdjacency()
{
  Output::print("TTetGenMeshLoader::CreateAdjacency() start");
  int n_trifaces = meshTetGenOut.numberoftrifaces;
  int n_tetrahedra = meshTetGenOut.numberoftetrahedra;
  int triface;
  int N_BoundComp = n_trifaces;

  // create the meshTrifaceHash vector (used in findTriFace)
  hashTriFaces();
   
  int FaceVertex[][3] = { {0, 1, 2} ,{0, 3, 1}, {2, 1, 3}, {0, 2, 3} };

  if(meshTetGenOut.adjtetlist) delete [] meshTetGenOut.adjtetlist;

  // adjtetlist[2*triface] and adjtetlist[2*triface+1] contains
  // indices of two neighboring tetrahedra
  // (or they are both equal to -1 if triface is on the boundary)
  
  // initialize the adjacency vector (size 2*number of faces) with -1
  meshTetGenOut.adjtetlist = new int [2*n_trifaces];
  for(int i=0;i<2*n_trifaces;++i)
    meshTetGenOut.adjtetlist[i] = -1;

  for(int i=0;i<n_tetrahedra;++i)
  {
    Output::print<4>("TTetGenMeshLoader::CreateAdjacency() tetrahedra ",i);
  
    for(int face=0;face<4;++face)
    {
      // find a triface with the same vertices
      triface = findTriFace(meshTetGenOut.tetrahedronlist[4*i+FaceVertex[face][0]],
        meshTetGenOut.tetrahedronlist[4*i+FaceVertex[face][1]],
        meshTetGenOut.tetrahedronlist[4*i+FaceVertex[face][2]]);

      assert (triface != -1);

      // if the element 2*triface is still = -1, replace with the index
      // of the current tetrahedra
      if(meshTetGenOut.adjtetlist[2*triface] == -1)
        meshTetGenOut.adjtetlist[2*triface] = i;
      else
      {
	// if the face has been already associated to a tetrahedra,
	// set 2*triface + 1 equal to the current one
        assert(meshTetGenOut.adjtetlist[2*triface+1] == -1);
        meshTetGenOut.adjtetlist[2*triface+1] = i;
        --N_BoundComp;
      }
    }
  }
  return N_BoundComp;
}

/// @todo This function is not used now (rewritten in Mesh.C)
int TTetGenMeshLoader::findTriFace(int a, int b, int c)
{
  int hash = a+b+c;
  int count ;
  int *trifacelist = meshTetGenOut.trifacelist;
  int triface, vertex, found;

  assert(meshTrifaceHash[hash] != nullptr);
  count = meshTrifaceHash[hash][0];

  for(int i=1;i<=count;++i)
  {
    triface = meshTrifaceHash[hash][i];
    found = 0;
    for(int j=0;j<3;++j)
    {
      vertex = trifacelist[3*triface+j];

      if(a == vertex || b == vertex || c == vertex)
        ++found;
    }
    if(found == 3)
      return triface;
  }
  return -1;
}


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
}


void TTetGenMeshLoader::GenerateMesh()
{
  std::string remove_extension;
  // remove the extension b/c load_poly need the file without extension
  size_t lastindex = meshFileName.find_last_of("."); 
  std::string rawname = meshFileName.substr(0, lastindex); 

  // checking first the file extension
  if(meshFileName.compare(rawname))
  {

    // set true for smesh files 
    plc = true;
    // load the .smesh file 
    // convert std::string to char* b/c load_poly(char*)
    meshTetGenIn.load_poly((char*)rawname.c_str());
    Output::print("NO MARKERS");
    this->Tetgen();

    this->hashTriFaces();
    this->nBoundaryComponents = this->CreateAdjacency();
  }
  else
  {

    ErrThrow("currently only .smesh files are supported");
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
  
  if(db["tetgen_steiner"].is(1) && plc)
  {
    options += "Y"; // Suppresses boundary facets/segments splitting.
  }
  options += "C"; // Checks the consistency of the final mesh.
  if(db["tetgen_quite"].is(0))
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
      tetrahedralize(&TetBeh, &meshTetGenIn, &meshTetGenOut);      
  }
  else
  {
    ErrThrow("parse_commandline failed !");    
  }
#else
  tetrahedralize((char*) options.c_str(), &meshTetGenIn, &meshTetGenOut);
#endif

  Output::print("TetGen - Mesh generated");
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
  }

  for(int i=0;i<meshTetGenOut.numberoftrifaces;++i)
  {
    a = meshTetGenOut.trifacelist[3*i  ];
    b = meshTetGenOut.trifacelist[3*i+1];
    c = meshTetGenOut.trifacelist[3*i+2];

    hash = a+b+c;

    if(meshTrifaceHash.at(hash) == NULL)
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

int TTetGenMeshLoader::CreateAdjacency()
{
  int n_trifaces = meshTetGenOut.numberoftrifaces;
  int n_tetrahedra = meshTetGenOut.numberoftetrahedra;
  int triface;
  int N_BoundComp = n_trifaces;

  // create the meshTrifaceHash vector (used in findTriFace)
  hashTriFaces();
   
  int FaceVertex[][3] = { {0, 1, 2} ,{0, 3, 1}, {2, 1, 3}, {0, 2, 3} };

  if(meshTetGenOut.adjtetlist) delete [] meshTetGenOut.adjtetlist;

  meshTetGenOut.adjtetlist = new int [2*n_trifaces];
  for(int i=0;i<2*n_trifaces;++i)
    meshTetGenOut.adjtetlist[i] = -1;

  for(int i=0;i<n_tetrahedra;++i)
  {
    for(int face=0;face<4;++face)
    {
      triface = findTriFace(meshTetGenOut.tetrahedronlist[4*i+FaceVertex[face][0]],
        meshTetGenOut.tetrahedronlist[4*i+FaceVertex[face][1]],
        meshTetGenOut.tetrahedronlist[4*i+FaceVertex[face][2]]);

      assert (triface != -1);

      if(meshTetGenOut.adjtetlist[2*triface] == -1)
        meshTetGenOut.adjtetlist[2*triface] = i;
      else
      {
        assert(meshTetGenOut.adjtetlist[2*triface+1] == -1);
        meshTetGenOut.adjtetlist[2*triface+1] = i;
        --N_BoundComp;
      }
    }
  }
  return N_BoundComp;
}

int TTetGenMeshLoader::findTriFace(int a, int b, int c)
{
  int hash = a+b+c;
  int count ;
  int *trifacelist = meshTetGenOut.trifacelist;
  int triface, vertex, found;

  assert(meshTrifaceHash[hash] != NULL);
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


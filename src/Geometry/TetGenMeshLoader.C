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
  db.set_name("db mesh generator");
  return db;
}

TTetGenMeshLoader::TTetGenMeshLoader(std::string filename, 
                             ParameterDatabase &parmoon_db)
 : meshFileName(filename), db(get_default_param_database()),
   plc(false), reconstruct(false), insertpoints(false)  
{
  this->db.merge(parmoon_db);
}

void TTetGenMeshLoader::Generate(TDomain& domain)
{
  //char *suffix = NULL;
  //suffix = strrchr(const_cast<char*>(meshFileName.c_str()), '.');
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
    
    this->buildBoundary(domain.BdParts, domain.N_BoundParts, 
                           domain.N_BoundComps, domain.StartBdCompID, domain.Interfaces);

    this->buildParMooNMDMesh(domain.CellTree, domain.N_RootCells,
                      domain.StartX, domain.StartY, domain.StartZ,
                      domain.BoundX, domain.BoundY, domain.BoundZ);
  }
  else
  {
    ErrThrow("currently only .smesh files are supported");
  }  
  
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_LE]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_Finest]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_Between]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_OCAF]->SetParam(std::addressof(domain));  
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

void TTetGenMeshLoader::buildBoundary(TBoundPart**& BdParts, int& N_BoundParts, int& N_BoundComps, 
                                      int*& StartBdCompID, int*& Interfaces)
{
  int n_trifaces = meshTetGenOut.numberoftrifaces;
  int *trifaces = meshTetGenOut.trifacelist;
  int counter=0;
  double p[3], a[3], b[3], n[3];

  this->hashTriFaces();
  N_BoundComps = CreateAdjacency();
  Output::print("N_BoundComps: ", N_BoundComps);

  N_BoundParts = 1;
  BdParts = new TBoundPart* [N_BoundParts];
  StartBdCompID = new int [N_BoundParts+1];

  StartBdCompID[0] = 0;
  StartBdCompID[1] = N_BoundComps;

  meshBoundComps.resize(N_BoundComps); 
  // missuse of trifacemarkerlist :(
  if(meshTetGenOut.trifacemarkerlist == NULL)
    meshTetGenOut.trifacemarkerlist = new int [n_trifaces];

  for(int i=0;i<n_trifaces;++i)
  {
    if(meshTetGenOut.adjtetlist[2*i+1] == -1 ||
      meshTetGenOut.adjtetlist[2*i  ] == -1)
    {
      // boundary face
      meshBoundComps.at(counter) = new TBdPlane(counter);

      p[0] = meshTetGenOut.pointlist[3*trifaces[3*i  ]  ];
      a[0] = meshTetGenOut.pointlist[3*trifaces[3*i+1]  ] - p[0];
      b[0] = meshTetGenOut.pointlist[3*trifaces[3*i+2]  ] - p[0];
      p[1] = meshTetGenOut.pointlist[3*trifaces[3*i  ]+1];
      a[1] = meshTetGenOut.pointlist[3*trifaces[3*i+1]+1] - p[1];
      b[1] = meshTetGenOut.pointlist[3*trifaces[3*i+2]+1] - p[1];
      p[2] = meshTetGenOut.pointlist[3*trifaces[3*i  ]+2];
      a[2] = meshTetGenOut.pointlist[3*trifaces[3*i+1]+2] - p[2];
      b[2] = meshTetGenOut.pointlist[3*trifaces[3*i+2]+2] - p[2];
      // normalize vector a
      double fac = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
      a[0] /= fac;
      a[1] /= fac;
      a[2] /= fac;
      // normalize vector b
      fac = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
      b[0] /= fac;
      b[1] /= fac;
      b[2] /= fac;
      // cross product of a and b
      n[0] = a[1]*b[2] - a[2]*b[1];
      n[1] = a[2]*b[0] - a[0]*b[2];
      n[2] = a[0]*b[1] - a[1]*b[0];
      // normalize n
      fac = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      n[0] /= fac;
      n[1] /= fac;
      n[2] /= fac;

      ((TBdPlane*) meshBoundComps[counter])->SetParams(p[0], p[1], p[2],
                                                       a[0], a[1], a[2],
                                                       n[0], n[1], n[2]);
      ++counter;

      meshTetGenOut.trifacemarkerlist[i] = counter;
    }
    else // not a boundary face
      meshTetGenOut.trifacemarkerlist[i] = 0;
  }
  
  assert(counter==N_BoundComps);
  BdParts[0] = new TBoundPart(N_BoundComps);
  for(int i=0; i<N_BoundComps; i++)
  {
    BdParts[0]->SetBdComp(i, meshBoundComps[i]);
  }
}

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

void TTetGenMeshLoader::buildParMooNMDMesh(TBaseCell**& CellTree, int& N_RootCells, 
                         double& StartX, double& StartY, double& StartZ, 
                         double& BoundX, double& BoundY, double& BoundZ)
{
  this->setVertices(StartX, StartY, StartZ, BoundX, BoundY, BoundZ);
  this->allocRootCells(CellTree, N_RootCells);  
  this->distributeJoints(CellTree);
}


void TTetGenMeshLoader::setVertices(double& StartX, double& StartY, double& StartZ, 
                         double& BoundX, double& BoundY, double& BoundZ)
{
  double x, y, z;
  double xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0;

  meshVertices.resize(meshTetGenOut.numberofpoints); 
  //= new TVertex* [mTetOut.numberofpoints];

  for(int i=0;i<meshTetGenOut.numberofpoints;++i)
  {
    x = meshTetGenOut.pointlist[3*i  ];
    y = meshTetGenOut.pointlist[3*i+1];
    z = meshTetGenOut.pointlist[3*i+2];

    if(i > 0)
    {
      if(xmin > x) xmin = x;
      if(xmax < x) xmax = x;

      if(ymin > y) ymin = y;
      if(ymax < y) ymax = y;

      if(zmin > z) zmin = z;
      if(zmax < z) zmax = z;
    }
    else
    {
      xmin = xmax = x;
      ymin = ymax = y;
      zmin = zmax = z;
    }

    meshVertices.at(i)=new TVertex (x, y, z);
  }

  StartX = xmin;
  BoundX = xmax - xmin;

  StartY = ymin;
  BoundY = ymax - ymin;

  StartZ = zmin;
  BoundZ = zmax - zmin;

  Output::print("number of vertices: ", meshTetGenOut.numberofpoints);
}

void TTetGenMeshLoader::allocRootCells(TBaseCell**& CellTree, int& N_RootCells)
{
  TMacroCell *Cell;
  TVertex *Vertex;
  
  N_RootCells = meshTetGenOut.numberoftetrahedra;

  CellTree = new TBaseCell* [N_RootCells];

  Output::print("number of tetrahedron attributes: ", meshTetGenOut.numberoftetrahedronattributes);
  
  for(int i=0;i<N_RootCells;++i)
  {
    Cell = new TMacroCell (TDatabase::RefDescDB[Tetrahedron], 0);
    CellTree[i] = Cell;

    if(meshTetGenOut.numberofcorners != 4)
    {
      ErrThrow("Wrong number of corners !");
    }

    for(int j=0;j<4;++j)
    {
      Vertex = meshVertices.at(meshTetGenOut.tetrahedronlist[4*i+j]);

      Cell->SetVertex(j, Vertex);
    }
    
    if(meshTetGenOut.numberoftetrahedronattributes == 1)
    {
      Cell->SetPhase_ID((int) meshTetGenOut.tetrahedronattributelist[i]);
    }
    else
    {
      ErrThrow("multiple attributes is not yet implemented");
    }
  }

  Output::print("number of tetrahedra: ", meshTetGenOut.numberoftetrahedra);
}

void TTetGenMeshLoader::distributeJoints(TBaseCell** CellTree)
{
  // size of meshjoint
  meshJoints.resize(meshTetGenOut.numberoftrifaces);
  // search face which belongs to current bdComp
  for(int i=0;i<meshTetGenOut.numberoftrifaces;++i)
  {
    // find element that contain this face
    if(meshTetGenOut.trifacemarkerlist[i] == 0)// inner joints
    {
      int left = meshTetGenOut.adjtetlist[2*i];
      int right = meshTetGenOut.adjtetlist[2*i+1];

      meshJoints.at(i) = new TJointEqN (CellTree[left], CellTree[right]);
    }
    else // boundary joints
    {
      Output::print<5>("boundary joint");
      
      int bdcomp = meshTetGenOut.trifacemarkerlist[i] - 1;

     TBoundComp3D* BoundComp = meshBoundComps[bdcomp];

      meshJoints.at(i)=new TBoundFace (BoundComp);
    }
  }

  Output::print("number of joints: ", meshTetGenOut.numberoftrifaces);
  
  const int *TmpFV, *TmpLen;
  int MaxLen, triface;
  TShapeDesc *ShapeDesc;

  for(int i=0;i<meshTetGenOut.numberoftetrahedra;++i)
  {
    ShapeDesc = CellTree[i]->GetShapeDesc();
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);

    for(int j=0;j<4;++j)
    {
      triface = findTriFace(meshTetGenOut.tetrahedronlist[4*i+TmpFV[j*MaxLen  ]],
                            meshTetGenOut.tetrahedronlist[4*i+TmpFV[j*MaxLen+1]],
                            meshTetGenOut.tetrahedronlist[4*i+TmpFV[j*MaxLen+2]]);

      assert (triface != -1);

      CellTree[i]->SetJoint(j, meshJoints.at(triface));

      // correct params
      if(meshJoints.at(triface)->GetType() == BoundaryFace)
      {

        TBoundFace *BoundFace = (TBoundFace*) meshJoints.at(triface);
        TBoundComp3D *BoundComp = BoundFace->GetBoundComp();
        double x, y, z, t, s;
        double param1[4], param2[4];

        for(int k=0;k<TmpLen[j];++k)
        {
          CellTree[i]->GetVertex(TmpFV[MaxLen*j+k])->GetCoords(x,y,z);
          BoundComp->GetTSofXYZ(x,y,z, t, s);

          param1[k] = t;
          param2[k] = s;
        }
        BoundFace->SetParameters(param1, param2);
      }
    }
  }
  // set map type
  for(int i=0;i<meshTetGenOut.numberoftrifaces;++i)
  {
    meshJoints[i]->SetMapType();
  }
}


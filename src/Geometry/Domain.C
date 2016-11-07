
#ifdef _MPI
#  include "mpi.h"
#include <MeshPartition.h>
#endif

#include <Database.h>
#include <Domain.h>
#include <JointEqN.h>
#include <MacroCell.h>
#include <MooNMD_Io.h>
#include <BdSphere.h>

#include <Quadrangle.h>

#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>

#ifdef __2D__
  #include <IsoBoundEdge.h>
  #include <IsoInterfaceJoint.h>
#endif

#ifdef __3D__
  #include <BoundFace.h>
  #include <BdPlane.h>
  #include <BdWall.h>
  #include <BdCylinder.h>
  #include <IsoInterfaceJoint3D.h>
  #include <InterfaceJoint3D.h>
  #include <IsoBoundFace.h>
  
  #include <TetGenMeshLoader.h>
  #include <BDEdge3D.h>
  #include <InnerEdge.h>
  #include <IsoEdge3D.h>
#endif

#ifdef __MORTAR__
  #include <BoundPoint.h>
  #include <MortarBaseJoint.h>
  #include <MortarJoint.h>
#endif


#include <string.h>
#include <vector>

extern "C"
{
  #include "triangle.h"
  void triangulate(char *, struct triangulateio *, struct triangulateio *,
                   struct triangulateio *);
}

ParameterDatabase get_default_domain_parameters()
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();

  db.add("refinement_n_initial_steps", (size_t)0,
         "This is the number of refinement steps before any computation "
         "starts. Usually the mesh is uniformly refined. In a multigrid "
         "program, this determines the number of uniform refinements until the "
         "finest mesh.", (size_t)0, (size_t)20);

  db.add("refinement_max_n_adaptive_steps", (size_t) 0,
         "A maximum number of adaptive refinement steps"
         "which may be applied to this domain."
         "THIS IS UNUSED AT THE MOMENT!",
         (size_t) 0, size_t (10));
  
   db.add("boundary_file", "Default_UnitSquare",
        "This is a file describing the boundary of the computational domain. "
        "You probably want to adjust this to be the path to some file which "
        "typically has the extension 'PRM'. See the documentation for GEO and "
        "PRM files.",
        {"Default_UnitSquare", "Default_UnitCube"}
       );
   
   db.add("geo_file", "Default_UnitSquare",
        "This files describes the computational mesh. You probably want to "
        "adjust this to be the path to some file which typically has the "
        "extension 'GEO' or 'xGEO'. See the documentation for GEO and PRM "
        "files.",
        {"UnitSquare", "TwoTriangles", "Default_UnitCube_Hexa", 
         "Default_UnitCube_Tetra","Default_UnitSquare"});
   
   db.add("mesh_tetgen_file", std::string("Wuerfel"),
         "This files describes the computational mesh. Typically this files"
         " has the extension 'mesh, 'smesh', 'node' or 'poly'. "
         " currently only the smesh files are supported");

  return db;
}


// Constructor
TDomain::TDomain(const ParameterDatabase& param_db) :
  Interfaces(nullptr),db(get_default_domain_parameters())
{
  RefLevel = 0;
  Output::print<4>("domain is initialized");
  db.merge(param_db, false);
  
  std::string geoname = db["geo_file"];
  std::string boundname = db["boundary_file"];
  std::string smesh = db["mesh_tetgen_file"];
  
  if( (geoname.substr(geoname.find_last_of(".") + 1) == "GEO" )
      && (boundname.substr(boundname.find_last_of(".") + 1) == "PRM" )
      )
  {
    this->Init(boundname.c_str(), geoname.c_str());
    Output::print<4>("GEO and PRM files are selected");
  }
  else if( (geoname.substr(geoname.find_last_of(".") + 1) == "mesh" )
	   && (boundname.substr(boundname.find_last_of(".") + 1) == "PRM" )
	   )
  {
    Output::print<4>("mesh and PRM files are selected");
    this->InitFromMesh(boundname.c_str(), geoname.c_str());
    
  }
#ifdef __3D__
  else if( (geoname.substr(geoname.find_last_of(".") + 1) == "mesh" ) )
  {
    // note: in 3D, we allow domain initialization withtou PRM
    Output::print<>("Initialize (3D) Domain using mesh file ",geoname.c_str());
    this->InitFromMesh(boundname.c_str(), geoname.c_str());
  }
  else if(smesh.substr(smesh.find_last_of(".")+1) == "smesh")
  {
    // initialize a TTetGenMeshLoader using the surface mesh and generate volume mesh
    TTetGenMeshLoader tetgen(smesh, db);   
    // convert the mesh into ParMooN Mesh object
    Mesh m(tetgen.meshTetGenOut);
    // generate vertices, cells and joint
    GenerateFromMesh(m);
    // test: write the mesh on a file
    m.writeToMesh("test_parmoon.mesh"); 
  }
#endif
  else 
  {
    // default cases for the tests
    this->Init(boundname.c_str(), geoname.c_str());
  }
}

//TODO This domain constructor, which is also responsible for read-in of the
// old database, is to be reomved soon.
TDomain::TDomain(char *ParamFile, const ParameterDatabase& param_db) :
  Interfaces(nullptr),db(get_default_domain_parameters())
{
  RefLevel = 0;
  
  db.merge(param_db, true);

  // This will be removed as soon as we got entirely rid of the
  // global database.
  /** set variables' value in TDatabase using ParamFile */
  this->ReadParam(ParamFile);
  
  std::string geoname = db["geo_file"];
  std::string boundname = db["boundary_file"];
  std::string smesh = db["mesh_tetgen_file"];
  
  if( (geoname.substr(geoname.find_last_of(".") + 1) == "GEO" ) 
      && (boundname.substr(boundname.find_last_of(".") + 1) == "PRM" ) 
      )
  {
    this->Init(boundname.c_str(), geoname.c_str());
    Output::print<4>("GEO and PRM files are selected");
    
  }
  else if( (geoname.substr(geoname.find_last_of(".") + 1) == "mesh" )
	   && (boundname.substr(boundname.find_last_of(".") + 1) == "PRM" )
	   )
  {
    Output::print<>("Initialize Domain using mesh file ",geoname.c_str(),
		    " and PRM file ",boundname.c_str());
    this->InitFromMesh(boundname.c_str(), geoname.c_str());
  }
    
#ifdef __3D__
  else if( (geoname.substr(geoname.find_last_of(".") + 1) == "mesh" ) )
  {
    // note: in 3D, we allow domain initialization withtou PRM
    Output::print<>("Initialize (3D) Domain using mesh file ",geoname.c_str());
    this->InitFromMesh(boundname.c_str(), geoname.c_str());
  }
  
  else if( (smesh.substr(smesh.find_last_of(".")+1) == "smesh")||
	   (smesh.substr(smesh.find_last_of(".")+1) == "mesh") )
  {
    // intialize mesh using tetegen mesh loader
    Output::print<>("Initialize Domain using the surface mesh file ", smesh.c_str());

    // initialize a TTetGenMeshLoader using the surface mesh and generate volume mesh
    TTetGenMeshLoader tetgen(smesh, db);   
    // convert the mesh into ParMooN Mesh object
    Mesh m(tetgen.meshTetGenOut);
    // generate vertices, cells and joint
    GenerateFromMesh(m);
    // test: write the mesh on a file
    m.writeToMesh("test2.mesh");
}
#endif
  else 
  {
    // default cases for the tests
    this->Init(boundname.c_str(), geoname.c_str());
  }
}

TDomain::~TDomain()
{
  delete [] StartBdCompID;
  if (Interfaces)
    delete [] Interfaces;
  for(int i = 0; i < N_BoundParts; ++i)
    delete BdParts[i];
  delete [] BdParts;
  delete [] CellTree;
}

// Methods
int TDomain::GetBdPartID(int BdCompID)
{
  int i;

  for (i=N_BoundParts-1;i>0;i--)
    if (BdCompID >= StartBdCompID[i]) break;

  return i;
}

int TDomain::GetLocalBdCompID(int BdCompID)
{
  int i;

  for (i=N_BoundParts-1;i>0;i--)
    if (BdCompID >= StartBdCompID[i]) break;

  return BdCompID - StartBdCompID[i];
}

#ifdef __MORTAR__
int TDomain::SetSubGridIDs(IntFunct2D *TestFunc)
{
  int i, j, N_, ID, CurrID;
  TVertex *Vert;

  for (i=0;i<N_RootCells;i++)
  {
    ID = 0;
    N_ = CellTree[i]->GetRefDesc()->GetN_OrigVertices();
    for(j=0;j<N_;j++)
    {
      Vert = CellTree[i]->GetVertex(j);
      CurrID = TestFunc(Vert->GetX(), Vert->GetY());

      if (CurrID > ID) ID = CurrID;
    }

    ((TMacroCell *) CellTree[i])->SetSubGridID(ID);
  }
  
  return 0;
}

int TDomain::GenMortarStructs()
{
  int i, j, k, N_1, N_2, ID1, ID2;
  TBaseCell *MeBase, *Me, *Neighb;
  TJoint *Joint;
  //TVector<int> *MortarAux;

  //MortarAux = new TVector<int>(30,30);
  std::vector<int> MortarAux;

  N_MortarFaces = 0;

  for (i=0;i<N_RootCells;i++)
    CellTree[i]->SetClipBoard(0);

  for (i=0;i<N_RootCells;i++)
  {
    Me = CellTree[i];
    N_1 = Me->GetRefDesc()->GetN_OrigEdges();
    for(j=0;j<N_1;j++)
      if (Neighb = Me->GetJoint(j)->GetNeighbour(CellTree[i]))
        if ((ID1 = Me->GetSubGridID()) != (ID2 = Neighb->GetSubGridID()))
        {
          N_2 = Neighb->GetRefDesc()->GetN_OrigEdges();
          for (k=0;k<N_2;k++)
            if (Neighb->GetJoint(k) == Me->GetJoint(j)) break;

          if (!(Me->GetClipBoard() & 1 << j))
          {
            N_MortarFaces++;
            Me->SetClipBoard(Me->GetClipBoard() | 1 << j);
            Neighb->SetClipBoard(Neighb->GetClipBoard() | 1 << k);

            if (ID1 < ID2)
              //MortarAux->AddElement(i);
              MortarAux.push_back(i);
            else
              //MortarAux->AddElement(-i-1);
              MortarAux.push_back(-i-1);

            //MortarAux->AddElement(j);
            //MortarAux->AddElement(k);
            MortarAux.push_back(j);
            MortarAux.push_back(k);
          }
        }
  } 

  MortarFaces = new TMortarFace[N_MortarFaces];

  for (i=0;i<N_MortarFaces;i++)
  {
    //if (MortarAux->GetElement(3*i) >= 0)
    if (MortarAux.at(3*i) >= 0)
    {
      //MortarFaces[i].Cell = CellTree[MortarAux->GetElement(3*i)];
      //MortarFaces[i].LocFaceNumber[0] = MortarAux->GetElement(3*i + 1);
      //MortarFaces[i].LocFaceNumber[1] = MortarAux->GetElement(3*i + 2);
      MortarFaces[i].Cell = CellTree[MortarAux.at(3*i)];
      MortarFaces[i].LocFaceNumber[0] = MortarAux.at(3*i + 1);
      MortarFaces[i].LocFaceNumber[1] = MortarAux.at(3*i + 2);
    }
    else
    {
      //MeBase = CellTree[-MortarAux->GetElement(3*i) - 1];
      //MortarFaces[i].LocFaceNumber[0] = MortarAux->GetElement(3*i + 2);
      //MortarFaces[i].LocFaceNumber[1] = MortarAux->GetElement(3*i + 1);
      MeBase = CellTree[-MortarAux.at(3*i) - 1];
      MortarFaces[i].LocFaceNumber[0] = MortarAux.at(3*i + 2);
      MortarFaces[i].LocFaceNumber[1] = MortarAux.at(3*i + 1);
      MortarFaces[i].Cell = MeBase->GetJoint(MortarFaces[i].LocFaceNumber[1])->
                            GetNeighbour(MeBase);
    }
  }

  for (i=0;i<N_MortarFaces;i++)
  {
    k = MortarFaces[i].LocFaceNumber[0];
    Me = MortarFaces[i].Cell;
    Joint = Me->GetJoint(k);
    Neighb = Joint->GetNeighbour(Me);

    delete Joint;

    Joint = new TMortarBaseJoint(Me, Neighb);
    Me->SetJoint(k, Joint);
    Neighb->SetJoint(MortarFaces[i].LocFaceNumber[1], Joint);
  }
  
  //delete MortarAux;
  
  return 0;
}

TCollection *TDomain::GetMortarColl(Iterators it, int level)
{
  TCollection *coll;
  int i, n_cells = 0, lev, info;
  const int *TmpEV;
  TBaseCell **cells, *CurrCell, *LastCell, *NewCell;
  TJoint *Joint;

  BeginMFace = new int[N_MortarFaces + 1];
  BeginMFace[0] = 0;
  
  for (i=0;i<N_MortarFaces;i++)
  {
    lev = level + (i << 8);
    TDatabase::IteratorDB[it]->Init(lev);
    while (TDatabase::IteratorDB[it]->Next(info)) n_cells++;
    BeginMFace[i+1] = n_cells;
  }

  cells = new TBaseCell*[n_cells];

  for (n_cells=i=0;i<N_MortarFaces;i++)
  {
    lev = level + (i << 8);
    TDatabase::IteratorDB[it]->Init(lev);
    LastCell = NULL;
    while ((CurrCell = TDatabase::IteratorDB[it]->Next(info)))
    {
      NewCell = new TGridCell(TDatabase::RefDescDB[S_Line], RefLevel);

      CurrCell->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV);
      NewCell->SetVertex(0, CurrCell->GetVertex(TmpEV[2*info]));
      NewCell->SetVertex(1, CurrCell->GetVertex(TmpEV[2*info + 1]));

      if (LastCell)
      {
        Joint = new TJointEqN(NewCell, LastCell);
        NewCell->SetJoint(0, Joint);
        LastCell->SetJoint(1, Joint);
      }
      else
        NewCell->SetJoint(0, new TJointEqN(NewCell));

      cells[n_cells++] = LastCell = NewCell;
    }
    NewCell->SetJoint(1, new TJointEqN(NewCell));
  }

  coll = new TCollection(n_cells, cells);

  #ifdef  _MPI 
  coll->SetN_OwnCells(N_OwnCells);
  #endif

 return coll;
}

int TDomain::InitMortarJoints(Iterators it, int level, TCollection *coll)
{
  int i, j, ncell, Ncell, lev, info;
  const int *TmpEV;
  bool Left;
  double X0, Y0, X1, Y1, X, Y;
  TBaseCell *CurrCell, *Cell1D;

  for (ncell=i=0;i<N_MortarFaces;i++)
  {
    lev = level + (i << 8);
    TDatabase::IteratorDB[it]->Init(lev);
    while ((CurrCell = TDatabase::IteratorDB[it]->Next(info)))
      if (CurrCell->GetParent())
        ((TMortarJoint *) CurrCell->GetJoint(info))->SetMEdgeInColl(ncell++);
      else
        ((TMortarBaseJoint *) CurrCell->GetJoint(info))->
                                SetMEdgeInColl(ncell++);
  }

  for (i=0;i<N_MortarFaces;i++)
  {
    lev = - level - (i << 8);
    TDatabase::IteratorDB[it]->Init(lev);
    ncell = BeginMFace[i];
    Ncell = BeginMFace[i+1] - 1;
    while ((CurrCell = TDatabase::IteratorDB[it]->Next(info)))
    {
      CurrCell->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV);
      X0 = CurrCell->GetVertex(TmpEV[2*info])->GetX();
      Y0 = CurrCell->GetVertex(TmpEV[2*info])->GetY();
      X1 = CurrCell->GetVertex(TmpEV[2*info+1])->GetX();
      Y1 = CurrCell->GetVertex(TmpEV[2*info+1])->GetY();

      do
      {
        Cell1D = coll->GetCell(ncell);
        X = Cell1D->GetVertex(0)->GetX();
        Y = Cell1D->GetVertex(0)->GetY();
        Left = (bool) ((X1-X)*(X1-X0) + (Y1-Y)*(Y1-Y0) < 0);

        if (Left) ncell++;

      } while (ncell < Ncell && Left);

      if (ncell > Ncell)
        ncell--;
      else
      {
        Cell1D = coll->GetCell(ncell);
        X = Cell1D->GetVertex(0)->GetX();
        Y = Cell1D->GetVertex(0)->GetY();
        if ((X1-X)*(X1-X0) + (Y1-Y)*(Y1-Y0) > 0) ncell--;
      }

      if (CurrCell->GetParent())
        ((TMortarJoint *) CurrCell->GetJoint(info))->
                            SetMEdgeInColl(-ncell - 1);
      else
        ((TMortarBaseJoint *) CurrCell->GetJoint(info))->
                                SetMEdgeInColl(-ncell - 1);
    }
  }

  return 0;
}
#endif // __MORTAR__

#ifdef __2D__
//void SaveTri(struct triangulateio &outt);

int TDomain::GenInitGrid()
{
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  int *VertsPerComp;
  int i, j, k, N_, N_BdComps, N_InitVerts = 0;
  int CurrComp = 0;
  int CurrPoint = 0, CurrSeg = 0, BdStartPoint;
  TBoundComp2D *CurrCompPtr;
  double T, X, Y, T_a, T_b, *Coordinates;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;
  int *PartMarker, *Triangles, *PointNeighb, maxEpV = 0;
  int a, b, len1, len2, aux, part, comp, Neib[2], CurrNeib;
  TJoint *Joint;
  TVertex **NewVertices;
  struct triangulateio inn,  outt;
  bool AllowEdgeRef = (bool) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
  std::ostringstream opts;

  VertsPerComp = new int[N_BoundComps];

  for (i=0;i<N_BoundParts;i++)
  {
    N_BdComps = BdParts[i]->GetN_BdComps();
    for (j=0;j<N_BdComps;j++)
      VertsPerComp[CurrComp++] = BdParts[i]->GetBdComp(j)->GetN_InitVerts();
  }

  for (i=0;i<N_BoundComps;i++)
    N_InitVerts += VertsPerComp[i] - 1;

  inn.numberofsegments = N_InitVerts;
  inn.segmentlist = new int[2*N_InitVerts];
  inn.segmentmarkerlist = new int[N_InitVerts];
  
  for (i=0;i<N_BoundParts;i++)
    if (Interfaces[i] < 0) N_InitVerts++;

  inn.numberofpoints = N_InitVerts;             
  inn.numberofpointattributes = 0;
  inn.numberofholes = N_Holes;
  inn.numberofregions = N_Regions;
  inn.regionlist = PointInRegion;

  inn.pointlist = new double[2*N_InitVerts + 2];
  inn.pointmarkerlist = new int[N_InitVerts + 1];
  inn.holelist = PointInHole;

  CurrComp = 1;

  for (i=0;i<N_BoundParts;i++)
  {
    N_BdComps = BdParts[i]->GetN_BdComps();
    BdStartPoint = CurrPoint;

    for (j=0;j<N_BdComps;j++)
    {
      CurrCompPtr = (TBoundComp2D*)(BdParts[i]->GetBdComp(j));
      switch (CurrCompPtr->GetType())
      {
        case Line:
          inn.segmentlist[2*CurrSeg  ] = CurrPoint;
          inn.segmentlist[2*CurrSeg+1] = CurrPoint + 1;
          inn.segmentmarkerlist[CurrSeg++] = CurrComp;

          CurrCompPtr->GetXYofT(0.0, X, Y);
          inn.pointlist[2*CurrPoint  ] = X;
          inn.pointlist[2*CurrPoint+1] = Y;
          inn.pointmarkerlist[CurrPoint++] = CurrComp;

          CurrCompPtr->GetXYofT(1.0, X, Y);
          inn.pointlist[2*CurrPoint  ] = X;
          inn.pointlist[2*CurrPoint+1] = Y;
          inn.pointmarkerlist[CurrPoint] = CurrComp;
        break;

        case Polygon:
          inn.segmentlist[2*CurrSeg  ] = CurrPoint;
          inn.segmentlist[2*CurrSeg+1] = CurrPoint + 1;
          inn.segmentmarkerlist[CurrSeg++] = CurrComp + 100000;

          CurrCompPtr->GetXYofT(0.0, X, Y);
          inn.pointlist[2*CurrPoint  ] = X;
          inn.pointlist[2*CurrPoint+1] = Y;
          inn.pointmarkerlist[CurrPoint++] = CurrComp;

          CurrCompPtr->GetXYofT(1.0, X, Y);
          inn.pointlist[2*CurrPoint  ] = X;
          inn.pointlist[2*CurrPoint+1] = Y;
          inn.pointmarkerlist[CurrPoint] = CurrComp;
        break;

        case Circle:
        case Spline2D:
          BdParts[i]->GetBdComp(j)->GenInitVerts(inn.pointlist, CurrPoint,
                        inn.segmentlist, CurrSeg);

          N_ = VertsPerComp[CurrComp-1];
          for (k=0;k<N_;k++)
          {
            if (AllowEdgeRef)
            {
              if (k) 
                inn.segmentmarkerlist[CurrSeg++] = CurrComp;
            }
            else
            {
              if (k) 
                inn.segmentmarkerlist[CurrSeg++] = 100000 + CurrComp;
            }

            inn.pointmarkerlist[CurrPoint++] = CurrComp;
          }
          CurrPoint--;
        break;

      default:
         cerr << "Not a 1D  or 2D BDcomp" << endl;
             exit(-1);
        break;	
	
      }

      CurrComp++;
    }

    if (Interfaces[i] > 0)
      inn.segmentlist[2*CurrSeg-1] = BdStartPoint;
    else
      if (ABS(inn.pointlist[2*CurrPoint] -
              inn.pointlist[2*BdStartPoint]) < 1e-4 &&
          ABS(inn.pointlist[2*CurrPoint + 1] -
              inn.pointlist[2*BdStartPoint + 1]) < 1e-4)
      {
        inn.segmentlist[2*CurrSeg-1] = BdStartPoint;
        inn.numberofpoints--;
      }
      else
        CurrPoint++;
  }

  outt.pointlist = NULL;        
  outt.pointmarkerlist = NULL; 
  outt.trianglelist = NULL;    
  outt.neighborlist = NULL;
  outt.segmentlist = NULL;
  outt.segmentmarkerlist = NULL;
  outt.edgelist = NULL;        
  outt.edgemarkerlist = NULL;  
  outt.triangleattributelist = NULL;

  opts << "pq" << TDatabase::ParamDB->MESHGEN_REF_QUALITY;
  opts << "zneQA" << ends;
  triangulate((char *)(opts.str().c_str()), &inn, &outt, (struct triangulateio *) NULL);
  //SaveTri(outt);

  N_RootCells = outt.numberoftriangles;

  // allocate auxillary fields
  Coordinates = outt.pointlist;
  Triangles = outt.trianglelist;
  PartMarker = new int[outt.numberofpoints];

  N_ = outt.numberofpoints;
  for (i=0;i<N_;i++)
    if (outt.pointmarkerlist[i])
      PartMarker[i] = GetBdPartID(outt.pointmarkerlist[i]-1);
    else
      PartMarker[i] = 0;

  // generate all vertices
  N_ = outt.numberofpoints;
  NewVertices = new TVertex*[N_];

  for (i=0;i<N_;i++)
  {
    X = Coordinates[2*i];
    Y = Coordinates[2*i+1];

    if ((CurrComp = outt.pointmarkerlist[i]) && AllowEdgeRef)
    {
      CurrComp--;
      part = GetBdPartID(CurrComp);
      comp = GetLocalBdCompID(CurrComp);
      if (BdParts[part]->GetBdComp(comp)->GetType() == Circle)
      {
        BdParts[part]->GetBdComp(comp)->GetTofXY(X, Y, T);
        BdParts[part]->GetBdComp(comp)->GetXYofT(T, X, Y);
      }
    }
        
    NewVertices[i] = new TVertex(X, Y);

    if (X > Xmax) Xmax = X;
    if (X < Xmin) Xmin = X;
    if (Y > Ymax) Ymax = Y;
    if (Y < Ymin) Ymin = Y;
  }

  // set bounding box
  StartX = Xmin;
  StartY = Ymin;
  BoundX = Xmax - Xmin;
  BoundY = Ymax - Ymin;

  // generate cells
  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle],
                                 RefLevel);

    CellTree[i]->SetVertex(0, NewVertices[outt.trianglelist[3*i    ]]);
    CellTree[i]->SetVertex(1, NewVertices[outt.trianglelist[3*i + 1]]);
    CellTree[i]->SetVertex(2, NewVertices[outt.trianglelist[3*i + 2]]);

    if(N_Regions)
    {
      ((TMacroCell *) CellTree[i])->SetSubGridID(
          (int)(outt.triangleattributelist[i]+0.05));
      CellTree[i]->SetReference_ID((int)(outt.triangleattributelist[i]+1.05));
    }
    else
      ((TMacroCell *) CellTree[i])->SetSubGridID(0);
  }

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  #ifdef __MORTAR__
    TDatabase::IteratorDB[It_Mortar1]->SetParam(this);
    TDatabase::IteratorDB[It_Mortar2]->SetParam(this);
  #endif 

  // search neighbours
  N_ = outt.numberofpoints;
  PointNeighb = new int[N_];

  memset(PointNeighb, 0, N_ * SizeOfInt);

  for (i=0;i<3*N_RootCells;i++)
    PointNeighb[Triangles[i]]++;

  for (i=0;i<N_;i++)
    if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];

  delete PointNeighb;
  PointNeighb = new int[++maxEpV * N_];

  memset(PointNeighb, 0, maxEpV * N_ * SizeOfInt);
  
  // first colomn contains the number of following elements
  for (i=0;i<3*N_RootCells;i++)
  {
    j = Triangles[i]*maxEpV;
    PointNeighb[j]++;
    PointNeighb[j + PointNeighb[j]] = i / 3;
  }

  // generate edges
  N_ = outt.numberofedges;
  for (i=0;i<N_;i++)
  {
    a = outt.edgelist[2*i];
    b = outt.edgelist[2*i+1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;

    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];

    for (j=1;j<=len1;j++)
    {
      aux = PointNeighb[a*maxEpV + j];

      for (k=1;k<=len2;k++)
        if (aux == PointNeighb[b*maxEpV + k])
        {
          Neib[CurrNeib++] = aux;
          break;
        }

      if (CurrNeib == 2) break;
    }

    if (outt.edgemarkerlist[i])
    {
      CurrComp = outt.edgemarkerlist[i] - 1;
      if (CurrComp >= 100000) CurrComp -= 100000;

      part = GetBdPartID(CurrComp);
      comp = GetLocalBdCompID(CurrComp);


     #ifdef _MPI
     if(rank==0)
     #endif
      {
      if (BdParts[part]->GetBdComp(comp)->GetTofXY(
            NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          BdParts[part]->GetBdComp(comp)->GetTofXY(
            NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
        cerr << "Warning in GenInitGrid: could not get parameter values!"
             << endl;
       }

      if (BdParts[part]->GetBdComp(comp)->GetType() != Line)
        if (ABS(T_a) < 1e-4 || ABS(T_a) > 0.9999 ||
            ABS(T_b) < 1e-4 || ABS(T_b) > 0.9999)
        {
          X = (NewVertices[a]->GetX() + NewVertices[b]->GetX()) / 2;
          Y = (NewVertices[a]->GetY() + NewVertices[b]->GetY()) / 2;

          BdParts[part]->GetBdComp(comp)->GetTofXY(X, Y, T);

          if ((T_a - T)*(T - T_b) < 0)
          {
            if (ABS(T_a) < 1e-4) T_a = 1.0;
            if (ABS(T_b) < 1e-4) T_b = 1.0;
            if (ABS(T_a) > 0.9999) T_a = 0.0;
            if (ABS(T_b) > 0.9999) T_b = 0.0;
          }
        }

      if (CurrNeib == 2)
        if(BdParts[part]->GetBdComp(comp)->IsFreeBoundary())
          Joint = new TIsoInterfaceJoint(BdParts[part]->GetBdComp(comp), 
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
        else
          Joint = new TInterfaceJoint(BdParts[part]->GetBdComp(comp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
      else
        if(BdParts[part]->GetBdComp(comp)->IsFreeBoundary())
          Joint = new TIsoBoundEdge(BdParts[part]->GetBdComp(comp), T_a, T_b);
        else
          Joint = new TBoundEdge(BdParts[part]->GetBdComp(comp), T_a, T_b);
    }
    else
    {
      if (CurrNeib != 2)
        cerr << "Warning in GenInitGrid: not enough neighbours!" << endl;

      Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
    }

    for (j=0;j<3;j++)
      if (Triangles[3*Neib[0]+j] == a) break;
        
    for (k=0;k<3;k++)
      if (Triangles[3*Neib[0]+k] == b) break;

    k = k*10 + j;

    switch (k)
    {
      case  1:
      case 10:
        aux = 0;
        break;
      case 12:
      case 21:
        aux = 1;
        break;
      case  2:
      case 20:
        aux = 2;
        break;
    }
    
    CellTree[Neib[0]]->SetJoint(aux, Joint);
    
    if (Neib[1] != -1)
    {
      for (j=0;j<3;j++)
        if (Triangles[3*Neib[1]+j] == a) break;
        
      for (k=0;k<3;k++)
        if (Triangles[3*Neib[1]+k] == b) break;

      k = k*10 + j;

      switch (k)
      {
        case  1:
        case 10:
          aux = 0;
          break;
        case 12:
        case 21:
          aux = 1;
          break;
        case  2:
        case 20:
          aux = 2;
          break;
      }
    
      CellTree[Neib[1]]->SetJoint(aux, Joint);
    }

    if (Joint->GetType() == InterfaceJoint ||
        Joint->GetType() == IsoInterfaceJoint)
      ((TInterfaceJoint *) Joint)->CheckOrientation();
  }
  
  delete NewVertices;
  delete PointNeighb;

  return 0;
}


void TDomain::Init(const char *PRM, const char *GEO)
{
  int Flag;

  if(PRM == nullptr || GEO == nullptr)
  {
    ErrThrow("No boundary description or initial mesh specified");
  }

  //start with treatment of the boundary description
  if (!strcmp(PRM, "Default_UnitSquare"))
  {//catch the only implemented default case - boundary of the unit square
    initializeDefaultUnitSquareBdry();
  }
  else
  {
    // non-default: read in from file
    //make an input file string from the file "PRM"
    std::ifstream bdryStream(PRM);
    if (!bdryStream)
    {
      ErrThrow("cannot open PRM file");
    }

    //do the actual read in
    ReadBdParam(bdryStream, Flag);

    //close the stream
    bdryStream.close();
  }

  // continue with initial mesh
  if (!strcmp(GEO, "InitGrid"))
  {
    GenInitGrid();
  }
  else if (!strcmp(GEO, "TwoTriangles"))
  {
    TwoTriangles();
  }
  else
  if (!strcmp(GEO, "TwoTrianglesRef"))
  {
    TwoTrianglesRef();
  }
  else if (!strcmp(GEO, "UnitSquare"))
  {
    UnitSquare();
  }
  else if (!strcmp(GEO, "UnitSquareRef"))
  {
    UnitSquareRef();
  }
  else if (!strcmp(GEO, "SquareInSquare"))
  {
    SquareInSquare();
  }
  else if (!strcmp(GEO, "UnitSquare_US22"))
  {
    UnitSquare_US22();
  }
  else if (!strcmp(GEO, "SquareInSquareRef"))
  {
    SquareInSquareRef();
  }
  else if (!strcmp(GEO, "PeriodicSquares"))
  {
    PeriodicSquares();
  }
  else if (!strcmp(GEO, "PeriodicSquaresLarge"))
  {
    PeriodicSquaresLarge();
  }
  else if (!strcmp(GEO, "PeriodicTrianglesLarge"))
  {
    PeriodicTrianglesLarge();
  }
  else if (!strcmp(GEO, "PeriodicRectangle_2_4"))
  {
    PeriodicRectangle_2_4();
  }
  else
  {

    // error message: .GEO file are no longer supported. They can be used to initialize a domain
    // only if the purpose is to covert the .GEO into a .mesh

    if (!strcmp(db["mesh_file"],"__nofile__"))
    {
      Output::print("** **************************************************************** **");
      Output::print("** ERROR: the (x)GEO files are no longer supported in ParMooN (2D): **");
      Output::print("** The file format .mesh should be used instead                     **");
      Output::print("**                                                                  **");
      Output::print("** If you use a standard geometry input (e.g. UnitSquare.*)         **");
      Output::print("** you can find the corresponding .mesh in the src/data/ directory. **");
      Output::print("** Otherwise, you can convert your .GEO file to a .mesh, using the  **");
      Output::print("** same PRM file for the boundary description, using the geo2mes2d. **");
      Output::print("**                                                                  **");
      Output::print("** To do this:                                                      **");
      Output::print("** (1) compile the conversion progam with 'make geo2mesh2d'         **");
      Output::print("** (2) create an input.dat with                                     **");
      Output::print("**     boundary_file: [your PRM file]                               **");
      Output::print("**     geo_file: [your GEO or xGEO file]                            **");
      Output::print("**     mesh_file: [the new mesh file]                               **");
      Output::print("**  (an example of suitable input file is available in src/data/    **");
      Output::print("** (3) run '[path/]geo2mesh2d input.dat'                            **");
      Output::print("** (4) replace 'geo_file: [the new mesh file]' in your input file   **");
      Output::print("**                                                                  **");
      Output::print("** Note: using .GEO and .PRM files is only allowed with the         **");
      Output::print("** purpose of converting .GEO -> .mesh                              **");      
      Output::print("** **************************************************************** **");
      exit(1);
    }      
    else 
    {
      Output::print("** Converting .GEO to .mesh **");    
    }
	
    // "GEO" interpreted as file path to GEO-file.
    // check if it actually is an .xGEO-file
    bool isxGEO = isExtendedGEO(GEO);

    // make an input file string from the file "GEO"
    std::ifstream bdryStream(GEO);
    if (!bdryStream)
    {
      ErrThrow("cannot open GEO file");
    }

    // and call the read-in method
    ReadGeo(bdryStream,isxGEO);
  }
}



#else // 3D
void TDomain::Init(const char *PRM, const char *GEO)
{
  int IsSandwich = 0;

  if(PRM == nullptr || GEO == nullptr)
  {
    ErrThrow("No boundary description or initial mesh specified");
  }

  // start with read in of boundary description

  if (!strcmp(PRM, "Default_UnitCube"))
  {
    // one implemented default case: unit cube
    IsSandwich = initializeDefaultCubeBdry();
  }
  else
  {
    // "PRM" interpreted as file path to PRM-file.
    // make an input file string from the file "PRM"
    std::ifstream bdryStream(PRM);
    if (!bdryStream)
    {
      ErrThrow("cannot open PRM file");
    }
    // otherwise: try to read in given .PRM file
    ReadBdParam(bdryStream, IsSandwich);
  }

  if(!IsSandwich)
  {
    if (!strcmp(GEO, "TestGrid3D"))
    {
      TestGrid3D();
    }
    else if (!strcmp(GEO, "Default_UnitCube_Hexa"))
    {
      initialize_cube_hexa_mesh();
    }
    else if (!strcmp(GEO, "Default_UnitCube_Tetra"))
    {
      initialize_cube_tetra_mesh();
    }
    else
    {
      // "GEO" interpreted as file path to GEO-file.
      // check if it actually is an .xGEO-file
      bool isxGEO = isExtendedGEO(GEO);
      // make an input file string from the file "GEO"
      std::ifstream bdryStream(GEO);
      if (!bdryStream)
      {
        ErrThrow("cannot open GEO file");
      }
      ReadGeo( bdryStream, isxGEO );
    }
  }
  else
  {
    // make an input file string from the file "GEO"
    std::ifstream bdryStream(GEO);
    if (!bdryStream)
    {
      ErrThrow("cannot open GEO file");
    }
    // then read in sandwich geo.
    ReadSandwichGeo(bdryStream);
  }
}
#endif // __2D__


// initialize a domain from a mesh and a boundary(PRM) file
void TDomain::InitFromMesh(std::string PRM, std::string MESHFILE)
{
#ifdef __2D__
  Output::print("TDomain:: InitFromMesh using ", PRM, " and ", MESHFILE);
  //make an input file string from the file "PRM"
  std::ifstream bdryStream(PRM);
  if (!bdryStream)
    {
      Output::print(" ** Error(TDomain::Init) cannot open PRM file ", PRM);
      exit(-1);
    }
  //do the actual read in
  int Flag;
  ReadBdParam(bdryStream, Flag);
  //close the stream
  bdryStream.close();

  // read mesh
  Mesh m(MESHFILE);
  m.setBoundary(PRM);
  unsigned int numberOfElements = m.triangle.size() + m.quad.size();
  unsigned int maxNVertexPerElem = 3;
  // make the ParMooN-grid
  if (m.quad.size()) {
    maxNVertexPerElem = 4;
  }
  if (numberOfElements==0) {
    Output::print(" ** Error(Domain::Init) the mesh has no elements");
    exit(-1);
  }

  // vertices data
  double *DCORVG;
  int *KNPR;
  unsigned int n_vertices = m.vertex.size();
  KNPR = new int[n_vertices];
  DCORVG =  new double[2*n_vertices];

  // fill the DCORVG array with GEO-like coordinates of two-dimensional points
  for (unsigned int i=0; i<n_vertices; i++) {

    // check if the vertex is on the boundary
    double localParam;
    int partID = m.boundary.isOnComponent(m.vertex[i].x,m.vertex[i].y,localParam);
    if (partID>=0) {
      DCORVG[2*i] = localParam;
      DCORVG[2*i+1] = 0.;
      KNPR[i] = partID+1;
    } else {
      DCORVG[2*i] = m.vertex[i].x;
      DCORVG[2*i+1] = m.vertex[i].y;
      KNPR[i] = 0;
    }  
  }
  int *KVERT,*ELEMSREF;
  KVERT = new int[maxNVertexPerElem * numberOfElements];
  ELEMSREF = new int[numberOfElements];

  // store triangles (+ 0 when using a mixed mesh)
  for (unsigned int i=0;i<m.triangle.size();i++)
  {
    for (unsigned int  j=0; j<3; j++) 
      KVERT[maxNVertexPerElem*i + j] = m.triangle[i].nodes[j];
    
    if (maxNVertexPerElem==4) KVERT[maxNVertexPerElem*i + 3] = 0;
    
    ELEMSREF[i] = m.triangle[i].reference;
    
  }

  // store quadrilaterals
  for (unsigned int i=0;i<m.quad.size();i++)
  {
    for (unsigned int  j=0; j<4; j++) 
      KVERT[maxNVertexPerElem* (m.triangle.size() + i) + j] = m.quad[i].nodes[j];

    ELEMSREF[m.triangle.size()+i] = m.quad[i].reference;
    
  }
  // N_RootCells is an internal Domain variable used in other functions
  N_RootCells = numberOfElements;
  
  MakeGrid(DCORVG,KVERT,KNPR,ELEMSREF,m.vertex.size(),maxNVertexPerElem);
  delete [] DCORVG;
  delete [] KVERT;
  delete [] KNPR;
  delete [] ELEMSREF;
#else

  int IsSandwich = 0;
  std::ifstream bdryStream(PRM);
  if (bdryStream)
  {
    // read .PRM if available (for example for sandwich grid)
    ReadBdParam(bdryStream, IsSandwich);
    if (IsSandwich)
    {
      Output::print(" WARNING: Sandwich grid not yet implemented with .mesh input file");
      Output::print("          Use a .GEO instead.");

      // =======================================
      // read the .GEO file
      // note: this is the same as the above 2D code
      
      // read mesh
      Mesh m(MESHFILE);
      m.setBoundary(PRM);
      unsigned int numberOfElements = m.triangle.size() + m.quad.size();
      unsigned int maxNVertexPerElem = 3;
      // make the ParMooN-grid
      if (m.quad.size()) {
	maxNVertexPerElem = 4;
      }
      if (numberOfElements==0) {
	Output::print(" ** Error(Domain::Init) the mesh has no elements");
	exit(-1);
      }
      
      // vertices data
      double *DCORVG;
      int *KNPR;
      unsigned int n_vertices = m.vertex.size();
      KNPR = new int[n_vertices];
      DCORVG =  new double[2*n_vertices];
      
      // fill the DCORVG array with GEO-like coordinates of two-dimensional points
      for (unsigned int i=0; i<n_vertices; i++) {
	
	// check if the vertex is on the boundary
	double localParam;
	int partID = m.boundary.isOnComponent(m.vertex[i].x,m.vertex[i].y,localParam);
	if (partID>=0) {
	  DCORVG[2*i] = localParam;
	  DCORVG[2*i+1] = 0.;
	  KNPR[i] = partID+1;
	} else {
	  DCORVG[2*i] = m.vertex[i].x;
	  DCORVG[2*i+1] = m.vertex[i].y;
	  KNPR[i] = 0;
	}  
      }
      int *KVERT,*ELEMSREF;
      KVERT = new int[maxNVertexPerElem * numberOfElements];
      ELEMSREF = new int[numberOfElements];
      
      // store triangles (+ 0 when using a mixed mesh)
      for (unsigned int i=0;i<m.triangle.size();i++)
	{
	  for (unsigned int  j=0; j<3; j++) 
	    KVERT[maxNVertexPerElem*i + j] = m.triangle[i].nodes[j];
	  
	  if (maxNVertexPerElem==4) KVERT[maxNVertexPerElem*i + 3] = 0;
	  
	  ELEMSREF[i] = m.triangle[i].reference;
	  
	}
      
      // store quadrilaterals
      for (unsigned int i=0;i<m.quad.size();i++)
	{
	  for (unsigned int  j=0; j<4; j++) 
	    KVERT[maxNVertexPerElem* (m.triangle.size() + i) + j] = m.quad[i].nodes[j];
	  
	  ELEMSREF[m.triangle.size()+i] = m.quad[i].reference;
	  
	}
      // N_RootCells is an internal Domain variable used in other functions
      N_RootCells = numberOfElements;
      // =======================================

      // =======================================
      // parameters specific to the sandwich case
      double DriftX = TDatabase::ParamDB->DRIFT_X;
      double DriftY = TDatabase::ParamDB->DRIFT_Y;
      double DriftZ = TDatabase::ParamDB->DRIFT_Z;
      int N_Layers = TDatabase::ParamDB->N_CELL_LAYERS+1;
      
      double *Lambda = new double[N_Layers];
      for(int i=0;i<N_Layers;i++)
      {
	Lambda[i] = i * (1.0/(N_Layers-1));
      }
      // =======================================
      
      MakeSandwichGrid(DCORVG, KVERT, KNPR,  m.vertex.size(),maxNVertexPerElem,
		       DriftX, DriftY, DriftZ, N_Layers, Lambda);
      //exit(1);
    }
    else
    {
      Output::print("** ERROR: I can only initialize 3D domain using PRM for sandwich grid");
      Output::print("** ERROR: check your input files." );
      exit(1);
    }
    
  } else {

    // otherwise: proceed without a PRM file
    // read mesh
    Mesh m(MESHFILE);
    GenerateFromMesh(m);
  }

#endif // __2D__
}

int TDomain::PS(const char *name, Iterators iterator, int arg)
{
  int BX, BY, info, i;
  double scale;
  std::ofstream dat(name);
  TBaseCell *CurrCell;

  if (!dat)
  {
    cerr << "cannot open '" << name << "' for output" << endl;
    return -1;
  }

  cout << "Generating postscript file " << name << endl;
  
  // scale = 5350 / BoundX;
  // if (7820 / BoundY < scale) scale = 7820 / BoundY;
  scale = 535 / BoundX;
  if (782 / BoundY < scale) scale = 782 / BoundY;

  BX = (int) (BoundX * scale + .5);
  BY = (int) (BoundY * scale + .5);

  dat << "%!PS-Adobe-2.0" << endl;
  dat << "%%Creator: MooN_MD (Volker Behns)" << endl;
  dat << "%%DocumentFonts: Helvetica" << endl;
  // dat << "%%BoundingBox: 300 300 " << 300+BX << " " << 300+BY << endl;
  dat << "%%BoundingBox: 25 25 " << 35+BX << " " << 35+BY << endl;
  dat << "%%Pages: 1" << endl;
  dat << "%%EndComments" << endl;
  dat << "%%EndProlog" << endl;
  dat << "%%Page: 1 1" << endl;
  dat << "/Helvetica findfont 14 scalefont setfont" << endl;
  dat << "/Cshow { dup stringwidth pop 2 div neg 0 rmoveto show } def" << endl;
  dat << "/M { moveto } def" << endl;
  dat << "/L { lineto } def" << endl;
  //dat << "0.10 0.10 scale" << endl;
  //dat << "10.0 setlinewidth" << endl;
  dat << "0.5 setlinewidth" << endl;

  TDatabase::IteratorDB[iterator]->Init(arg);

  // loop over all cells
  i = 0;
  while((CurrCell = TDatabase::IteratorDB[iterator]->Next(info)))
  {
    CurrCell->SetClipBoard(i);
    CurrCell->PS(dat, scale, StartX, StartY);
    i++;
  }

  dat << "stroke" << endl;
  dat << "showpage" << endl;
  dat << "%%Trailer" << endl;
  dat << "%%Pages: 1" << endl;

  return 0;
}

int TDomain::PS(const char *name, TCollection *Coll)
{
  int BX, BY;
  double scale;
  std::ofstream dat(name);
  TBaseCell *CurrCell;
  int i, N_;

  if (!dat)
  {
    cerr << "cannot open '" << name << "' for output" << endl;
    return -1;
  }

  Output::info("PS","Generating postscript file ", name);

  N_ = Coll->GetN_Cells();
  
  // scale = 5350 / BoundX;
  // if (7820 / BoundY < scale) scale = 7820 / BoundY;
  scale = 535 / BoundX;
  if (782 / BoundY < scale) scale = 782 / BoundY;

  BX = (int) (BoundX * scale + .5);
  BY = (int) (BoundY * scale + .5);

  dat << "%!PS-Adobe-2.0" << endl;
  dat << "%%Creator: MooN_MD (Volker Behns)" << endl;
  dat << "%%DocumentFonts: Helvetica" << endl;
  // dat << "%%BoundingBox: 300 300 " << 300+BX << " " << 300+BY << endl;
  dat << "%%BoundingBox: 25 25 " << 35+BX << " " << 35+BY << endl;
  dat << "%%Pages: 1" << endl;
  dat << "%%EndComments" << endl;
  dat << "%%EndProlog" << endl;
  dat << "%%Page: 1 1" << endl;
  dat << "/Helvetica findfont 14 scalefont setfont" << endl;
  dat << "/Cshow { dup stringwidth pop 2 div neg 0 rmoveto show } def" << endl;
  dat << "/M { moveto } def" << endl;
  dat << "/L { lineto } def" << endl;
  // dat << "0.10 0.10 scale" << endl;
  // dat << "10.0 setlinewidth" << endl;
  dat << "0.5 setlinewidth" << endl;

  // loop over all cells
  for(i=0;i<N_;i++)
  {
    CurrCell = Coll->GetCell(i);
    CurrCell->SetClipBoard(i);
    CurrCell->PS(dat, scale, StartX, StartY);
  }

  dat << "stroke" << endl;
  dat << "showpage" << endl;
  dat << "%%Trailer" << endl;
  dat << "%%Pages: 1" << endl;

  return 0;
}

int TDomain::MD_raw(const char *name, Iterators iterator, int arg)
{
  int N_Tri = 0, N_Quad = 0, info;
  std::ofstream dat(name);
  TBaseCell *CurrCell;
  int OutValue[6];

  if (!dat)
  {
    cerr << "cannot open '" << name << "' for output" << endl;
    return -1;
  }

  TDatabase::IteratorDB[iterator]->Init(arg);

  // loop over all cells
  while ((CurrCell = TDatabase::IteratorDB[iterator]->Next(info)))
    if (CurrCell->GetRefDesc()->GetShapeDesc()->GetType() == Triangle)
      N_Tri++;
    else
      N_Quad++;

  // write heads of files
  OutValue[0] = 2;           // dimension
  OutValue[1] = 2;           // number of object classes
  OutValue[2] = Triangle;
  OutValue[3] = N_Tri;
  OutValue[4] = Quadrangle;
  OutValue[5] = N_Quad;

#ifdef __COMPAQ__
  SwapIntArray(OutValue, 6);
#endif
  dat.write((const char *) OutValue, 6*SizeOfInt);

  // fill data in files
  if (N_Tri)
  {
    TDatabase::IteratorDB[iterator]->Init(arg);

    // loop over all cells
    while( (CurrCell = TDatabase::IteratorDB[iterator]->Next(info)))
      if (CurrCell->GetRefDesc()->GetShapeDesc()->GetType() == Triangle)
        CurrCell->MD_raw(dat);
  }

  if (N_Quad)
  {
    TDatabase::IteratorDB[iterator]->Init(arg);

    // loop over all cells
    while ((CurrCell = TDatabase::IteratorDB[iterator]->Next(info)))
      if (CurrCell->GetRefDesc()->GetShapeDesc()->GetType() == Quadrangle)
        CurrCell->MD_raw(dat);
  }

  dat.close();

  return 0;
}

int TDomain::Refine()
{
  TBaseCell *CurrCell;
  int info;

  RefLevel++;

  TDatabase::IteratorDB[It_Finest]->Init(0);

  // loop over all cells
  while( (CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
    CurrCell->Refine(RefLevel);

  return 0;
}

int TDomain::RegRefineAll()
{
  TBaseCell *CurrCell;
  int info;

  RefLevel++;

  TDatabase::IteratorDB[It_Finest]->Init(0);

  // loop over all cells
  while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
  {
    CurrCell->SetRegRefine();
    CurrCell->Refine(RefLevel);
  }

  return 0;
}

int TDomain::RegRefineSub(int ID)
{
  TBaseCell *CurrCell;
  int info;

  RefLevel++;

  TDatabase::IteratorDB[It_Finest]->Init(0);

  // loop over all cells
  while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
    if (CurrCell->GetSubGridID() == ID)
    {
      CurrCell->SetRegRefine();
      CurrCell->Refine(RefLevel);
    }

  return 0;
}

int TDomain::RefineByIndicator(DoubleFunct2D *Indicator)
{
  TBaseCell *CurrCell;
  TVertex *vert;
  int j, k, info;
  int Inner, Outer;
  double x,y,val;

  TDatabase::IteratorDB[It_Finest]->Init(0);

  // loop over all cells
  while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
  {
    Inner=0;
    Outer=0;
    k=CurrCell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      vert=CurrCell->GetVertex(j);
      x=vert->GetX();
      y=vert->GetY();
      Indicator(x,y,&val);
      if(val<=1e-10) Inner++;
      if(val>=-1e-10) Outer++;
    }
    if((Inner>0) && (Outer>0)) 
    {
      // there vertices on both sides
      CurrCell->SetRegRefine();
    }
  }

  Refine();
  return Gen1RegGrid();
}

int TDomain::RefineByErrorEstimator(TCollection *Collection,
                                    double *eta_K, double eta_max,
                                    double eta, bool ConfClosure)
{
  TBaseCell *cell;
  int m, N_,changed, refine_strategy,max_cell_level,it=0,max_it=100;
  double reftol, min_fraction_to_change, coarsetol, decrease_reftol_factor;
  double sum_local_errors, increase_coarsetol_factor;
  double fraction_of_error, min_changed;

  // the following are recommendations
  // all variables in capital letters can be defined in readin.dat
  refine_strategy=0;
  reftol = 0.5;
  coarsetol =0.0;  
  min_fraction_to_change = 0.1;
  decrease_reftol_factor = 0.8;
  increase_coarsetol_factor = 1.1;
  fraction_of_error = 0.25;
  max_cell_level = 100000;

  // the actual values are read from readin.dat

  refine_strategy=TDatabase::ParamDB->REFINE_STRATEGY;  // refinement strategy
  reftol=TDatabase::ParamDB->REFTOL; 
                                  // tolerance for refinement

  coarsetol=TDatabase::ParamDB->COARSETOL; 
                                  // tolerance for derefinement

  if (coarsetol>=reftol)
    coarsetol=0.001*reftol;
  // minimal fraction of cells to refine
  min_fraction_to_change=TDatabase::ParamDB->MIN_FRACTION_TO_CHANGE; 

  // decrease reftol if necessary by this
  decrease_reftol_factor=TDatabase::ParamDB->DECREASE_REFTOL_FACTOR;

  // increase coarsetol if necessary by this
  increase_coarsetol_factor=TDatabase::ParamDB->INCREASE_COARSETOL_FACTOR; 

  // fraction of global error    
  fraction_of_error=TDatabase::ParamDB->FRACTION_OF_ERROR;

  // maximal geo level of a cell
  max_cell_level=TDatabase::ParamDB->MAX_CELL_LEVEL;

  if (refine_strategy==0)         // compare with the maximal local error
  {
    reftol = reftol*eta_max;
    reftol *= reftol;           // since local estimates are stored as square
    coarsetol = coarsetol*eta_max;
    coarsetol = coarsetol*coarsetol;
    N_=Collection->GetN_Cells();       // # cells
    min_changed=N_*min_fraction_to_change; // minimal number of cells to change
    changed=-1;
    while ((changed< min_changed)&&(it<max_it))
    {
      changed=0;
      for(m=0;m<N_;m++)              // loop over all cells
      {
        cell=Collection->GetCell(m);
        // mark cell for refinement
        if ((eta_K[m]>=reftol)&&(cell->GetGeoLevel()<max_cell_level))
        {
          cell->SetRegRefine();    
          changed++;
        }
        if (eta_K[m]<=coarsetol)   // mark cell for coarsening
        {
          ;
        }
      }
      if (changed< min_changed)      // not enough cells marked,
                                     // change tolerances
      {
        reftol*=decrease_reftol_factor;
        coarsetol *=increase_coarsetol_factor;
      }
      Output::print<2>("total ", N_, " changed ",  changed);
      it++;
    }
    if (ConfClosure)
    {
      Output::print<3>(" before ");
      MakeConfClosure();
      Output::print<3>(" after ");
      
      return 0;
    }
    else
    {  
      Refine();                          // refine marked cells
      return Gen1RegGrid();              // make grid 1--regular 
    }
  }

  if (refine_strategy==1)       // mark smalles set of cells whose
                                // sum of local errors
  {                             // is prescribed fraction of global error
    reftol = 0.9*eta_max;       // set refinement tolerance
    reftol *= reftol;           // since local estimates are stored as square
    coarsetol = coarsetol*eta_max;
    coarsetol = coarsetol*coarsetol;
    N_=Collection->GetN_Cells();       // # cells
    min_changed=N_*min_fraction_to_change;   // minimal number of cells
                                             // to be changed
    changed=-1;
    sum_local_errors = 0.0;            // sum of local errors of marked cells
    while ((changed< min_changed)&&(sum_local_errors <=
            fraction_of_error * eta)&&(it<max_it))
    {
      changed=0;
      for(m=0;m<N_;m++)              // loop over all cells
      {
        cell=Collection->GetCell(m);
        // mark cell for refinement
        if ((eta_K[m]>=reftol)&&(cell->GetGeoLevel()<max_cell_level))
        {
          cell->SetRegRefine();    
          sum_local_errors += eta_K[m];
          changed++;
        }
        if (eta_K[m]<=coarsetol)   // mark cell for derefinement
        {
          ;
        }
      }
      Output::stat<2>("total ", N_ , " changed ",  changed, " global error "
                  , eta, " error of refined cells ", sqrt(sum_local_errors));
      if ((changed< min_changed)||(sqrt(sum_local_errors) <=
           fraction_of_error * eta))
      {                  // criteria not fulfilled, change tolerances
        reftol*=decrease_reftol_factor;
        coarsetol *=increase_coarsetol_factor;
        sum_local_errors = 0.0;
      }
      it++;
    }
    if (ConfClosure)
    {
      MakeConfClosure();
      return 0;
    }
    else
    {
      Refine();
      return Gen1RegGrid();
    }
  }
  return 0;
}

#ifdef __2D__
int TDomain::MakeConfClosure()
{
  TBaseCell *CurrCell, *parent;
  int i, info, MaxLevel, clip;
  Refinements type;

  const char filename[] = "before_closure.ps";

  // delete existing closures
  TDatabase::IteratorDB[It_Finest]->Init(0);
  while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
    if ((parent = CurrCell->GetParent()))
    {
      type = parent->GetRefDesc()->GetType();
      if (type >= TriBis0 && type <= Quad2Conf3)
      {
        // check whether a child is marked for refinement
        // TriReg means regular refinement for triangles and quadrangles
        type = CurrCell->GetRefDesc()->GetType();
        CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info);
        if (type < CurrCell->GetRefDesc()->GetType())
          type = TriReg;

        if (parent->GetN_Children() == 3)
        {
          CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info);
          if (type < CurrCell->GetRefDesc()->GetType())
            type = TriReg;
        }

        parent->Derefine();

        if (type == TriReg)
          parent->SetRegRefine();
      }
    }
    

  // generate a new 1-regular grid
  Refine();
  Gen1RegGrid();

  // initialize clipboards
  TDatabase::IteratorDB[It_Finest]->Init(0);
  while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
    CurrCell->SetClipBoard(0);

  // look for elements which have to be refined
  TDatabase::IteratorDB[It_Finest]->Init(0);
  MaxLevel = TDatabase::IteratorDB[It_Finest]->GetMaxLevel();
  PS(filename,It_Finest,0);
  OutPut("MaxLevel " << MaxLevel << endl);
  
  for (i=MaxLevel;i>=0;i--)
  {
      // get iterator on level i
      TDatabase::IteratorDB[It_EQ]->Init(i);
      // loop over the mesh cells on level i
      while ((CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)))
	  // if current cell does not possess children, compute conforming closure
	  if (!CurrCell->ExistChildren())
	  {
	      CurrCell->MakeConfClosure();
	  }
      // get iterator on finest level
      TDatabase::IteratorDB[It_Finest]->Init(0);
      // loop over mesh cells on finest level
      while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
      {
	  // get clip board
	  clip = CurrCell->GetClipBoard();
	  //cout << i << " clip " << clip << endl;
	  // compute refinement rule
	  switch (clip)
	  {
	      case 145: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + TriBis0]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 146: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + TriBis1]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 148: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + TriBis2]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 273: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + Quad1Conf2]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 274: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + Quad1Conf3]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 276: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + Quad1Conf0]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 280: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + Quad1Conf1]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 291: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + Quad2Conf1]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 293: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + QuadBis0]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 294: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + Quad2Conf2]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 297: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + Quad2Conf3]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 298: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + QuadBis1]);
                  CurrCell->Refine(RefLevel);
                  break;
	      case 300: CurrCell->SetRefDesc(TDatabase::RefDescDB[
						 N_SHAPES + Quad2Conf0]);
                  CurrCell->Refine(RefLevel);
                  break;
	      default: if (clip > 512)
	      {
		  CurrCell->SetRegRefine();
		  CurrCell->Refine(RefLevel);
	      }
	  }
    }
  }

  return 0;
}

#else

int TDomain::CloseGrid(int level)
{
  TBaseCell *CurrCell, *Cell, *Child, *Neigh;
  int i, j, k, clip;
  int N_Edges, N_Children;
  int info, LocEdge, CurrLocEdge, LocFace=0, NeighLocFace;
  int MapType;
  int first_bis, second_bis, neigh_first_bis, neigh_second_bis, neigh_type;
  TJoint *Joint, *LastJoint;
  const int *TmpEF, *TmpEF2;
  int TmpEFMaxLen, TmpEF2MaxLen;
  bool RefineRegular;
  Refinements NeighFaceRef, MyFaceRef;
  int N_ToRefine = 0;

  // Reset ClipBoards
  TDatabase::IteratorDB[It_EQ]->Init(level);
  k=0;
  while( (CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)))
  {
    N_Children = CurrCell->GetRefDesc()->GetN_Children();

    // Check if CurrCell is refined irregularly but has a son which has to be refined
    if(2 <= N_Children && N_Children < 8) // CurrCell is refined irregularly
    {
      RefineRegular = false;

      for(i=0; i<N_Children; ++i)
      {
        Child = CurrCell->GetChild(i);
        // Check if Child is refined directly
        
        //if(Child->GetRefDesc()->GetType() != NoRef) // TODO
        if(Child->GetClipBoard() > 0)
          RefineRegular = true;

        // Check if Child contains an edge that is refined by a neighbour
        if(!RefineRegular)
        {         
          // TODO
        }
        
        if(RefineRegular)
        {
          CurrCell->SetRegRefine();
          break;
        }
        else
        CurrCell->SetNoRefinement();
      }
    }
    
    // Initialize Clipboard
    if(CurrCell->GetRefDesc()->GetType() >= TetraReg 
       && CurrCell->GetRefDesc()->GetType() <= TetraReg2)
      CurrCell->SetClipBoard(63);
    else
      CurrCell->SetClipBoard(-1);
  }

  /*
   * Set Clipboard to 0 if element is unrefined but contains an edge that is refined by an other tetrahedron
   */
  k = -1;
  TDatabase::IteratorDB[It_EQ]->Init(level);
  while ( (Cell = TDatabase::IteratorDB[It_EQ]->Next(info)))
    {
      Cell->GetShapeDesc()->GetEdgeFace(TmpEF, TmpEFMaxLen);
      N_Edges = Cell->GetN_Edges();

      k++;
      if(Cell->GetClipBoard() == 63)
  {
    for (i=0;i<N_Edges;i++)
      {
        if(Cell->GetClipBoard() & (pow(2, i) == 0))
    continue;

        LocEdge = i;
        for(j=0; j<TmpEFMaxLen; ++j)
    {
      LastJoint = Cell->GetJoint(TmpEF[2*LocEdge+j]);

      if(!(CurrCell = LastJoint->GetNeighbour(Cell)))
        continue;

      CurrLocEdge = LastJoint->GetNeighbourEdgeIndex(Cell, LocEdge);

      while(CurrCell != Cell && CurrCell)
        {
          // Check if Edge is refined in CurrCell
          if(CurrCell->GetClipBoard() == -1)
      {
        CurrCell->SetClipBoard(0);
        N_ToRefine++;
      }

          // Get next joint which contains this edge
          CurrCell->GetShapeDesc()->GetEdgeFace(TmpEF2, TmpEF2MaxLen);
          if(CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]) == LastJoint)
      Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge+1]);
          else
      Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]);

          // Get new element and the index of our edge in this element
          CurrLocEdge = Joint->GetNeighbourEdgeIndex(CurrCell, CurrLocEdge);
          CurrCell = Joint->GetNeighbour(CurrCell);
          LastJoint = Joint;
        }
    } // j=0..N_JointsPerEdge
      } // i=0..N_edges
  }
    }

  while(N_ToRefine > 0)
    {
      TDatabase::IteratorDB[It_EQ]->Init(level);
      k=0;
      while ( (CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)) )
  {
    if(CurrCell->GetClipBoard() == 0)
      {
        //std::cout << "Scanning Cell " << k << "\n";
        N_ToRefine += CurrCell->MakeConfClosure() - 1;
      }
    k++;
  }
    }
      
  k=-1;
  TDatabase::IteratorDB[It_EQ]->Init(level);
  while ((CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)))
    {
      k++;
      clip = CurrCell->GetClipBoard();

      switch(clip)
  {
    /*
     * Unique Refinements
     */
    // NoRef
  case 0:
    std::cerr << "This should not happen\n";
    break;
  case -1:
    CurrCell->SetNoRefinement();
    break;
    // Bisections
  case 1:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis0]);
    break;
  case 2:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis1]);
    break;
  case 4:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis2]);
    break;
  case 8:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis3]);
    break;
  case 16:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis4]);
    break;
  case 32:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis5]);
    break;
    // Diagonal Bisections
  case 33:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis05]);
    break;
  case 10:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis13]);
    break;
  case 20:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis24]);
    break;
    // QuadX
  case 7:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraQuad0]);
    break;
  case 25:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraQuad1]);
    break;
  case 50:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraQuad2]);
    break;
  case 44:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraQuad3]);
    break;
    // Reg
  case 63:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraReg]);
    break;
  case 3:
  case 5:
  case 9:
  case 17:
  case 6:
  case 18:
  case 34:
  case 12:
  case 36:
  case 24:
  case 40:
  case 48:
    // Determine Local Face Index of face which is refined by double bisection
    switch(clip)
      {
      case 3: LocFace = 0; break;
      case 5: LocFace = 0; break;
      case 9: LocFace = 1; break;
      case 17: LocFace = 1; break;
      case 6: LocFace = 0; break;
      case 18: LocFace = 2; break;
      case 34: LocFace = 2; break;
      case 12: LocFace = 3; break;
      case 36: LocFace = 3; break;
      case 24: LocFace = 1; break;
      case 40: LocFace = 3; break;
      case 48: LocFace = 2; break;
      }

    // Check if Neighbour is already refined
    Joint = CurrCell->GetJoint(LocFace);
    Neigh = Joint->GetNeighbour(CurrCell);

    if(Neigh && Neigh->ExistChildren())
      {
        // Find Local Face Index at Neigh
        for(NeighLocFace=0; NeighLocFace<Neigh->GetN_Faces(); ++NeighLocFace)
    if(Joint == Neigh->GetJoint(NeighLocFace))
      break;

        if(NeighLocFace == Neigh->GetN_Faces())
    {
      std::cerr << "Face was not found at Neighbour\n";

    }

        NeighFaceRef = Neigh->GetRefDesc()->GetFaceRef(NeighLocFace);

        MapType = Joint->GetMapType();

        first_bis = (NeighFaceRef-TriBis01) / 2;
        second_bis = (NeighFaceRef-TriBis01) % 2;
        if(first_bis <= second_bis) second_bis++;

        neigh_first_bis = ((2-first_bis) + MapType) % 3;
        neigh_second_bis = ((2-second_bis) + MapType) % 3;

        neigh_type = 2*neigh_first_bis + neigh_second_bis;
        if(neigh_second_bis > neigh_first_bis) neigh_type--;

        MyFaceRef = Refinements(TriBis01 + neigh_type);

        switch(clip)
    {
      /*
       * Non-unique refinements
       */
      // Bis0X
    case 3:
      if(MyFaceRef == TriBis01)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis01]);
      else if(MyFaceRef == TriBis10)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis10]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 5:
      if(MyFaceRef == TriBis02)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis02]);
      else if(MyFaceRef == TriBis20)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis20]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 9:
      if(MyFaceRef == TriBis20)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis03]);
      else if(MyFaceRef == TriBis02)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis30]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 17:
      if(MyFaceRef == TriBis21)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis04]);
      else if(MyFaceRef == TriBis12)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis40]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
      // Bis 1X
    case 6:
      if(MyFaceRef == TriBis12)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis12]);
      else if(MyFaceRef == TriBis21)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis21]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 18:
      if(MyFaceRef == TriBis01)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis14]);
      else if(MyFaceRef == TriBis10)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis41]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 34:
      if(MyFaceRef == TriBis02)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis15]);
      else if(MyFaceRef == TriBis20)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis51]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
      // Bis2X
    case 12:
      if(MyFaceRef == TriBis02)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis23]);
      else if(MyFaceRef == TriBis20)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis32]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 36:
      if(MyFaceRef == TriBis01)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis25]);
      else if(MyFaceRef == TriBis10)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis52]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
      // Bis 3X
    case 24:
      if(MyFaceRef == TriBis01)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis34]);
      else if(MyFaceRef == TriBis10)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis43]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 40:
      if(MyFaceRef == TriBis21)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis35]);
      else if(MyFaceRef == TriBis12)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis53]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
      // Bis 4X
    case 48:
      if(MyFaceRef == TriBis12)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis45]);
      else if(MyFaceRef == TriBis21)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis54]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    }
      }
    else
      {
        switch(clip)
    {
      /*
       * Non-unique refinements
       */
      // Bis0X
    case 3:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis01]);
      break;
    case 5:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis02]);
      break;
    case 9:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis03]);
      break;
    case 17:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis04]);
      break;
      // Bis 1X
    case 6:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis12]);
      break;
    case 18:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis14]);
      break;
    case 34:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis15]);
      break;
      // Bis2X
    case 12:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis23]);
      break;
    case 36:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis25]);
      break;
      // Bis 3X
    case 24:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis34]);
      break;
    case 40:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis35]);
      break;
      // Bis 4X
    case 48:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis45]);
      break;
    }
      }


    break;
  default:
    std::cout << "Undefined Reference Description: " << CurrCell->GetClipBoard() << "\n";
    break;
  }

      CurrCell->Refine(level+1);
    }
  return 0;
}

/**
 * This function is used to set the right refinement descriptions to each cell such that the grid is conforming afterwards.
 * After refinement of one ore more selected Cells this method has to be called.
 *
 * Max Winkler (23.02.2012)
 */
int TDomain::MakeConfClosure()
{
  int m;
  int MaxLevel;
  
  // Generate 1-regular grid
  Gen1RegGrid();

  TDatabase::IteratorDB[It_Finest]->Init(0);
  MaxLevel = TDatabase::IteratorDB[It_Finest]->GetMaxLevel();

  for(m=MaxLevel; m>=0; --m)
    {
      CloseGrid(m);
    }
  return 0;
}
#endif

#ifdef __2D__
int TDomain::Gen1RegGrid()
{
  int MaxLevel, CurrLevel, info;
  TBaseCell *CurrCell;

  TDatabase::IteratorDB[It_Finest]->Init(0);
  MaxLevel = TDatabase::IteratorDB[It_Finest]->GetMaxLevel();

  for (CurrLevel=MaxLevel;CurrLevel>0;CurrLevel--)
  {
    TDatabase::IteratorDB[It_EQ]->Init(CurrLevel);

    while ((CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)))
      if (!CurrCell->ExistChildren())
        CurrCell->Gen1RegGrid();
  }
  return 0;
}
#else
int TDomain::Gen1RegGrid()
{
  TBaseCell* CurrCell;
  int info;

  TDatabase::IteratorDB[It_Finest]->Init(0);
  while((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)) )   
    CurrCell->Gen1RegGrid();

  return 0;
}
#endif

int TDomain::ConvertQuadToTri(int type)
{
  TBaseCell *CurrCell;
  Shapes CellType;
  int info;

  if (type)
  {
    RefLevel++;

    TDatabase::IteratorDB[It_Finest]->Init(0);

    // loop over all cells
    while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
    {
      CellType = CurrCell->GetType();
      if (CellType == Parallelogram || CellType == Quadrangle ||
          CellType == Rectangle)
      {
        if (type == 1)
        {
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES +
                               QuadToTri0]);
        }
        else
        {
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES +
                               QuadToTri1]);
        }
        CurrCell->Refine(RefLevel);
      }
    }
  }

  return 0;
}

/*
 Extract a subcollection from a collection object:
*/
TCollection *TDomain::GetCollection(TCollection *coll, int reference)
{
  Output::info<2>("GetCollection", " Domain::GetCollection with reference: ", reference);
  int n_cells;
  TBaseCell **cells, *CurrCell;
  TCollection *subcoll;

  n_cells = 0;
  for (int i=0;i<coll->GetN_Cells();i++)
  {
    CurrCell = coll->GetCell(i);
    if (CurrCell->GetReference_ID()==reference)
    {
      n_cells++;
    }
  }

  // fill array of cells
  cells = new TBaseCell*[n_cells];
  int j = 0;
  for (int i=0;i<coll->GetN_Cells();i++)
  {
    CurrCell = coll->GetCell(i);
    if (CurrCell->GetReference_ID()==reference)
    {
      cells[j] = CurrCell;
      j++;
    }
  }
  Output::print<2>("TDomain::GetCollection() creating collection, n_cells = ", n_cells);
  // create collection from an array of cells
  subcoll = new TCollection(n_cells, cells);
  
  return subcoll;
}


/*
 Extract a collection from a Domain object:
 For a given Iterator
   - count how many cells it contain
   - add cells to collection
*/
TCollection *TDomain::GetCollection(Iterators it, int level) const
{
  TCollection *coll;
  int i, n_cells, info;
  TBaseCell **cells, *CurrCell;

  // initialize the iterator
  // note:
  // enum Iterators {It_EQ, It_LE, It_Finest, It_EQLevel, It_LELevel,
  //            It_Between, It_OCAF, It_Mortar1, It_Mortar2};
  TDatabase::IteratorDB[it]->Init(level);
  n_cells=0;
  // if the pointer to the next item is not empty increase number of cells
  // info: set to the current level
  while (TDatabase::IteratorDB[it]->Next(info)) n_cells++;

  // fill array of cells
  cells=new TBaseCell*[n_cells];
  TDatabase::IteratorDB[it]->Init(level);
  i=0;
  while (((CurrCell = TDatabase::IteratorDB[it]->Next(info))))
  {
    cells[i]=CurrCell;
    cells[i]->SetCellIndex(i);
    i++;
  }

  // create collection from an array of cells
  coll=new TCollection(n_cells, cells);

  #ifdef  _MPI 
  coll->SetN_OwnCells(N_OwnCells);
  #endif

  return coll;
}

TCollection *TDomain::GetCollection(Iterators it, int level, int ID) const
{
  if(ID == -4711) 
    return this->GetCollection(it, level);
  
  int info;
  // initialize the iterator
  // note:
  // enum Iterators {It_EQ, It_LE, It_Finest, It_EQLevel, It_LELevel,
  //            It_Between, It_OCAF, It_Mortar1, It_Mortar2};
  TDatabase::IteratorDB[it]->Init(level);
  
  //cout << " count cells " << endl;
  int n_cells = 0;
  // loop over all cells, check their reference id
  // info: set to the current level
  while(TBaseCell *CurrCell = TDatabase::IteratorDB[it]->Next(info))
  {
    //cout << " id: " << CurrCell->GetPhase_ID() << endl;
    if (CurrCell->GetReference_ID() == ID)
      n_cells++;
  }
  
  //cout << " fill cells " << endl;
  // fill array of cells
  TBaseCell **cells = new TBaseCell*[n_cells];
  int i=0;
  TDatabase::IteratorDB[it]->Init(level);
  while(TBaseCell *CurrCell = TDatabase::IteratorDB[it]->Next(info))
  {
    if (CurrCell->GetReference_ID() == ID)
    {
      cells[i] = CurrCell;
      cells[i]->SetCellIndex(i);
      i++;
    }
  }
  
  // create collection from an array of cells
  TCollection *coll = new TCollection(n_cells, cells);
  
  #ifdef  _MPI 
  coll->SetN_OwnCells(N_OwnCells);
  #endif
  
  return coll;
}


#ifdef  _MPI 
TCollection *TDomain::GetOwnCollection(Iterators it, int level, int ID) const
{
  TCollection *coll;
  int i, n_cells, info;
  TBaseCell **cells, *CurrCell;


  TDatabase::IteratorDB[it]->Init(level);
  n_cells=0;
  while ((CurrCell = TDatabase::IteratorDB[it]->Next(info)))
  {
   if(ID == CurrCell->GetSubDomainNo() )
    n_cells++;
  }

  cells=new TBaseCell*[n_cells];
  TDatabase::IteratorDB[it]->Init(level);
  i=0;
  while ((CurrCell = TDatabase::IteratorDB[it]->Next(info)))
  {
   if(ID == CurrCell->GetSubDomainNo())
    {
     cells[i]=CurrCell;
     i++;
    }
  }

  coll = new TCollection(n_cells, cells);
  coll->SetN_OwnCells(N_OwnCells);

  return coll;
}
#endif

void TDomain::GetSortedCollection(TCollection* &Coll, int* &Indices)
{
  int ActiveLevel=0;
  int ActiveRootCell=0;
  // int N_RootCells=N_InitVCells;
  int Level=RefLevel;
  TBaseCell *ActiveCell=nullptr, **Cells;
  int i,*N_Children, *CurrentChild, *CurrentIndex;

  N_Children=new int[RefLevel+1];
  CurrentChild=new int[RefLevel+1];
  Indices=new int[RefLevel+1];
  CurrentIndex=new int[RefLevel+1];

  for(i=0;i<=RefLevel;i++)
  {
    N_Children[i]=0;
    CurrentChild[i]=0;
    Indices[i]=0;
  }

  do {

    if (ActiveLevel)
      do
      {
        ActiveCell = ActiveCell->GetParent();
        N_Children[ActiveLevel] = 0;
        CurrentChild[ActiveLevel] = 0;
        ActiveLevel--;
      } while (ActiveLevel && 
                  N_Children[ActiveLevel] == CurrentChild[ActiveLevel]);
  
    if (!Level  || (!ActiveLevel &&
          N_Children[ActiveLevel] == CurrentChild[ActiveLevel]))
    {
      if (ActiveRootCell < N_RootCells)
      {
        ActiveCell = CellTree[ActiveRootCell];
        ActiveRootCell++;
        N_Children[ActiveLevel] = ActiveCell->GetN_Children();
        CurrentChild[ActiveLevel] = 0;
      }
      else
      { ActiveCell=NULL;}
     }
    while (N_Children[ActiveLevel] && ActiveLevel != Level && ActiveCell)
    {
      ActiveCell = ActiveCell->GetChild(CurrentChild[ActiveLevel]);
      CurrentChild[ActiveLevel]++;
      ++ActiveLevel;
      N_Children[ActiveLevel] = ActiveCell->GetN_Children();
    }

    if (ActiveCell) 
    { 
      // cout << "new Cell at level: " << ActiveLevel << endl;
      Indices[ActiveLevel]++;
    }
  
  } while(ActiveCell);

  CurrentIndex[0]=0;
  for(i=1;i<=RefLevel;i++)
  {
    Indices[i]+=Indices[i-1];
    CurrentIndex[i]=Indices[i-1];
  }

  Cells=new TBaseCell* [Indices[RefLevel]];

  ActiveLevel=0;
  ActiveRootCell=0;

  for(i=0;i<=RefLevel;i++)
  {
    N_Children[i]=0;
    CurrentChild[i]=0;
  }

  do {
    if (ActiveLevel)
      do
      {
        ActiveCell = ActiveCell->GetParent();
        N_Children[ActiveLevel] = 0;
        CurrentChild[ActiveLevel] = 0;
        ActiveLevel--;
      } while (ActiveLevel && 
                  N_Children[ActiveLevel] == CurrentChild[ActiveLevel]);
  
    if (!Level  || (!ActiveLevel &&
          N_Children[ActiveLevel] == CurrentChild[ActiveLevel]))
     {
      if (ActiveRootCell < N_RootCells)
      {
        ActiveCell = CellTree[ActiveRootCell];
        ActiveRootCell++;
        N_Children[ActiveLevel] = ActiveCell->GetN_Children();
        CurrentChild[ActiveLevel] = 0;
      }
      else
      { ActiveCell=NULL;}
     }
    while (N_Children[ActiveLevel] && ActiveLevel != Level && ActiveCell)
    {
      ActiveCell = ActiveCell->GetChild(CurrentChild[ActiveLevel]);
      CurrentChild[ActiveLevel]++;
      ++ActiveLevel;
      N_Children[ActiveLevel] = ActiveCell->GetN_Children();
    }

    if (ActiveCell) 
    { 
      // cout << "new Cell at level: " << ActiveLevel << endl;
      Cells[CurrentIndex[ActiveLevel]]=ActiveCell;
      CurrentIndex[ActiveLevel]++;
    }
  } while(ActiveCell);

  Coll=new TCollection(Indices[RefLevel], Cells);

  #ifdef  _MPI 
  Coll->SetN_OwnCells(N_OwnCells);
  #endif
  
  delete [] N_Children;
  delete [] CurrentChild;
  delete [] CurrentIndex;

}

/*
void TDomain::Refine1Reg(int MinLevel, int MaxLevel)
{
  TCollection *coll;
  TBaseCell *cell;
  int *indices, maxlevel;
  int i,j,k,N_;

  RefLevel++; // to remove later
  GetSortedCollection(coll, indices);

  N_=indices[RefLevel];
  // cout << "number of cells: " << N_ << endl;
  for(i=0;i<N_;i++)
  {
    cell=coll->GetCell(i);
    if(cell->IsToRefine())
      cell->SetClipBoard(Refinement);
    else
      cell->SetClipBoard(NoRefinement);
  }
  coll->~TCollection();

  if(RefLevel==3)
  {
  CellTree[0]->GetChild(1)->GetChild(0)->SetClipBoard(Remove);
  CellTree[0]->GetChild(1)->GetChild(1)->SetClipBoard(Remove);
  CellTree[0]->GetChild(1)->GetChild(2)->SetClipBoard(Remove);
  CellTree[0]->GetChild(1)->GetChild(3)->SetClipBoard(Remove);

  // CellTree[1]->GetChild(0)->GetChild(2)->SetClipBoard(Refinement);
  }

  for(k=RefLevel;k>=0;k--)
  {
    // cout << "Level " << k << endl;

    coll=GetCollection(It_EQ, k);
    N_=coll->GetN_Cells();

    if(k>MaxLevel)
    {
      for(i=0;i<N_;i++)
      {
        cell=coll->GetCell(i);
        cell->SetClipBoard(Remove);
      }
    }

    if(k==MaxLevel)
    {
      for(i=0;i<N_;i++)
      {
        cell=coll->GetCell(i);
        if(cell->GetClipBoard()==Refinement)
          cell->SetClipBoard(NoRefinement);
      }
    }

    if(k<MinLevel)
    {
      for(i=0;i<N_;i++)
      {
        cell=coll->GetCell(i);
        cell->SetClipBoard(Refinement);
      }
    }

    for(i=0;i<N_;i++)
    {
      cell=coll->GetCell(i);
      // cout << cell->GetClipBoard() << endl;
      cell->Gen1RegMarks();
      //cout << i << endl;
      cell->RefDeref();
    }
    coll->~TCollection();
  }

}
*/

/*
void TDomain::Refine1Reg(int MinLevel, int MaxLevel)
{
  TCollection *coll;
  TBaseCell *cell;
  int *indices, maxlevel;
  int i,j,k,N_;

  for(k=MaxLevel;k>=MinLevel;k--)
  {
    // cout << "Level " << k << endl;

    coll=GetCollection(It_EQ, k);
    N_=coll->GetN_Cells();

    for(i=0;i<N_;i++)
    {
      cell=coll->GetCell(i);
      // cout << cell->GetClipBoard() << endl;
      cell->Gen1RegMarks();
      // cout << i << endl;
      cell->RefDeref();
    }

    coll->~TCollection();
    delete coll;
  }
  RefLevel++;
}
*/

/*
void EvaluateMarks(TCollection *coll)
{
  TBaseCell *cell;
  int i,N_, j,k, all;

  N_=coll->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    cell=coll->GetCell(i);
    if((cell->GetClipBoard()==Refinement) && (!cell->ExistChildren()))
    {
      // cout << "checking 1-regularity to coarser cells" << endl;
      cell->Gen1RegMarks();
    }

    if(cell->ExistChildren())
    {
      k=cell->GetN_Children();
      all=0;
      for(j=0;j<k;j++)
        if(cell->GetChild(j)->GetClipBoard()==Remove) all++;
      if(all==k) // all children are marks for removing
      {
        cell->SetClipBoard(NoRefinement);
        cell->Gen1RegMarks();
        // cout << "checking 1-regularity to finer cells" << endl;
      }
    }
  }
}
*/

/*
void CoarseGrid(TCollection *coll)
{
  TBaseCell *cell;
  int i,N_;

  N_=coll->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    cell=coll->GetCell(i);
    if(cell->ExistChildren() && cell->GetClipBoard()==NoRefinement)
    {
      // cout << "derefine" << endl;
      cell->Derefine();
    }
  }
}
*/

/*
void RefineGrid(TCollection *coll)
{
  TBaseCell *cell;
  int i,N_;

  N_=coll->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    cell=coll->GetCell(i);
    if(!cell->ExistChildren() && cell->GetClipBoard()==Refinement)
    {
      // cout << "refine" << endl;
      cell->SetRegRefine();
      cell->Refine(0);
    }
  }
}
*/

/*
void TDomain::Refine1Reg(int MinLevel, int MaxLevel)
{
  TCollection *coll;
  int k,N_;

  for(k=MaxLevel;k>=0;k--)
  {
    coll=GetCollection(It_EQ, k);
    EvaluateMarks(coll);
    delete coll;
  }

  for(k=0;k<=MaxLevel;k++)
  {
    coll=GetCollection(It_EQ, k);
    N_=coll->GetN_Cells();
    if(k<MaxLevel) CoarseGrid(coll);
    RefineGrid(coll);
    coll->~TCollection();
    delete coll;
  }
}
*/

void TDomain::DeRefine()
{
  CellTree[0]->Derefine();
}

#ifdef __3D__
static void Sort(TVertex **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, *rr, len;
  TVertex *Mid, *Temp;

  len=(int)(2*log((double)length)/log(2.0)+2.0);
  rr=new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete[] rr;
}

static int GetIndex(TVertex **Array, int Length, TVertex *Element)
{
  int l=0, r=Length, m=(r+l)/2;
  TVertex *Mid;

  Mid=Array[m];
  while(Mid!=Element)
  {
    if(Mid>Element)
    {
      l=m;
    }
    else
    {
      r=m;
    }
    m=(r+l)/2;
    Mid=Array[m];
  }
  return m;
}

int TDomain::Grape(const char *name, TCollection *coll)
{
  int Format, Version, N_Elements, N_Vertices, N_LocVertices;
  int Header[5];
  int i,j,k,l,m, N_;
  int *Type;
  int *VertexNumbers, *NumberVertex;
  double *Coords;
  TBaseCell *cell;
  TVertex **Vertices, *Last, *Current;

  std::ofstream dat(name);
  if(!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }

  Format = 1;
  Version = 0;
  N_Elements = coll->GetN_Cells();

  Type = new int[N_Elements];

  N_LocVertices = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
    switch(cell->GetType())
    {
      case Tetrahedron:
        Type[i] = 0;
      break;

      case Hexahedron:
      case Brick:
        Type[i] = 1;
      break;
      default:
        ErrMsg("TDomain::Grape only in 3D");
        exit(1);
    }
  }

  Vertices = new TVertex*[N_LocVertices];
  N_ = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = coll->GetCell(i);
    k = cell->GetN_Vertices();
    for(j=0;j<k;j++)
      Vertices[N_++] = cell->GetVertex(j);
  }

  // sort the Vertices array
  Sort(Vertices, N_);

  Last = NULL;
  N_Vertices = 0;
  for(i=0;i<N_LocVertices;i++)
    if ((Current = Vertices[i]) != Last)
    {
      N_Vertices++;
      Last = Current;
    }

  Header[0] = Format;
  Header[1] = Version;
  Header[2] = N_Elements;
  Header[3] = N_Vertices;
  Header[4] = N_LocVertices;

  dat.write((char *)Header,SizeOfInt*5);
  dat.write((char *)Type,SizeOfInt*N_Elements);
  delete Type;

  Coords = new double[3*N_Vertices];
  VertexNumbers = new int[N_LocVertices];
  NumberVertex = new int[N_LocVertices];

  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1],Coords[N_+2]);
      // cout << N_/2 << "  "  << Coords[N_] << "       ";
      // cout << Coords[N_+1] << "       " << Coords[N_+2] << endl;
      k++;
      N_ += 3;
      Last=Current;
    }
    NumberVertex[i]=k;
  }

  m=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = coll->GetCell(i);
    k = cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      // cout << VertexNumbers[m] << endl;
      m++;
    } // endfor j
  } //endfor i

  dat.write((char *)VertexNumbers,sizeof(int)*N_LocVertices);
  dat.write((char *)Coords,SizeOfDouble*3*N_Vertices);
  
  delete Coords;
  delete NumberVertex;
  delete VertexNumbers;
  delete Vertices;

  return 0;
}

/** make boundary parameter consistent */
void TDomain::MakeBdParamsConsistent(TCollection *coll)
{
  int i,j,k;
  int N_Cells, N_Joints;
  double Param1[4], Param2[4];
  TBaseCell *cell;
  TJoint *joint;
  TBoundFace *boundface;
  TBoundComp3D *bdcomp;
  const int *TmpFV, *TmpLen;
  int MaxLen;
  double x,y,z;

  N_Cells = coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();
    cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == BoundaryFace ||
                joint->GetType() == IsoBoundFace)
      {
        boundface = (TBoundFace *)joint;
        bdcomp = boundface->GetBoundComp();

        if(bdcomp->GetType() == Plane)
        {
          for(k=0;k<TmpLen[j];k++)
          {
            cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(x,y,z);
            bdcomp->GetTSofXYZ(x,y,z,Param1[k], Param2[k]);
          } // endfor k
          boundface->SetParameters(Param1, Param2);
        } // endif Plane

        if(bdcomp->GetType() == Wall) 
        {
          if( ((TBdWall*)bdcomp)->GetBdComp2D()->GetType() == Line  )
          {
            for(k=0;k<TmpLen[j];k++)
            {
              cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(x,y,z);
              bdcomp->GetTSofXYZ(x,y,z,Param1[k], Param2[k]);
            } // endfor k
            boundface->SetParameters(Param1, Param2);
          } // endif Line
        } // endif Wall

      } // endif BoundaryFace
    } // endfor j
  } // endfor i
} // MakeBdParamsConsistent

#endif

#ifdef __2D__
/** set boundary parameters and cell shape according to
    possibly moving vertices */
//  TDatabase::IteratorDB[It_Finest]->Init(0);
//  MaxLevel = TDatabase::IteratorDB[It_Finest]->GetMaxLevel();
//  
//  for (CurrLevel=MaxLevel;CurrLevel>0;CurrLevel--)
//  {
//    TDatabase::IteratorDB[It_EQ]->Init(CurrLevel);
//
//    while ((CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)))
//      if (!CurrCell->ExistChildren())
//        CurrCell->Gen1RegGrid();
//  }
void TDomain::CorrectParametersAndShapes()
{
  int MaxLevel, CurrLevel;
  int info, j, N_Edges;
  TVertex **Vertices;
  TJoint *joint;
  TGridCell *ActiveCell;

  TDatabase::IteratorDB[It_Finest]->Init(0);
  MaxLevel = TDatabase::IteratorDB[It_Finest]->GetMaxLevel();
  
  for (CurrLevel=MaxLevel;CurrLevel>=0;CurrLevel--)
  {
    TDatabase::IteratorDB[It_EQ]->Init(CurrLevel);

    while ((ActiveCell = (TGridCell *)(TDatabase::IteratorDB[It_EQ]->Next(info))))
    {
      Vertices = ActiveCell->GetVertices();
      N_Edges = ActiveCell->GetN_Edges();
      for(j=0;j<N_Edges;j++)
      {
        joint = ActiveCell->GetJoint(j);
        if(joint->GetType() == BoundaryEdge)
        {
          ((TBoundEdge *)joint)->UpdateParameters(Vertices[j],
                            Vertices[(j+1)%N_Edges]);
        }
        if(joint->GetType() == InterfaceJoint)
        {
          ((TInterfaceJoint *)joint)->UpdateParameters(Vertices[j],
                            Vertices[(j+1)%N_Edges]);
        }
      }
    
      // if Cell is a quadrilateral
      if(N_Edges == 4)
      {
        if(ActiveCell->ExistChildren())
          ActiveCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES+((TQuadrangle *)
               TDatabase::RefDescDB[Quadrangle]->GetShapeDesc())->
               CheckQuad(Vertices)]);
        else
          ActiveCell->SetRefDesc(TDatabase::RefDescDB[((TQuadrangle *)
               TDatabase::RefDescDB[Quadrangle]->GetShapeDesc())->
               CheckQuad(Vertices)]);
      } // endif
    } // endwhile
  } // endfor
} // CorrectParametersAndShapes
//  TGridCell *Cell, *ActiveCell;
//  int i,j,k;
//  int N_Edges;
//  TJoint *joint;
//  TVertex **Vertices;
//  int ActiveLevel, ActiveRootCell;
//  int *N_Children, *CurrentChild;
//  int MaxN_Levels;
//
//  ActiveLevel = 0;
//  ActiveRootCell = 0;
//  MaxN_Levels = TDatabase::IteratorDB[It_Finest]->GetMaxLevel();
//  MaxN_Levels = 3;
//  N_Children = new int[MaxN_Levels];
//  CurrentChild = new int[MaxN_Levels];
//  for(k=0;k<MaxN_Levels;k++)
//  {
//    N_Children[k] = 0;
//    CurrentChild[k] = 0;
//  }
//
//  do
//  {
//    if (ActiveLevel)
//      do
//      {
//        ActiveCell = (TGridCell *)(ActiveCell->GetParent());
//        N_Children[ActiveLevel] = 0;
//        CurrentChild[ActiveLevel--] = 0;
//      } while (ActiveLevel &&
//          N_Children[ActiveLevel] == CurrentChild[ActiveLevel]);
//
//    if (!ActiveLevel &&
//          N_Children[ActiveLevel] == CurrentChild[ActiveLevel])
//    {
//      if (ActiveRootCell < N_RootCells)
//      {
//        ActiveCell = (TGridCell *)(CellTree[ActiveRootCell++]);
//        N_Children[ActiveLevel] = ActiveCell->GetN_Children();
//        CurrentChild[ActiveLevel] = 0;
//      }
//      else
//        return;
//    }
//
//    while(N_Children[ActiveLevel])
//    {
//      ActiveCell = (TGridCell *)(ActiveCell->GetChild(CurrentChild[ActiveLevel]++));
//      N_Children[++ActiveLevel] = ActiveCell->GetN_Children();
//    }
//
//    for(k=0;k<MaxN_Levels;k++)
//      cout << "xyz" << k << " " << N_Children[k] << " " << CurrentChild[k] << endl;
//
//    Vertices = ActiveCell->GetVertices();
//    N_Edges = ActiveCell->GetN_Edges();
//    for(j=0;j<N_Edges;j++)
//    {
//      joint = ActiveCell->GetJoint(j);
//      if(joint->GetType() == BoundaryEdge)
//      {
//        ((TBoundEdge *)joint)->UpdateParameters(Vertices[j],
//                          Vertices[(j+1)%N_Edges]);
//      }
//      if(joint->GetType() == InterfaceJoint)
//      {
//        ((TInterfaceJoint *)joint)->UpdateParameters(Vertices[j],
//                          Vertices[(j+1)%N_Edges]);
//      }
//    }
//
//    // if Cell is a quadrilateral
//    if(N_Edges == 4)
//    {
//      if(ActiveCell->ExistChildren())
//        ActiveCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES+((TQuadrangle *)
//             TDatabase::RefDescDB[Quadrangle]->GetShapeDesc())->
//             CheckQuad(Vertices)]);
//      else
//        ActiveCell->SetRefDesc(TDatabase::RefDescDB[((TQuadrangle *)
//             TDatabase::RefDescDB[Quadrangle]->GetShapeDesc())->
//             CheckQuad(Vertices)]);
//    }
//
//  } while(ActiveCell);
#endif

/*
#include <stdio.h>

void SaveTri(struct triangulateio &outt)
{
  FILE *f;
  int i;

  f = fopen ("ch.node", "wt");
  fprintf(f, "%d   2   0  1\n", outt.numberofpoints);
  for (i=0; i<=outt.numberofpoints-1; i++)
    fprintf(f,"%6d  %10.4f  %10.4f  %6d\n",i, outt.pointlist[2*i],
                   outt.pointlist[2*i+1], outt.pointmarkerlist[i]);
  fclose(f);

  f = fopen ("ch.ele", "wt");
  fprintf(f, "%d   3\n", outt.numberoftriangles);
  for (i=0; i<=outt.numberoftriangles-1; i++)
    fprintf(f,"%6d  %6d  %6d  %6d\n", i, outt.trianglelist[3*i], 
             outt.trianglelist[3*i+1], outt.trianglelist[3*i+2]);
  fclose(f);
 
  f = fopen ("ch.neigh", "wt");
  fprintf(f, "%d   3\n", outt.numberoftriangles);
  for (i=0; i<=outt.numberoftriangles-1; i++)
    fprintf(f,"%6d  %6d  %6d  %6d\n", i, outt.neighborlist[3*i], 
             outt.neighborlist[3*i+1], outt.neighborlist[3*i+2]);
  fclose(f);
 
  f = fopen ("ch.edge", "wt");
  fprintf(f, "%d   1\n", outt.numberofedges);
  for (i=0; i<=outt.numberofedges-1; i++)
    fprintf(f,"%6d  %6d  %6d  %6d\n", i, outt.edgelist[2*i], 
             outt.edgelist[2*i+1], outt.edgemarkerlist[i]);
  fclose(f);
}
*/


int TDomain::RefineallxDirection()
{
	TBaseCell *CurrCell;
	int info;
	static int s = 0;  // changed from 0 to 1, if refine in y direction  

	RefLevel ++;
	TDatabase::IteratorDB[It_Finest]->Init(0);

	//loop over all cells

	while((CurrCell=TDatabase::IteratorDB[It_Finest]->Next(info)))
	{
	  switch(CurrCell->GetType())
	  {
		case Quadrangle:
	  	case Parallelogram:
		case Rectangle:
			CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES+QuadBis0+s]);
			CurrCell->Refine(RefLevel);
		break;
	  default:
            cerr << "RefineallxDirection only for Quad cells" << endl;
             exit(-1);
          break;	
	  }
	}

	s = 1;
	
	return 0;
}

bool TDomain::isExtendedGEO(const char* GEO)
  {
      bool isXgeo{false};
      // check if input file is an extended geo file (.xGEO)
      int nn=0;
      while (GEO[nn] != 0)
      {
        ++nn;
      }

      //check if we found the correct place in the char arary
      if(GEO[nn-3] != 'G' || GEO[nn-2] != 'E' || GEO[nn-1] != 'O')
      {
        ErrThrow("Incorrect read-in of .(x)GEO-filename! (Make sure the "
            "filename ends on '.GEO' or '.xGEO')" );
      }

      if (GEO[nn-4]=='x')
      {
	Output::print<2>(" *** reading xGEO file (with physical references) ***");
        isXgeo = true;
      }
      return isXgeo;
  }



#ifdef __3D__
// int TDomain::Tetgen(const char *FileName)
// {
//   TTetGenMeshLoader MeshLoader(FileName);
//   MeshLoader.Generate(this);
// 
//   // initialize iterators
//   TDatabase::IteratorDB[It_EQ]->SetParam(this);
//   TDatabase::IteratorDB[It_LE]->SetParam(this);
//   TDatabase::IteratorDB[It_Finest]->SetParam(this);
//   TDatabase::IteratorDB[It_Between]->SetParam(this);
//   TDatabase::IteratorDB[It_OCAF]->SetParam(this);
// 
//   return 0;
// }

/** @brief added by sashi */
int TDomain::GenerateEdgeInfo()
{
  int i, ii, j, k, m, I, n_cells, info, N_RootVertices, *PointNeighb, MaxCpV = 0, EMaxLen;
  int level=0, N, NVert_All, *NumberVertex, *VertexNumbers, N_Edges, N_FaceEdges, N_FaceVert;
  int *RootVer_Loc, *NeibRootVer_Loc, a, b, len1, len2, cell_a, N_Neibs, *Neibs, *Neib_EdgeNo;
  int NeibN_Edges, *VertBeginIndex, Neiba, Neibb, BoundEdgeMarker[MAXN_EDGES], N_Faces, MaxLen;
  const int *EdgeVertex, *TmpFV, *TmpLen,  *TmpFE, *ETmpLen;;
 
  TBaseCell **cells, *CurrCell, *NeibCell, **NeibCells;
  Iterators it=It_Finest;
  TVertex **Vertices_All, *Last, *Current;
  TShapeDesc *ShapeDesc, *NeibShapeDesc;
  TEdge *edge, *Neibedge;
  TJoint *joint;

#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
#endif

  NVert_All=0;
  n_cells=0;
  TDatabase::IteratorDB[it]->Init(level);
  while((CurrCell=TDatabase::IteratorDB[it]->Next(info)))
  {
    NVert_All +=CurrCell->GetN_Vertices();
    n_cells++;
  }

  cells=new TBaseCell*[n_cells];
  VertBeginIndex = new int[n_cells+1];
  NumberVertex = new int[NVert_All];
  VertexNumbers = new int[NVert_All];
  Vertices_All = new TVertex*[NVert_All];

  NVert_All=0;
  n_cells=0;
  VertBeginIndex[0] = 0;
  TDatabase::IteratorDB[it]->Init(level);
  while((CurrCell = TDatabase::IteratorDB[it]->Next(info)))
  {
    cells[n_cells]=CurrCell;
    N=CurrCell->GetN_Vertices();

    for(i=0;i<N;i++)
    {
      Vertices_All[NVert_All] = CurrCell->GetVertex(i);
      NVert_All++;
    }

    n_cells++;
    VertBeginIndex[n_cells] = NVert_All;

  } // while
// NVert_All--;
  // sort the Vertices array
  Sort(Vertices_All, NVert_All);
// NVert_All++;
  Last=NULL;
  N_RootVertices=-1;
  for(i=0;i<NVert_All;i++)
   {
    if((Current=Vertices_All[i])!=Last)
    {
      N_RootVertices++;
      Last=Current;
    }
    NumberVertex[i]=N_RootVertices;
   }
  N_RootVertices++;

//   cout << "N_RootVertices " << N_RootVertices << " NVert_All " << NVert_All << endl;

  m=0;
  for(ii=0;ii<n_cells;ii++)
   {
    CurrCell = cells[ii];
    N=CurrCell->GetN_Vertices();
    for(i=0;i<N;i++)
    {
      Current= CurrCell->GetVertex(i);
      I=GetIndex(Vertices_All, NVert_All, Current);
      VertexNumbers[m]=NumberVertex[I];
//       cout << VertexNumbers[m] << endl;
      m++;
    } // endfor j
   } //while

  /** find max No. cells met any vertex in the collection **/
//    VertBound = new int[N_RootVertices];
   PointNeighb = new int[N_RootVertices];
   memset(PointNeighb, 0, N_RootVertices*SizeOfInt);
//    memset(VertBound, 0, N_RootVertices*SizeOfInt);

  for(i=0;i<NVert_All;i++)
    PointNeighb[VertexNumbers[i]]++;

  for(i=0;i<N_RootVertices;i++)
    if (PointNeighb[i] > MaxCpV) MaxCpV = PointNeighb[i];

   //printf("Max number of cells per vertex %d \n", MaxCpV);
   /** No. of edge neibs will not exceed MaxCpV */ 
   Neibs = new int[MaxCpV];
   Neib_EdgeNo = new int[MaxCpV];
   NeibCells = new TBaseCell*[MaxCpV];

  /** PointNeighb's first column contains number of neib cells met with each vertex **/
  /** further columns contain the cell numbers met with this vertex **/
   MaxCpV++;

   delete [] PointNeighb;
   PointNeighb = new int[N_RootVertices*MaxCpV];
   memset(PointNeighb, 0, (N_RootVertices*MaxCpV)*SizeOfInt);

  NVert_All=0;
  for(ii=0;ii<n_cells;ii++)
   {
    CurrCell = cells[ii];
    N=CurrCell->GetN_Vertices();
//     CurrCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

    for(i=0;i<N;i++)
     {
      I = VertexNumbers[NVert_All] *MaxCpV;
      PointNeighb[I]++;
      PointNeighb[I + PointNeighb[I] ] = ii;
      NVert_All ++;
     }

   } // for(ii=0;ii<n_cells;ii++)

// exit(0);

//   for(i=0;i<N_RootVertices;i++)
//     if ( VertBound[i]==0) printf("Vert %d VertBound %d \n", i, VertBound[i]);
//     //printf("Number of cells in vertex %d \n",  PointNeighb[i*MaxCpV]);
// 
// exit(0);
//   NVert_All=0;
  for(ii=0;ii<n_cells;ii++)
   {
    CurrCell = cells[ii];
    N=CurrCell->GetN_Vertices();
    N_Edges=CurrCell->GetN_Edges();
    N_Faces=CurrCell->GetN_Joints();
    
    ShapeDesc= CurrCell->GetShapeDesc();
    ShapeDesc->GetEdgeVertex(EdgeVertex);
    ShapeDesc->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    
    RootVer_Loc = VertexNumbers+VertBeginIndex[ii];
    
    
    /** set the marker for Edges */
    for(i=0;i<MAXN_EDGES;i++)
     BoundEdgeMarker[i] = 0;;
    
    /** set bound and isobound Edges */
    for(i=0;i<N_Faces;i++)
     {
      joint=CurrCell->GetJoint(i);

      if(joint->GetType() == BoundaryFace)
       {         
        N_FaceEdges = ETmpLen[i];
        for(j=0;j<N_FaceEdges;j++) 
         {
          /** not modified if it is isobound edge */
          if (BoundEdgeMarker[TmpFE[i*EMaxLen+j] ] == 0)
           BoundEdgeMarker[TmpFE[i*EMaxLen+j] ] = 1;
         }
         
         /** set all vertices on the bound face as bound vert */
         N_FaceVert = TmpLen[i];        
         for(j=0;j<N_FaceVert;j++) 
          CurrCell->GetVertex(TmpFV[i*MaxLen + j])->SetAsBoundVert();
          
       } // if(joint->GetType() == Boundar
      else if(joint->GetType() == IsoBoundFace)
       {
        N_FaceEdges = ETmpLen[i];
        for(j=0;j<N_FaceEdges;j++)
         BoundEdgeMarker[TmpFE[i*EMaxLen+j] ] = 2;

        /** set all vertices on the bound face as bound vert */
        N_FaceVert = TmpLen[i];        
        for(j=0;j<N_FaceVert;j++) 
         CurrCell->GetVertex(TmpFV[i*MaxLen + j])->SetAsBoundVert();
       } // 
     } // for(i=0;i<N_    
    
    for(i=0;i<N_Edges;i++)
     {
      N_Neibs = 0;
      edge = CurrCell->GetEdge(i);
      
      if(edge==NULL) // edge not yet created
       {
        a = RootVer_Loc[EdgeVertex[2*i]];
        b = RootVer_Loc[EdgeVertex[2*i+1]];

        len1 = PointNeighb[a*MaxCpV]; // No. cells met with vert a
        len2 = PointNeighb[b*MaxCpV]; // No. cells met with vert b

        /** find No. cells having this edge */
        for(j=1;j<=len1;j++)
         {
          cell_a = PointNeighb[a*MaxCpV + j];
// 	  if(rank==0 && cell_a==16)
// 	  cout  <endl;
//            cout <<    " cell_a " << cell_a<<endl;

          for(k=1;k<=len2;k++)
           {

            if(cell_a == PointNeighb[b*MaxCpV +k] )
             {  
              Neibs[N_Neibs] = cell_a;
              N_Neibs ++;
              break;
             }
           } // for(k=1;k<=
         } // for(j=1;j<

        /** find the local index of the edge in all cells */
        for(j=0;j<N_Neibs;j++)
         {

          if(Neibs[j]==ii)
           {
             Neib_EdgeNo[j] = i;
           }
          else 
           {
            NeibCell = cells[Neibs[j]];
            NeibShapeDesc= NeibCell->GetShapeDesc();
            NeibN_Edges=NeibCell->GetN_Edges();
            NeibShapeDesc->GetEdgeVertex(EdgeVertex);
            NeibRootVer_Loc = VertexNumbers+VertBeginIndex[ Neibs[j]];
            //cout <<    "  Neib " << Neibs[j]<<endl;

            for(k=0;k<NeibN_Edges;k++)
             {
              Neibedge = NeibCell->GetEdge(k);

              if(Neibedge==NULL) // edge not yet created
               {
                 Neiba = NeibRootVer_Loc[EdgeVertex[2*k]];
                 Neibb = NeibRootVer_Loc[EdgeVertex[2*k+1]];

                 if( (Neiba==a &&  Neibb==b) || (Neiba==b &&  Neibb==a) )
                  {
                   //cout <<    "  Neib edge found "  <<endl;
                   Neib_EdgeNo[j] = k;
                   break;
                  }
               }
             } //for(k=0;k<NeibN_Edges;k

            if(k==NeibN_Edges)
             {
              cout << "Error could not find edge index in the Neib cell !" <<endl;
              exit(0);
             }
            }
         } // for(j=0;j<N_Neib

       /**generate the edge and set in all cells */
       for(j=0;j<N_Neibs;j++)
        {
         NeibCells[j] = cells[Neibs[j]];
         NeibCells[j]->SetEdge(Neib_EdgeNo[j], edge);
        }

//   	  if(rank==0 && ii==35 && i==1)      
// 	  {
//            for(j=0;j<N_Neibs;j++)
// 	     printf("rank %d j %d NeibCells %d BoundEdgeMarker %d \n", rank, j, Neibs[j], BoundEdgeMarker[i]); 
// 	     
// 	  }
//         
        
        
       if(BoundEdgeMarker[i]==2 )
        {
         edge = new TIsoEdge3D(N_Neibs, NeibCells);
        }
       else if(BoundEdgeMarker[i]==1)
        {
         edge = new TBDEdge3D(N_Neibs, NeibCells);
        }
       else
        {
         edge = new TInnerEdge(N_Neibs, NeibCells);
        }                
       } //if(edge==NULL)
//      else
//       {
//          cout <<  ii <<  " Edge already set " << i  <<endl;
// // exit(0);
//       }

       for(j=0;j<N_Neibs;j++)
         NeibCells[j]->SetEdge(Neib_EdgeNo[j], edge);

     } // for(i=0;i<N_Edges;i++)
   } //  for(ii=0;ii<n_cel

   delete [] Neibs;
   delete [] Neib_EdgeNo;
   delete [] NeibCells;
   delete [] cells;
   delete [] VertBeginIndex;
   delete [] NumberVertex;
   delete [] VertexNumbers;
   delete [] Vertices_All;
//    delete [] VertBound;
   delete [] PointNeighb;

#ifdef _MPI
   if(rank==0)
#endif
    Output::info<5>("Domain.C","3D Mesh Edges Generated ");

  return 0;
}


#endif

size_t TDomain::get_n_initial_refinement_steps() const
{
  return db["refinement_n_initial_steps"];
}

size_t TDomain::get_max_n_adaptive_steps() const
{
  return db["refinement_max_n_adaptive_steps"];
}

void TDomain::print_info(std::string name) const
{
#ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &my_rank);
  MPI_Comm_size(TDatabase::ParamDB->Comm, &size);

  int sbuf_own = N_OwnCells;
  int sbuf_halo = N_RootCells - N_OwnCells;

  std::vector<int> ns_own_cells(size,0);
  std::vector<int> ns_halo_cells(size,0);

  {
    MPI_Gather(
      &sbuf_own, 1, MPI_INT,            //send
      &ns_own_cells.at(0), 1, MPI_INT,  //receive
      0, MPI_COMM_WORLD);               //control
    MPI_Gather(
      &sbuf_halo, 1, MPI_INT,            //send
      &ns_halo_cells.at(0), 1, MPI_INT,  //receive
      0, MPI_COMM_WORLD);               //control
  }
  if(my_rank == 0)
  {
    Output::stat("Domain", name);
    size_t sum_cells_total = 0;
    for(int i =0; i < size ;++i)
    {
      Output::dash("Process", i, "\t n_own_cells: ", ns_own_cells.at(i),
                    "\t n_halo_cells: ", ns_halo_cells.at(i));
      sum_cells_total += ns_own_cells.at(i);
    }
    Output::dash("Total number of cells: ", sum_cells_total);

  }

#else
  Output::stat("Domain ", name);
  Output::dash("No domain statistics printout in non-MPI case so far.");
#endif
}


void determine_n_refinement_steps_multigrid(
  const std::string& multigrid_type,
  int n_multigrid_levels,
  int n_initial_refinement_steps,
  int& n_ref_before, int& n_ref_after)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank=0;
#endif
  //Check if input multigrid type makes sense.
  if(multigrid_type!=std::string("standard") && multigrid_type!=std::string("mdml"))
    ErrThrow("Unknown multigrid type ", multigrid_type,". Use 'standard' or 'mdml'.");
  //Check if input numbers make sense
  if(n_multigrid_levels > n_initial_refinement_steps + 1)
    ErrThrow("You requested ", n_multigrid_levels, " geometric levels for "
        "multigrid, but only ", n_initial_refinement_steps + 1, " grid levels "
        "will be present after refining ", n_initial_refinement_steps, " times. "
        " Choose a smaller multigrid_n_levels.");

  //split the number of refinement steps
  n_ref_after =  n_multigrid_levels - 1;
  n_ref_before =  n_initial_refinement_steps - n_ref_after;

  //print information
  if(my_rank == 0)
  {
    Output::info("REFINEMENT FOR MULTIGRID","You are using ",multigrid_type,
                 " multigrid with ",n_multigrid_levels," geometric levels. ");
#ifdef _MPI
    Output::dash("Number of refinement steps before domain partitioning: ", n_ref_before);
    Output::dash("Number of refinement steps after domain partitioning: ", n_ref_after);
    Output::dash("Total number of refinement steps: ", n_initial_refinement_steps);
#else
    Output::dash("Number of refinement steps before picking grids: ", n_ref_before);
    Output::dash("Number of refinement steps after picking grids: ", n_ref_after);
    Output::dash("Total number of refinement steps: ", n_initial_refinement_steps);
#endif
  }
}

std::list<TCollection* > TDomain::refine_and_get_hierarchy_of_collections(
    const ParameterDatabase& parmoon_db
#ifdef _MPI
    , int& maxSubDomainPerDof
#endif
    )
{
#ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  // Split the number of refinement steps - see doc of the method
  // TODO Removing this strange construction is a TODO
  int n_ref_before = this->get_n_initial_refinement_steps();
  int n_ref_after = 0;
  if(parmoon_db.contains("preconditioner")
     && parmoon_db["preconditioner"].is(std::string("multigrid"))
     && parmoon_db.contains("solver_type")
     && parmoon_db["solver_type"].is(std::string("iterative")))
  {
    determine_n_refinement_steps_multigrid(
      parmoon_db["multigrid_type"],
      parmoon_db["multigrid_n_levels"],
      this->get_n_initial_refinement_steps(),
      n_ref_before, n_ref_after);
  }

  for(int i = 0; i < n_ref_before; i++)
    this->RegRefineAll();

#ifdef _MPI
  // Partition the by now finest grid using Metis and distribute among processes.

  // 1st step: Analyse interfaces and create edge objects,.
  this->GenerateEdgeInfo();

  // 2nd step: Call the mesh partitioning.

  int maxCellsPerVertex;
  //do the actual partitioning, and examine the return value
  if ( Partition_Mesh3D(MPI_COMM_WORLD, this, maxCellsPerVertex) == 1)
  {
    ErrThrow("Partitioning did not succeed.");
  }

  // 3rd step: Generate edge info anew
  //(since distributing changed the domain).
  this->GenerateEdgeInfo();

  // calculate largest possible number of processes which share one dof
  maxSubDomainPerDof = MIN(maxCellsPerVertex, size);

#endif

  std::list<TCollection* > gridCollections;
  gridCollections.push_front(this->GetCollection(It_Finest, 0));

  //this is only relevant for multigrid
  for(int level=0; level <  n_ref_after; ++level)
  {
    this->RegRefineAll();
#ifdef _MPI
    this->GenerateEdgeInfo();  // has to be called anew after every refinement step
    Domain_Crop(MPI_COMM_WORLD, this); // remove unwanted cells in the halo after refinement
#endif
    // Grab collection.
    gridCollections.push_front(this->GetCollection(It_Finest, 0));
  }

  return gridCollections;
}


#ifdef __3D__


void TDomain::GenerateFromMesh(Mesh& m)
{
  // create inner faces (if needed)
  m.createInnerFaces();
  
  /** 
      compute number boundary faces:
      This step at the moment is necessary to initialize some useful arrays
      (e.g. boundary faces)
  */
  m.computeNumberOfBoundaryFaces();
  
  // build the boundary parts
  this->buildBoundary(m);

  // build the mesh
  this->buildParMooNMesh(m);
  
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);  
}

// note this function needs:
// list of triangles (faces)
// number of boundary faces
// adjacency (1,-1 for boundary faces)
// list of points
///@todo write part of the function inside the Mesh class, here set only the BdParts
void TDomain::buildBoundary(Mesh& m)
{
  // we consider only 1 boundary part
  ///@todo do we need to take into account the general case?
  this->N_BoundParts = 1;
  // vector of boundary parts
  this->BdParts = new TBoundPart* [this->N_BoundParts];

  // the number of boundary components is computed from the adjacency of
  // tetgen mesh. See TTetGenMeshLoaden::CreateAdjacency()
  this->N_BoundComps = m.n_boundary_faces;
  meshBoundComps.resize(this->N_BoundComps); 
  Output::print("TDomain::buildBoundary() - N_BoundComps: ", this->N_BoundComps);
  this->BdParts[0] = new TBoundPart(this->N_BoundComps);
  
  // StartBdCompID[i] gives the iindex of the element in BdParts
  // where the BdComp i starts. The last element is equal to to N_BoundParts
  this->StartBdCompID = new int [this->N_BoundParts+1];

  // we now consider only a boundary components,
  // i.e. StartBdCompID[0] = 0, StartBdCompID[1] = N_BoundComps
  ///@todo check how we can extend this to multiple boundary markers
  this->StartBdCompID[0] = 0; 
  this->StartBdCompID[1] = this->N_BoundComps;


  // in order to describe the boundary, each boundary triangle is considered
  // to be a boundary component of the (single) boundary part.
  // Hence,
  // (*) first, we check if the triangle is on the boundary
  // (*) second, given the vertices, we define a TBdPlane describing the triangle
  int n_trifaces = m.triangle.size();

  int counter=0; // count the number of boundary faces found

  if(m.boundaryFacesMarker.size() == 0)
  {
    m.boundaryFacesMarker.resize(n_trifaces);
  }
  for(int i=0;i<n_trifaces;++i)
  {
    if(m.faceToTetra[i][0] == -1 || m.faceToTetra[i][1] == -1)
    {
      // the face is on the boundary
      meshBoundComps.at(counter) = new TBdPlane(counter,1),m.triangle[i].reference);

      double p[3], a[3], b[3], n[3];
      ///@attention the lists of nodes starts from 1 (not from 0)
      p[0] = m.vertex[ m.triangle[i].nodes[0] - 1].x;
      a[0] = m.vertex[ m.triangle[i].nodes[1] - 1].x - p[0];
      b[0] = m.vertex[ m.triangle[i].nodes[2] - 1].x - p[0];

      p[1] = m.vertex[ m.triangle[i].nodes[0] - 1].y;
      a[1] = m.vertex[ m.triangle[i].nodes[1] - 1].y - p[1];
      b[1] = m.vertex[ m.triangle[i].nodes[2] - 1].y - p[1];

      p[2] = m.vertex[ m.triangle[i].nodes[0] - 1].z;
      a[2] = m.vertex[ m.triangle[i].nodes[1] - 1].z - p[2];
      b[2] = m.vertex[ m.triangle[i].nodes[2] - 1].z - p[2];
      
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
      // the id of the boundary face is taken as the attribute of the triangle
      // note: later, trifacemarkerlist=0 is used to identify inner faces (see below)
      m.boundaryFacesMarker[i] = counter;
      Output::print<4>("TDomain::buildBoundary() Triangle ",i,
		       " boundaryMarker = ",m.boundaryFacesMarker[i],
		       " attribute = ",m.triangle[i].reference);
    }
    else // not a boundary face
      m.boundaryFacesMarker[i] = 0;
  }

  // check if the number of boundary components is consistent
  assert(counter==N_BoundComps);
  // initialize ParMooN BdParts and the BdComp
  for(int i=0; i<this->N_BoundComps; i++)
    this->BdParts[0]->SetBdComp(i, meshBoundComps[i]);
  

}


void TDomain::buildParMooNMesh(Mesh& m)
{
  // create a vector<TVertex*> from the list of pointa
  this->setVertices(m);
  // allocate memory for CellTree and set the vertices
  this->allocRootCells(m);  
  // create joints and set joints type
  this->distributeJoints(m);
}

/**
 * @brief set the coordinates of the vertices
 *
 * This functions reads the vertices coordinates from meshTetGenOut
 * and creates TVertex*, stored in the vector meshVertices
 * Moreover, bounds for the domain are computed
 */
void TDomain::setVertices(Mesh& m)
{
 
  double xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0;
  meshVertices.resize(m.vertex.size()); 
  Output::print("TDomain::setVertices() total number of vertices: ",  meshVertices.size());

  for(unsigned int i=0;i<meshVertices.size();++i)
  {
    double x = m.vertex[i].x;
    double y = m.vertex[i].y;
    double z = m.vertex[i].z;

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

}

void TDomain::allocRootCells(Mesh& m)
{
  TVertex *Vertex;
  TMacroCell *Cell;
  // set the number of total cells and allocate the cell array
  this->N_RootCells = m.tetra.size();
  Output::print<3>("TDomain::allocRootCells() number of tetrahedra: ",
		   this->N_RootCells);

  this->CellTree = new TBaseCell* [this->N_RootCells];
  // set the vertices of each cell, reading from meshVertices
  for(int i=0;i<this->N_RootCells;++i)
  {
    Cell = new TMacroCell (TDatabase::RefDescDB[Tetrahedron], 0);
    this->CellTree[i] = Cell;
    
    for(int j=0;j<4;++j)
    {
      Vertex = meshVertices.at(m.tetra[i].nodes[j]-1);
      Cell->SetVertex(j, Vertex);
    }
    Cell->SetPhase_ID((int) m.tetra[i].reference);
    
  }

}

///@attention the functions hashTriFaces(), CreateAdjacency() must have been called before
void TDomain::distributeJoints(Mesh& m)
{
  // allocate the joints array
  meshJoints.resize(m.triangle.size());
  Output::print("number of joints: ", meshJoints.size());

  // set the faces for each joint
  for(unsigned int i=0; i<meshJoints.size(); ++i)
  {
    // find element that contain this face
    if(m.boundaryFacesMarker[i] == 0)// inner joints
    {
      int left = m.faceToTetra[i][0];
      int right = m.faceToTetra[i][1];
      Output::print<4>("TDomain::distributeJoints() Joint ", i ," Tetra: ", left, " and ", right);
      meshJoints.at(i) = new TJointEqN (CellTree[left], CellTree[right]);
    }
    else // boundary joints
    {
      Output::print<5>("Joint ", i, " is a boundary joint");
      int bdcomp = m.boundaryFacesMarker[i] - 1;
      TBoundComp3D* BoundComp = meshBoundComps[bdcomp];
      meshJoints.at(i)=new TBoundFace (BoundComp);
    }
  }

  
  for(unsigned int i=0;i<m.tetra.size();++i)
  {
    const int *TmpFV, *TmpLen;
    int MaxLen;
    
    TShapeDesc *ShapeDesc = CellTree[i]->GetShapeDesc();
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);

    for(int j=0;j<4;++j)
    {
      ///@attention make sure that the vector meshTrifaceHash has been created 
      int triface = m.findTriFace(m.tetra[i].nodes[TmpFV[j*MaxLen  ]],
				  m.tetra[i].nodes[TmpFV[j*MaxLen+1]],
				  m.tetra[i].nodes[TmpFV[j*MaxLen+2]]);

      // check that a face has been found
      assert (triface != -1);

      CellTree[i]->SetJoint(j, meshJoints.at(triface));

      // correct parameters for boundary faces
      if(meshJoints.at(triface)->GetType() == BoundaryFace)
      {

        TBoundFace *BoundFace = (TBoundFace*) meshJoints.at(triface);
        TBoundComp3D *BoundComp = BoundFace->GetBoundComp();
        double param1[4], param2[4];

        for(int k=0;k<TmpLen[j];++k)
        {
	  double x, y, z, t, s;
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
  for(unsigned int i=0;i<meshJoints.size();++i)
  {
    meshJoints[i]->SetMapType();
  }

}

void TDomain::GenerateFromTetgen(TTetGenMeshLoader& tgml)
{
  // check file extensions, call tetgen and generate the volume mesh,
  // the has list and the adjacency between faces
  // note: the mesh is contained in the  (tetgenio) meshTetGenOut object
  tgml.GenerateMesh();

  /*
  // build the boundary parts
  this->buildBoundary(tgml);

  // build the mesh
  this->buildParMooNMesh(tgml);
  
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);  
  */
}

// =======================================
// note this function needs:
// list of triangles (faces)
// number of boundary faces
// adjacency (1,-1 for boundary faces)
// list of points
// Can we use push_back to fill the meshBoundComps vector?
void TDomain::buildBoundary(TTetGenMeshLoader& tgml)
{
  // we consider only 1 boundary part
  ///@todo do we need to take into account the general case?
  this->N_BoundParts = 1;
  // vector of boundary parts
  this->BdParts = new TBoundPart* [this->N_BoundParts];

  // the number of boundary components is computed from the adjacency of
  // tetgen mesh. See TTetGenMeshLoaden::CreateAdjacency()
  this->N_BoundComps = tgml.nBoundaryComponents;
  meshBoundComps.resize(this->N_BoundComps); 
  Output::print("TDomain::buildBoundary() - N_BoundComps: ", this->N_BoundComps);
  Output::print("here");
  this->BdParts[0] = new TBoundPart(this->N_BoundComps);
  Output::print("here");
  // StartBdCompID[i] gives the iindex of the element in BdParts
  // where the BdComp i starts. The last element is equal to to N_BoundParts
  this->StartBdCompID = new int [this->N_BoundParts+1];

  // we now consider only a boundary components,
  // i.e. StartBdCompID[0] = 0, StartBdCompID[1] = N_BoundComps
  ///@todo check how we can extend this to multiple boundary markers
  this->StartBdCompID[0] = 0; 
  this->StartBdCompID[1] = this->N_BoundComps;


  // in order to describe the boundary, each boundary triangle is considered
  // to be a boundary component of the (single) boundary part.
  // Hence,
  // (*) first, we check if the triangle is on the boundary
  // (*) second, given the vertices, we define a TBdPlane describing the triangle
  int n_trifaces = tgml.meshTetGenOut.numberoftrifaces;
  int *trifaces = tgml.meshTetGenOut.trifacelist;
  int counter=0; // count the number of boundary faces found

  ///@attention misuse of trifacemarkerlist
  //  0 = boundary faces, i>0 = i-th boundary face
  ///@todo this should be changed in order to allow physical boundary markers
  if(tgml.meshTetGenOut.trifacemarkerlist == NULL)
    tgml.meshTetGenOut.trifacemarkerlist = new int [n_trifaces];

  for(int i=0;i<n_trifaces;++i)
  {
    if(tgml.meshTetGenOut.adjtetlist[2*i+1] == -1 ||
      tgml.meshTetGenOut.adjtetlist[2*i  ] == -1)
    {
      // the face is on the boundary
      meshBoundComps.at(counter) = new TBdPlane(counter);
      double p[3], a[3], b[3], n[3];
      p[0] = tgml.meshTetGenOut.pointlist[3*trifaces[3*i  ]  ];
      a[0] = tgml.meshTetGenOut.pointlist[3*trifaces[3*i+1]  ] - p[0];
      b[0] = tgml.meshTetGenOut.pointlist[3*trifaces[3*i+2]  ] - p[0];
      p[1] = tgml.meshTetGenOut.pointlist[3*trifaces[3*i  ]+1];
      a[1] = tgml.meshTetGenOut.pointlist[3*trifaces[3*i+1]+1] - p[1];
      b[1] = tgml.meshTetGenOut.pointlist[3*trifaces[3*i+2]+1] - p[1];
      p[2] = tgml.meshTetGenOut.pointlist[3*trifaces[3*i  ]+2];
      a[2] = tgml.meshTetGenOut.pointlist[3*trifaces[3*i+1]+2] - p[2];
      b[2] = tgml.meshTetGenOut.pointlist[3*trifaces[3*i+2]+2] - p[2];
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
      Output::print(" counter:",counter);
      // the id of the boundary face is taken as the attribute of the triangle
      // note: later, trifacemarkerlist=0 is used to identify inner faces (see below)
      tgml.meshTetGenOut.trifacemarkerlist[i] = counter;
    }
    else // not a boundary face
      tgml.meshTetGenOut.trifacemarkerlist[i] = 0;
  }

  // check if the number of boundary components is consistent
  assert(counter==N_BoundComps);
  // initialize ParMooN BdParts and the BdComp
  for(int i=0; i<this->N_BoundComps; i++)
    this->BdParts[0]->SetBdComp(i, meshBoundComps[i]);
  

}

void TDomain::buildParMooNMesh(TTetGenMeshLoader& tgml)
{
  // create a vector<TVertex*> from the list of pointa
  this->setVertices(tgml);
  // allocate memory for CellTree and set the vertices
  this->allocRootCells(tgml);  
  // create joints and set joints type
  this->distributeJoints(tgml);
}

/**
 * @brief set the coordinates of the vertices
 *
 * This functions reads the vertices coordinates from meshTetGenOut
 * and creates TVertex*, stored in the vector meshVertices
 * Moreover, bounds for the domain are computed
 * 
 * @todo should this function be split in two parts:
 * e.g. setVertices and computeDomainBounds?
 */
void TDomain::setVertices(TTetGenMeshLoader& tgml)
{
  double x, y, z;
  double xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0;
  meshVertices.resize(tgml.meshTetGenOut.numberofpoints); 
  
  for(int i=0;i<tgml.meshTetGenOut.numberofpoints;++i)
  {
    x = tgml.meshTetGenOut.pointlist[3*i  ];
    y = tgml.meshTetGenOut.pointlist[3*i+1];
    z = tgml.meshTetGenOut.pointlist[3*i+2];

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

  Output::print("number of vertices: ", tgml.meshTetGenOut.numberofpoints);
}

void TDomain::allocRootCells(TTetGenMeshLoader& tgml)
{
  TVertex *Vertex;
  TMacroCell *Cell;
  // set the number of total cells and allocate the cell array
  this->N_RootCells = tgml.meshTetGenOut.numberoftetrahedra;
  Output::print<3>("TDomain::allocRootCells() number of tetrahedra: ",
		   tgml.meshTetGenOut.numberoftetrahedra);

  this->CellTree = new TBaseCell* [this->N_RootCells];
  // set the vertices of each cell, reading from meshVertices
  for(int i=0;i<this->N_RootCells;++i)
  {
    Cell = new TMacroCell (TDatabase::RefDescDB[Tetrahedron], 0);
    this->CellTree[i] = Cell;
    
    for(int j=0;j<4;++j)
    {
      Vertex = meshVertices.at(tgml.meshTetGenOut.tetrahedronlist[4*i+j]);
      Cell->SetVertex(j, Vertex);
    }
    Cell->SetPhase_ID((int) tgml.meshTetGenOut.tetrahedronattributelist[i]);
    
  }

}

///@attention the functions hashTriFaces(), CreateAdjacency() must have been called before
void TDomain::distributeJoints(TTetGenMeshLoader& tgml)
{
  // size of meshjoint
  meshJoints.resize(tgml.meshTetGenOut.numberoftrifaces);
  // search face which belongs to current bdComp
  for(int i=0;i<tgml.meshTetGenOut.numberoftrifaces;++i)
  {
    // find element that contain this face
    if(tgml.meshTetGenOut.trifacemarkerlist[i] == 0)// inner joints
    {
      int left = tgml.meshTetGenOut.adjtetlist[2*i];
      int right = tgml.meshTetGenOut.adjtetlist[2*i+1];
      
      meshJoints.at(i) = new TJointEqN (CellTree[left], CellTree[right]);
    }
    else // boundary joints
    {
      Output::print<5>("boundary joint");
      int bdcomp = tgml.meshTetGenOut.trifacemarkerlist[i] - 1;
      TBoundComp3D* BoundComp = meshBoundComps[bdcomp];
      meshJoints.at(i)=new TBoundFace (BoundComp);
    }
  }

  Output::print("number of joints: ", tgml.meshTetGenOut.numberoftrifaces);
  
  for(int i=0;i<tgml.meshTetGenOut.numberoftetrahedra;++i)
  {
    const int *TmpFV, *TmpLen;
    int MaxLen;
    
    TShapeDesc *ShapeDesc = CellTree[i]->GetShapeDesc();
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);

    for(int j=0;j<4;++j)
    {
      ///@attention make sure that the vector meshTrifaceHash has been created 
      int triface =
	tgml.findTriFace(tgml.meshTetGenOut.tetrahedronlist[4*i+TmpFV[j*MaxLen  ]],
			 tgml.meshTetGenOut.tetrahedronlist[4*i+TmpFV[j*MaxLen+1]],
			 tgml.meshTetGenOut.tetrahedronlist[4*i+TmpFV[j*MaxLen+2]]);

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
  for(int i=0;i<tgml.meshTetGenOut.numberoftrifaces;++i)
  {
    meshJoints[i]->SetMapType();
  }
}

#endif

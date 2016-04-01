// =======================================================================
// @(#)ReadGeo.C        1.8 11/15/99
//
// Purpose:     read geometry file and
//              generate coarse grid
//
// Author:      Volker Behns  20.01.98
// Version:     1.0
//
// =======================================================================

#include <Database.h>
#include <Domain.h>
#include <MacroCell.h>
#include <JointEqN.h>
#include <PeriodicJoint.h>
#include <Quadrangle.h>
#include <MooNMD_Io.h>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>

#include <InnerInterfaceJoint.h>
#ifdef __2D__
  #include <IsoBoundEdge.h>
  #include <IsoInterfaceJoint.h>
#else
  #include <BoundFace.h>
  #include <IsoBoundFace.h>
  #include <BdWall.h>
  #include <Hexahedron.h>
  #include <BdPlane.h>
  #include <IsoInterfaceJoint3D.h>
  #include <BdNoPRM.h>
#endif

int TDomain::ReadGeo(std::istream& dat, bool readXgeo)
{
  char line[100];
  int i, j, N_Vertices, NVpF, NVE, NBCT;
  double *DCORVG;
  int *KVERT, *KNPR;
  #ifdef __3D__
    int NBF, *BoundFaces, *FaceParam;
    int *InterfaceParam, N_Interfaces;
  #endif

  // physical references
  int *ELEMSREF;
  
  dat.getline (line, 99);
  dat.getline (line, 99);
  
  // determine dimensions for creating arrays
  dat >> N_RootCells >> N_Vertices >> NVpF >> NVE >> NBCT;
  dat.getline (line, 99);
  dat.getline (line, 99);

  // allocate auxillary fields
  #ifdef __2D__
    DCORVG =  new double[2*N_Vertices];
  #else
    DCORVG =  new double[3*N_Vertices];
  #endif
  KVERT = new int[NVE*N_RootCells];
  KNPR = new int[N_Vertices];
  
  // additional array of physical properties of elements (only if readXgeo)
  ELEMSREF = new int[N_RootCells];

  
  // read fields
  for (i=0;i<N_Vertices;i++)
  {
    #ifdef __2D__
      dat >> DCORVG[2*i] >> DCORVG[2*i + 1];
    #else
      dat >> DCORVG[3*i] >> DCORVG[3*i + 1] >> DCORVG[3*i + 2];
    #endif
    dat.getline (line, 99);
  }

  dat.getline (line, 99);
   
  for (i=0;i<N_RootCells;i++)
  {
    for (j=0;j<NVE;j++)
      dat >> KVERT[NVE*i + j];
    if(readXgeo)
      dat >> ELEMSREF[i];
    else
      ELEMSREF[i] = 0;
    dat.getline (line, 99);
  }

  dat.getline (line, 99);
  
  for (i=0;i<N_Vertices;i++)
    dat >> KNPR[i];

  #ifdef __3D__
    dat.getline (line, 99);
    dat.getline (line, 99);

    NBF = NBCT;
    BoundFaces = new int[NBF*NVpF];
    FaceParam = new int[4*NBF];

    for (i=0;i<NBF;i++)
    {
      for (j=0;j<NVpF;j++)
        dat >> BoundFaces[i*NVpF+j];
      dat.getline (line, 99);
    }

    dat.getline (line, 99);

    for (i=0;i<NBF;i++)
    {
	dat >> FaceParam[i*4] >> FaceParam[i*4+1] >>
	    FaceParam[i*4+2] >>  FaceParam[i*4+3];
	dat.getline (line, 99);
    }

    dat.getline (line, 99);
    N_Interfaces = -12345;
    if (dat.eof())
      N_Interfaces = 0;
    else
      dat >> N_Interfaces;

    if(N_Interfaces == -12345)
      N_Interfaces = 0;

    dat.getline (line, 99);

    if(N_Interfaces)
    {
      InterfaceParam = new int[N_Interfaces*6];

      dat.getline (line, 99);
      for(i=0;i<N_Interfaces;i++)
      {
        dat >> InterfaceParam[i*6]   >> InterfaceParam[i*6+1] >>
               InterfaceParam[i*6+2] >> InterfaceParam[i*6+3] >>
               InterfaceParam[i*6+4] >> InterfaceParam[i*6+5];
        dat.getline (line, 99);
      } // endfor i
    } // endif N_Interfaces
    else
      InterfaceParam = NULL;

  #endif

  #ifdef __2D__
    MakeGrid(DCORVG, KVERT, KNPR, ELEMSREF, N_Vertices, NVE);
  #else
    MakeGrid(DCORVG, KVERT, KNPR, ELEMSREF, N_Vertices, NVE,
             BoundFaces, FaceParam, NBF, NVpF,
             InterfaceParam, N_Interfaces);
  delete [] BoundFaces;
  delete [] FaceParam;
  #endif

  delete [] DCORVG;
  delete [] KVERT;
  delete [] KNPR;
  delete [] ELEMSREF;

  return 0;
}

#ifdef __2D__
// extended version
// Alfonso, 2011
int TDomain::MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                      int N_Vertices, int NVE)
{
  int a, b, j, k, l, comp, Part, NeighborID, N_Edges,maxElpV = 0;
  int aux1, aux2, aux3;
  double T_a, T_b, T, X, Y;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;
  int *KVEL;
  TVertex **NewVertices, *LocVerts[4];
  TJoint **KMT, *Joint;
  TBaseCell *JNeib1, *JNeib2;
  Shapes CellType;

  double bd_parameter_a, bd_parameter_b;

  if(TDatabase::ParamDB->SC_VERBOSE>1)
    cout << " Domain::MakeGrid() creating 2D grid " << endl;
 
  // generate vertices, edges and cells
  // search neighbours
  KVEL = new int[N_Vertices];
  
  memset(KVEL, 0, N_Vertices * SizeOfInt);
  
  //KVEL(j) = how many times vertices j appears in the list
  for (int i=0;i<NVE*N_RootCells;i++)
  {
    if (KVERT[i])
    {
      KVEL[KVERT[i]-1]++;
    }
  }
  
  // get the maximum
  for (int i=0;i<N_Vertices;i++)
    if (KVEL[i] > maxElpV) maxElpV = KVEL[i];
  
  delete [] KVEL;

  //KVEL = new int[++maxElpV * N_Vertices];
  maxElpV = maxElpV+1;
  KVEL = new int[maxElpV * N_Vertices];
  memset(KVEL, 0, maxElpV * N_Vertices * SizeOfInt);
  // first column contains the number of following elements
  for (int i=0;i<NVE*N_RootCells;i++)
  {
    if (KVERT[i])
    {
      j = (KVERT[i] - 1)*maxElpV;
      KVEL[j]++;
      KVEL[j + KVEL[j]] = i / NVE;
    }
  }


  // generate vertices
  //cout << " Domain::MakeGrid() generate " << N_Vertices << " vertices " << endl;
  NewVertices = new TVertex*[N_Vertices];
  
  for (int i=0;i<N_Vertices;i++)
  {
    //cout << " i = " << i << endl;
    if (KNPR[i])
    {
      // vertices on a boundary (inner or outer) joint
      // described by bd parametrization of 
      // part BdParts[ KNPR[i] ] [ (int) DCORVG[2*i]] ;

      T = DCORVG[2*i];
      comp = (int) T;
      if (GetLastLocalComp(KNPR[i]-1) == comp - 1)
      {
        comp--;
        T = 1.0;
      }
      else
        T -= comp;

      BdParts[KNPR[i] - 1]->GetXYofT(comp, T, X, Y);

      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;
  
      NewVertices[i] = new TVertex(X, Y);
      //cout << " vertex: " << i << " comp: " << comp << " X,Y= " << X << "," << Y << endl;
      //OutPut("bd " << i << " "<< X << " " << Y << endl);
    }
    else
    { 
      // inner vertex (described by coordinates)
      NewVertices[i] = new TVertex(DCORVG[2*i], DCORVG[2*i+1]);
    } 
  } // for (i=0;i<N_Vertices;i++) {

  // set bounding box
  StartX = Xmin;
  StartY = Ymin;
  BoundX = Xmax - Xmin;
  BoundY = Ymax - Ymin;
  //StartZ = 0; // this is 2D
  //BoundZ = 0; // this is 2D
  
  // create the CellTree and set references
  CellTree = new TBaseCell*[N_RootCells];
  for (int i=0;i<N_RootCells;i++)
  {
    CellType = Quadrangle;
    if (NVE == 3)
      CellType = Triangle;
    else
      if (!KVERT[NVE*i + 3]) CellType = Triangle;

    if (CellType == Quadrangle)
    {
      LocVerts[0] = NewVertices[KVERT[NVE*i    ]-1];
      LocVerts[1] = NewVertices[KVERT[NVE*i + 1]-1];
      LocVerts[2] = NewVertices[KVERT[NVE*i + 2]-1];
      LocVerts[3] = NewVertices[KVERT[NVE*i + 3]-1];

      CellType = ((TQuadrangle *) TDatabase::RefDescDB[Quadrangle]->
                 GetShapeDesc())->CheckQuad(LocVerts);

      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[CellType],
                                   RefLevel);

      CellTree[i]->SetVertex(0, LocVerts[0]);
      CellTree[i]->SetVertex(1, LocVerts[1]);
      CellTree[i]->SetVertex(2, LocVerts[2]);
      CellTree[i]->SetVertex(3, LocVerts[3]);

    }
    else
    {
      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[
                                   Triangle], RefLevel);
      CellTree[i]->SetVertex(0, NewVertices[KVERT[NVE*i    ]-1]);
      CellTree[i]->SetVertex(1, NewVertices[KVERT[NVE*i + 1]-1]);
      CellTree[i]->SetVertex(2, NewVertices[KVERT[NVE*i + 2]-1]);
    }
    // ReferenceID: used to specify physical/geometrical reference
    CellTree[i]->SetReference_ID(ELEMSREF[i]);
    // set the index in the cell list
    CellTree[i]->SetCellIndex(i);
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

  // generate edges
  KMT = new TJoint*[N_RootCells*4];
  for (int i=0; i<N_RootCells*4; i++)
    KMT[i] = NULL;

  for (int i=0;i<N_RootCells;i++)
  {
    N_Edges = CellTree[i]->GetN_Edges();

    for (j=0;j<N_Edges;j++)
    {
      a = KVERT[NVE*i + j] - 1;
      b = KVERT[NVE*i + (j+1) % N_Edges] - 1;
      // Part: -1 for inner, 0,1,2,... ofr bdpart 0,1,2,...
      Part = KNPR[a] - 1; // Bdpart
      // find neighbor ID
      NeighborID = -1;
      aux1 = KVEL[a*maxElpV];
      aux2 = KVEL[b*maxElpV];

      for (k=1;k<=aux1;k++)
      {
        aux3 = KVEL[a*maxElpV + k];
        if (aux3 == i) continue;

        for (l=1;l<=aux2;l++)
          if (aux3 == KVEL[b*maxElpV + l])
          {
            NeighborID = aux3;
            break;
          }
        if (NeighborID >= 0) break;
      }
      
      //OutPut(NeighborID << endl);
      if (NeighborID > i)
      {
        if( (KNPR[a]) && (KNPR[b]) )
        {
          if( (KNPR[a] == KNPR[b]) && (Interfaces[KNPR[a]-1] < 0) )
          {
            KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
            // bd.comp ID of the point is given as integer part of first coordinate
            comp = (int) DCORVG[2*a]; 
            Part = KNPR[a] - 1;
            bd_parameter_a = DCORVG[2*a] - comp;
            // correction for the endpoint of nonclosed interfaces
            if(BdParts[Part]->GetN_BdComps() == comp)
            {
              bd_parameter_a = 1.0;
              comp--;
            }
      
            // check if vertex belongs to boundary component
            if(BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[a]->GetX(), 
                                                        NewVertices[a]->GetY(),
                                                        bd_parameter_a))
            {
              // if not: 
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);

              int check_bd_node = 
                BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[a]->GetX(),
                                                         NewVertices[a]->GetY(),
                                                         bd_parameter_a);
              if (!check_bd_node)
              {
                cout << " MakeGrid(), line " << __LINE__ << ": Error: vertex " 
                     << a << " does not belong to a bd. component " << endl;
              }
            }

            comp = (int) DCORVG[2*b];
            Part = KNPR[b] - 1;

            if(comp > GetLastLocalComp(Part) ) 
              comp = GetLastLocalComp(Part);

            bd_parameter_b = DCORVG[2*b] - comp;
            if(BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[b]->GetX(), 
                                                        NewVertices[b]->GetY(),
                                                        bd_parameter_b))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);
              
              int check_bd_node = BdParts[Part]->GetBdComp(comp)->GetTofXY(
                NewVertices[b]->GetX(), NewVertices[b]->GetY(), bd_parameter_b);
              if (!check_bd_node)
              {
                cout << " MakeGrid(), line " << __LINE__ << ": Error: vertex " 
                     << b << " does not belong to a bd. component " << endl;
              }
            }
      
            T_a = bd_parameter_a;
            T_b = bd_parameter_b;
      
            if (BdParts[Part]->GetBdComp(comp)->GetType() != Line)
            {
              if (ABS(T_a) < 1e-4 || ABS(T_a) > 0.9999 ||
                  ABS(T_b) < 1e-4 || ABS(T_b) > 0.9999)
              {
                X = (NewVertices[a]->GetX() + NewVertices[b]->GetX()) / 2.;
                Y = (NewVertices[a]->GetY() + NewVertices[b]->GetY()) / 2.;

                BdParts[Part]->GetBdComp(comp)->GetTofXY(X, Y, T);

                if ((T_a - T)*(T - T_b) < 0)
                {
                  cout << __FILE__ << " line: " << __LINE__ 
                       << " changing boundary parameters of vertices " 
                       << a << " and " << b << endl;

                  if (ABS(T_a) < 1e-4)
                    T_a = 1.0;
                  else
                    if (ABS(T_a) > 0.9999)
                      T_a = 0.0;

                  if (ABS(T_b) < 1e-4)
                    T_b = 1.0;
                  else
                    if (ABS(T_b) > 0.9999)
                      T_b = 0.0;
                }
              }
            }

            if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
              KMT[i*4 + j] = new TIsoInterfaceJoint(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b,CellTree[i], CellTree[NeighborID]);
            else
              KMT[i*4 + j] = new TInterfaceJoint(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b,CellTree[i], CellTree[NeighborID]);

            CellTree[i]->SetJoint(j, KMT[i*4 + j]);
            ((TInterfaceJoint *)KMT[i*4 + j])->CheckOrientation();
          } // two vertices on the same Interfaces
          else
          {
            // two vertices on the boundary
            //KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
            
            // check if neighbor belongs to different subdomain:
            if (CellTree[i]->GetReference_ID() 
              == CellTree[NeighborID]->GetReference_ID())
            {
              KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
            }
            else
            {
              //KMT[i*4+j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
              KMT[i*4+j] = new TInnerInterfaceJoint(CellTree[i], 
                                                    CellTree[NeighborID]);
              ((TInnerInterfaceJoint*)KMT[i*4 + j])->SetParams(
                 NewVertices[a]->GetX(),
                 NewVertices[a]->GetY(),
                 NewVertices[b]->GetX()-NewVertices[a]->GetX(),
                 NewVertices[b]->GetY()-NewVertices[a]->GetY());
              ((TInnerInterfaceJoint*)KMT[i*4 + j])->SetIndexInNeighbor(
                 CellTree[i],j);
            }
          }
        }
        else
        {
          // at least one vertex is inside the domain
          // check if neighbor belongs to different subdomain:
          if(CellTree[i]->GetReference_ID() 
             == CellTree[NeighborID]->GetReference_ID())
          {
            KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
          }
          else
          {
            //KMT[i*4+j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
            KMT[i*4+j] = new TInnerInterfaceJoint(CellTree[i], 
                                                  CellTree[NeighborID]);
            ((TInnerInterfaceJoint*)KMT[i*4 + j])->SetParams(
                   NewVertices[a]->GetX(),
                   NewVertices[a]->GetY(),
                   NewVertices[b]->GetX()-NewVertices[a]->GetX(),
                   NewVertices[b]->GetY()-NewVertices[a]->GetY());
            ((TInnerInterfaceJoint*)KMT[i*4 + j])->SetIndexInNeighbor(
                   CellTree[i],j);
          }
        }
        CellTree[i]->SetJoint(j, KMT[i*4 + j]);
      }
      else 
      {
        if (NeighborID == -1)
        {
          if (Interfaces[KNPR[a]-1] < 0)
          {
            comp = (int) DCORVG[2*b];
            T_b = DCORVG[2*b] - comp;
            Part = KNPR[b] - 1;
            
            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[a]->GetX(), 
                                                         NewVertices[a]->GetY(), 
                                                         T_a))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);
        
              BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[a]->GetX(), 
                                                       NewVertices[a]->GetY(), 
                                                       T_a);
            }
          }
          else if (Interfaces[KNPR[b]-1] < 0)
          {
            comp = (int) DCORVG[2*a];
            T_a = DCORVG[2*a] - comp;
            
            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[b]->GetX(), 
                                                         NewVertices[b]->GetY(), 
                                                         T_b))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);
              
              BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[b]->GetX(), 
                                                       NewVertices[b]->GetY(), 
                                                       T_b);
            }
          }
          else
          {
            //cout << "a=" << a << endl;
            comp = (int) DCORVG[2*a];
            T_a = DCORVG[2*a] - comp;
            T_b = DCORVG[2*b] - comp;
            //cout << "comp: " << comp << endl;
          }
          //cout << "Part = " << Part << endl;
          //cout << " test = " << BdParts[Part]->GetBdComp(comp)->IsFreeBoundary() << endl;
          if (T_b < T_a) 
            T_b = 1.;
          if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
            Joint = new TIsoBoundEdge(BdParts[Part]->GetBdComp(comp), T_a, T_b);
          else 
          {
            Joint = new TBoundEdge(BdParts[Part]->GetBdComp(comp), T_a, T_b);
          }
          CellTree[i]->SetJoint(j, Joint);
        }
        else
        {
          // joint already created (from the neighboring cell)
          TJoint *nJoint;
          for (int jj=0;jj<CellTree[NeighborID]->GetN_Edges();jj++)
          {
            nJoint = CellTree[NeighborID]->GetJoint(jj);
            JNeib1 = nJoint->GetNeighbour(0);
            JNeib2 = nJoint->GetNeighbour(1);
            if((JNeib1 == CellTree[NeighborID] && JNeib2 == CellTree[i]) 
               || (JNeib1 == CellTree[i] && JNeib2 == CellTree[NeighborID])) 
              break;
          }
          CellTree[i]->SetJoint(j, nJoint);
          
          if (CellTree[i]->GetReference_ID()
            != CellTree[NeighborID]->GetReference_ID())
          {
            ((TInnerInterfaceJoint*)nJoint)->SetIndexInNeighbor(CellTree[i],j);
          }
        }
      }
    }
  } //for (i=0;i<N_RootCells;i++) {

  // free memory
  delete [] KVEL;
  delete [] NewVertices;
  delete [] KMT;
  
  return 0;
}

int TDomain::MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int N_Vertices, 
                      int NVE)
{
   int a, b, i, j, k, l, comp, Part, Neib, N_E, maxElpV = 0;
  int aux1, aux2, aux3;
  double T_a, T_b, T, X, Y;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;
  int *KVEL;
  TVertex **NewVertices, *LocVerts[4];
  TJoint **KMT, *Joint;
  TBaseCell *JNeib1, *JNeib2;
  Shapes CellType;
 
  // generate vertices, edges and cells
  // search neighbours
  KVEL = new int[N_Vertices];

  memset(KVEL, 0, N_Vertices * SizeOfInt);
 
  for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i]) KVEL[KVERT[i]-1]++;

   for (i=0;i<N_Vertices;i++)
    if (KVEL[i] > maxElpV) maxElpV = KVEL[i];

  delete [] KVEL;
  KVEL = new int[++maxElpV * N_Vertices];

  memset(KVEL, 0, maxElpV * N_Vertices * SizeOfInt);
  
  // first column contains the number of following elements
   for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i])
    {
      j = (KVERT[i] - 1)*maxElpV;
      KVEL[j]++;
      KVEL[j + KVEL[j]] = i / NVE;
    }

  // generate vertices
    NewVertices = new TVertex*[N_Vertices];
  
  for (i=0;i<N_Vertices;i++)
    if (KNPR[i])
    {
      T = DCORVG[2*i];
      comp = (int) T;
      if (GetLastLocalComp(KNPR[i]-1) == comp - 1)
      {
        comp--;
        T = 1.0;
      }
      else
        T -= comp;
      
      BdParts[KNPR[i] - 1]->GetXYofT(comp, T, X, Y);

      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;

      NewVertices[i] = new TVertex(X, Y);
      //OutPut("bd " << i << " "<< X << " " << Y << endl);
    }
    else
      NewVertices[i] = new TVertex(DCORVG[2*i], DCORVG[2*i+1]);

   // set bounding box
  StartX = Xmin;
  StartY = Ymin;
  BoundX = Xmax - Xmin;
  BoundY = Ymax - Ymin;

  // generate cells
  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellType = Quadrangle;
    if (NVE == 3)
      CellType = Triangle;
    else
      if (!KVERT[NVE*i + 3]) CellType = Triangle;

    if (CellType == Quadrangle)
    {
      LocVerts[0] = NewVertices[KVERT[NVE*i    ]-1];
      LocVerts[1] = NewVertices[KVERT[NVE*i + 1]-1];
      LocVerts[2] = NewVertices[KVERT[NVE*i + 2]-1];
      LocVerts[3] = NewVertices[KVERT[NVE*i + 3]-1];

      CellType = ((TQuadrangle *) TDatabase::RefDescDB[Quadrangle]->
                 GetShapeDesc())->CheckQuad(LocVerts);

      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[CellType],
                                   RefLevel);

      CellTree[i]->SetVertex(0, LocVerts[0]);
      CellTree[i]->SetVertex(1, LocVerts[1]);
      CellTree[i]->SetVertex(2, LocVerts[2]);
      CellTree[i]->SetVertex(3, LocVerts[3]);

    }
    else
    {
      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[
                                   Triangle], RefLevel);
      CellTree[i]->SetVertex(0, NewVertices[KVERT[NVE*i    ]-1]);
      CellTree[i]->SetVertex(1, NewVertices[KVERT[NVE*i + 1]-1]);
      CellTree[i]->SetVertex(2, NewVertices[KVERT[NVE*i + 2]-1]);
    }
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

  // generate edges
  KMT = new TJoint*[N_RootCells*4];
  for (i=0;i<N_RootCells*4;i++)
    KMT[i] = NULL;
  
  for (i=0;i<N_RootCells;i++)
  {
    switch (CellTree[i]->GetType())
    {
      case Triangle: 
           N_E = 3;
           break;
      case Parallelogram: 
      case Quadrangle: 
      case Rectangle: 
           N_E = 4;
           break;
      default:
           cerr << "Error in ReadGeo" << endl;
           return -1;
    }

    //OutPut("edges " << i << " " << N_RootCells << endl);
    //cout << i;
   for (j=0;j<N_E;j++)
    {
      a = KVERT[NVE*i + j] - 1;
      b = KVERT[NVE*i + (j+1) % N_E] - 1;
      Part = KNPR[a] - 1;
      Neib = -1;

      aux1 = KVEL[a*maxElpV];
      aux2 = KVEL[b*maxElpV];

      for (k=1;k<=aux1;k++)
      {
        aux3 = KVEL[a*maxElpV + k];
        if (aux3 == i) continue;

        for (l=1;l<=aux2;l++)
          if (aux3 == KVEL[b*maxElpV + l])
          {
            Neib = aux3;
            break;
          }
        if (Neib >= 0) break;
      }
      //OutPut(Neib << endl);
      if (Neib > i)
      {
        if( (KNPR[a]) && (KNPR[b]) )
        {
          if( (KNPR[a] == KNPR[b]) && (Interfaces[KNPR[a]-1] < 0) )
          {
            KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[Neib]);
            comp = (int) DCORVG[2*a];
            Part = KNPR[a] - 1;
            T_a = DCORVG[2*a] - comp;

            // correction for the endpoint of nonclosed interfaces
            if(BdParts[Part]->GetN_BdComps() == comp)
            {
              T_a = 1.0;
              comp--;
            }

            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(
                  NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);

              BdParts[Part]->GetBdComp(comp)->GetTofXY(
                NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a);
            }

            comp = (int) DCORVG[2*b];
            Part = KNPR[b] - 1;

            // correction needed for not closed interfaces
            if(comp > GetLastLocalComp(Part) ) comp = GetLastLocalComp(Part);

            T_b = DCORVG[2*b] - comp;
            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(
                  NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);

              BdParts[Part]->GetBdComp(comp)->GetTofXY(
                NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_a);
            }

            if (BdParts[Part]->GetBdComp(comp)->GetType() != Line)
              if (ABS(T_a) < 1e-4 || ABS(T_a) > 0.9999 ||
                  ABS(T_b) < 1e-4 || ABS(T_b) > 0.9999)
              {
                X = (NewVertices[a]->GetX() + NewVertices[b]->GetX()) / 2;
                Y = (NewVertices[a]->GetY() + NewVertices[b]->GetY()) / 2;

                BdParts[Part]->GetBdComp(comp)->GetTofXY(X, Y, T);

                if ((T_a - T)*(T - T_b) < 0)
                {
                  if (ABS(T_a) < 1e-4)
                    T_a = 1.0;
                  else
                    if (ABS(T_a) > 0.9999)
                      T_a = 0.0;

                  if (ABS(T_b) < 1e-4)
                    T_b = 1.0;
                  else
                    if (ABS(T_b) > 0.9999)
                      T_b = 0.0;
                }
              }

            if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
              KMT[i*4 + j] = new TIsoInterfaceJoint(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b,CellTree[i], CellTree[Neib]);
            else
              KMT[i*4 + j] = new TInterfaceJoint(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b,CellTree[i], CellTree[Neib]);

            CellTree[i]->SetJoint(j, KMT[i*4 + j]);
            ((TInterfaceJoint *)KMT[i*4 + j])->CheckOrientation();
          }
          else
          {
            KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[Neib]);
          }
        }
        else
        {
          KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[Neib]);
        }
        CellTree[i]->SetJoint(j, KMT[i*4 + j]);
      }
      else
        if (Neib == -1)
        {
          if (Interfaces[KNPR[a]-1] < 0)
          {
	      comp = (int) DCORVG[2*b];
            T_b = DCORVG[2*b] - comp;
            Part = KNPR[b] - 1;

            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(
                  NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);

              BdParts[Part]->GetBdComp(comp)->GetTofXY(
                NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a);
            }
          }
          else
            if (Interfaces[KNPR[b]-1] < 0)
            {
		comp = (int) DCORVG[2*a];
              T_a = DCORVG[2*a] - comp;

              if (BdParts[Part]->GetBdComp(comp)->GetTofXY(
                    NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
              {
                if (comp)
                  comp--;
                else
                  comp = GetLastLocalComp(Part);

                BdParts[Part]->GetBdComp(comp)->GetTofXY(
                  NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b);
              }
            }
            else
            {
              comp = (int) DCORVG[2*a];
	      T_a = DCORVG[2*a] - comp;
              T_b = DCORVG[2*b] - comp;
            }

          if (T_b < T_a) T_b = 1.;
          if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
            Joint = new TIsoBoundEdge(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b);
          else
            Joint = new TBoundEdge(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b);
          CellTree[i]->SetJoint(j, Joint);
        }
        else
        {
          for (k=0;k<4;k++)
          {
            Joint = KMT[Neib*4 + k];
            if (Joint)
            {
              aux1 = Neib*4 + k;
              JNeib1 = KMT[aux1]->GetNeighbour(0);
              JNeib2 = KMT[aux1]->GetNeighbour(1);

              if ((JNeib1 == CellTree[Neib] && JNeib2 == CellTree[i]) ||
                  (JNeib1 == CellTree[i] && JNeib2 == CellTree[Neib]) ){ break;}
            }
          }

          CellTree[i]->SetJoint(j, KMT[Neib*4 + k]);
        }
    }
   //OutPut("edges done " << i << " " << N_RootCells << endl);

  }

   // free memory
  delete [] KVEL;
  delete [] NewVertices;
  delete [] KMT;
  
  return 0;
}

#else
int TDomain::ReadSandwichGeo(std::istream& dat)
{
  char line[100];
  int i, j, N_Vertices, NVpF, NVE, NBCT;
  double *DCORVG;
  int *KVERT, *KNPR;
  double DriftX, DriftY, DriftZ, *Lambda;
  int N_Layers, grid_type;

  grid_type = TDatabase::ParamDB->GRID_TYPE;

  dat.getline (line, 99);
  dat.getline (line, 99);

  // determine dimensions for creating arrays
  dat >> N_RootCells >> N_Vertices >> NVpF >> NVE >> NBCT;
  dat.getline (line, 99);
  dat.getline (line, 99);

  // allocate auxillary fields
  DCORVG =  new double[2*N_Vertices];
  KVERT = new int[NVE*N_RootCells];
  KNPR = new int[N_Vertices];
  // read fields
  for (i=0;i<N_Vertices;i++)
  {
    dat >> DCORVG[2*i] >> DCORVG[2*i + 1];
    dat.getline (line, 99);
  }

  dat.getline (line, 99);

  for (i=0;i<N_RootCells;i++)
  {
    for (j=0;j<NVE;j++)
      dat >> KVERT[NVE*i + j];
    dat.getline (line, 99);
  }

  dat.getline (line, 99);

  for (i=0;i<N_Vertices;i++)
    dat >> KNPR[i];

  DriftX = TDatabase::ParamDB->DRIFT_X;
  DriftY = TDatabase::ParamDB->DRIFT_Y;
  DriftZ = TDatabase::ParamDB->DRIFT_Z;

  if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 1234)
   {
    TDatabase::ParamDB->N_CELL_LAYERS = 3;
    N_Layers = 4;
    Lambda = new double[N_Layers];
    Lambda[0] = 0.0;
    Lambda[1] = 7.66/20.32;
    Lambda[2] = 15.32/20.32;
    Lambda[3] = 1.0;
   }
/*  else if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 1356)
   {
    N_Layers = TDatabase::ParamDB->N_CELL_LAYERS+1;     
    Lambda = new double[N_Layers];     
    double tmp [] = { 0, 0.0119047619, 0.02380952381, 0.03571428571, 0.04761904762,
                 0.05952380952, 0.07142857143, 0.08333333333, 0.09523809524, 0.10714285714,
                 0.11904761905, 0.13333333333, 0.14761904762, 0.16666666667, 0.18571428571, 
                 0.20952380952, 0.23333333333, 0.2619047619, 0.29523809524, 0.33333333333, 
                 0.38095238095, 0.42857142857, 0.47619047619, 0.52380952381, 0.57142857143, 
                 0.61904761905, 0.66666666667, 0.71428571429, 0.7619047619, 0.80952380952, 
                 0.85714285714, 0.90476190476, 0.95238095238, 1  };
                 
    for(i=0;i<N_Layers; i++)     
     Lambda[i] = tmp[i];
                 
//      MPI_Finalize();
//      exit(0);
   }  */   
  else if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 180)
  {
    // turbulent channel flow for Re_tau = 180, IMPORTANT !!!
     N_Layers = TDatabase::ParamDB->N_CELL_LAYERS+1;
     Lambda = new double[N_Layers];
     for(i=0;i<N_Layers;i++)
      {
      if (grid_type == 0)
       Lambda[i] = 1 + tanh(2.75*(2.0*i/(N_Layers-1) -1))/tanh(2.75);
      else
       Lambda[i] = 1-cos(i*Pi/(N_Layers-1));

       OutPut("z coordinate " << Lambda[i] <<endl);
       Lambda[i] /= 2.0;
      }
     Lambda[0] = 0;
     Lambda[N_Layers-1] = 1;
   }
  else if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 6)
   {
    // for Channel_Carolina.h
    N_Layers = TDatabase::ParamDB->N_CELL_LAYERS+1;
    Lambda = new double[N_Layers];
    for(i=0;i<N_Layers;i++)
     {
      Lambda[i] = tanh(4.50*i/(N_Layers-1))/tanh(4.50);
      OutPut("z coordinate " << Lambda[i] <<endl);
     }
    Lambda[0] = 0;
    Lambda[N_Layers-1] = 1;     
   }
  else
   {
    N_Layers = TDatabase::ParamDB->N_CELL_LAYERS+1;
    Lambda = new double[N_Layers];
    for(i=0;i<N_Layers;i++)
     {
      Lambda[i] = i * (1.0/(N_Layers-1));
     //OutPut("lambda " << i << " " << Lambda[i] << endl);
     }

  }




  // for WIND TUNNEL
  /*if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 1506)
  {
      N_Layers = 19;
      Lambda = new double[N_Layers];
      for(i=0;i<N_Layers;i++)
      {
	  Lambda[i] = i*0.01;
	  OutPut("wind tunnel: z coordinate " << Lambda[i] <<endl);
      }
      }*/

  MakeSandwichGrid(DCORVG, KVERT, KNPR, N_Vertices, NVE,
                   DriftX, DriftY, DriftZ, N_Layers, Lambda);

  delete DCORVG;
  delete KVERT;
  delete KNPR;

  return 0;
}

// this makes a 3D coarse grid
int TDomain::MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                      int N_Vertices, int NVE, int *BoundFaces,
                      int *FaceParam, int NBF, int NVpF,
                      int *InterfaceParam, int N_Interfaces)
{
  int i, j, k, l, m, maxElpV = 0, comp, auxi, auxj;
  double S, T, X, Y, Z;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;
  double Zmin = 1e10, Zmax = -1e10;
  int *KVEL, *NewJoints_aux1, *NewJoints_aux2, NewJoints_aux3;
  int aux1, aux2, aux3;
  TVertex **NewVertices, *CurrVert, *Vert[4];
  std::vector<TJoint *> NewJoints;
  Shapes CellType;
  TBaseCell *CurrCell, *CurrCell_aux;
  int maxlen, maxlen_aux, N_Faces, N_Verts, N_Faces_aux;
  const int *TmpFV, *TmpFV_aux, *TmpLen, *TmpLen_aux;

  TBoundFace *bdface;
  TBoundComp3D *bdcomp;
  double Param1[4], Param2[4];
  JointType CurrJointType;
  TInterfaceJoint3D *IFJoint;
  TIsoInterfaceJoint3D *IIJoint;

  // generate vertices, faces and cells
  // search neighbours
  KVEL = new int[N_Vertices];

  memset(KVEL, 0, N_Vertices * SizeOfInt);

  for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i]) KVEL[KVERT[i]-1]++;

  for (i=0;i<N_Vertices;i++)
    if (KVEL[i] > maxElpV) maxElpV = KVEL[i];

  delete[] KVEL;
  KVEL = new int[++maxElpV * N_Vertices];

  memset(KVEL, 0, maxElpV * N_Vertices * SizeOfInt);
  
  // first colomn contains the number of following elements
  for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i])
    {
      j = (KVERT[i] - 1)*maxElpV;
      KVEL[j]++;
      KVEL[j + KVEL[j]] = i / NVE;
    }

  // generate vertices
  NewVertices = new TVertex*[N_Vertices];

  for (i=0;i<N_Vertices;i++)
    if (KNPR[i])
    {
	// parametrisation of the boundary is given
	if (KNPR[i]>0)
	{
	    // component of boundary
	    comp = (int) DCORVG[3*i];
	    // value of first parameter
	    T = DCORVG[3*i + 1];
	    // value of second parameter
	    S = DCORVG[3*i + 2];
	    // get coordinates from parameters
	    BdParts[KNPR[i] - 1]->GetXYZofTS(comp, T, S, X, Y, Z);
	}
	else
	{
	    X = DCORVG[3*i];
	    Y = DCORVG[3*i+1];
	    Z = DCORVG[3*i+2];
	    //cout << "read " << X << " " << Y << " " << Z << endl;
	}
      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;
      if (Z > Zmax) Zmax = Z;
      if (Z < Zmin) Zmin = Z;

      NewVertices[i] = new TVertex(X, Y, Z);
    }
    else
      NewVertices[i] = new TVertex(DCORVG[3*i], DCORVG[3*i + 1],
                                   DCORVG[3*i + 2]);

  // set bounding box
  StartX = Xmin;
  StartY = Ymin;
  StartZ = Zmin;
  BoundX = Xmax - Xmin;
  BoundY = Ymax - Ymin;
  BoundZ = Zmax - Zmin;

  // generate cells
  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellType = Hexahedron;
    if (NVE == 4)
      CellType = Tetrahedron;
    else
      if (!KVERT[NVE*i + 7]) CellType = Tetrahedron;

    if (CellType == Hexahedron)
    {
      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Hexahedron], RefLevel);

      auxi = NVE*i;
      for (j=0;j<8;j++)
        CellTree[i]->SetVertex(j, NewVertices[KVERT[auxi++]-1]);

      // check if Hexahedron is even a brick and if yes, change the refinement
      // descriptor accordingly. This allows us to use an affine reference 
      // mapping instead of the (more expensive) trilinear one.
      CellType = ((THexahedron *) TDatabase::RefDescDB[Hexahedron]->
                 GetShapeDesc())->CheckHexa(((TGridCell *)(CellTree[i]))
                        ->GetVertices());
      if(CellType != Hexahedron)
      {
        int ret = CellTree[i]->SetRefDesc(TDatabase::RefDescDB[CellType]);
        if(ret==-1)
        {
          ErrMsg("setRefDesc(" << CellType << ") was not sucessfull");
          exit(0);
        }
      }
    }
    else
    {
      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron],
                                   RefLevel);

      auxi = NVE*i;
      for (j=0;j<4;j++)
        CellTree[i]->SetVertex(j, NewVertices[KVERT[auxi++]-1]);
    }
  }

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  // generate faces
  NewJoints.resize(N_RootCells*6, nullptr);

  //OutPut("NBF a " << NBF << endl);
  // first generate boundary joints
  // NBF -- number of faces on boundaries
  for (i=0;i<NBF;i++)
  {
    auxi = i*4;
    //cout << FaceParam[auxi + 2] <<  endl;
    bdcomp = BdParts[FaceParam[auxi + 2] - 1]->GetBdComp(FaceParam[auxi+3] - 1);
    //cout << "free " << bdcomp->IsFreeBoundary()  << endl;
    if(bdcomp->IsFreeBoundary())
    {
      NewJoints.at(FaceParam[auxi]*6 + FaceParam[auxi + 1] - 7) =
        new TIsoBoundFace(bdcomp);
    }
    else
    {
      NewJoints.at(FaceParam[auxi]*6 + FaceParam[auxi + 1] - 7) =
        new TBoundFace(bdcomp);
    }
  }
//   OutPut("NBF b " << NBF << endl);

  // generate interface joints
  for (i=0;i<N_Interfaces;i++)
  {
    auxi = i*6;

    bdcomp = BdParts[InterfaceParam[auxi + 4] - 1]->
                GetBdComp(InterfaceParam[auxi + 5] - 1);
    if(bdcomp->IsFreeBoundary())
    {
      NewJoints.at(InterfaceParam[auxi    ]*6 + InterfaceParam[auxi + 1] - 7) =
      NewJoints.at(InterfaceParam[auxi + 2]*6 + InterfaceParam[auxi + 3] - 7) =
        new TIsoInterfaceJoint3D(bdcomp, Param1, Param2,
                CellTree[InterfaceParam[auxi]-1],
                CellTree[InterfaceParam[auxi+2]-1]);
    }
    else
    {
      NewJoints.at(InterfaceParam[auxi    ]*6 + InterfaceParam[auxi + 1] - 7) =
      NewJoints.at(InterfaceParam[auxi + 2]*6 + InterfaceParam[auxi + 3] - 7) =
        new TInterfaceJoint3D(bdcomp, Param1, Param2,
                CellTree[InterfaceParam[auxi]-1],
                CellTree[InterfaceParam[auxi+2]-1]);
    }
  }

  // new generate EqN-joints
  NewJoints_aux1 = new int[maxElpV];
  NewJoints_aux2 = new int[maxElpV];

  for (i=0;i<N_RootCells;i++)
  {
    CurrCell = CellTree[i];
    N_Faces = CurrCell->GetN_Faces();
    CurrCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, maxlen);

    auxi = i*6;
    for (j=0;j<N_Faces;j++)
    {
      if (!NewJoints.at(auxi++))
      {
        N_Verts = TmpLen[j];
        aux2 = (KVERT[NVE*i + TmpFV[maxlen*j]] - 1)*maxElpV;
        aux1 = KVEL[aux2++];
        for (k=0;k<aux1;k++)
          NewJoints_aux1[k] = KVEL[aux2++];

        for (k=1;k<N_Verts;k++)
        {
          aux3 = (KVERT[NVE*i + TmpFV[maxlen*j+k]] - 1)*maxElpV;
          aux2 = KVEL[aux3++];
          for (l=0;l<aux2;l++)
            NewJoints_aux2[l] = KVEL[aux3++];

          for (aux3=l=0;l<aux1;l++)
          {
            NewJoints_aux3 = NewJoints_aux1[l];
            for (m=0;m<aux2;m++)
              if (NewJoints_aux3 == NewJoints_aux2[m])
              {
                NewJoints_aux1[aux3++] = NewJoints_aux3;
                break;
              }
          }

          aux1 = aux3;
        }

        if (aux1 == 2)
        {
          if (i == NewJoints_aux1[0])
          {
            aux1 = i;
            aux2 = NewJoints_aux1[1];
          }
          else
          {
            aux1 = i;
            aux2 = NewJoints_aux1[0];
          }
          
          auxj = j * maxlen;
          for (k=0;k<N_Verts;k++)
            Vert[k] = CurrCell->GetVertex(TmpFV[auxj++]);

          CurrCell_aux = CellTree[aux2];

          N_Faces_aux = CurrCell_aux->GetN_Faces();
          CurrCell_aux->GetShapeDesc()->GetFaceVertex(TmpFV_aux,
                          TmpLen_aux, maxlen_aux);

          for (k=0;k<N_Faces_aux;k++)
            if (N_Verts == TmpLen_aux[k])
            {
              auxj = k * maxlen_aux;
              CurrVert = CurrCell_aux->GetVertex(TmpFV_aux[auxj]);
              for (m=l=0;l<N_Verts;l++)
                if (Vert[l] == CurrVert) break;

              if (l != N_Verts)
                for (m=1;m<N_Verts;m++)
                  if (Vert[(N_Verts + l - m) % N_Verts] != CurrCell_aux->
                        GetVertex(TmpFV_aux[auxj + m]))
                    break;

              if (m == N_Verts) break;
            }

          if (k == N_Faces_aux)
          {
            cerr << "Error in ReadGeo: could not find local face" << endl;
            exit (-1);
          }

          NewJoints.at(aux1*6 + j) = NewJoints.at(aux2*6 + k) = new
            TJointEqN(CellTree[aux2], CellTree[aux1]);
          NewJoints.at(aux1*6 + j)->SetMapType(l);
        }
        else
        {
          cerr << "Error in ReadGeo: no element on face" << endl;
          exit (-1);
        }
      } // endif (!NewJoints[auxi++])
    } // endfor j

    // copy joints to cells
    auxi = i*6;
    for (j=0;j<N_Faces;j++)
    {
      CurrCell->SetJoint(j, NewJoints.at(auxi));
      //if(NewJoints[auxi]->GetType() == InterfaceJoint3D)
      //  cout << auxi << " cell: " << i << " joint: " << j << endl;
      auxi++;
    }
  }

  for(i=0;i<N_RootCells;i++)
  {
    CurrCell = CellTree[i];
    N_Faces = CurrCell->GetN_Faces();
    for(j=0;j<N_Faces;j++)
    {
      if(CurrCell->GetJoint(j)->GetType() == InterfaceJoint3D)
       ((TInterfaceJoint3D *)(CurrCell->GetJoint(j)))->CheckOrientation();
    } // endfor j
  } // endfor i

  for(i=0;i<N_RootCells;i++)
  {
      //cout << "Cell number: " << i << endl;
    CurrCell = CellTree[i];
    N_Faces = CurrCell->GetN_Faces();
    CurrCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, maxlen);

    for(j=0;j<N_Faces;j++)
    {
      CurrJointType = CurrCell->GetJoint(j)->GetType();
      if(CurrJointType == BoundaryFace ||
         CurrJointType == IsoBoundFace ||
         CurrJointType == InterfaceJoint3D ||
         CurrJointType == IsoInterfaceJoint3D)
      {
        if(CurrJointType == BoundaryFace ||
           CurrJointType == IsoBoundFace)
        {
          bdface = (TBoundFace*)(CurrCell->GetJoint(j));
          bdcomp = bdface->GetBoundComp();
        }

        if(CurrJointType == InterfaceJoint3D ||
           CurrJointType == IsoInterfaceJoint3D )
        {
          bdcomp = ((TInterfaceJoint3D*)(CurrCell->GetJoint(j)))->GetBoundComp();
          OutPut("Interface ");
          OutPut(((TInterfaceJoint3D*)(CurrCell->GetJoint(j)))
                        ->CheckInside(CurrCell) << endl);
        }

        //cout << "Face number: " << j << " vertices:" << endl;
        for(k=0;k<TmpLen[j];k++)
        {
          //cout << TmpFV[maxlen*j+k] << ": ";
          CurrCell->GetVertex(TmpFV[maxlen*j+k])->GetCoords(X, Y, Z);
          //cout << "X: " << X << " Y: " << Y << " Z: " << Z;
          bdcomp->GetTSofXYZ(X, Y, Z, T, S);
          //cout << " T: " << T << " S: " << S << endl;
          Param1[k] = T;
          Param2[k] = S;
        } // endfor k


        if(CurrJointType == BoundaryFace ||
           CurrJointType == IsoBoundFace)
        {
          bdface->SetParameters(Param1, Param2);
        }
        else
        {
          ((TInterfaceJoint3D*)(CurrCell->GetJoint(j)))
                ->SetMapType();
          ((TInterfaceJoint3D*)(CurrCell->GetJoint(j)))
                ->SetParameters(Param1, Param2);
        }
      } // endif BoundaryFace
    } // endfor j
  } // endfor i

  for(int i=0; i<N_RootCells; i++)
  {
    CurrCell = CellTree[i];
    CurrCell->SetReference_ID(ELEMSREF[i]);
  }
  
  // free memory
  delete[] KVEL;
  delete[] NewVertices;
  delete[] NewJoints_aux1;
  delete[] NewJoints_aux2;

  return 0;
}



int TDomain::MakeSandwichGrid(double *DCORVG, int *KVERT, int *KNPR,
                              int N_Vertices, int NVE,
                              double DriftX, double DriftY, double DriftZ,
                              int N_Layers, double *Lambda)
{
  int a, b, i, j, k, l, comp, Part, Neib, N_E, maxElpV = 0;
  int aux1, aux2, aux3;
  double T_a, T_b, T, X, Y, Z, Xa, Xb, Ya, Yb, Za, Zb, S;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;
  double Zmin = 0;
  int *KVEL;
  TVertex **NewVertices, *LocVerts[8], *vert;
  TJoint **KMT, *Joint;
  TBaseCell *JNeib1, *JNeib2, *CurrCell;
  Shapes CellType;
  TBdWall *BdWall;
  double x, y, z;
  double Param1[4], Param2[4];
  TBdPlane *Top, *Bottom;
  int SortOrder, *OrderArray;
  TVertex *vert0, *vert1, *vert2;
  int k0, k1, k2, k3;
  int i0, i1, i2;
  int *KMTupper, *KMTlower;

  double x0, x1, x2, x3;
  double y0, y1, y2, y3;
  double z0, z1, z2, z3;
  double det;

  // generate vertices, edges and cells
  // search neighbours
  KVEL = new int[N_Vertices];

  memset(KVEL, 0, N_Vertices * SizeOfInt);

  for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i]) KVEL[KVERT[i]-1]++;

  for (i=0;i<N_Vertices;i++)
    if (KVEL[i] > maxElpV) maxElpV = KVEL[i];

  delete KVEL;

  maxElpV++;
  KVEL = new int[maxElpV * N_Vertices];

  memset(KVEL, 0, maxElpV * N_Vertices * SizeOfInt);
  
  // first colomn contains the number of following elements
  for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i])
    {
      j = (KVERT[i] - 1)*maxElpV;
      KVEL[j]++;
      KVEL[j + KVEL[j]] = i / NVE;
    }

  // generate vertices
  NewVertices = new TVertex*[N_Vertices*N_Layers];
  
  for (i=0;i<N_Vertices;i++)
    if (KNPR[i])
    {
      T = DCORVG[2*i];
      comp = (int) T;
      if (GetLastLocalComp(KNPR[i]-1) == comp - 1)
      {
        comp--;
        T = 1.0;
      }
      else
        T -= comp;

      // cout << T << " " << comp << endl;
      
      BdWall = (TBdWall *)(BdParts[KNPR[i] - 1]->GetBdComp(comp));
      BdWall->GetBdComp2D()->GetXYofT(T, X, Y);
      BdWall->SetParams(DriftX, DriftY, DriftZ);

      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;

      for(j=0;j<N_Layers;j++)
      {
        x = X + Lambda[j] * DriftX;
        y = Y + Lambda[j] * DriftY;
        z = Lambda[j] * DriftZ;
        NewVertices[j*N_Vertices+i] = new TVertex(x, y, z);
      }
    }
    else
    {
      for(j=0;j<N_Layers;j++)
      {
        x = DCORVG[2*i] + Lambda[j] * DriftX;
        y = DCORVG[2*i+1] + Lambda[j] * DriftY;
        z = Lambda[j] * DriftZ;
        NewVertices[j*N_Vertices+i] = new TVertex(x, y, z);
      }
    }

  // for(i=0;i<N_Vertices*N_Layers;i++)
  //   cout << i << NewVertices[i] << endl;

  // set bounding box
  StartX = Xmin;
  StartY = Ymin;
  StartZ = Zmin;
  BoundX = Xmax - Xmin + DriftX;
  BoundY = Ymax - Ymin + DriftY;
  BoundZ = DriftZ-Zmin;

  switch(NVE)
  {
    case 4: // quadrilaterals in 2D mesh
      // generate cells
      CellTree = new TBaseCell*[N_RootCells*N_Layers];
    
      for (i=0;i<N_RootCells;i++)
      {
        CellType = Quadrangle;
        if (NVE == 3)
          CellType = Triangle;
        else
          if (!KVERT[NVE*i + 3]) CellType = Triangle;
    
        if (CellType == Quadrangle)
        {
          for(j=1;j<N_Layers;j++)
          {
            k = (j-1)*N_Vertices + KVERT[NVE*i    ]-1;
            LocVerts[0] = NewVertices[k];
            k += N_Vertices;
            LocVerts[4] = NewVertices[k];
    
            k = (j-1)*N_Vertices + KVERT[NVE*i + 1]-1;
            LocVerts[1] = NewVertices[k];
            k += N_Vertices;
            LocVerts[5] = NewVertices[k];
    
            k = (j-1)*N_Vertices + KVERT[NVE*i + 2]-1;
            LocVerts[2] = NewVertices[k];
            k += N_Vertices;
            LocVerts[6] = NewVertices[k];
    
            k = (j-1)*N_Vertices + KVERT[NVE*i + 3]-1;
            LocVerts[3] = NewVertices[k];
            k += N_Vertices;
            LocVerts[7] = NewVertices[k];
    
            CellType = ((THexahedron *) TDatabase::RefDescDB[Hexahedron]->
                       GetShapeDesc())->CheckHexa(LocVerts);
    
            CellTree[(j-1)*N_RootCells + i] =
                    new TMacroCell(TDatabase::RefDescDB[CellType], RefLevel);
    
            for(k=0;k<8;k++) 
            {
              CellTree[(j-1)*N_RootCells + i]->SetVertex(k, LocVerts[k]);
              // cout << (j-1)*N_RootCells + i << ": ";
              // cout << LocVerts[k]->GetX() << " ";
              // cout << LocVerts[k]->GetY() << " ";
              // cout << LocVerts[k]->GetZ() << endl;
            }
          } // endfor N_Layers
        }
        else
        {
          Error("Mixed meshes with triangles AND quadrilaterals are not allowed" << endl);
          exit(-1);
        }
      } // endfor N_RootCells
    
      // generate edges
      KMT = new TJoint*[N_RootCells*4*N_Layers];
//       memset(KMT, 0, N_RootCells*4*N_Layers * SizeOfInt); //changed by sashi on 7 Nov 2012
      
      j=N_RootCells*4*N_Layers;
      for(i=0;i<j;i++)
        KMT[i] = NULL;
    
      Bottom = new TBdPlane(1000);
      Bottom->SetParams(0,0,Zmin, 1,0,0, 0,0,-1);
      Top = new TBdPlane(1001);
      Top->SetParams(0,0,DriftZ-Zmin, 1,0,0, 0,0,1); 
    
      for (i=0;i<N_RootCells;i++)
      {
        switch (CellTree[i]->GetType())
        {
          case Tetrahedron: 
               N_E = 3;
               break;
          case Hexahedron: 
          case Brick: 
               N_E = 4;
               break;
          default:
               cerr << "Error in ReadGeo" << endl;
               return -1;
        }
    
        for (j=0;j<N_E;j++)
        {
          a = KVERT[NVE*i + j] - 1;
          b = KVERT[NVE*i + (j+1) % N_E] - 1;
          Part = KNPR[a] - 1;
          Neib = -1;
    
          aux1 = KVEL[a*maxElpV];
          aux2 = KVEL[b*maxElpV];
    
          for (k=1;k<=aux1;k++)
          {
            aux3 = KVEL[a*maxElpV + k];
            if (aux3 == i) continue;
    
            for (l=1;l<=aux2;l++)
              if (aux3 == KVEL[b*maxElpV + l])
              {
                Neib = aux3;
                break;
              }
            if (Neib >= 0) break;
          }
          
          // cout << Neib << "  " << i << endl;
    
          if (Neib > i)
          {
            if( (KNPR[a]) && (KNPR[b]) )
            {
              if( (KNPR[a] == KNPR[b]) && (Interfaces[KNPR[a]-1] < 0) )
              {
                comp = (int) DCORVG[2*a];
                T_a = DCORVG[2*a] - comp;
                T_b = DCORVG[2*b] - comp;

                BdParts[Part]->GetBdComp(comp)->GetXYZofTS(T_a, Lambda[0], Xa, Ya, Za);
                BdParts[Part]->GetBdComp(comp)->GetXYZofTS(T_b, Lambda[0], Xb, Yb, Zb);

                X = 0.5*(Xa+Xb);
                Y = 0.5*(Ya+Yb);
                Z = 0.5*(Za+Zb);
                BdParts[Part]->GetBdComp(comp)->GetTSofXYZ(X, Y, Z, T, S);

                if ((T_a - T)*(T - T_b) < 0)
                {
                  if (ABS(T_a) < 1e-4)
                    T_a = 1.0;
                  else
                    if (ABS(T_a) > 0.9999)
                      T_a = 0.0;

                  if (ABS(T_b) < 1e-4)
                    T_b = 1.0;
                  else
                    if (ABS(T_b) > 0.9999)
                      T_b = 0.0;
                }

                for(k=0;k<N_Layers-1;k++)
                {
                  if(j == N_E-1)
                  {
                    Param1[0] = T_b;
                    Param1[1] = T_a;
                    Param1[2] = T_a;
                    Param1[3] = T_b;
    
                    Param2[0] = Lambda[k];
                    Param2[1] = Lambda[k];
                    Param2[2] = Lambda[k+1];
                    Param2[3] = Lambda[k+1];
                  }
                  else
                  {
                    Param1[0] = T_a;
                    Param1[1] = T_a;
                    Param1[2] = T_b;
                    Param1[3] = T_b;
    
                    Param2[0] = Lambda[k];
                    Param2[1] = Lambda[k+1];
                    Param2[2] = Lambda[k+1];
                    Param2[3] = Lambda[k];
                  }
                  // usual bpoundary joint
                  Joint = new TInterfaceJoint3D(BdParts[Part]->GetBdComp(comp),
                            Param1, Param2, CellTree[k*N_RootCells + i],
                            CellTree[k*N_RootCells + Neib]);
                  KMT[k*4*N_RootCells + i*4 + j] = Joint;
                  // cout << "hier" << T_a << "  " << T_b << endl;
                  // cout << k*N_RootCells + i << " " << j+1 << endl;
                } // endfor k
              }
              else
              {
                for(k=0;k<N_Layers-1;k++)
                {
                  KMT[k*4*N_RootCells + i*4 + j] =
                      new TJointEqN(CellTree[k*N_RootCells + i],
                                    CellTree[k*N_RootCells + Neib]);
                }
              }
            }
            else
            {
              for(k=0;k<N_Layers-1;k++)
              {
                KMT[k*4*N_RootCells + i*4 + j] =
                    new TJointEqN(CellTree[k*N_RootCells + i],
                                  CellTree[k*N_RootCells + Neib]);
              }
              // cout << "KMT2: " << i*4+j << endl;
            }
    
            for(k=0;k<N_Layers-1;k++)
            {
              CellTree[k*N_RootCells + i]->SetJoint(j+1,
                            KMT[k*4*N_RootCells + i*4 + j]);
            }
          }
          else
          {
            if (Neib == -1)
            {
              if (Interfaces[KNPR[a]-1] < 0)
              {
                Error("Error in ReadGeo, line " << __LINE__ << endl);
                exit(-1);
              }
              else
                if (Interfaces[KNPR[b]-1] < 0)
                {
                  Error("Error in ReadGeo, line " << __LINE__ << endl);
                  exit(-1);
                }
                else
                {
                  comp = (int) DCORVG[2*a];
                  T_a = DCORVG[2*a] - comp;
                  T_b = DCORVG[2*b] - comp;
                }
    
              
              if (T_b < T_a)
              {
                T_b = 1.;
                // cout << T_a << " ... " << T_b << endl;
              }
    
              if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
              {
                Error("Error in ReadGeo, line " << __LINE__ << endl);
                exit(-1);
              }
              else
              {
                for(k=0;k<N_Layers-1;k++)
                {
                  if(j == N_E-1)
                  {
                    Param1[0] = T_b;
                    Param1[1] = T_a;
                    Param1[2] = T_a;
                    Param1[3] = T_b;
    
                    Param2[0] = Lambda[k];
                    Param2[1] = Lambda[k];
                    Param2[2] = Lambda[k+1];
                    Param2[3] = Lambda[k+1];
                  }
                  else
                  {
                    Param1[0] = T_a;
                    Param1[1] = T_a;
                    Param1[2] = T_b;
                    Param1[3] = T_b;
    
                    Param2[0] = Lambda[k];
                    Param2[1] = Lambda[k+1];
                    Param2[2] = Lambda[k+1];
                    Param2[3] = Lambda[k];
                  }
                  /*
                  // test for periodic boundary conditions
                  if(TDatabase::ParamDB->P5 == 1234567)
                  {
                    OutPut("cell: " << i << " " << j << endl);
                    if(! (CellTree[k*N_RootCells + i]->GetJoint(j+1) ) )
                    {
                      Joint = new TPeriodicJoint(CellTree[k*N_RootCells+i],
                                            CellTree[k*N_RootCells+i+3]);
                      CellTree[k*N_RootCells+i]->SetJoint(j+1, Joint);
                      CellTree[k*N_RootCells+i+3]->SetJoint(j+1, Joint);
                      Joint->SetMapType(1);
                      OutPut("making periodic joint" << endl);
                    }
                  }
                  else
                  */
                  {
                    // usual bpoundary joint
                    Joint = new TBoundFace(BdParts[Part]->GetBdComp(comp),
                              Param1, Param2);
                    CellTree[k*N_RootCells + i]->SetJoint(j+1, Joint);
                    // cout << "hier" << T_a << "  " << T_b << endl;
                    // cout << k*N_RootCells + i << " " << j+1 << endl;
                  }
                }
              }
    
            }
            else
            {
              for (k=0;k<4;k++)
              {
                Joint = KMT[Neib*4 + k];
                if (Joint != NULL)
                {
                  aux1 = Neib*4 + k;
                  JNeib1 = KMT[aux1]->GetNeighbour(0);
                  JNeib2 = KMT[aux1]->GetNeighbour(1);
    
                  if ( (JNeib1 == CellTree[Neib] && JNeib2 == CellTree[i]) ||
                       (JNeib1 == CellTree[i] && JNeib2 == CellTree[Neib])) break;
                }
              }
    
              for(l=0;l<N_Layers-1;l++)
              {
                CellTree[l*N_RootCells + i]->SetJoint(j+1,
                              KMT[l*4*N_RootCells + Neib*4 + k]);
              }
              // cout << "dort" << endl;
            }
          }
        } // endfor N_E
    
        // create top and bottom joints
        Param1[0] = NewVertices[KVERT[NVE*i    ]-1]->GetX();
        Param1[1] = NewVertices[KVERT[NVE*i + 1]-1]->GetX();
        Param1[2] = NewVertices[KVERT[NVE*i + 2]-1]->GetX();
        Param1[3] = NewVertices[KVERT[NVE*i + 3]-1]->GetX();
        Param2[0] = NewVertices[KVERT[NVE*i    ]-1]->GetY();
        Param2[1] = NewVertices[KVERT[NVE*i + 1]-1]->GetY();
        Param2[2] = NewVertices[KVERT[NVE*i + 2]-1]->GetY();
        Param2[3] = NewVertices[KVERT[NVE*i + 3]-1]->GetY();
        Joint = new TBoundFace(Bottom, Param1, Param2);
        CellTree[i]->SetJoint(0, Joint);
    
        Param1[0] = NewVertices[KVERT[NVE*i    ]-1]->GetX()+DriftX;
        Param1[3] = NewVertices[KVERT[NVE*i + 1]-1]->GetX()+DriftX;
        Param1[2] = NewVertices[KVERT[NVE*i + 2]-1]->GetX()+DriftX;
        Param1[1] = NewVertices[KVERT[NVE*i + 3]-1]->GetX()+DriftX;
        Param2[0] = -(NewVertices[KVERT[NVE*i    ]-1]->GetY()+DriftY);
        Param2[3] = -(NewVertices[KVERT[NVE*i + 1]-1]->GetY()+DriftY);
        Param2[2] = -(NewVertices[KVERT[NVE*i + 2]-1]->GetY()+DriftY);
        Param2[1] = -(NewVertices[KVERT[NVE*i + 3]-1]->GetY()+DriftY);
        Joint = new TBoundFace(Top, Param1, Param2);
        CellTree[(N_Layers-2)*N_RootCells + i]->SetJoint(5, Joint);
        
        // create horizontal joints
        for(k=1;k<N_Layers-1;k++)
        {
          Joint = new TJointEqN(CellTree[(k-1)*N_RootCells+i],
                                CellTree[ k   *N_RootCells+i]);
          CellTree[(k-1)*N_RootCells+i]->SetJoint(5, Joint);
          CellTree[ k   *N_RootCells+i]->SetJoint(0, Joint);
        }
      }
    
      N_RootCells *= (N_Layers-1);
      delete [] KMT; // [] added by sashi
    break;

    case 3: // triangles in 2d mesh
      // generate cells
      CellTree = new TBaseCell*[3*N_RootCells*N_Layers];
      OrderArray = new int[N_RootCells];
      KMTupper = new int[N_RootCells*N_Layers*3];
      KMTlower = new int[N_RootCells*N_Layers*3];
    
      for (i=0;i<N_RootCells;i++)
      {
        // sort bottom vertices according to pointer
        vert0 = NewVertices[KVERT[NVE*i    ]-1];
        vert1 = NewVertices[KVERT[NVE*i + 1]-1];
        vert2 = NewVertices[KVERT[NVE*i + 2]-1];

        if( (vert0 < vert1) && (vert1 < vert2) )
          SortOrder = 0; // 0-1-2
        if( (vert0 < vert2) && (vert2 < vert1) )
          SortOrder = 1; // 0-2-1
        if( (vert1 < vert0) && (vert0 < vert2) )
          SortOrder = 2; // 1-0-2
        if( (vert1 < vert2) && (vert2 < vert0) )
          SortOrder = 3;  // 1-2-0
        if( (vert2 < vert0) && (vert0 < vert1) )
          SortOrder = 4; // 2-0-1
        if( (vert2 < vert1) && (vert1 < vert0) )
          SortOrder = 5; // 2-1-0

        OrderArray[i] = SortOrder;
        // cout << "SortOrder: " << SortOrder << endl;
        for(j=1;j<N_Layers;j++)
        {
          switch(SortOrder)
          {
            case 0:
              k = (j-1)*N_Vertices + KVERT[NVE*i    ]-1;
              LocVerts[0] = NewVertices[k];
              k += N_Vertices;
              LocVerts[3] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i + 1]-1;
              LocVerts[1] = NewVertices[k];
              k += N_Vertices;
              LocVerts[4] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i + 2]-1;
              LocVerts[2] = NewVertices[k];
              k += N_Vertices;
              LocVerts[5] = NewVertices[k];

              k0 = 0; k1 = 1; k2 = 2; k3 = 3;

              KMTupper[((j-1)*N_RootCells+i)*3 + 0] = 2*32 + 1;
              KMTupper[((j-1)*N_RootCells+i)*3 + 1] = 1*32 + 2;
              KMTupper[((j-1)*N_RootCells+i)*3 + 2] = 2*32 + 3;

              KMTlower[((j-1)*N_RootCells+i)*3 + 0] = 1*32 + 1;
              KMTlower[((j-1)*N_RootCells+i)*3 + 1] = 0*32 + 2;
              KMTlower[((j-1)*N_RootCells+i)*3 + 2] = 0*32 + 3;
            break;

            case 1:
              k = (j-1)*N_Vertices + KVERT[NVE*i + 0]-1;
              LocVerts[0] = NewVertices[k];
              k += N_Vertices;
              LocVerts[3] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i + 2]-1;
              LocVerts[1] = NewVertices[k];
              k += N_Vertices;
              LocVerts[4] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i + 1]-1;
              LocVerts[2] = NewVertices[k];
              k += N_Vertices;
              LocVerts[5] = NewVertices[k];

              k0 = 1; k1 = 1; k2 = 3; k3 = 2;

              KMTupper[((j-1)*N_RootCells+i)*3 + 2] = 2*32 + 3;
              KMTupper[((j-1)*N_RootCells+i)*3 + 1] = 1*32 + 2;
              KMTupper[((j-1)*N_RootCells+i)*3 + 0] = 2*32 + 1;

              KMTlower[((j-1)*N_RootCells+i)*3 + 2] = 1*32 + 0;
              KMTlower[((j-1)*N_RootCells+i)*3 + 1] = 0*32 + 3;
              KMTlower[((j-1)*N_RootCells+i)*3 + 0] = 0*32 + 2;
            break;

            case 2:
              k = (j-1)*N_Vertices + KVERT[NVE*i + 1]-1;
              LocVerts[0] = NewVertices[k];
              k += N_Vertices;
              LocVerts[3] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i    ]-1;
              LocVerts[1] = NewVertices[k];
              k += N_Vertices;
              LocVerts[4] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i + 2]-1;
              LocVerts[2] = NewVertices[k];
              k += N_Vertices;
              LocVerts[5] = NewVertices[k];

              k0 = 1; k1 = 0; k2 = 3; k3 = 2;

              KMTupper[((j-1)*N_RootCells+i)*3 + 0] = 2*32 + 3;
              KMTupper[((j-1)*N_RootCells+i)*3 + 2] = 1*32 + 2;
              KMTupper[((j-1)*N_RootCells+i)*3 + 1] = 2*32 + 1;

              KMTlower[((j-1)*N_RootCells+i)*3 + 0] = 1*32 + 0;
              KMTlower[((j-1)*N_RootCells+i)*3 + 2] = 0*32 + 3;
              KMTlower[((j-1)*N_RootCells+i)*3 + 1] = 0*32 + 2;
            break;

            case 3:
              k = (j-1)*N_Vertices + KVERT[NVE*i + 1]-1;
              LocVerts[0] = NewVertices[k];
              k += N_Vertices;
              LocVerts[3] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i + 2]-1;
              LocVerts[1] = NewVertices[k];
              k += N_Vertices;
              LocVerts[4] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i    ]-1;
              LocVerts[2] = NewVertices[k];
              k += N_Vertices;
              LocVerts[5] = NewVertices[k];

              k0 = 0; k1 = 1; k2 = 2; k3 = 3;

              KMTupper[((j-1)*N_RootCells+i)*3 + 1] = 2*32 + 1;
              KMTupper[((j-1)*N_RootCells+i)*3 + 2] = 1*32 + 2;
              KMTupper[((j-1)*N_RootCells+i)*3 + 0] = 2*32 + 3;

              KMTlower[((j-1)*N_RootCells+i)*3 + 1] = 1*32 + 1;
              KMTlower[((j-1)*N_RootCells+i)*3 + 2] = 0*32 + 2;
              KMTlower[((j-1)*N_RootCells+i)*3 + 0] = 0*32 + 3;
            break;

            case 4:
              k = (j-1)*N_Vertices + KVERT[NVE*i + 2]-1;
              LocVerts[0] = NewVertices[k];
              k += N_Vertices;
              LocVerts[3] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i    ]-1;
              LocVerts[1] = NewVertices[k];
              k += N_Vertices;
              LocVerts[4] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i + 1]-1;
              LocVerts[2] = NewVertices[k];
              k += N_Vertices;
              LocVerts[5] = NewVertices[k];

              k0 = 0; k1 = 1; k2 = 2; k3 = 3;

              KMTupper[((j-1)*N_RootCells+i)*3 + 2] = 2*32 + 1;
              KMTupper[((j-1)*N_RootCells+i)*3 + 0] = 1*32 + 2;
              KMTupper[((j-1)*N_RootCells+i)*3 + 1] = 2*32 + 3;

              KMTlower[((j-1)*N_RootCells+i)*3 + 2] = 1*32 + 1;
              KMTlower[((j-1)*N_RootCells+i)*3 + 0] = 0*32 + 2;
              KMTlower[((j-1)*N_RootCells+i)*3 + 1] = 0*32 + 3;
            break;

            case 5:
              k = (j-1)*N_Vertices + KVERT[NVE*i + 2]-1;
              LocVerts[0] = NewVertices[k];
              k += N_Vertices;
              LocVerts[3] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i + 1]-1;
              LocVerts[1] = NewVertices[k];
              k += N_Vertices;
              LocVerts[4] = NewVertices[k];
      
              k = (j-1)*N_Vertices + KVERT[NVE*i    ]-1;
              LocVerts[2] = NewVertices[k];
              k += N_Vertices;
              LocVerts[5] = NewVertices[k];

              k0 = 1; k1 = 0; k2 = 3; k3 = 2;

              KMTupper[((j-1)*N_RootCells+i)*3 + 1] = 2*32 + 3;
              KMTupper[((j-1)*N_RootCells+i)*3 + 0] = 1*32 + 2;
              KMTupper[((j-1)*N_RootCells+i)*3 + 2] = 2*32 + 1;

              KMTlower[((j-1)*N_RootCells+i)*3 + 1] = 1*32 + 0;
              KMTlower[((j-1)*N_RootCells+i)*3 + 0] = 0*32 + 3;
              KMTlower[((j-1)*N_RootCells+i)*3 + 2] = 0*32 + 2;
            break;

          }
  
          if(SortOrder == 0 || SortOrder == 3 || SortOrder == 4)
          {
            k = ((j-1)*N_RootCells +i)*3;
            CellTree[k] = 
                 new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);
            CurrCell = CellTree[k];
            CurrCell->SetVertex(0, LocVerts[0]);
            CurrCell->SetVertex(1, LocVerts[1]);
            CurrCell->SetVertex(k2, LocVerts[2]);
            CurrCell->SetVertex(k3, LocVerts[5]);
            CurrCell->SetClipBoard(k);
  
            /*
            CurrCell->GetVertex(0)->GetCoords(x0, y0, z0);
            CurrCell->GetVertex(1)->GetCoords(x1, y1, z1);
            CurrCell->GetVertex(2)->GetCoords(x2, y2, z2);
            CurrCell->GetVertex(3)->GetCoords(x3, y3, z3);
  
            x1 -= x0; y1 -= y0; z1 -= z0;
            x2 -= x0; y2 -= y0; z2 -= z0;
            x3 -= x0; y3 -= y0; z3 -= z0;
  
            det =   x1*y2*z3 + y1*z2*x3 + z1*x2*y3
                 -( z1*y2*x3 + y1*x2*z3 + x1*z2*y3);
            cout << "cell nr: " << k << " det: " << det << endl,
            */
  
            k++;
            CellTree[k] = 
                 new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);
            CurrCell = CellTree[k];
            CurrCell->SetVertex(0, LocVerts[0]);
            CurrCell->SetVertex(1, LocVerts[1]);
            CurrCell->SetVertex(k2, LocVerts[5]);
            CurrCell->SetVertex(k3, LocVerts[4]);
            CurrCell->SetClipBoard(k);
  
            /*
            CurrCell->GetVertex(0)->GetCoords(x0, y0, z0);
            CurrCell->GetVertex(1)->GetCoords(x1, y1, z1);
            CurrCell->GetVertex(2)->GetCoords(x2, y2, z2);
            CurrCell->GetVertex(3)->GetCoords(x3, y3, z3);
  
            x1 -= x0; y1 -= y0; z1 -= z0;
            x2 -= x0; y2 -= y0; z2 -= z0;
            x3 -= x0; y3 -= y0; z3 -= z0;
  
            det =   x1*y2*z3 + y1*z2*x3 + z1*x2*y3
                 -( z1*y2*x3 + y1*x2*z3 + x1*z2*y3);
            cout << "cell nr: " << k << " det: " << det << endl,
            */
  
            k++;
            CellTree[k] = 
                 new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);
            CurrCell = CellTree[k];
            CurrCell->SetVertex(0, LocVerts[0]);
            CurrCell->SetVertex(1, LocVerts[4]);
            CurrCell->SetVertex(k2, LocVerts[5]);
            CurrCell->SetVertex(k3, LocVerts[3]);
            CurrCell->SetClipBoard(k);
  
            /*
            CurrCell->GetVertex(0)->GetCoords(x0, y0, z0);
            CurrCell->GetVertex(1)->GetCoords(x1, y1, z1);
            CurrCell->GetVertex(2)->GetCoords(x2, y2, z2);
            CurrCell->GetVertex(3)->GetCoords(x3, y3, z3);
  
            x1 -= x0; y1 -= y0; z1 -= z0;
            x2 -= x0; y2 -= y0; z2 -= z0;
            x3 -= x0; y3 -= y0; z3 -= z0;
  
            det =   x1*y2*z3 + y1*z2*x3 + z1*x2*y3
                 -( z1*y2*x3 + y1*x2*z3 + x1*z2*y3);
            cout << "cell nr: " << k << " det: " << det << endl;
            */
  
            Joint = new TJointEqN(CellTree[((j-1)*N_RootCells +i)*3 + 0],
                                  CellTree[((j-1)*N_RootCells +i)*3 + 1]);
            CellTree[((j-1)*N_RootCells +i)*3 + 1]->SetJoint(0, Joint);
            CellTree[((j-1)*N_RootCells +i)*3 + 0]->SetJoint(1, Joint);

            Joint = new TJointEqN(CellTree[((j-1)*N_RootCells +i)*3 + 1],
                                  CellTree[((j-1)*N_RootCells +i)*3 + 2]);
            CellTree[((j-1)*N_RootCells +i)*3 + 1]->SetJoint(3, Joint);
            CellTree[((j-1)*N_RootCells +i)*3 + 2]->SetJoint(0, Joint);
          }
          else
          {
            k = ((j-1)*N_RootCells +i)*3;
            CellTree[k] = 
                 new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);
            CurrCell = CellTree[k];
            CurrCell->SetVertex(0, LocVerts[1]);
            CurrCell->SetVertex(1, LocVerts[0]);
            CurrCell->SetVertex(2, LocVerts[2]);
            CurrCell->SetVertex(3, LocVerts[5]);
            CurrCell->SetClipBoard(k);
  
            /*
            CurrCell->GetVertex(0)->GetCoords(x0, y0, z0);
            CurrCell->GetVertex(1)->GetCoords(x1, y1, z1);
            CurrCell->GetVertex(2)->GetCoords(x2, y2, z2);
            CurrCell->GetVertex(3)->GetCoords(x3, y3, z3);
  
            x1 -= x0; y1 -= y0; z1 -= z0;
            x2 -= x0; y2 -= y0; z2 -= z0;
            x3 -= x0; y3 -= y0; z3 -= z0;
  
            det =   x1*y2*z3 + y1*z2*x3 + z1*x2*y3
                 -( z1*y2*x3 + y1*x2*z3 + x1*z2*y3);
            cout << "cell nr: " << k << " det: " << det << endl,
            */
  
            k++;
            CellTree[k] = 
                 new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);
            CurrCell = CellTree[k];
            CurrCell->SetVertex(0, LocVerts[0]);
            CurrCell->SetVertex(1, LocVerts[1]);
            CurrCell->SetVertex(2, LocVerts[4]);
            CurrCell->SetVertex(3, LocVerts[5]);
            CurrCell->SetClipBoard(k);
  
            /*
            CurrCell->GetVertex(0)->GetCoords(x0, y0, z0);
            CurrCell->GetVertex(1)->GetCoords(x1, y1, z1);
            CurrCell->GetVertex(2)->GetCoords(x2, y2, z2);
            CurrCell->GetVertex(3)->GetCoords(x3, y3, z3);
  
            x1 -= x0; y1 -= y0; z1 -= z0;
            x2 -= x0; y2 -= y0; z2 -= z0;
            x3 -= x0; y3 -= y0; z3 -= z0;
  
            det =   x1*y2*z3 + y1*z2*x3 + z1*x2*y3
                 -( z1*y2*x3 + y1*x2*z3 + x1*z2*y3);
            cout << "cell nr: " << k << " det: " << det << endl,
            */
  
            k++;
            CellTree[k] = 
                 new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);
            CurrCell = CellTree[k];
            CurrCell->SetVertex(0, LocVerts[0]);
            CurrCell->SetVertex(1, LocVerts[5]);
            CurrCell->SetVertex(2, LocVerts[4]);
            CurrCell->SetVertex(3, LocVerts[3]);
            CurrCell->SetClipBoard(k);
  
            /*
            CurrCell->GetVertex(0)->GetCoords(x0, y0, z0);
            CurrCell->GetVertex(1)->GetCoords(x1, y1, z1);
            CurrCell->GetVertex(2)->GetCoords(x2, y2, z2);
            CurrCell->GetVertex(3)->GetCoords(x3, y3, z3);
  
            x1 -= x0; y1 -= y0; z1 -= z0;
            x2 -= x0; y2 -= y0; z2 -= z0;
            x3 -= x0; y3 -= y0; z3 -= z0;
  
            det =   x1*y2*z3 + y1*z2*x3 + z1*x2*y3
                 -( z1*y2*x3 + y1*x2*z3 + x1*z2*y3);
            cout << "cell nr: " << k << " det: " << det << endl;
            */
  
            Joint = new TJointEqN(CellTree[((j-1)*N_RootCells +i)*3 + 0],
                                  CellTree[((j-1)*N_RootCells +i)*3 + 1]);
            CellTree[((j-1)*N_RootCells +i)*3 + 1]->SetJoint(1, Joint);
            CellTree[((j-1)*N_RootCells +i)*3 + 0]->SetJoint(1, Joint);

            Joint = new TJointEqN(CellTree[((j-1)*N_RootCells +i)*3 + 1],
                                  CellTree[((j-1)*N_RootCells +i)*3 + 2]);
            CellTree[((j-1)*N_RootCells +i)*3 + 1]->SetJoint(3, Joint);
            CellTree[((j-1)*N_RootCells +i)*3 + 2]->SetJoint(0, Joint);
          }

        } // endfor N_Layers
      } // endfor N_RootCells
    
      Bottom = new TBdPlane(1000);
      Bottom->SetParams(0,0,0, 1,0,0, 0,0,-1);
      Top = new TBdPlane(1001);
      Top->SetParams(0,0,DriftZ, 1,0,0, 0,0,1); 
    
      N_E = 3;
      for (i=0;i<N_RootCells;i++)
      {
        SortOrder = OrderArray[i];
        // cout << "i: " << i << " " << SortOrder << endl;
        for (j=0;j<N_E;j++)
        {
          a = KVERT[NVE*i + j] - 1;
          b = KVERT[NVE*i + (j+1) % N_E] - 1;
          Part = KNPR[a] - 1;
          Neib = -1;
    
          aux1 = KVEL[a*maxElpV];
          aux2 = KVEL[b*maxElpV];
    
          for (k=1;k<=aux1;k++)
          {
            aux3 = KVEL[a*maxElpV + k];
            if (aux3 == i) continue;
    
            for (l=1;l<=aux2;l++)
              if (aux3 == KVEL[b*maxElpV + l])
              {
                Neib = aux3;
                break;
              }
            if (Neib >= 0) break;
          }

          // on which local joint of Neib is i
          if(Neib>=0)
          {
            vert = NewVertices[KVERT[3*i+(j+1)%3]];
            for(l=0;l<3;l++)
              if(vert == NewVertices[KVERT[3*Neib+l]]) break;
          }
          
          // cout << "Neib: " << Neib << " i: " << i << endl;
    
          if (Neib > i)
          {
            if( (KNPR[a]) && (KNPR[b]) )
            {
              if( (KNPR[a] == KNPR[b]) && (Interfaces[KNPR[a]-1] < 0) )
              {
                Error("Error in ReadGeo, line " << __LINE__ << endl);
                exit(-1);
              }
              else
              {
                // set neighbours
                // cout << "hier: ";
                // cout << i << " " << j << " " << Neib << " " << l << endl;
                // cout << "upper: " << KMTupper[3*i+j] << " " << KMTupper[3*Neib+l] << endl; 
                k0 = KMTupper[3*i+j] % 32;
                k1 = (KMTupper[3*i+j] / 32) + 3*i;
                k2 = KMTupper[3*Neib+l] % 32;
                k3 = (KMTupper[3*Neib+l] / 32) + 3*Neib;
                // cout << k1 << " " << k0 << " ::: " << k3 << " " << k2 << endl;
                for(k=0;k<N_Layers-1;k++)
                {
                  Joint = new TJointEqN(CellTree[k1],CellTree[k3]);
                  CellTree[k1]->SetJoint(k0, Joint);
                  CellTree[k3]->SetJoint(k2, Joint);
                  k1 += 3*N_RootCells;
                  k3 += 3*N_RootCells;
                } // endfor k

                // cout << "lower: " << KMTlower[3*i+j] << " " << KMTlower[3*Neib+l] << endl; 
                k0 = KMTlower[3*i+j] % 32;
                k1 = (KMTlower[3*i+j] / 32) + 3*i;
                k2 = KMTlower[3*Neib+l] % 32;
                k3 = (KMTlower[3*Neib+l] / 32) + 3*Neib;
                // cout << k1 << " " << k0 << " ::: " << k3 << " " << k2 << endl;
                for(k=0;k<N_Layers-1;k++)
                {
                  Joint = new TJointEqN(CellTree[k1],CellTree[k3]);
                  CellTree[k1]->SetJoint(k0, Joint);
                  CellTree[k3]->SetJoint(k2, Joint);
                  k1 += 3*N_RootCells;
                  k3 += 3*N_RootCells;
                } // endfor k
              }
            }
            else
            {
              // set neighbours
              // cout << "KMT2: ";
              // cout << i << " " << j << " " << Neib << " " << l << endl;
              // cout << "upper: " << KMTupper[3*i+j] << " " << KMTupper[3*Neib+l] << endl; 
              k0 = KMTupper[3*i+j] % 32;
              k1 = (KMTupper[3*i+j] / 32) + 3*i;
              k2 = KMTupper[3*Neib+l] % 32;
              k3 = (KMTupper[3*Neib+l] / 32) + 3*Neib;
              // cout << k1 << " " << k0 << " ::: " << k3 << " " << k2 << endl;
              for(k=0;k<N_Layers-1;k++)
              {
                Joint = new TJointEqN(CellTree[k1],CellTree[k3]);
                CellTree[k1]->SetJoint(k0, Joint);
                CellTree[k3]->SetJoint(k2, Joint);
                k1 += 3*N_RootCells;
                k3 += 3*N_RootCells;
              } // endfor k

              // cout << "lower: " << KMTlower[3*i+j] << " " << KMTlower[3*Neib+l] << endl; 
              k0 = KMTlower[3*i+j] % 32;
              k1 = (KMTlower[3*i+j] / 32) + 3*i;
              k2 = KMTlower[3*Neib+l] % 32;
              k3 = (KMTlower[3*Neib+l] / 32) + 3*Neib;
              // cout << k1 << " " << k0 << " ::: " << k3 << " " << k2 << endl;
              for(k=0;k<N_Layers-1;k++)
              {
                Joint = new TJointEqN(CellTree[k1],CellTree[k3]);
                CellTree[k1]->SetJoint(k0, Joint);
                CellTree[k3]->SetJoint(k2, Joint);
                k1 += 3*N_RootCells;
                k3 += 3*N_RootCells;
              } // endfor k
            }
          }
          else
          {
            if (Neib == -1)
            {
              if (Interfaces[KNPR[a]-1] < 0 || Interfaces[KNPR[b]-1] < 0)
              {
                Error("Error in ReadGeo, line " << __LINE__ << endl);
                exit(-1);
              }
              else
              {
                comp = (int) DCORVG[2*a];
                T_a = DCORVG[2*a] - comp;
                T_b = DCORVG[2*b] - comp;
              }
              
              if (T_b < T_a)
              {
                T_b = 1.;
                // cout << "AT ";
              }
              // cout << T_a << " ... " << T_b << endl;

              if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
              {
                Error("Error in ReadGeo, line " << __LINE__ << endl);
                exit(-1);
              }
              else
              {
                // upper cells
                // cout << KMTupper[3*i+j] << endl; 
                k0 = KMTupper[3*i+j] % 32;
                k1 = (KMTupper[3*i+j] / 32) + 3*i;
                // cout << "B:" << k1 << " " << k0 << endl;
                CurrCell = CellTree[k1];
                switch(k0)
                {
                  case 0:
                    Error("This case should not appear! " << __LINE__ << endl);
                    exit(-1);
                  break;

                  case 1:
                    Param1[0] = T_a;
                    Param1[1] = T_a;
                    Param1[2] = T_b;
                    i0 = 0; i1 = 1; i2 = 1;
                  break;

                  case 2:
                    switch(SortOrder)
                    {
                      case 0:
                      case 3:
                      case 4:
                        Param1[0] = T_b;
                        Param1[1] = T_a;
                        Param1[2] = T_a;
                        i0 = 1; i1 = 0; i2 = 1;
                      break;

                      case 1:
                      case 2:
                      case 5:
                        Param1[0] = T_b;
                        Param1[1] = T_b;
                        Param1[2] = T_a;
                        i0 = 1; i1 = 0; i2 = 1;
                      break;
                    }
                  break;

                  case 3:
                    Param1[0] = T_b;
                    Param1[1] = T_a;
                    Param1[2] = T_b;
                    i0 = 0; i1 = 1; i2 = 1;
                  break;
                }
                                        
                for(k=0;k<N_Layers-1;k++)
                {
                  Param2[0] = Lambda[k+i0];
                  Param2[1] = Lambda[k+i1];
                  Param2[2] = Lambda[k+i2];
                  Joint = new TBoundFace(BdParts[Part]->GetBdComp(comp),
                            Param1, Param2);
                  CellTree[k*N_RootCells*3 + k1]->SetJoint(k0, Joint);
                  // cout << "upper: " << k*N_RootCells*3 + k1 << " joint: " << k0 << endl;
                  // cout << Param1[0] << " " << Param2[0] << endl;
                  // cout << Param1[1] << " " << Param2[1] << endl;
                  // cout << Param1[2] << " " << Param2[2] << endl;
                } // endfor k

                // lower cells
                // cout << KMTlower[3*i+j] << endl; 
                k2 = KMTlower[3*i+j] % 32;
                k3 = (KMTlower[3*i+j] / 32) + 3*i;
                // cout << "B:" << k3 << " " << k2 << endl;
                CurrCell = CellTree[k3];
                switch(k2)
                {
                  case 0:
                    switch(SortOrder)
                    {
                      case 0:
                      case 3:
                      case 4:
                        Error("This case should not appear! " << __LINE__ << endl);
                        exit(-1);
                      break;

                      case 1:
                      case 2:
                      case 5:
                        Param1[0] = T_b;
                        Param1[1] = T_a;
                        Param1[2] = T_a;
                        i0 = 0; i1 = 0; i2 = 1;
                      break;
                    }
                  break;

                  case 1:
                    switch(SortOrder)
                    {
                      case 0:
                      case 3:
                      case 4:
                        Param1[0] = T_a;
                        Param1[1] = T_b;
                        Param1[2] = T_b;
                        i0 = 0; i1 = 1; i2 = 0;
                      break;

                      case 1:
                      case 2:
                      case 5:
                        Error("This case should not appear! " << __LINE__ << endl);
                        exit(-1);
                      break;
                    }
                  break;

                  case 2:
                    Param1[0] = T_b;
                    Param1[1] = T_a;
                    Param1[2] = T_b;
                    i0 = 0; i1 = 0; i2 = 1;
                  break;

                  case 3:
                    Param1[0] = T_b;
                    Param1[1] = T_a;
                    Param1[2] = T_a;
                    i0 = 0; i1 = 0; i2 = 1;
                  break;
                }
                                        
                for(k=0;k<N_Layers-1;k++)
                {
                  Param2[0] = Lambda[k+i0];
                  Param2[1] = Lambda[k+i1];
                  Param2[2] = Lambda[k+i2];
                  Joint = new TBoundFace(BdParts[Part]->GetBdComp(comp),
                            Param1, Param2);
                  CellTree[k*N_RootCells*3 + k3]->SetJoint(k2, Joint);
                  // cout << "lower: " << k*N_RootCells*3 + k3 << " joint: " << k2 << endl;
                  // cout << Param1[0] << " " << Param2[0] << endl;
                  // cout << Param1[1] << " " << Param2[1] << endl;
                  // cout << Param1[2] << " " << Param2[2] << endl;
                } // endfor k
              } // endif FreeBoundary
            } // endif Neib == -1
            // else
            // {
            //   cout << "already done" << endl;
            // }
          } // endif
        } // endfor N_E
    
        // create top and bottom joints
        CurrCell = CellTree[i*3];
        CurrCell->GetVertex(0)->GetCoords(x0, y0, z0);
        CurrCell->GetVertex(1)->GetCoords(x1, y1, z1);
        CurrCell->GetVertex(2)->GetCoords(x2, y2, z2);
        Bottom->GetTSofXYZ(x0, y0, z0, Param1[0], Param2[0]);
        Bottom->GetTSofXYZ(x1, y1, z1, Param1[1], Param2[1]);
        Bottom->GetTSofXYZ(x2, y2, z2, Param1[2], Param2[2]);
        Joint = new TBoundFace(Bottom, Param1, Param2);
        CellTree[i*3]->SetJoint(0, Joint);
    
        CurrCell = CellTree[((N_Layers-2)*N_RootCells + i)*3+2];
        CurrCell->GetVertex(2)->GetCoords(x0, y0, z0);
        CurrCell->GetVertex(1)->GetCoords(x1, y1, z1);
        CurrCell->GetVertex(3)->GetCoords(x2, y2, z2);
        Top->GetTSofXYZ(x0, y0, z0, Param1[0], Param2[0]);
        Top->GetTSofXYZ(x1, y1, z1, Param1[1], Param2[1]);
        Top->GetTSofXYZ(x2, y2, z2, Param1[2], Param2[2]);
        Joint = new TBoundFace(Top, Param1, Param2);
        CellTree[((N_Layers-2)*N_RootCells + i)*3+2]->SetJoint(2, Joint);
        
        // create horizontal joints
        for(k=1;k<N_Layers-1;k++)
        {
          Joint = new TJointEqN(CellTree[((k-1)*N_RootCells+i)*3+2],
                                CellTree[( k   *N_RootCells+i)*3+0]);
          CellTree[((k-1)*N_RootCells+i)*3+2]->SetJoint(2, Joint);
          CellTree[( k   *N_RootCells+i)*3+0]->SetJoint(0, Joint);
        }
      }
    
      N_RootCells *= (N_Layers-1)*3;
      delete KMTupper;
      delete KMTlower;
    break;

    default:
      Error("wrong NVE! Only 3 and 4 are allowed!" << endl);
      exit(-1);
  } // endswitch(NVE)

//   cout << "N_RootCells: " << N_RootCells << endl;
  for(i=0;i<N_RootCells;i++)
    CellTree[i]->SetClipBoard(i);

  for(i=0;i<N_RootCells;i++)
  {
    CurrCell = CellTree[i];
    k = CurrCell->GetN_Joints();
    for(j=0;j<k;j++)
    {
      // cout << i << " --  " << j;
      // cout << " " << (int)(CurrCell->GetJoint(j)) << endl;
      if( CurrCell->GetJoint(j) != NULL )
        CurrCell->GetJoint(j)->SetMapType();
    }

    // k = CurrCell->GetN_Vertices();
    // for(j=0;j<k;j++)
    //   cout << i << " " << j << " " << CurrCell->GetVertex(j) << endl;
  }

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

/*
  for(i=0;i<N_RootCells;i++)
  {
    CurrCell = CellTree[i];
    cout << "cell number: " << CurrCell->GetClipBoard() << endl;
    k = CurrCell->GetN_Joints();
    for(j=0;j<k;j++)
    {
      cout << "joint: " << j << "  ";
      if(CurrCell->GetJoint(j))
      {
        if(CurrCell->GetJoint(j)->GetNeighbour(CurrCell))
        {
          cout << "N" << 
          CurrCell->GetJoint(j)->GetNeighbour(CurrCell)->GetClipBoard()
          << endl;
        }
        else
        {
          cout << "T" << CurrCell->GetJoint(j)->GetType() << endl;
        }
      }
      else
      {
        cout << endl;
      }
    }
  }
*/

/*
  for(i=0;i<N_RootCells;i++)
  {
    cout << "cell: " << i << endl;
    CurrCell = CellTree[i];
    for(j=0;j<4;j++)
    {
      cout << "joint: " << j << endl;
      if(CurrCell->GetJoint(j)->GetType() == BoundaryFace ||
                CurrCell->GetJoint(j)->GetType() == IsoBoundFace)
      {
        cout << "on boundary" << endl;
        switch(j)
        {
          case 0:
            k0 = 0; k1 = 1; k2 = 2;
          break;

          case 1:
            k0 = 0; k1 = 3; k2 = 1;
          break;

          case 2:
            k0 = 2; k1 = 1; k2 = 3;
          break;

          case 3:
            k0 = 0; k1 = 2; k2 = 3;
          break;
        } // endswitch
        cout << CurrCell->GetVertex(k0) << endl;
        CurrCell->GetVertex(k0)->GetCoords(x0, y0, z0);
        ((TBoundFace *)(CurrCell->GetJoint(j)))->GetParameters(Param1, Param2);
        ((TBoundFace *)(CurrCell->GetJoint(j)))->GetBoundComp()
                ->GetXYZofTS(Param1[0], Param2[0], x1, y1, z1);
        cout << x1 << " " << y1 << " " << z1 << endl;
        if(fabs(x0-x1)+fabs(y0-y1)+fabs(z0-z1) > 1e-8) 
          cout << "error0" << endl;

        cout << CurrCell->GetVertex(k1) << endl;
        CurrCell->GetVertex(k1)->GetCoords(x0, y0, z0);
        ((TBoundFace *)(CurrCell->GetJoint(j)))->GetParameters(Param1, Param2);
        ((TBoundFace *)(CurrCell->GetJoint(j)))->GetBoundComp()
                ->GetXYZofTS(Param1[1], Param2[1], x1, y1, z1);
        cout << x1 << " " << y1 << " " << z1 << endl;
        if(fabs(x0-x1)+fabs(y0-y1)+fabs(z0-z1) > 1e-8) 
          cout << "error1" << endl;

        cout << CurrCell->GetVertex(k2) << endl;
        CurrCell->GetVertex(k2)->GetCoords(x0, y0, z0);
        ((TBoundFace *)(CurrCell->GetJoint(j)))->GetParameters(Param1, Param2);
        ((TBoundFace *)(CurrCell->GetJoint(j)))->GetBoundComp()
                ->GetXYZofTS(Param1[2], Param2[2], x1, y1, z1);
        cout << x1 << " " << y1 << " " << z1 << endl;
        if(fabs(x0-x1)+fabs(y0-y1)+fabs(z0-z1) > 1e-8) 
          cout << "error2" << endl;
      } // endif
    } // endfor j
  } // endfor i
*/

  // free memory
  delete KVEL;

  delete NewVertices;

  // exit(-1);
  
  return 0;
}

#endif // __2D__

bool TDomain::checkIfxGEO(const char* GEO)
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
        if(TDatabase::ParamDB->SC_VERBOSE>1)
          cout << " *** reading xGEO file (with physical references) ***" << endl;
        isXgeo = true;
      }
      return isXgeo;
  }


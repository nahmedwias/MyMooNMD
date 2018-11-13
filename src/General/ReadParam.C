// =======================================================================
// @(#)ReadParam.C        1.28 06/27/00
//
// Purpose:     read parameter file and
//              read boundary paramerization
//.
// Author:      Volker Behns  22.07.97
// Version:     1.0
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <BdCircle.h>
#include <BdLine.h>
#include <BdSpline.h>
#include <BdPolygon.h>
#include <BdNonUniformSpline.h>
#include <Domain.h>
#include <Database.h>
#include <Joint.h>
#include <MacroCell.h>
#include <MooNMD_Io.h>
#include <fstream>
#include <string.h>

#ifdef __3D__
#include <BdNoPRM.h>
#include <BdPlane.h>
#include <BdWall.h>
#include <BdSphere.h>
#endif
#include <stdlib.h>

#ifdef __2D__
void TDomain::ReadBdParam(std::istream& dat)
#else
void TDomain::ReadBdParam(std::istream& dat, bool& sandwich_flag)
#endif
{
#ifdef _MPI
  int rank; // out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  char line[100];
  int i, j, CurrentBdPart, N_BdComp, CompType, CompID = 0, N_Spls;
#ifdef __2D__
  TBoundComp2D *BdComp;
#else
  TBoundComp3D *BdComp=nullptr;
  TBoundComp2D *BdComp2D=nullptr;
  sandwich_flag = false;
#endif



  // determine dimensions for creating arrays

  // get the number of boundaries (inner and outer) of the domain
  dat.getline (line, 99);
  dat >> N_BoundParts;
  dat.getline (line, 99);

  BdParts = new TBoundPart*[N_BoundParts];
  Interfaces = new int[N_BoundParts];
  N_BoundComps = 0;
  // StartBdCompID: index in the bdpart list where the new BdPart starts
  StartBdCompID = new int[N_BoundParts + 1];
  StartBdCompID[0] = 0; // set the first to 0

  for (i=0;i<N_BoundParts;i++)
  {
    dat.getline (line, 99); // IBCT
    dat >> CurrentBdPart;
    dat.getline (line, 99);
    Interfaces[i] = CurrentBdPart;
    CurrentBdPart = ABS(CurrentBdPart); // it can be negative (for orientation)
    if (i+1 != CurrentBdPart)
    {
#ifdef _MPI
      if(rank==0)
#endif
        cerr << "Error: different number of boundary part" << endl;
      cerr << "CurrentBdPart " << i << "  " << CurrentBdPart << endl;
      exit(-1);
    }

    // get number of components f the bdpart i
    dat.getline (line, 99);
    dat >> N_BdComp;
    dat.getline (line, 99);

    BdParts[i] = new TBoundPart(N_BdComp);
    N_BoundComps += N_BdComp;
    SetStartBdCompID(N_BoundComps, i+1); // set the start ID of the next bdcomp

    dat.getline (line, 99);
    for (j=0;j<N_BdComp;j++)
    {
      dat >> CompType >> N_Spls;  //ITYP NSPLINE NPAR
      dat.getline (line, 99);

#ifdef __2D__
      // 2D types: Line (1), Circle (2), Spline (3), Poygon (4), NonUnif Spline (5)
      switch (abs(CompType))
      {
        case 1: BdComp = new TBdLine(CompID++);
                break;
        case 2: BdComp = new TBdCircle(CompID++);
                break;
        case 3: BdComp = new TBdSpline(CompID++, N_Spls);
                break;
        case 4: BdComp = new TBdPolygon(CompID++, N_Spls);
                break;
        case 5: BdComp = new TBdNonUniformSpline(CompID++, N_Spls);
                break;
        default:
#ifdef _MPI
                if(rank==0)
#endif
                  OutPut("ReadParam.C: Boundary type not implemented" << endl);
                exit(-1);
      }
#else
      // 3D types: Line (1), Circle (2), Spline (3), Poygon (4), NonUnif Spline (5),
      //           Plane (10), Sphere (11)
      switch (abs(CompType))
      {
        case 1: BdComp2D = new TBdLine(CompID++);
                break;
        case 2: BdComp2D = new TBdCircle(CompID++);
                break;
        case 3: BdComp2D = new TBdSpline(CompID++, N_Spls);
                break;
        case 4: BdComp2D = new TBdPolygon(CompID++, N_Spls);
                break;
        case 5: BdComp2D = new TBdNonUniformSpline(CompID++, N_Spls);
                break;
        case 10: BdComp = new TBdPlane(CompID++);
                 break;
        case 11: BdComp = new TBdSphere(CompID++);
                 break;
        case 4711: BdComp = new TBdNoPRM(CompID++); // create grid without PRM file
                   break;
        default:
#ifdef _MPI
                   if(rank==0)
#endif
                     OutPut("ReadParam.C: Boundary type (3D) not implemented" << endl);
                   exit(-1);
      }

      if(abs(CompType)<10)
      {
        BdComp = new TBdWall(CompID-1, BdComp2D);
        sandwich_flag = true;
      }
#endif // 3D
      BdParts[i]->SetBdComp(j, BdComp);

      if(CompType<0)
      {
        BdComp->SetFreeBoundaryStatus(true);
        cout <<i<< " ReadBdParam : " << j << endl;
      }
    }
  }

  dat.getline (line, 99);
  for (i=0;i<N_BoundParts;i++)
  {
    N_BdComp = BdParts[i]->GetN_BdComps();
    for (j=0;j<N_BdComp;j++)
    {
      BdParts[i]->GetBdComp(j)->ReadIn(dat);
    }
  }

  // read HOLES (if any)
  dat.getline (line, 99);
  N_Holes = -12345;
  if (dat.eof())
    N_Holes = 0;
  else
    dat >> N_Holes;

  if(N_Holes == -12345)
    N_Holes = 0;

  dat.getline (line, 99);

  if (N_Holes)
  {
    // coordinates of a point in a hole
    PointInHole = new double[2*N_Holes];

    dat.getline (line, 99);
    for (i=0;i<N_Holes;i++)
    {
      dat >> PointInHole[2*i] >> PointInHole[2*i+1];
      dat.getline (line, 99);
    }
  }
  else
    PointInHole = nullptr;

  dat.getline (line, 99);
  N_Regions = -12345;
  if (dat.eof())
    N_Regions = 0;
  else
    dat >> N_Regions;

  if(N_Regions == -12345)
    N_Regions = 0;

  dat.getline (line, 99);

  // read REGIONS (if any)
  if (N_Regions)
  {
    PointInRegion = new double[4*N_Regions];

    dat.getline (line, 99);
    for (i=0;i<N_Regions;i++)
    {
      dat >> PointInRegion[4*i] >> PointInRegion[4*i+1];
      PointInRegion[4*i+2] = i;
      PointInRegion[4*i+3] = 10000;
      dat.getline (line, 99);
    }
  }
  else
    PointInRegion = nullptr;
}

int TDomain::ReadMapFile(char *MapFile, TDatabase *Database)
{
#ifdef _MPI
  int rank; // out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  char line[100];
  int i, N_Cells, ID;
  std::ifstream dat(MapFile);

  if (!dat)
  {
#ifdef _MPI
    if(rank==0)
#endif
      cerr << "cannot open '" << MapFile << "' for input" << endl;
    exit(-1);
  }

  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);

  // get number of elements
  dat >> N_Cells;
  dat.getline (line, 99);

  if (N_Cells != N_RootCells)
  {
#ifdef _MPI
    if(rank==0)
#endif
      cerr << "Error in ReadMapFile: wrong number of elements!!!" << endl;
    exit(-1);
  }

  dat.getline (line, 99);
  for (i=0;i<N_Cells;i++)
  {
    dat >> ID;
    ((TMacroCell *) CellTree[i])->SetSubGridID(ID);
    dat.getline (line, 99);
  }

  dat.close();

  return 0;
}

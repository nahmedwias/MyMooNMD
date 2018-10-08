// =======================================================================
// @(#)NSE3DMGLevel4.C        1.10 07/03/00
//
// Class:       TNSE3DMGLevel4
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system
//              type 2 (A, B1, B2, B1T, B2T)
//
// Author:      Volker John 25.08.1999
//
// History:     24.08.1999 start of implementation
//              25.08.1999 CellVanka
//              25.08.1999 step length control
//              25.08.1999 ExactSolve
//              26.08.1999 NodalVanka
//
// =======================================================================

#include <NSE_MGLevel4.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <OldSolver.h>
#ifdef __2D__
  #include <FESpace2D.h>
  #include <FEDatabase2D.h>
#endif
#ifdef __3D__
  #include <FESpace3D.h>
  #include <FEDatabase3D.h>
#endif

#include <stdlib.h>
#include <string.h>

#include <LinAlg.h>
#include <ItMethod.h>
#include <FgmresIte.h>
#include <DirectSolver.h>

#include <LocalProjection.h>

/** constructor */
#ifdef __2D__
  TNSE_MGLevel4::TNSE_MGLevel4(int level,
                               TSquareMatrix2D *a11, TSquareMatrix2D *a12,
                               TSquareMatrix2D *a21, TSquareMatrix2D *a22,
                               TMatrix2D *b1, TMatrix2D *b2,
                               TMatrix2D *b1t, TMatrix2D *b2t,
                               double *f1, double *u1,
                               int n_aux, double *al, int velocity_space,
                               int pressure_space, TCollection *Coll,
			       int *dw)
#endif
#ifdef __3D__
  TNSE_MGLevel4::TNSE_MGLevel4(int level,
                               TSquareMatrix3D *a11, TSquareMatrix3D *a12,
                               TSquareMatrix3D *a13, TSquareMatrix3D *a21,
                               TSquareMatrix3D *a22, TSquareMatrix3D *a23,
                               TSquareMatrix3D *a31, TSquareMatrix3D *a32,
                               TSquareMatrix3D *a33,
                               TMatrix3D *b1, TMatrix3D *b2, TMatrix3D *b3,
                               TMatrix3D *b1t, TMatrix3D *b2t, TMatrix3D *b3t,
                               double *f1, double *u1,
                               int n_aux, double *al, int velocity_space,
                               int pressure_space, TCollection *Coll,
                               int *dw)

#endif
  : TNSE_MGLevel(level, f1, u1, n_aux, al,
                 velocity_space, pressure_space, Coll),
    vankaSmoother_(*this, b1->GetN_Columns(),level)
{
  int i;
  double *aux;

  A11 = a11;
  A12 = a12;
  A21 = a21;
  A22 = a22;
  StructureA = &(A11->GetStructure());
  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
#ifdef __3D__
  A13 = a13;
  A23 = a23;
  A31 = a31;
  A32 = a32;
  A33 = a33;
  A13Entries = A13->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();
#endif

  B1T = b1t;
  B2T = b2t;
  StructureBT = &(B1T->GetStructure());
  BTRowPtr = StructureBT->GetRowPtr();
  BTKCol = StructureBT->GetKCol();
  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  B1 = b1;
  B2 = b2;
  StructureB = &(B1->GetStructure());
  BRowPtr = StructureB->GetRowPtr();
  BKCol = StructureB->GetKCol();
  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

#ifdef __3D__
  B3T = b3t;
  B3TEntries = B3T->GetEntries();
  B3 = b3;
  B3Entries = B3->GetEntries();
#endif

  C = nullptr;
  StructureC = nullptr;
  CRowPtr = nullptr;
  CKCol = nullptr;
  CEntries = nullptr;

#ifdef __2D__
  USpace = A11->GetFESpace2D();
  PSpace = B1->GetTestSpace2D();
#endif
#ifdef __3D__
  USpace = A11->GetFESpace3D();
  PSpace = B1->GetTestSpace3D();
#endif

  N_Active = USpace->GetActiveBound();
  HangingNodeBound = USpace->GetHangingBound();
  N_Dirichlet = USpace->GetN_Dirichlet();

  N_UDOF = USpace->GetN_DegreesOfFreedom();
  N_PDOF = PSpace->GetN_DegreesOfFreedom();

  N_DOF = GEO_DIM*N_UDOF+N_PDOF;

  U1 = u1;
  U2 = u1 + N_UDOF;
  P  = u1 + GEO_DIM*N_UDOF;

  Rhs1 = f1;
  Rhs2 = f1 + N_UDOF;
  RhsP = f1 + GEO_DIM*N_UDOF;

#ifdef __3D__
  U3 = u1 + 2*N_UDOF;
  Rhs3 = f1 + 2*N_UDOF;
#endif

  N_Aux = n_aux;
  Aux = new double* [N_Aux];
  aux = new double[N_Aux*N_DOF];
  for(i=0;i<N_Aux;i++)
    Aux[i] = aux+i*N_DOF;

  Type = 4;

  alpha = al[0];
  downwind = dw;
}

/** constructor with matrix C */
#ifdef __2D__
  TNSE_MGLevel4::TNSE_MGLevel4(int level,
                               TSquareMatrix2D *a11, TSquareMatrix2D *a12,
                               TSquareMatrix2D *a21, TSquareMatrix2D *a22,
                               TMatrix2D *b1, TMatrix2D *b2,
                               TMatrix2D *b1t, TMatrix2D *b2t,
                               TMatrix2D *c,
                               double *f1, double *u1,
                               int n_aux, double *al, int velocity_space,
                               int pressure_space, TCollection *Coll,
                   int *dw)
#endif
#ifdef __3D__
  TNSE_MGLevel4::TNSE_MGLevel4(int level,
                               TSquareMatrix3D *a11, TSquareMatrix3D *a12,
                               TSquareMatrix3D *a13, TSquareMatrix3D *a21,
                               TSquareMatrix3D *a22, TSquareMatrix3D *a23,
                               TSquareMatrix3D *a31, TSquareMatrix3D *a32,
                               TSquareMatrix3D *a33,
                               TMatrix3D *b1, TMatrix3D *b2, TMatrix3D *b3,
                               TMatrix3D *b1t, TMatrix3D *b2t, TMatrix3D *b3t,
                               TMatrix3D *c,
                               double *f1, double *u1,
                               int n_aux, double *al, int velocity_space,
                               int pressure_space, TCollection *Coll,
                               int *dw)
#endif
  : TNSE_MGLevel(level, f1, u1, n_aux, al,
                 velocity_space, pressure_space, Coll),
    vankaSmoother_(*this, b1->GetN_Columns(),level) 
{
  int i;
  double *aux;

  A11 = a11;
  A12 = a12;
  A21 = a21;
  A22 = a22;
  StructureA = &(A11->GetStructure());
  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
#ifdef __3D__
  A13 = a13;
  A23 = a23;
  A31 = a31;
  A32 = a32;
  A33 = a33;
  A13Entries = A13->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();
#endif

  B1T = b1t;
  B2T = b2t;
  StructureBT = &(B1T->GetStructure());
  BTRowPtr = StructureBT->GetRowPtr();
  BTKCol = StructureBT->GetKCol();
  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  B1 = b1;
  B2 = b2;
  StructureB = &(B1->GetStructure());
  BRowPtr = StructureB->GetRowPtr();
  BKCol = StructureB->GetKCol();
  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

#ifdef __3D__
  B3T = b3t;
  B3TEntries = B3T->GetEntries();
  B3 = b3;
  B3Entries = B3->GetEntries();
#endif

  C = c;
  StructureC = &(C->GetStructure());
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

#ifdef __2D__
  USpace = A11->GetFESpace2D();
  PSpace = B1->GetTestSpace2D();
#endif
#ifdef __3D__
  USpace = A11->GetFESpace3D();
  PSpace = B1->GetTestSpace3D();
#endif

  N_Active = USpace->GetActiveBound();
  HangingNodeBound = USpace->GetHangingBound();
  N_Dirichlet = USpace->GetN_Dirichlet();

  N_UDOF = USpace->GetN_DegreesOfFreedom();
  N_PDOF = PSpace->GetN_DegreesOfFreedom();

  N_DOF = GEO_DIM*N_UDOF+N_PDOF;

  U1 = u1;
  U2 = u1 + N_UDOF;
  P  = u1 + GEO_DIM*N_UDOF;

  Rhs1 = f1;
  Rhs2 = f1 + N_UDOF;
  RhsP = f1 + GEO_DIM*N_UDOF;

#ifdef __3D__
  U3 = u1 + 2*N_UDOF;
  Rhs3 = f1 + 2*N_UDOF;
#endif

  N_Aux = n_aux;
  Aux = new double* [N_Aux];
  aux = new double[N_Aux*N_DOF];
  for(i=0;i<N_Aux;i++)
    Aux[i] = aux+i*N_DOF;

  Type = 4;

  alpha = al[0];
  downwind = dw;
}

/** destructor */
TNSE_MGLevel4::~TNSE_MGLevel4()
{
} // ~TNSE3DMGLevel4

/** calculate defect */
void TNSE_MGLevel4::Defect(double *u1, double *f1,  double *d1, double &res)
{

#ifdef __2D__
  // compute defect
  //if(C)
  //  CoupledDefect(A11,A12,A21,A22,B1,B2,B1T,B2T,C,u1,f1,d1);
  //else
    CoupledDefect(A11,A12,A21,A22,B1,B2,B1T,B2T,u1,f1,d1);
  // project defect into L20
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector2D(d1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif
#ifdef __3D__
  // compute defect
  // if(C)
  //   CoupledDefect(A11,A12,A13,A21,A22,A23,A31,A32,A33,
  //               B1,B2,B3,B1T,B2T,B3T,C,u1,f1,d1);
  // else
    CoupledDefect(A11,A12,A13,A21,A22,A23,A31,A32,A33,
                B1,B2,B3,B1T,B2T,B3T,u1,f1,d1);
  // project defect into L20
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector3D(d1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif

  memset(d1+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
  memset(d1+N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#ifdef __3D__
  memset(d1+2*N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#endif
  // compute residual
  res = sqrt(Ddot(N_DOF,d1,d1));
}

/** correct Dirichlet and hanging nodes */
void TNSE_MGLevel4::CorrectNodes(double *u1)
{
  int i,j,k, index;
  double s, t;

  double *u2 = u1+N_UDOF;
#ifdef __3D__
  double u;
  double *u3 = u2+N_UDOF;
#endif

  // set Dirichlet nodes
  memset(u1+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
  memset(u2+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memset(u3+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
#endif

  // set hanging nodes
  j = ARowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = 0;
    t = 0;
#ifdef __3D__
    u = 0;
#endif
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      if(index != i)
      {
        s -= A11Entries[j] * u1[index];
        t -= A22Entries[j] * u2[index];
#ifdef __3D__
        u -= A22Entries[j] * u3[index];
#endif
      }
    } // endfor j
    u1[i] = s;
    u2[i] = t;
  } // endfor i
}

/** cellwise Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel4::CellVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother,int N_Levels)
{
#ifdef __2D__
  const int RhsDim =  3*MaxN_BaseFunctions2D;
  TFE2D *UEle, *PEle;
  TSquareMatrix2D *sqmatrix[1];
#endif
#ifdef __3D__
  const int RhsDim =  4*MaxN_BaseFunctions3D;
  TFE3D *UEle, *PEle;
  TSquareMatrix3D *sqmatrix[1];
#endif
  int j,k,l,m, N_Cells, ii;
  int j1, j2, k1, k2;
  double value, value1, value2;
  double value11,value12,value21,value22;
  double *uold;
  TCollection *Coll;
  double System[RhsDim*RhsDim];
  double Rhs[RhsDim], sol[RhsDim];
  int *UGlobalNumbers, *UBeginIndex, *UDOFs, UDOF, N_U;
  int *PGlobalNumbers, *PBeginIndex, *PDOFs, PDOF, N_P;
  int N_LocalDOF;
  int begin, end, ActiveBound;
  double damp = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
  TBaseCell *Cell;
  double *u2, *p, *rhs2, *rhsp;
  TItMethod *itmethod = nullptr;
  int LargestDirectSolve = TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  MatVecProc *MatVect=MatVectFull;
  DefectProc *Defect=DefectFull;
  TSquareMatrix **matrix= (TSquareMatrix **)sqmatrix;
#ifdef __3D__
  int j3, k3;
  double value3;
  double value13, value23,value31,value32,value33;
  double *u3, *rhs3;
#endif
 
  TDatabase::ParamDB->INTERNAL_LOCAL_DOF = -1;
#ifdef __2D__
  sqmatrix[0] = (TSquareMatrix2D *)System;
#endif
#ifdef __3D__
  sqmatrix[0] = (TSquareMatrix3D *)System;
#endif

  if(VankaColl)
    Coll = VankaColl;
  else
    Coll = USpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  UGlobalNumbers = USpace->GetGlobalNumbers();
  UBeginIndex = USpace->GetBeginIndex();
  ActiveBound = USpace->GetActiveBound();

  PGlobalNumbers = PSpace->GetGlobalNumbers();
  PBeginIndex = PSpace->GetBeginIndex();

  // set pointers
  u2 = u1 + N_UDOF;
#ifdef __3D__
  u3 = u2 + N_UDOF;
#endif
  p  = u1 + GEO_DIM*N_UDOF;

  rhs2 = rhs1 + N_UDOF;
#ifdef __3D__
  rhs3 = rhs2 + N_UDOF;
#endif
  rhsp = rhs1 + GEO_DIM*N_UDOF;

  // set Dirichlet values
  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#endif

  SetHangingNodes(u1);

  // old values
  uold = aux;
  //double *pold  = uold+GEO_DIM*N_UDOF;

  // save current solution on 'old' vectors
  memcpy(uold, u1, N_DOF*SizeOfDouble);

  // start of cell loop
  for(ii=0;ii<N_Cells;ii++)
  {
    //ii = downwind[i];
    Cell = Coll->GetCell(ii);
//    OutPut(i << downwind[i] << endl);
#ifdef __2D__
    UEle = TFEDatabase2D::GetFE2D(USpace->GetFE2D(ii, Cell));
    PEle = TFEDatabase2D::GetFE2D(PSpace->GetFE2D(ii, Cell));
#endif
#ifdef __3D__
    UEle = TFEDatabase3D::GetFE3D(USpace->GetFE3D(ii, Cell));
    PEle = TFEDatabase3D::GetFE3D(PSpace->GetFE3D(ii, Cell));
#endif

    // get local number of dof
    N_U = UEle->GetN_DOF();
    N_P = PEle->GetN_DOF();
    N_LocalDOF = GEO_DIM*N_U+N_P;

    // reset local systems
    /*if (N_LocalDOF > RhsDim)
    {
      OutPut(
    "TNSE_MGLevel4::CellVanka - Not enough memory in array Rhs!!!"
    << endl << "available " << RhsDim << " needed " <<
    N_LocalDOF << endl);
      exit(4711);
      }*/
    // reset local systems
    memset(System, 0, SizeOfDouble*N_LocalDOF*N_LocalDOF);
    if (N_LocalDOF > LargestDirectSolve)
    {
      // size of local system has changed
      if (N_LocalDOF != TDatabase::ParamDB->INTERNAL_LOCAL_DOF)
      {
        // itmethod exists already
        if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0)
          delete itmethod;
        // allocate new itmethod
        itmethod = new TFgmresIte(MatVect, Defect, nullptr, 0, N_LocalDOF, 1);
        TDatabase::ParamDB->INTERNAL_LOCAL_DOF = N_LocalDOF;
      }
    }

    UDOFs = UGlobalNumbers+UBeginIndex[ii];
    PDOFs = PGlobalNumbers+PBeginIndex[ii];

    // fill local matrix
    for(j=0;j<N_U;j++)
    {
      j1 = j;
      j2 = j+N_U;
#ifdef __3D__
      j3 = j2+N_U;
#endif
      UDOF = UDOFs[j];

      // A block
      begin = ARowPtr[UDOF];
      end = ARowPtr[UDOF+1];

      Rhs[j1] = rhs1[UDOF];
      Rhs[j2] = rhs2[UDOF];
#ifdef __3D__
      Rhs[j3] = rhs3[UDOF];
#endif

      for(k=begin;k<end;k++)
      {
        l = AKCol[k];

        value11 = A11Entries[k];
        value12 = A12Entries[k];
        value21 = A21Entries[k];
        value22 = A22Entries[k];
#ifdef __3D__
        value13 = A13Entries[k];
        value23 = A23Entries[k];
        value31 = A31Entries[k];
        value32 = A32Entries[k];
        value33 = A33Entries[k];
#endif

#ifdef __2D__
        if (UDOF>=ActiveBound) // Dirichlet node
          value12 = value21 = 0;

        Rhs[j1] -= value11*u1[l]+value12*u2[l];
        Rhs[j2] -= value21*u1[l]+value22*u2[l];
#endif
#ifdef __3D__
        if (UDOF>=ActiveBound) // Dirichlet node
          value12 = value13 = value21 = value23 = value31 = value32 = 0;

        Rhs[j1] -= value11*u1[l]+value12*u2[l]+value13*u3[l];
        Rhs[j2] -= value21*u1[l]+value22*u2[l]+value23*u3[l];
        Rhs[j3] -= value31*u1[l]+value32*u2[l]+value33*u3[l];
#endif

        for(m=0;m<N_U;m++)
          if(UDOFs[m]==l)
          {
            // column belongs to local system
            k1 = m*N_LocalDOF;
            k2 = (m+N_U)*N_LocalDOF;
            System[k1+j1] = value11;
            System[k2+j1] = value12;
            System[k1+j2] = value21;
            System[k2+j2] = value22;
#ifdef __3D__
            k3 = (m+2*N_U)*N_LocalDOF;
            System[k3+j1] = value13;
            System[k3+j2] = value23;
            System[k1+j3] = value31;
            System[k2+j3] = value32;
            System[k3+j3] = value33;
#endif
            break;
          }
      } // endfor k

      if(UDOF<ActiveBound)  // active dof
      {
        // transpose(B) block for non-Dirichlet nodes
        begin = BTRowPtr[UDOF];
        end = BTRowPtr[UDOF+1];

        for(k=begin;k<end;k++)
        {
          l = BTKCol[k];
          value1 = B1TEntries[k];
          value2 = B2TEntries[k];
#ifdef __3D__
          value3 = B3TEntries[k];
#endif
          value = p[l];
          Rhs[j1] -= value1*value;
          Rhs[j2] -= value2*value;
#ifdef __3D__
          Rhs[j3] -= value3*value;
#endif

          for(m=0;m<N_P;m++)
            if(PDOFs[m]==l)
            {
              // column belongs to local system
              k1 = (m+GEO_DIM*N_U)*N_LocalDOF;
              System[k1+j1] = value1;
              System[k1+j2] = value2;
#ifdef __3D__
              System[k1+j3] = value3;
#endif
              break;
            }

        } // endfor k
      } // endif UDOF<ActiveBound
    } // endfor j

    for(j=0;j<N_P;j++)
    {
      j1 = j+GEO_DIM*N_U;
      PDOF = PDOFs[j];
      begin = BRowPtr[PDOF];
      end = BRowPtr[PDOF+1];
      Rhs[j1] = rhsp[PDOF];

      for(k=begin;k<end;k++)
      {
        l=BKCol[k];
        value1 = B1Entries[k];
        value2 = B2Entries[k];
#ifdef __3D__
        value3 = B3Entries[k];
#endif
        Rhs[j1] -= value1*u1[l];
        Rhs[j1] -= value2*u2[l];
#ifdef __3D__
        Rhs[j1] -= value3*u3[l];
#endif

        for(m=0;m<N_U;m++)
          if(UDOFs[m]==l)
          {
            // column belongs to local system
            k1 = m;
            k2 = m+N_U;
#ifdef __3D__
            k3 = k2 + N_U;
#endif
            System[k1*N_LocalDOF+j1] = value1;
            System[k2*N_LocalDOF+j1] = value2;
#ifdef __3D__
            System[k3*N_LocalDOF+j1] = value3;
#endif
            break;
          }
      } // endfor k
    } // endfor j

    if(C) // handle matrix C if present
    {
      for(j=0;j<N_P;j++)
      {
        j1 = j+GEO_DIM*N_U;
        PDOF = PDOFs[j];
        begin = CRowPtr[PDOF];
        end = CRowPtr[PDOF+1];
        for(k=begin;k<end;k++)
        {
          l = CKCol[k];
          value = -CEntries[k]; // minus is right sign
          Rhs[j1] -= value*p[l];

          for(m=0;m<N_P;m++)
            if(PDOFs[m] == l)
            {
              // column belongs to local system
              System[(m+GEO_DIM*N_U)*N_LocalDOF + j1] = value;
              break;
            } // endif
        } // endfor k
      } // endfor j
    } // endif C


    // solve local system
    if (smoother==1 && !C) // no diagonal Vanka for matrix C
    {
      // diagonal Vanka
#ifdef __2D__
      SolveDiagonalVanka2D(System, Rhs, N_U, N_P, N_LocalDOF);
#endif
#ifdef __3D__
      SolveDiagonalVanka3D(System, Rhs, N_U, N_P, N_LocalDOF);
#endif
    }
    else
    {
      // full Vanka
      if (N_LocalDOF > LargestDirectSolve)
      {
        memset(sol,0,N_LocalDOF*SizeOfDouble);
        itmethod->Iterate(matrix,nullptr,sol,Rhs);
        memcpy(Rhs, sol, N_LocalDOF*SizeOfDouble);
      }
      else
      {
        SolveLinearSystemLapack(System, Rhs, N_LocalDOF, N_LocalDOF);
      }
    }
#ifdef __3D__
    j1 = 2*N_U;
#endif
    for(j=0;j<N_U;j++)
    {
      l = UDOFs[j];
      u1[l] += damp*Rhs[j];
      u2[l] += damp*Rhs[j+N_U];
#ifdef __3D__
      u3[l] += damp*Rhs[j+j1];
#endif  
    }

    j1 = GEO_DIM*N_U;
    for(j=0;j<N_P;j++)
    {
      l = PDOFs[j];
      p[l] += damp*Rhs[j+j1];
    }
  } // endfor loop over cells

  // apply damping
  if (fabs(1-alpha)>1e-3)
    for(j=0;j<N_DOF;j++)
       u1[j] = uold[j]+alpha*(u1[j]-uold[j]);

  // set Dirichlet values
  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#endif

  SetHangingNodes(u1);

  if(TDatabase::ParamDB->INTERNAL_LOCAL_DOF > 0)
  {
    delete itmethod;
  }
/*  OutPut("test downwind started" << endl);
  for (i=0;i<N_Cells;i++)
  {
     if (test[i] != 1)
     {
        OutPut("test component " << i << " is "<< test[i] << endl);
     }
  }
  OutPut("test downwind done" << endl);
  delete test;
*/
} // end Vanka

/** nodal Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel4::NodalVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother, int N_Levels)
{
#ifdef __2D__
  const int MaxN_LocalU = 2*MaxN_BaseFunctions2D;
  const int SystemRhs = 3*MaxN_BaseFunctions2D;
  TSquareMatrix2D *sqmatrix[1];
#endif
#ifdef __3D__
  const int MaxN_LocalU = 4*MaxN_BaseFunctions3D;
  const int SystemRhs = 8*MaxN_BaseFunctions3D;
  TSquareMatrix3D *sqmatrix[1];
#endif
  int i,j,k,l,m;
  int j1, j2, j4, k1, k2;
  double value, value1, value2;
  double value11,value12,value21,value22;
  double *uold;
  double System[SystemRhs*SystemRhs];
  double Rhs[SystemRhs], sol[SystemRhs];
  int N_LocalDOF;
  int begin, end, HangingBound;
  int UDOFs[MaxN_LocalU], UDOF, N_U, N_UGEO;
  double *u2, *p, *rhs2, *rhsp;
  TItMethod *itmethod = nullptr;
  double damp = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
  int LargestDirectSolve = TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
  MatVecProc *MatVect=MatVectFull;
  DefectProc *Defect=DefectFull;
  TSquareMatrix **matrix= (TSquareMatrix **)sqmatrix;
#ifdef __3D__
  int j3, k3;
  double value3;
  double value13, value23,value31,value32,value33;
  double *u3, *rhs3;
#endif

  TDatabase::ParamDB->INTERNAL_LOCAL_DOF = -1;
  
#ifdef __2D__
  sqmatrix[0] = (TSquareMatrix2D *)System;
#endif
#ifdef __3D__
  sqmatrix[0] = (TSquareMatrix3D *)System;
#endif

  HangingBound = USpace->GetHangingBound();

  // set pointers
  u2 = u1 + N_UDOF;
#ifdef __3D__
  u3 = u2 + N_UDOF;
#endif
  p  = u1 + GEO_DIM*N_UDOF;

  rhs2 = rhs1 + N_UDOF;
#ifdef __3D__
  rhs3 = rhs2 + N_UDOF;
#endif
  rhsp = rhs1 + GEO_DIM*N_UDOF;

  // set Dirichlet values
  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#endif

  // old values
  uold = aux;
  //double *pold = uold+GEO_DIM*N_UDOF;

  // save current solution on 'old' vectors
  memcpy(uold, u1, N_DOF*SizeOfDouble);

  // start of loop over pressure dofs
  for(i=0;i<N_PDOF;i++)
  {
    N_U = 0;
    // go through row i of B1 and B2
    begin = BRowPtr[i];
    end = BRowPtr[i+1];
    value = rhsp[i];      // rhs of this pressure value
   for(k=begin;k<end;k++)
    { 
      l=BKCol[k];
      UDOFs[N_U] = l;
      value1 = B1Entries[k];
      value2 = B2Entries[k];
#ifdef __3D__
      value3 = B3Entries[k];
#endif
      j1 = GEO_DIM*N_U;
     System[j1] = value1;  // save values for local B
     System[j1+1] = value2;
#ifdef __2D__
      value -= value1*u1[l]+value2*u2[l]; // update rhs
#endif
#ifdef __3D__
      System[j1+2] = value3;
      value -= value1*u1[l]+value2*u2[l]+value3*u3[l]; // update rhs
#endif
      N_U++;           // count # velo dof connected to the pressure dof
    }                    // row done
    /*if (N_U>=MaxN_LocalU)
    {
      OutPut("TNSE_MGLevel4::NodalVanka - N_U too large !!! " << N_U
         << " " << MaxN_LocalU << __FILE__ << endl);
      exit(4711);
      }*/
#ifdef __3D__
    int N_U2 = 2 * N_U;
#endif
    N_UGEO = GEO_DIM * N_U;
    N_LocalDOF = N_UGEO +1;
    /*if (N_LocalDOF > SystemRhs)
    {
      OutPut(
         "TNSE_MGLevel4::NodalVanka - Not enough memory in array Rhs !!!"
         << endl << "available " << SystemRhs << " needed " <<
         N_LocalDOF << endl);
      exit(4711);
      }*/
    if (N_LocalDOF > LargestDirectSolve)
    {
      // size of local system has changed
      if (N_LocalDOF != TDatabase::ParamDB->INTERNAL_LOCAL_DOF)
      {
        // itmethod exists already
        if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0)
        {
          delete itmethod;
        }
        // allocate new itmethod
        itmethod = new TFgmresIte(MatVect, Defect, nullptr, 0, N_LocalDOF, 1);
        TDatabase::ParamDB->INTERNAL_LOCAL_DOF = N_LocalDOF;
      }
    }
    memset(System+N_UGEO, 0, SizeOfDouble*(N_LocalDOF*N_LocalDOF-N_UGEO));

    Rhs[N_LocalDOF-1] = value;  // set rhs

    for (k=0;k<N_U;k++)         // copy local B to the right place
    {
      j4 = GEO_DIM*k;
      System[k*N_LocalDOF+N_UGEO]=System[j4];
      System[(k+N_U)*N_LocalDOF+N_UGEO]=System[j4+1];
#ifdef __3D__
      System[(k+N_U2)*N_LocalDOF+N_UGEO]=System[j4+2];
#endif
    }
    // reset first part of array System
    memset(System, 0, SizeOfDouble*N_UGEO);

    // fill local matrix
    for(j=0;j<N_U;j++)
    {
      j1 = j;
      j2 = j+N_U;
#ifdef __3D__
      j3 = j2+N_U;
#endif
      UDOF = UDOFs[j];

      // A block
      begin = ARowPtr[UDOF];
      end = ARowPtr[UDOF+1];

      Rhs[j1] = rhs1[UDOF];
      Rhs[j2] = rhs2[UDOF];
#ifdef __3D__
      Rhs[j3] = rhs3[UDOF];
#endif

      for(k=begin;k<end;k++)
      {
        l = AKCol[k];
        value11 = A11Entries[k];
        value12 = A12Entries[k];
        value21 = A21Entries[k];
        value22 = A22Entries[k];
#ifdef __3D__
        value13 = A13Entries[k];
        value23 = A23Entries[k];
        value31 = A31Entries[k];
        value32 = A32Entries[k];
        value33 = A33Entries[k];
#endif

#ifdef __2D__
        if (UDOF>=HangingBound) // Dirichlet node
          value21 = value12 = 0;

        Rhs[j1] -= value11*u1[l]+value12*u2[l];
        Rhs[j2] -= value21*u1[l]+value22*u2[l];
#endif
#ifdef __3D__
        if (UDOF>=HangingBound) // Dirichlet node
          value12 = value13 = value21 = value23 = value31 = value32 = 0;

        Rhs[j1] -= value11*u1[l]+value12*u2[l]+value13*u3[l];
        Rhs[j2] -= value21*u1[l]+value22*u2[l]+value23*u3[l];
        Rhs[j3] -= value31*u1[l]+value32*u2[l]+value33*u3[l];
#endif

        for(m=0;m<N_U;m++)
          if(UDOFs[m]==l)
          {
            k1 = m*N_LocalDOF;
            k2 = (m+N_U)*N_LocalDOF;

            System[k1+j1] = value11;
            System[k2+j1] = value12;
            System[k1+j2] = value21;
            System[k2+j2] = value22;
#ifdef __3D__
            k3 = (m+2*N_U)*N_LocalDOF;
            System[k3+j1] = value13;
            System[k3+j2] = value23;
            System[k1+j3] = value31;
            System[k2+j3] = value32;
            System[k3+j3] = value33;
#endif
            break;
          }
      } // endfor k

      if(UDOF<HangingBound)  // active dof
      {
        // transpose(B) block for non-Dirichlet nodes
        begin = BTRowPtr[UDOF];
        end = BTRowPtr[UDOF+1];

        for(k=begin;k<end;k++)
        {
          l = BTKCol[k];
          value1 = B1TEntries[k];
          value2 = B2TEntries[k];
#ifdef __3D__
          value3 = B3TEntries[k];
#endif
          value = p[l];
          Rhs[j1] -= value1*value;
          Rhs[j2] -= value2*value;
#ifdef __3D__
          Rhs[j3] -= value3*value;
#endif

          if(i==l)
          {
        j4 = N_UGEO*N_LocalDOF;
            System[j4+j1] = value1;
            System[j4+j2] = value2;
#ifdef __3D__
            System[j4+j3] = value3;
#endif
          }
        } // endfor k
      } // endif UDOF<HangingBound
    } // endfor j

    if(C)
    {
      // fill C block if present
      begin = CRowPtr[i];
      end = CRowPtr[i+1];
      for(k=begin;k<end;k++)
      {
        l = CKCol[k];
        value = -CEntries[k]; // minus is right sign
        Rhs[N_LocalDOF-1] -= value*p[l];
        if(l==i) // main diagonal
          System[N_LocalDOF*N_LocalDOF-1] = value;
      } // endfor k
    } // endif C

    // solve local system
    if (smoother==3 && !C) // no diagonal Vanka for matrix C
    {
#ifdef __2D__
      // diagonal Vanka
      SolveDiagonalVanka2D(System,  Rhs, N_U, 1, N_LocalDOF);
#endif
#ifdef __3D__
      // diagonal Vanka
      SolveDiagonalVanka3D(System,  Rhs, N_U, 1, N_LocalDOF);
#endif
    }
    else
    {
       // full Vanka
      if (N_LocalDOF > LargestDirectSolve)
      {
        memset(sol,0,N_LocalDOF*SizeOfDouble);
        itmethod->Iterate(matrix,nullptr,sol,Rhs);
        memcpy(Rhs, sol, N_LocalDOF*SizeOfDouble);
      }
      else
      {
        //for (m=0;m<N_LocalDOF;m++)
        SolveLinearSystemLapack(System, Rhs, N_LocalDOF, N_LocalDOF);
      }
    }

    // update dof
    for(j=0;j<N_U;j++)
    {
      l = UDOFs[j];
      u1[l] += damp*Rhs[j];
      u2[l] += damp*Rhs[j+N_U];
#ifdef __3D__
      u3[l] += damp*Rhs[j+N_U2];
#endif  
    }
    p[i] += damp*Rhs[N_UGEO];
  } // endfor loop over pressure nodes

  // apply damping
  for(j=0;j<N_DOF;j++)
    u1[j] = uold[j]+alpha*(u1[j]-uold[j]);

  // set Dirichlet values
  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound,
         N_Dirichlet*SizeOfDouble);
#endif
  // itmethod exists
  if ( TDatabase::ParamDB->INTERNAL_LOCAL_DOF >0)
  {
    delete itmethod;
  }
} // end Vanka

/** solve exact on this level */
void TNSE_MGLevel4::SolveExact(double *u1, double *rhs1)
{
  double *a, *b;
  int i,j,k,l,index, N_DOF2 = N_DOF*N_DOF;
  int N_UDOFGEO = GEO_DIM *N_UDOF ;
  double value;
#ifdef __3D__
  int N_UDOF2 = 2*N_UDOF;
#endif

  a = new double[N_DOF2];
  b = new double[N_DOF];

  memset(a, 0, N_DOF2*SizeOfDouble);

  j = ARowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      index = AKCol[j];
      a[index*N_DOF+i] =  A11Entries[j];
      a[(index+N_UDOF)*N_DOF+i] = A12Entries[j];
      a[index*N_DOF+i+N_UDOF] = A21Entries[j];
      a[(index+N_UDOF)*N_DOF+i+N_UDOF] = A22Entries[j];
#ifdef __3D__
      a[(index+N_UDOF2)*N_DOF+i] = A13Entries[j];
      a[(index+N_UDOF2)*N_DOF+i+N_UDOF] = A23Entries[j];
      a[index*N_DOF+i+N_UDOF2] = A31Entries[j];
      a[(index+N_UDOF)*N_DOF+i+N_UDOF2] = A32Entries[j];
      a[(index+N_UDOF2)*N_DOF+i+N_UDOF2] = A33Entries[j];
#endif
    }
  } // endfor i

  // Dirichlet and hanging nodes
  for(i=N_Active;i<N_UDOF;i++)
  {
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      a[i+N_DOF*index] = A11Entries[j];
      a[(index+N_UDOF)*N_DOF+i+N_UDOF] = A22Entries[j];
#ifdef __3D__
      a[(index+N_UDOF2)*N_DOF+i+N_UDOF2] = A33Entries[j];
#endif
    }
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      l = (index+GEO_DIM*N_UDOF)*N_DOF;
      a[l+i] = B1TEntries[j];
      a[l+i+N_UDOF] = B2TEntries[j];
#ifdef __3D__
      a[l+i+2*N_UDOF] = B3TEntries[j];
#endif
    } // endfor j
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      a[(index)*N_DOF+i+N_UDOFGEO] = B1Entries[j];
      a[(index+N_UDOF)*N_DOF+i+N_UDOFGEO] = B2Entries[j];
#ifdef __3D__
      a[(index+N_UDOF2)*N_DOF+i+N_UDOFGEO] = B3Entries[j];
#endif
    } // endfor j
  } // endfor i

  if(C)
  {
    j = CRowPtr[0];
    for(i=0;i<N_PDOF;i++)
    {
      k = CRowPtr[i+1];
      for(;j<k;j++)
      {
        index = CKCol[j];
        value = -CEntries[j]; // minus is right sign
        a[(index+GEO_DIM*N_UDOF)*N_DOF + i+GEO_DIM*N_UDOF] = value;
      }
    } // endfor i
  } // endif C

  // condition for pressure, fix first value
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    for (i=0;i<N_DOF;i++)
      a[i*N_DOF+N_UDOFGEO] = 0;
    a[(N_UDOFGEO)*N_DOF+N_UDOFGEO] = 1;
  }

  // copy into local data
  memcpy(b, rhs1, N_DOF*SizeOfDouble);
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    b[N_UDOFGEO] = 0;

  /* for (i=0;i<N_DOF;i++)
      for (j=0;j<N_DOF;j++)
      OutPut("a("<<i+1<<","<<j+1<<") = " << a[i*N_DOF+j] << endl);*/
 SolveLinearSystemLapack(a, b, N_DOF, N_DOF);
	  

   
  // copy from local data
  memcpy(u1, b, N_DOF*SizeOfDouble);

  delete a;
  delete b;
}
/** solve exact on this level */
void TNSE_MGLevel4::SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag)
{
 
#ifdef __3D__
    if (umfpack_flag==-1)
    {
	DirectSolver_old(A11, A12, A13, A21, A22, A23, A31, A32, A33,
		     B1T, B2T, B3T, B1, B2, B3, rhs1, u1, 4);
	umfpack_flag = 0;
    }
    DirectSolver_old(A11, A12, A13, A21, A22, A23, A31, A32, A33,
		 B1T, B2T, B3T, B1, B2, B3, rhs1, u1, umfpack_flag);
    umfpack_flag = 1;
#endif
#ifdef __2D__

    // free memory from previous time
    if (umfpack_flag==-1)
    {
	DirectSolver_old(A11, A12, A21, A22, B1T, B2T, B1, B2, rhs1, u1, 4);
	umfpack_flag = 0;
    }
	
    //OutPut("TNSE_MGLevel4::SolveExactUMFPACK: Are we here?" << endl);
    DirectSolver_old(A11, A12, A21, A22, B1T, B2T, B1, B2, rhs1, u1, 
                     umfpack_flag);
    umfpack_flag = 1;
#endif
}

/** step length control for smoother */
double TNSE_MGLevel4::StepLengthControl (double *u1, double *u1old,
                                         double *def1,
                                         int N_Parameters,
                                         double *Parameters)
{
  double *x,*y,omega,numerator,nominator;
  int i;

  // allocate auxiliary array
  x = new double[2*N_DOF];
  y = x+N_DOF;

  // save current update in array x
  for (i=0;i<N_DOF;i++)
    x[i] = u1[i]-u1old[i];
  memset(y,0,N_DOF*SizeOfDouble);

  // compute matrix times update
#ifdef __2D__
//  if(C)
//    CoupledMatVect(A11,A12,A21,A22,B1,B2,B1T,B2T,C,x,y);
//  else
    CoupledMatVect(A11,A12,A21,A22,B1,B2,B1T,B2T,x,y);
#endif
#ifdef __3D__
  // if(C)
  //   CoupledMatVect(A11,A12,A13,A21,A22,A23,A31,A32,A33,
  //                  B1,B2,B3,B1T,B2T,B3T,C,x,y);
  // else
  CoupledMatVect(A11,A12,A13,A21,A22,A23,A31,A32,A33,
                 B1,B2,B3,B1T,B2T,B3T,x,y);
#endif

  numerator = Ddot(N_DOF,def1,y);
  nominator = Ddot(N_DOF,y,y);

  if (nominator > 0)
    omega = numerator/nominator;
  else
    {
      if(N_Parameters>0)
        omega = alpha;
      else
        omega = 0.5;
      OutPut("MESSAGE : Step length control failed. Set omega = " << omega<< endl);
    }
  if (fabs(omega)<0.0001)
    {
       if(N_Parameters>0)
        omega = Parameters[0];
      else
        omega = 0.9;
    }
  delete x;
  Output::print<2>("step length control ", omega);
  return(omega);
}


/** Braess--Sarazin smoother  */
void TNSE_MGLevel4::BraessSarazin(double *u1, double *rhs1,
                                  double *aux, int N_Parameters,
                                  double *Parameters,int N_Levels)
{
  double *sol;
  int j;

  // array for update
  sol = aux;

  // initialize array for update
  memset(sol, 0, (N_DOF)*SizeOfDouble);

  // set pointer to rhs and solution
  j = 0;

  // call the algebraic solver
  // the last input parameter is only a dummy
#ifdef __2D__
  if(C)
  {
    OutPut("Braess-Sarazin smoother not implemented for C block !!!" << endl);
    exit(-1);
  }
  else
    OldSolver(A11, A12, A21, A22, B1, B2, B1T, B2T, rhs1, sol, j);
#endif
#ifdef __3D__
  OutPut("Braess-Sarazin smoother not implemented !!!" << endl);
  exit(4711);
#endif

  // update and apply damping
  for(j=0;j<N_DOF;j++)
    u1[j] += alpha*sol[j];

}
/** print all matrices and oth right hand sides */
void TNSE_MGLevel4::PrintAll()
{
}

/*********************
 * Things belonging to new Vanka implementation.
 *********************/

//! @brief Initialize the smoother_ object.
void TNSE_MGLevel4::initializeSmoother(){
  // Call the vanka object to initialize itself with this level's information.
  vankaSmoother_.initialize();
}

/*! @brief Perform one step of Vanka smoothing on the level. Call only after initializer finished.
 *  @param[in,out] currentSolution The current global solution to be updated by the smoother.
 *  @param[in] currentRHS The global right hand side.
 *  @param[in,out] storeOldSolution Auxilliary, memory array to be used internally for the damping.
 */
void TNSE_MGLevel4::applySmoother(double *currentSolution, const double* const currentRHS,
    double *storeOldSolution){

  // save current solution in 'old' vectors
  memcpy(storeOldSolution, currentSolution, N_DOF*SizeOfDouble);

  //Puts up and solve all local systems and update currentSolution.
  vankaSmoother_.solveLocalSystems(currentSolution, currentRHS);

  // apply damping
  for(int j=0;j<N_DOF;j++)
    currentSolution[j] = storeOldSolution[j]+alpha*(currentSolution[j]-storeOldSolution[j]);
  // set Dirichlet values
  memcpy(currentSolution+HangingNodeBound, currentRHS+HangingNodeBound,
     N_Dirichlet*SizeOfDouble);
  memcpy(currentSolution + N_UDOF + HangingNodeBound, currentRHS+N_UDOF+HangingNodeBound,
     N_Dirichlet*SizeOfDouble);
}

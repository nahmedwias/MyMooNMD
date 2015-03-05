// =======================================================================
// LocalProjection.C
//
// Purpose:   routines for local projection stabilization
//
// Author:    Gunar Matthies  2007/03/06
//
// =======================================================================

#include <Database.h>
#include <MooNMD_Io.h>

#ifdef __2D__
  #include <FEDatabase2D.h>
  #include <FEFunction2D.h>
  #include <NodalFunctional2D.h>
#else  
  #include <FEDatabase3D.h>
  #include <FEFunction3D.h>
  #include <NodalFunctional3D.h> 
#endif

#include <LinAlg.h>

#include <string.h>
#include <stdlib.h>

#ifdef __2D__
FE2D GetElement2D(TBaseCell *cell, int CoarseOrder)
{
  FE2D ele = (FE2D)0;
  Shapes shapetype;

  shapetype = cell->GetType();
  switch(shapetype)
  {
    // regularly refined quadrilateral
    case Quadrangle:
    case Parallelogram:
    case Rectangle:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_Q0_2D_Q_M;
        break;

        case 1:
          ele = D_P1_2D_Q_M;
        break;

        case 2:
          ele = D_P2_2D_Q_M;
        break;

        case 3:
          ele = D_P3_2D_Q_M;
        break;

        case 4:
          ele = D_P4_2D_Q_M;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_Q00_2D_Q_M;
          }
          else
          {
            OutPut("CoarseOrder: " << CoarseOrder << endl);
            OutPut("Projection space is defined up to order 4" << endl);
            exit(-1);
          }
      } // end switch CoarseOrder
    break; // end regularly refined quadrilateral

    case Triangle:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_P0_2D_T_A;
        break;

        case 1:
          ele = D_P1_2D_T_A;
        break;

        case 2:
          ele = D_P2_2D_T_A;
        break;

        case 3:
          ele = D_P3_2D_T_A;
        break;

        case 4:
          ele = D_P4_2D_T_A;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_P00_2D_T_A;
          }
          else
          {
            OutPut("CoarseOrder: " << CoarseOrder << endl);
            OutPut("Projection space is defined up to order 4" << endl);
            exit(-1);
          }
      } // end switch CoarseOrder
    break;
    default:
      OutPut("Invalid shape" << endl);
      exit(-1);
  } // end switch reftype
  return ele;
}

/** Navier--Stokes type 2 (NSTYPE==2) with C*/
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  int N_Active;
  double *CEntries;
  int *CRowPtr, *CKCol;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;

  N_Active = A->GetActiveBound();

  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u1[index];
      t -= value * u2[index];
    }
    r1[i] = s;
    r2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s -= value1 * u1[index] + value2 * u2[index];
    } // endfor j
    r3[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    }
    r3[i] -= s;
  } // endfor i
}

void Defect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], B[0], B[1], B[2], B[3], B[4], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF, TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  int N_Active;
  double *CEntries;
  int *CRowPtr, *CKCol;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];

  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
      t += value * u2[index];
    }
    v1[i] = s;
    v2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s += value1 * u1[index] + value2 * u2[index];
    } // endfor j
    q[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    v1[i] += s;
    v2[i] += t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    } // endfor j
    q[i] += s;
  } // endfor i

  return;
}

void MatVect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], B[2], B[3], B[4], x, y);
  return;
}

// stabilisation of full gradient (velocity or pressure)
void UltraLocalProjection(void* A, boolean ForPressure)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  int OrderDiff;
  double lpcoeff, lpexponent;

  if(!(TDatabase::ParamDB->LP_FULL_GRADIENT) && !(ForPressure))
  {
    OutPut("Local projection stabilization is implemented only for full gradient!" << endl);
    exit(-1);
  }

  if(ForPressure)
  {
    lpcoeff = -(TDatabase::ParamDB->LP_PRESSURE_COEFF);
    lpexponent = TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
  }
  else
  {
    lpcoeff = TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
    lpexponent = TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  }

  if(ForPressure)
  {
    fespace = (TFESpace2D*)(((TMatrix2D *)A)->GetStructure()->GetTestSpace());
    ActiveBound = -1;
    RowPtr = ((TMatrix2D *)A)->GetRowPtr();
    KCol = ((TMatrix2D *)A)->GetKCol();
    Entries = ((TMatrix2D *)A)->GetEntries();
    // cout << "for pressure" << endl;
  }
  else
  {
    fespace = ((TSquareMatrix2D *)A)->GetFESpace();
    ActiveBound = fespace->GetActiveBound();
    RowPtr = ((TSquareMatrix2D *)A)->GetRowPtr();
    KCol = ((TSquareMatrix2D *)A)->GetKCol();
    Entries = ((TSquareMatrix2D *)A)->GetEntries();
//     cout << "not for pressure" << endl;
  }

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*( ChildValueX[k]*ChildValueX[l]
                                     +ChildValueY[k]*ChildValueY[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ]
                                           +H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if((dof<ActiveBound) || (ForPressure))
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // 

double UltraLocalError(TFEFunction2D *uh, DoubleFunct2D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values;
  double error, locerror;
  double exactval[4];

  error = 0.0;

  fespace = uh->GetFESpace2D();
  Values = uh->GetValues();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    // only two right-hand sides (x and y derivative)
    memset(H, 0, N_CoarseDOF*2*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];

      // calculate gradient of discrete uh in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        val = Values[DOF[k]];
        valx += ChildValueX[k]*val;
        valy += ChildValueY[k]*val;
      }

      // get gradient of exact u
      ExactU(X[j], Y[j], exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k*2  ] += val*valx;
        H[k*2+1] += val*valy;
      } // end for k

      // grad-grad term
      locerror += w*(valx*valx + valy*valy);

    } // end for j
    memcpy(P, H, N_CoarseDOF*2*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2, 2);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2  ]*H[i2*2  ]
                                       +H[i1*2+1]*H[i2*2+1]);
    locerror += s;

    // grad-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2*2  ] * H[i2*2  ];
      s += P[i2*2+1] * H[i2*2+1];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1*2  ] * H[i1*2  ];
      s += P[i1*2+1] * H[i1*2+1];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalError

void AddStreamlineTerm(TSquareMatrix2D* A, TFEFunction2D *uh1,
                       TFEFunction2D *uh2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  double *Values1, *Values2;
  double BValue[MaxN_BaseFunctions2D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      valx = 0.0;
      valy = 0.0;
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
      }
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm

void AddStreamlineTermPWConst(TSquareMatrix2D* A, TFEFunction2D *uh1,
                              TFEFunction2D *uh2,
                              double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  double *Values1, *Values2;
  double BValue[MaxN_BaseFunctions2D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    // calculate pw constant approximation of velocity field (uh1, uh2)
    val  = 0.0;
    valx = 0.0;
    valy = 0.0;
    for(j=0;j<N_Points;j++)
    {
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        val  += w*ChildValue[k]*1;
        valx += w*ChildValue[k]*Values1[l];
        valy += w*ChildValue[k]*Values2[l];
      }
    }
    valx /= val;
    valy /= val;

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      /*
      valx = 0.0;
      valy = 0.0;
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
      }
      */
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTermPWConst

void AddDivergenceTerm(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixA11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s22;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;

  fespace = A11->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();
  EntriesA11 = A11->GetEntries();
  EntriesA12 = A12->GetEntries();
  EntriesA21 = A21->GetEntries();
  EntriesA22 = A22->GetEntries();
  // cout << "" << endl;


  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    memset(LocMatrixA11, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA12, 0, N_DOF*N_DOF*SizeOfDouble);
    memset(LocMatrixA22, 0, N_DOF*N_DOF*SizeOfDouble);
    // since A21=transpose(A12) A21 is not explicitly needed

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixA11[k*N_DOF+l] += w*ChildValueX[k]*ChildValueX[l];
          LocMatrixA12[k*N_DOF+l] += w*ChildValueX[k]*ChildValueY[l];
          LocMatrixA22[k*N_DOF+l] += w*ChildValueY[k]*ChildValueY[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ];
            s12 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m+N_DOF];
            s22 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF];
          }
        LocMatrixA11[l*N_DOF+m] += s11;
        LocMatrixA12[l*N_DOF+m] += s12;
        LocMatrixA22[l*N_DOF+m] += s22;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s12 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m+N_DOF];
          s22 += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s12 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l      ];
          s22 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrixA11[l*N_DOF+m] -= s11;
        LocMatrixA12[l*N_DOF+m] -= s12;
        LocMatrixA22[l*N_DOF+m] -= s22;
//         LocMatrixA21[m*N_DOF+l] = LocMatrixA12[l*N_DOF+m];
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesA11[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA11[l*N_DOF+m];
              EntriesA12[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA12[l*N_DOF+m];
              // since A21=transpose(A12)
              EntriesA21[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA12[m*N_DOF+l];
              EntriesA22[p] += lpcoeff*pow(hK,lpexponent)*LocMatrixA22[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddDivergenceTerm

/** Navier--Stokes type 4 (NSTYPE==4) with C*/
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    TMatrix *C,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2,value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries;
  double *CEntries;
  int *CRowPtr, *CKCol;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  N_UDOF = A11->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  N_Active = A11->GetActiveBound();
  j = ARowPtr[0];
  // real dof
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A21Entries[j];
      value3 = A22Entries[j];
      s += value * u1[index] + value1 * u2[index];
      t += value2* u1[index] + value3 * u2[index];
    }
    v1[i] = s;
    v2[i] = t;
  } // endfor i
  // Dirichlet and hanging nodes
  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      s += value * u1[index];
      t += value1 * u2[index];
    }
    v1[i] = s;
    v2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s += value1 * u1[index] + value2 * u2[index];
    } // endfor j
    q[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    v1[i] += s;
    v2[i] += t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    } // endfor j
    q[i] += s;
  } // endfor i

  return;
}

void MatVect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], B[4], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   TMatrix *C,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries;
  int N_Active;
  double *CEntries;
  int *CRowPtr, *CKCol;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  N_UDOF = A11->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;

  N_Active = A11->GetActiveBound();

  j = ARowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = v1[i];
    t = v2[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A21Entries[j];
      value3 = A22Entries[j];
      s -= value * u1[index] + value1 * u2[index];
      t -= value2* u1[index] + value3 * u2[index];
    }
    r1[i] = s;
    r2[i] = t;
  } // endfor i

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      s -= value * u1[index];
      t -= value1 * u2[index];
    }
    r1[i] = s;
    r2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s -= value1 * u1[index] + value2 * u2[index];
    } // endfor j
    r3[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    }
    r3[i] -= s;
  } // endfor i
}

void Defect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], B[4], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

double UltraLocalErrorDivergence(TFEFunction2D *uh1, TFEFunction2D *uh2,
                       DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values1, *Values2;
  double error, locerror;
  double exactval1[4], exactval2[4];
  double div;

  error = 0.0;

  fespace = uh1->GetFESpace2D();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    // only one right-hand side (div (u-uh))
    memset(H, 0, N_CoarseDOF*1*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];

      // calculate x derivative of uh1 and y derivative of uh2 in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValueX[k]*Values1[l];
        valy += ChildValueY[k]*Values2[l];
      }

      // get x derivative of u1 and y derivative of u2
      ExactU1(X[j], Y[j], exactval1);
      ExactU2(X[j], Y[j], exactval2);

      valx -= exactval1[1];
      valy -= exactval2[2];

      div = valx+valy;

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k  ] += val*div;
      } // end for k

      // id-id term
      locerror += w*div*div;

    } // end for j
    memcpy(P, H, N_CoarseDOF*1*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 1, 1);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*H[i1  ]*H[i2  ];
    locerror += s;

    // id-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2  ] * H[i2  ];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1  ] * H[i1  ];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalError

double UltraLocalErrorStreamline(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values;
  double *BValues1, *BValues2;
  double error, locerror;
  double exactval[4];
  double valb1, valb2, StreamlineDerivative;

  error = 0.0;

  fespace = uh->GetFESpace2D();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    Values = uh->GetValues();

    BValues1 = b1->GetValues();
    BValues2 = b2->GetValues();
    
    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    // only one right-hand side (b.grad(ui))
    memset(H, 0, N_CoarseDOF*1*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];

      // calculate derivative of uh and components of b in this quadrature point
      valx = 0.0;
      valy = 0.0;
      valb1 = 0.0;
      valb2 = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valb1 += ChildValue[k]*BValues1[l];
        valb2 += ChildValue[k]*BValues2[l];
        valx += ChildValueX[k]*Values[l];
        valy += ChildValueY[k]*Values[l];
      }

      // get gradient of exact u
      ExactU(X[j], Y[j], exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      StreamlineDerivative = valb1*valx + valb2*valy;

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k  ] += val*StreamlineDerivative;
      } // end for k

      // id-id term
      locerror += w*StreamlineDerivative*StreamlineDerivative;

    } // end for j
    memcpy(P, H, N_CoarseDOF*1*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 1, 1);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*H[i1  ]*H[i2  ];
    locerror += s;

    // id-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2  ] * H[i2  ];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1  ] * H[i1  ];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalError

double UltraLocalErrorStreamlinePWConst(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values;
  double *BValues1, *BValues2;
  double error, locerror;
  double exactval[4];
  double valb1, valb2, StreamlineDerivative;

  error = 0.0;

  fespace = uh->GetFESpace2D();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    DOF = GlobalNumbers + BeginIndex[i];

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);

    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    Values = uh->GetValues();

    BValues1 = b1->GetValues();
    BValues2 = b2->GetValues();
    
    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);
    ChildValues  = TFEDatabase2D::GetOrigElementValues(BF_ID, D00);

    // only one right-hand side (b.grad(ui))
    memset(H, 0, N_CoarseDOF*1*SizeOfDouble);

    // calculate pw constant approximation of velocity field (uh1, uh2)
    val  = 0.0;
    valx = 0.0;
    valy = 0.0;
    for(j=0;j<N_Points;j++)
    {
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        val  += w*ChildValue[k]*1;
        valx += w*ChildValue[k]*BValues1[l];
        valy += w*ChildValue[k]*BValues2[l];
      }
    }
    valb1 = valx/val;
    valb2 = valy/val;

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];

      // calculate derivative of uh and components of b in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValueX[k]*Values[l];
        valy += ChildValueY[k]*Values[l];
      }

      // get gradient of exact u
      ExactU(X[j], Y[j], exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      StreamlineDerivative = valb1*valx + valb2*valy;

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k  ] += val*StreamlineDerivative;
      } // end for k

      // id-id term
      locerror += w*StreamlineDerivative*StreamlineDerivative;

    } // end for j
    memcpy(P, H, N_CoarseDOF*1*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 1, 1);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*H[i1  ]*H[i2  ];
    locerror += s;

    // id-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2  ] * H[i2  ];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1  ] * H[i1  ];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalErrorPWConst


// Adaptive post processing basd on Friedhelm Schiweck talk at MAFELAP 09 - Sashi
void AdaptivePostProcess(TFEFunction2D *FeFunction, double *PostSol, bool DirichletBC)
{
  int i, j, k, l, CoarseOrder, N_Points, N_U;
  int N_Cells, N_CoarseDOF, N_DOF;
  int N_UsedElements;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int OrderDiff;

  double **CoarseValues, *CoarseValue, *LPS_sol;
  double val, w, *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double rhs[MaxN_BaseFunctions2D];
  double Values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double PCValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double *Value, *CurrentValue, sol;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double *W, maxbubble=-1E8, minbubble=1E8;


  bool SecondDer[1] = { FALSE };

  TFESpace2D *fespace;
  TCollection *coll;
  TFE2D *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  FE2D CurrEleID, UsedElements[1];
  TNodalFunctional2D *nf;

  fespace = FeFunction->GetFESpace2D();
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();
  LPS_sol = FeFunction->GetValues();
  N_U = FeFunction->GetLength();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  W = new double[N_U];
  memset(W, 0, N_U*SizeOfDouble);


  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);
    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;
    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    N_UsedElements = 1;
    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();
    N_CoarseDOF = CoarseBF->GetDimension();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements, 
                           coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     BF->GetDerivatives(D00, xi[j], eta[j], Values[j]);


    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    memset(rhs, 0, N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
     {
      CoarseValue = CoarseValues[j];
      Value = Values[j];
      w = AbsDetjk[j]*weights[j];

//    find the lps solution at this quadrature point
      sol=0.;
      for(k=0;k<N_DOF;k++)
       sol+= LPS_sol[DOF[k]]*Value[k];

      for(k=0;k<N_CoarseDOF;k++)
       {
        val = w*CoarseValue[k];
        rhs[k] += sol*val;

        for(l=0;l<N_CoarseDOF;l++)
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];

       } // end for k
     } // for(j=0;j<N_Points;j++)

//     for(j=0;j<N_CoarseDOF;j++)
//       for(k=0;k<N_CoarseDOF;k++)
//         cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
//     cout  << endl;
//     for(j=0;j<N_CoarseDOF;j++)
//         cout << j << "  " << rhs[j] << endl;

    SolveLinearSystemLapack(G, rhs, N_CoarseDOF, N_CoarseDOF);

// interpolate the discont solution to the original LPS space

    nf = CurrentElement->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);

    for(j=0;j<N_Points;j++)
     CoarseBF->GetDerivatives(D00, xi[j], eta[j], PCValues[j]);

    memset(PointValues, 0, N_Points*SizeOfDouble);

    for(j=0;j<N_Points;j++)
     for(k=0;k<N_CoarseDOF;k++)
      PointValues[j] +=  rhs[k]*PCValues[j][k];

    nf->GetAllFunctionals(coll, cell, PointValues,
                          FunctionalValues);

    for(j=0;j<N_DOF;j++)
     {
      PostSol[DOF[j]] += FunctionalValues[j];
      W[DOF[j]] += 1.;
      if(j==N_DOF-1)
        {
//          cout<< j<< " PostSol " << PostSol[DOF[j]] << endl;
          if (maxbubble< PostSol[DOF[j]]) maxbubble = PostSol[DOF[j]];
          if (minbubble > PostSol[DOF[j]]) minbubble = PostSol[DOF[j]];
        }
     }
   } // for(i=0;i<N_Cells;i++)

  for(i=0;i<N_U;i++)
   PostSol[i] /=W[i];

//   cout<<" maxbubble " << maxbubble  <<" minbubble " << minbubble  << endl;
//     exit(0);
}


//Fefunction can be a different FEspace - Sashi
void AddALEStreamlineLPS(TSquareMatrix2D* A, int N_FeFunct, TFEFunction2D **FeFunct,
                         double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int *U_GlobalNumbers, *U_BeginIndex, *U_DOF;  
  int *W_GlobalNumbers, *W_BeginIndex, *W_DOF;  
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF, N_UDOF, N_WDOF;
  TCollection *Coll;
  TFESpace2D *fespace, *U_fespace, *W_fespace;
  FE2D CurrEleID, U_CurrEleID, W_CurrEleID, UsedElements[2];
  int N_UsedElements;
  TFE2D *CurrentElement, *CoarseElement, *U_CurrentElement, *W_CurrentElement;
  TBaseFunct2D *BF, *CoarseBF, *U_BF, *W_BF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  BaseFunct2D  U_BF_ID, W_BF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  TFEFunction2D *uh1, *uh2, *wh1, *wh2;  
  
  bool SecondDer[2] = { FALSE, FALSE };
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
//   double **ChildValues, *ChildValue;
  double **U_Values, *U_Value; 
  double **W_Values, *W_Value;
  
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  double *Values1, *Values2, *Values3, *Values4;
  double BValue[MaxN_BaseFunctions2D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  if(N_FeFunct==2)
   {
    uh1 = FeFunct[0];
    uh2 = FeFunct[1];      
   }
  else if(N_FeFunct==4)
   {
    uh1 = FeFunct[0];
    uh2 = FeFunct[1];   
    wh1 = FeFunct[2];
    wh2 = FeFunct[3]; 
   }
  else
   {
    cout << "N_FeFunct must be 2 or 4 " <<endl;
    exit(0);    
   }
   
  U_fespace = uh1->GetFESpace2D();
  U_BeginIndex = U_fespace->GetBeginIndex();
  U_GlobalNumbers = U_fespace->GetGlobalNumbers(); 
  Values1 = uh1->GetValues();
  Values2 = uh2->GetValues();
    
  if(N_FeFunct==4) 
   {
    W_fespace = wh1->GetFESpace2D();  
    W_BeginIndex = W_fespace->GetBeginIndex();
    W_GlobalNumbers = W_fespace->GetGlobalNumbers(); 
    Values3 = wh1->GetValues();
    Values4 = wh2->GetValues();       
   }
  
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE2D(i, cell);
    CurrentElement = TFEDatabase2D::GetFE2D(CurrEleID);
    BF = CurrentElement->GetBaseFunct2D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    U_CurrEleID = U_fespace->GetFE2D(i, cell);
    U_CurrentElement = TFEDatabase2D::GetFE2D(U_CurrEleID);
    U_BF = U_CurrentElement->GetBaseFunct2D();
    U_BF_ID = U_BF->GetID();
    N_UDOF = U_BF->GetDimension();    
    
    W_CurrEleID = W_fespace->GetFE2D(i, cell);
    W_CurrentElement = TFEDatabase2D::GetFE2D(W_CurrEleID);
    W_BF = W_CurrentElement->GetBaseFunct2D();
    W_BF_ID = W_BF->GetID();
    N_WDOF = W_BF->GetDimension();   
    
    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase2D::GetFE2D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct2D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];    
    U_DOF = U_GlobalNumbers + U_BeginIndex[i];  
    W_DOF = W_GlobalNumbers + W_BeginIndex[i];  

    PCValues = TFEDatabase2D::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = TFEDatabase2D::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = TFEDatabase2D::GetOrigElementValues(BF_ID, D01);

    U_Values  = TFEDatabase2D::GetOrigElementValues(U_BF_ID, D00);
    W_Values  = TFEDatabase2D::GetOrigElementValues(W_BF_ID, D00);    

    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);
    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      U_Value  = U_Values[j];
      W_Value  = W_Values[j];      
      
      w = AbsDetjk[j]*weights[j];
      valx = 0.0;
      valy = 0.0;     
      // compute components of uh in j
      for(k=0;k<N_UDOF;k++)
      {
        l = U_DOF[k];
        valx += U_Value[k]*Values1[l];
        valy += U_Value[k]*Values2[l];
      }
      // sub mesh velo (uh-wh)
      for(k=0;k<N_WDOF;k++)
      {
        l = W_DOF[k];
        valx -= W_Value[k]*Values3[l];
        valy -= W_Value[k]*Values4[l];
      }

      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm
#else


FE3D GetElement3D(TBaseCell *cell, int CoarseOrder)
{
  FE3D ele = (FE3D)0;
  Shapes shapetype;

  shapetype = cell->GetType();
  switch(shapetype)
  {
    // regularly refined Hexahedron
    case Hexahedron:
    case Brick:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_Q0_3D_H_M;
        break;

        case 1:
          ele = D_P1_3D_H_M;
        break;

        case 2:
          ele = D_P2_3D_H_M;
        break;

        case 3:
          ele = D_P3_3D_H_M;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_Q00_3D_H_M;
          }
          else
          {
            OutPut("CoarseOrder: " << CoarseOrder << endl);
            OutPut("Projection space is defined up to order 3" << endl);
            exit(-1);
          }
      } // end switch CoarseOrder
    break; // end regularly refined quadrilateral

    case Tetrahedron:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_P0_3D_T_A;
        break;

        case 1:
          ele = D_P1_3D_T_A;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_P00_3D_T_A;
          }
          else
          {
            OutPut("CoarseOrder: " << CoarseOrder << endl);
            OutPut("Projection space is defined up to order 1" << endl);
            exit(-1);
          }
      } // end switch CoarseOrder
    break;
    default:
      OutPut("Invalid shape" << endl);
      exit(-1);
  } // end switch reftype
  return ele;
}

// ADDED ON 17.06.2011 BY SASHI
void AddStreamlineTerm(TSquareMatrix3D* A, TFEFunction3D *uh1,
                       TFEFunction3D *uh2, TFEFunction3D *uh3,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  int N_UsedElements, N_Points;  
  int i1, i2;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;

  TCollection *Coll;
  TFESpace3D *fespace;
  FE3D CurrEleID, UsedElements[2];
  TFE3D *CurrentElement, *CoarseElement;
  TBaseFunct3D *BF, *CoarseBF;
  BaseFunct3D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };

  double *xi, *eta, *zeta, *weights;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[2*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[2*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;  
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy, valz;
  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s, hK;
  double *Entries;
  double *Values1, *Values2, *Values3;
  double BValue[MaxN_BaseFunctions3D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE3D(i, cell);
    CurrentElement = TFEDatabase3D::GetFE3D(CurrEleID);

    BF = CurrentElement->GetBaseFunct3D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement3D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase3D::GetFE3D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct3D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase3D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();
    Values3 = uh3->GetValues();
    
    PCValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    ChildValuesX = TFEDatabase3D::GetOrigElementValues(BF_ID, D100);
    ChildValuesY = TFEDatabase3D::GetOrigElementValues(BF_ID, D010);
    ChildValuesZ = TFEDatabase3D::GetOrigElementValues(BF_ID, D001);   
    ChildValues  = TFEDatabase3D::GetOrigElementValues(BF_ID, D000);

    memset(H, 0, N_CoarseDOF*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];     
      ChildValue  = ChildValues[j];
      w = AbsDetjk[j]*weights[j];
      valx = 0.0;
      valy = 0.0;
      valz = 0.0;
      
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
        valz += ChildValue[k]*Values3[l];        
      }
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k] + valz*ChildValueZ[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm


// stabilisation of full gradient (scalar)
void UltraLocalProjection(TSquareMatrix3D* A, 
                          double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  int N_UsedElements, N_Points;  
  int i1, i2;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;

  TCollection *Coll;
  TFESpace3D *fespace;
  FE3D CurrEleID, UsedElements[2];
  TFE3D *CurrentElement, *CoarseElement;
  TBaseFunct3D *BF, *CoarseBF;
  BaseFunct3D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { FALSE, FALSE };

  double *xi, *eta, *zeta, *weights;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue; 
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;  
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy, valz;
  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s, hK;
  double *Entries;
  double BValue[MaxN_BaseFunctions3D];

  fespace = A->GetFESpace();
  ActiveBound = fespace->GetActiveBound();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();
  // cout << "" << endl;

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    CurrEleID = fespace->GetFE3D(i, cell);
    CurrentElement = TFEDatabase3D::GetFE3D(CurrEleID);

    BF = CurrentElement->GetBaseFunct3D();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement3D(cell, CoarseOrder);

    // approx (index 1) and proj (index 0) space
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    CoarseElement = TFEDatabase3D::GetFE3D(UsedElements[0]);
    CoarseBF = CoarseElement->GetBaseFunct3D();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    TFEDatabase3D::GetOrig(N_UsedElements, UsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && fabs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*SizeOfDouble);
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    DOF = GlobalNumbers + BeginIndex[i];
    
    PCValues = TFEDatabase3D::GetOrigElementValues(CoarseBF_ID, D000);

    ChildValuesX = TFEDatabase3D::GetOrigElementValues(BF_ID, D100);
    ChildValuesY = TFEDatabase3D::GetOrigElementValues(BF_ID, D010);
    ChildValuesZ = TFEDatabase3D::GetOrigElementValues(BF_ID, D001);   
 
    memset(H, 0, N_CoarseDOF*3*N_DOF*SizeOfDouble);

    memset(LocMatrix, 0, N_DOF*N_DOF*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];     
      w = AbsDetjk[j]*weights[j];

      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*3*N_DOF+l      ] += val*ChildValueX[l];
          H[k*3*N_DOF+l+N_DOF] += val*ChildValueY[l];
          H[k*3*N_DOF+l+2*N_DOF] += val*ChildValueZ[l];  
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*(  ChildValueX[k]*ChildValueX[l]
                                     + ChildValueY[k]*ChildValueY[l]
                                     + ChildValueZ[k]*ChildValueZ[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*3*N_DOF*SizeOfDouble);

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 3*N_DOF, 3*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */
    
    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*3*N_DOF+l      ]*H[i2*3*N_DOF+m      ]
                                           +H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF]
                                           +H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*3*N_DOF+l      ] * H[i2*3*N_DOF+m      ];
          s += P[i2*3*N_DOF+l+N_DOF] * H[i2*3*N_DOF+m+N_DOF];
          s += P[i2*3*N_DOF+l+2*N_DOF] * H[i2*3*N_DOF+m+2*N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*3*N_DOF+m      ] * H[i1*3*N_DOF+l      ];
          s += P[i1*3*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l+N_DOF];
          s += P[i1*3*N_DOF+m+2*N_DOF] * H[i1*3*N_DOF+l+2*N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm



#endif
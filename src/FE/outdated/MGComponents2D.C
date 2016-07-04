// =======================================================================
// MGComponents2D.C
//
// Purpose:     components for multigrid in 2d
//
// Author:      Gunar Matthies          27.01.1999
//              Volker John             27.10.1999  
//
//              parallel methods  (Sashikumaar Ganesan) 13.10.2009
// =======================================================================

#ifdef _MPI
#  include "mpi.h"
#endif

#include <LinAlg.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <NodalFunctional2D.h>
#include <MooNMD_Io.h>

#include <stdlib.h>
#include <string.h>

#ifdef __2D__

#define AT(i,j) (a[(j)*LDA+(i)])
#define A(i,j) (a[(i)*LDA+(j)])

/* solve diagonal Vanka system */
void SolveDiagonalVanka2D(double *a, double *b, int N_U, int N_P, int LDA)
// Arguments:
//    a         double array which contains the system matrix
//              A(i,j) = A[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_U       number of velocity unknowns
//    N_P       number of presure unknowns
//    LDA       leading dimension of matrix a
{
  int i,j, row;
  int N_Eqn;
  double pp, dp, ffp, tmp;
  double Ai[MaxN_BaseFunctions2D];

  N_Eqn = 2*N_U+N_P;
  row = 2*N_U;

  /* for(ii=0;ii<N_Eqn;ii++)
  {
    cout << ii;
    for(jj=0;jj<N_Eqn;jj++)
      cout << setw(8) << A(ii,jj);
    cout << endl;
  }
  cout << endl;*/


  if(N_P==1)
  {
    for(i=0;i<N_U;i++)
      Ai[i] = 1/AT(i,i);
  
    dp = 0;
    ffp = b[row];
    for(i=0,j=N_U; i<N_U; i++,j++)
    {
      tmp = Ai[i];
      dp  -= tmp * (AT(row, i) * AT(i, row) + AT(row, j) * AT(j, row));
      ffp -= tmp * (AT(row, i) * b[i] + AT(row, j) * b[j]);
    }
  
    pp = ffp / dp;
    b[row] = pp;
  
    for(i=0, j=N_U; i<N_U; i++, j++)
    {
      tmp = Ai[i];
      b[i] = tmp * (b[i] - AT(i, row) * pp);
      b[j] = tmp * (b[j] - AT(j, row) * pp);
    }
  }
  else
    SolveLinearSystemLapack(a, b, N_Eqn, N_Eqn);

}

/** prolongate */
void Prolongate(const TFESpace2D *CoarseSpace, 
        const TFESpace2D *FineSpace, double *CoarseFunction, 
        double *FineFunction, double *aux)

{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs;// N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int Index;
  double *entry;

  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
//  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  //cout << "N_FineCells: " << N_FineCells << endl;
  //cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_FineDOFs);
  memset(FineFunction, 0, SizeOfDouble*N_FineDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Coarse;l++)
        Val[l] = CoarseFunction[CoarseDOF[l]];

      BaseFunctions->ChangeBF(CoarseColl, parent, Val);

      for(j=0;j<N_Children;j++)
      {
        // cout << "child: " << j << endl;
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do prolongation
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetProlongationMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];

        for(k=0;k<N_Fine;k++)
        {
          s = 0;
          entry = QQ+k*MaxN_BaseFunctions2D;
          for(l=0;l<N_Coarse;l++)
          {
            // s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
            s += entry[l] * Val[l];
            // cout << k << " " << l << " " << entry[l] << endl;
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val2);

        for(k=0;k<N_Fine;k++)
        {
          Index = FineDOF[k];
          FineFunction[Index] += Val2[k];
          aux[Index] += 1;
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      // do prolongation
      QQ = TFEDatabase2D::GetProlongationMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Coarse;l++)
        Val[l] = CoarseFunction[CoarseDOF[l]];

      BaseFunctions->ChangeBF(CoarseColl, cell, Val);

      for(k=0;k<N_Fine;k++)
      {
        s = 0;
        for(l=0;l<N_Coarse;l++)
        {
          s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase2D::GetBaseFunct2D(FineBF)
                      ->ChangeBF(FineColl, cell, Val2);
      for(k=0;k<N_Fine;k++)
      {
        Index = FineDOF[k];
        FineFunction[Index] += Val2[k];
        aux[Index] += 1;
      }
    } // endelse
  } // endfor i

  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] /= aux[i];
  }
}

void Prolongate(const TFESpace2D *CoarseSpace, const TFESpace2D *FineSpace,
        int N_Functions,
        double *CoarseFunction, double *FineFunction, double *aux)

{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int Index;
  double *entry;
  int CoarseOffset, FineOffset, IFunct;

  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  // cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_FineDOFs);
  memset(FineFunction, 0, SizeOfDouble*N_Functions*N_FineDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }


  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do prolongation
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetProlongationMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          CoarseOffset = IFunct*N_CoarseDOFs;
          FineOffset = IFunct*N_FineDOFs;

          for(l=0;l<N_Coarse;l++)
            Val[l] = CoarseFunction[CoarseOffset+CoarseDOF[l]];

          BaseFunctions->ChangeBF(CoarseColl, parent, Val);

          for(k=0;k<N_Fine;k++)
          {
            s = 0;
            entry = QQ+k*MaxN_BaseFunctions2D;
            for(l=0;l<N_Coarse;l++)
            {
              // s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
              s += entry[l] * Val[l];
            } // endfor l
            Val2[k] = s;
          } // endfor k

          TFEDatabase2D::GetBaseFunct2D(FineBF)
                          ->ChangeBF(FineColl, cell, Val2);

          for(k=0;k<N_Fine;k++)
          {
            Index = FineDOF[k];
            FineFunction[FineOffset+Index] += Val2[k];
            aux[Index] += 1;
          } // endfor k
        } // endfor IFunct
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      // do prolongation
      QQ = TFEDatabase2D::GetProlongationMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        CoarseOffset = IFunct*N_CoarseDOFs;
        FineOffset = IFunct*N_FineDOFs;

        for(l=0;l<N_Coarse;l++)
          Val[l] = CoarseFunction[CoarseOffset + CoarseDOF[l]];

        BaseFunctions->ChangeBF(CoarseColl, cell, Val);

        for(k=0;k<N_Fine;k++)
        {
          s = 0;
          for(l=0;l<N_Coarse;l++)
          {
            s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val2);

        for(k=0;k<N_Fine;k++)
        {
          Index = FineDOF[k];
          FineFunction[FineOffset + Index] += Val2[k];
          aux[Index] += 1;
        } // endfor k
      } // endfor IFunct
    } // endelse
  } // endfor i

  Dscal(N_FineDOFs, 1.0/N_Functions, aux);
  for(IFunct=0;IFunct<N_Functions;IFunct++)
  {
    FineOffset = IFunct*N_FineDOFs;
    for(i=0;i<N_FineDOFs;i++)
    {
      FineFunction[FineOffset + i] /= aux[i];
    }
  }
}

/** defect restriction from level+1 to level */
void DefectRestriction(const TFESpace2D *CoarseSpace,
        const TFESpace2D *FineSpace, double *CoarseFunction,
        double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int *DOF, Index;
  
  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  //cout << "N_CoarseCells: " << N_CoarseCells << endl;
  
  //for parallel, use the aux array from the input
#ifndef _MPI
  memset(aux, 0, SizeOfDouble*N_FineDOFs);
#endif    
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
#ifndef _MPI
    DOF = FineGlobalNumbers+FineBeginIndex[i];
    FineId = FineSpace->GetFE2D(i, cell);
    FineElement = TFEDatabase2D::GetFE2D(FineId);
    FineBF = FineElement->GetBaseFunct2D_ID();
    N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();
    for(j=0;j<N_Fine;j++)
      aux[DOF[j]] += 1;
#endif    
  }

  // modify fine function values, will be repaired at end
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] /= aux[i];
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetProlongationMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions2D+k] * Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                        ->ChangeBF(CoarseColl, parent, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          Index = CoarseDOF[k];
          CoarseFunction[Index] += Val2[k];
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;

      // do restriction
/*
      cout << "CoarseId: " << CoarseId << endl;
      cout << "Ref: " << Ref << endl;
      cout << "FineId: " << FineId << endl;
      cout << "j: " << j << endl;
      cout << endl;
*/
      QQ = TFEDatabase2D::GetProlongationMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[l*MaxN_BaseFunctions2D+k]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                      ->ChangeBF(CoarseColl, cell, Val2);

      for(k=0;k<N_Coarse;k++)
        CoarseFunction[CoarseDOF[k]] += Val2[k];
    } // endelse
  } // endfor i

  // repair fine function values since they are modified at beginning
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] *= aux[i];
  }
}

/** defect restriction from level+1 to level */
void DefectRestriction(const TFESpace2D *CoarseSpace, 
                       const TFESpace2D *FineSpace,
        int N_Functions,
        double *CoarseFunction, double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int *DOF, Index;
  int FineOffset, CoarseOffset, IFunct;

  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  // cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_FineDOFs);
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs*N_Functions);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);

    DOF = FineGlobalNumbers+FineBeginIndex[i];
    FineId = FineSpace->GetFE2D(i, cell);
    FineElement = TFEDatabase2D::GetFE2D(FineId);
    FineBF = FineElement->GetBaseFunct2D_ID();
    N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();
    for(j=0;j<N_Fine;j++)
      aux[DOF[j]] += 1;
  }

  // modify fine function values, will be repaired at end
  for(IFunct=0;IFunct<N_Functions;IFunct++)
  {
    FineOffset = IFunct*N_FineDOFs;
    for(i=0;i<N_FineDOFs;i++)
    {
      FineFunction[FineOffset + i] /= aux[i];
    }
  } // endfor IFunct

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetProlongationMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          FineOffset = IFunct*N_FineDOFs;
          CoarseOffset = IFunct*N_CoarseDOFs;

          for(l=0;l<N_Fine;l++)
            Val[l] = FineFunction[FineOffset + FineDOF[l]];

          TFEDatabase2D::GetBaseFunct2D(FineBF)
                            ->ChangeBF(FineColl, cell, Val);

          for(k=0;k<N_Coarse;k++)
          {
            s = 0;
            for(l=0;l<N_Fine;l++)
            {
              s += QQ[l*MaxN_BaseFunctions2D+k] * Val[l];
            } // endfor l
            Val2[k] = s;
          } // endfor k

          TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                          ->ChangeBF(CoarseColl, parent, Val2);

          for(k=0;k<N_Coarse;k++)
          {
            Index = CoarseDOF[k];
            CoarseFunction[CoarseOffset + Index] += Val2[k];
          } // endfor k
        } // endfor IFunct
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;

      // do restriction
/*
      cout << "CoarseId: " << CoarseId << endl;
      cout << "Ref: " << Ref << endl;
      cout << "FineId: " << FineId << endl;
      cout << "j: " << j << endl;
      cout << endl;
*/
      QQ = TFEDatabase2D::GetProlongationMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        FineOffset = IFunct*N_FineDOFs;
        CoarseOffset = IFunct*N_CoarseDOFs;

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineOffset + FineDOF[l]];

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions2D+k]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                        ->ChangeBF(CoarseColl, cell, Val2);

        for(k=0;k<N_Coarse;k++)
          CoarseFunction[CoarseOffset + CoarseDOF[k]] += Val2[k];
      } // endfor IFunct
    } // endelse
  } // endfor i

  // repair fine function values since they are modified at beginning
  for(IFunct=0;IFunct<N_Functions;IFunct++)
  {
    FineOffset = IFunct*N_FineDOFs;
    for(i=0;i<N_FineDOFs;i++)
    {
      FineFunction[FineOffset + i] *= aux[i];
    }
  } // endfor IFunct
}

/** function restriction from level+1 to level */
void RestrictFunction(const TFESpace2D *CoarseSpace, 
    const TFESpace2D *FineSpace,
    double *CoarseFunction, double *FineFunction,
    double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_CoarseDOFs; // N_FineDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  
  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
//  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  // cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_CoarseDOFs);
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    // cout << "i= " << i << "    ";
    // cout << "k= " << k << endl;
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, MaxN_BaseFunctions2D*SizeOfDouble);

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetRestrictionMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions2D+l] * Val[l];
          } // endfor l
          Val2[k] += s;
        } // endfor k
      } // endfor j

      TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                      ->ChangeBF(CoarseColl, parent, Val2);

      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
        aux[l] += 1;
        CoarseFunction[l] += Val2[k];
      } // endfor k
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;

      // do restriction
/*
      cout << "CoarseId: " << CoarseId << endl;
      cout << "Ref: " << Ref << endl;
      cout << "FineId: " << FineId << endl;
      cout << "j: " << j << endl;
      cout << endl;
*/
      QQ = TFEDatabase2D::GetRestrictionMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      TFEDatabase2D::GetBaseFunct2D(FineBF)
                      ->ChangeBF(FineColl, cell, Val);

      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                      ->ChangeBF(CoarseColl, cell, Val2);

      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
        CoarseFunction[l] += Val2[k];
        aux[l] += 1;
      } // endfor k
    } // endelse
  } // endfor i

  for(i=0;i<N_CoarseDOFs;i++)
    CoarseFunction[i] /= aux[i];

/*
  for(i=0;i<N_CoarseDOFs;i++)
    cout << "CoarseFunction[" << i << "]: " << CoarseFunction[i] << endl;

  for(i=0;i<N_FineDOFs;i++)
    cout << "FineFunction[" << i << "]: " << FineFunction[i] << endl;
*/

} // RestrictFunction

/** function restriction from level+1 to level */
void RestrictFunction(const TFESpace2D *CoarseSpace, 
                      const TFESpace2D *FineSpace,
    int N_Functions,
    double *CoarseFunction, double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE2D CoarseId, FineId;
  TFE2D *CoarseElement, *FineElement;
  BaseFunct2D CoarseBF, FineBF;
  TBaseFunct2D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double *Val2;
  int FineOffset, CoarseOffset, Offset, IFunct;

  // begin code
  CoarseColl = CoarseSpace->GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  CoarseBeginIndex = CoarseSpace->GetBeginIndex();
  CoarseGlobalNumbers = CoarseSpace->GetGlobalNumbers();
  N_CoarseDOFs = CoarseSpace->GetN_DegreesOfFreedom();
  
  FineColl = FineSpace->GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  FineBeginIndex = FineSpace->GetBeginIndex();
  FineGlobalNumbers = FineSpace->GetGlobalNumbers();
  N_FineDOFs = FineSpace->GetN_DegreesOfFreedom();

  // cout << "N_FineCells: " << N_FineCells << endl;
  // cout << "N_CoarseCells: " << N_CoarseCells << endl;

  memset(aux, 0, SizeOfDouble*N_CoarseDOFs);
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs*N_Functions);
  Val2 = new double[MaxN_BaseFunctions2D*N_Functions];

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    // cout << "i= " << i << "    ";
    // cout << "k= " << k << endl;
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, parent);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, N_Functions*MaxN_BaseFunctions2D*SizeOfDouble);

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE2D(FineNumber, cell);
        FineElement = TFEDatabase2D::GetFE2D(FineId);
        FineBF = FineElement->GetBaseFunct2D_ID();
        N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase2D::GetRestrictionMatrix2D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          FineOffset = IFunct*N_FineDOFs;
          for(l=0;l<N_Fine;l++)
            Val[l] = FineFunction[FineOffset + FineDOF[l]];

          TFEDatabase2D::GetBaseFunct2D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

          Offset = IFunct*N_Coarse;
          for(k=0;k<N_Coarse;k++)
          {
            s = 0;
            for(l=0;l<N_Fine;l++)
            {
              s += QQ[k*MaxN_BaseFunctions2D+l] * Val[l];
            } // endfor l
            Val2[Offset + k] += s;
          } // endfor k
        } // endfor IFunct
      } // endfor j

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        CoarseOffset = IFunct*N_CoarseDOFs;
        Offset = IFunct*N_Coarse;

        TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                        ->ChangeBF(CoarseColl, parent, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          l=CoarseDOF[k];
          aux[l] += 1;
          CoarseFunction[CoarseOffset + l] += Val2[Offset + k];
        } // endfor k
      } // endfor IFunct
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE2D(i, cell);
      FineElement = TFEDatabase2D::GetFE2D(FineId);
      FineBF = FineElement->GetBaseFunct2D_ID();
      N_Fine = TFEDatabase2D::GetBaseFunct2D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE2D(CoarseNumber, cell);

      CoarseElement = TFEDatabase2D::GetFE2D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct2D_ID();
      BaseFunctions = TFEDatabase2D::GetBaseFunct2D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;

      // do restriction
/*
      cout << "CoarseId: " << CoarseId << endl;
      cout << "Ref: " << Ref << endl;
      cout << "FineId: " << FineId << endl;
      cout << "j: " << j << endl;
      cout << endl;
*/
      QQ = TFEDatabase2D::GetRestrictionMatrix2D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        FineOffset = IFunct*N_FineDOFs;
        CoarseOffset = IFunct*N_CoarseDOFs;

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineOffset + FineDOF[l]];

        TFEDatabase2D::GetBaseFunct2D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase2D::GetBaseFunct2D(CoarseBF)
                        ->ChangeBF(CoarseColl, cell, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          l=CoarseDOF[k];
          CoarseFunction[CoarseOffset + l] += Val2[k];
          aux[l] += 1;
        } // endfor k
      } // endfor IFunct

    } // endelse
  } // endfor i

  Dscal(N_CoarseDOFs, 1.0/N_Functions, aux);
  for(IFunct=0;IFunct<N_Functions;IFunct++)
  {
    CoarseOffset = IFunct*N_CoarseDOFs;
    for(i=0;i<N_CoarseDOFs;i++)
      CoarseFunction[CoarseOffset + i] /= aux[i];
  } // endfor IFunct

  delete Val2;

} // RestrictFunction



/** Navier--Stokes type 1 (NSTYPE==1) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

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

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t;
        v2[index] += value2 * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
  return;
}
void MatVect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

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
      //if (i>=N_Active)
//	  OutPut(i << " " << index << " matv " << value << endl);
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

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
      }
    } // endfor j
    r3[i] = s;
  } // endfor i
}
void Defect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], B[0], B[1], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}
 
/** Navier--Stokes type 2 (NSTYPE==2) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T,
        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  int N_Active;

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

  return;
}
void MatVect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], B[2], B[3], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T,
        double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  int N_Active;

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
}

void Defect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], B[0], B[1], B[2], B[3], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();



  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

/** Navier--Stokes type 3 (NSTYPE==3) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *A12Entries, *A21Entries, *A22Entries;
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

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t;
        v2[index] += value2 * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
  return;
}
void MatVect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], B[0], B[1], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12,
                   TSquareMatrix *A21, TSquareMatrix *A22, 
                   TMatrix *B1, TMatrix *B2,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *A12Entries, *A21Entries, *A22Entries;
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

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
      }
    } // endfor j
    r3[i] = s;
  } // endfor i
}
void Defect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], A[1], A[2], A[3], B[0], B[1], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

/** Navier--Stokes type 4 (NSTYPE==4) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21, 
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2,value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries;
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

  return;
}
void MatVect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries;
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
}
void Defect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}


/** Navier--Stokes type 14 (NSTYPE==14) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
/** with equal order interpolation */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21, 
                    TSquareMatrix *A22, TSquareMatrix *C, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2,value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol, *CRowPtr, *CKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries, *CEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

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
      index = CKCol[j];
      value1 = CEntries[j];
      s += value1 * p[index];
    } // endfor j
    q[i] += s;
  } // endfor i

  return;
}

void MatVect_EquOrd_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  // A[4] is the pressure-pressure matrix  
  CoupledMatVect(A[0], A[1], A[2], A[3], A[4], B[0], B[1], B[2], B[3], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TSquareMatrix *C, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol, *CRowPtr, *CKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries, *CEntries;
  int N_Active;

  // Aij and C have the same structure
  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

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

  // velocity-velocity block, active dofs
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

  // velocity-velocity block, non-active dofs, like Dirichlet
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

  // divergence matrix with velocity
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

  // gradient matrix with pressure
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

  // pressure-pressure block
  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      index = CKCol[j];
      value1 = CEntries[j];
      value = p[index];
      s += value1 * value;
    }
    r3[i] -= s;
   } // endfor i
}

void Defect_EquOrd_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;
  // A[4] is the pressure-pressure matrix
  CoupledDefect(A[0], A[1], A[2], A[3], A[4], B[0], B[1], B[2], B[3], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

void MatVect_Raviart_Thomas_NSE4(TSquareMatrix **A, TMatrix **B, double *x, 
    double *y)
{
  CoupledMatVect(A[0],B[0],x,y);
  return;
}

/**  Darcy Raviart-Thomas
 ( A B' )
 ( B 0  )
 block A: flux x flux
*/
void CoupledMatVect(TSquareMatrix *A, TMatrix *B, double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value;
  double *u, *p;
  double *v, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *BEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();

  BEntries = B->GetEntries();
  
  N_UDOF = A->GetN_Rows();
  N_PDOF = B->GetN_Rows();

  u = x;
  p = u+N_UDOF;

  v = y;
  q = v+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u[index];
    }
    v[i] = s;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      s += value * u[index];

      if(index<N_Active)
      {
        t = p[i];
        v[index] += value * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B, 
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value;
  double *u, *p;
  double *v, *q;
  double *r1, *r2;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *BEntries;
  int N_Active;
  
  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();

  BEntries = B->GetEntries();
  
  N_UDOF = A->GetN_Rows();
  N_PDOF = B->GetN_Rows();

  u = x;
  p  = u+N_UDOF;

  v = b;
  q  = v+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  
  // ( r1 ) = ( v ) _ ( A B' ) ( u )
  // ( r2 )   ( q )   ( B 0  ) ( p )
  
  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = v[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u[index];
    }
    r1[i] = s;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      s -= value * u[index];

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value * t;
      }
    } // endfor j
    r2[i] = s;
  } // endfor i
}

void Defect_Raviart_Thomas_NSE4(TSquareMatrix **A, TMatrix **B, 
                                double *x, double *b, double *r)
{
  CoupledDefect(A[0], B[0], x, b, r);
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    int N_UDOF,N_PDOF;
    N_UDOF = A[0]->GetN_Rows();
    N_PDOF = B[0]->GetN_Rows();
    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  }
}




/** Convection-diffusion problem with VMM */
void MatVectCD_VMM(TSquareMatrix *A, TMatrix *B, TMatrix *C, TSquareMatrix *D, 
                   double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, value, value1;
  double *u1, *p;
  double *v1, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *CRowPtr, *DRowPtr, *CKCol, *DKCol;
  double *AEntries, *BEntries, *CEntries, *DEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();
  BEntries = B->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  DRowPtr = D->GetRowPtr();
  DKCol = D->GetKCol();
  DEntries = D->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = D->GetN_Rows();

  u1 = x;
  p  = u1+N_UDOF;

  v1 = y;
  q  = v1+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  // block A 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
    }
    v1[i] = s;
  } // endfor i

  // block B 
  j = BRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0; 
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      s += value * p[index];
    }
    v1[i] += s;
  } // endfor i

  // block C
  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      index = CKCol[j];
      value1 = CEntries[j];
      s += value1 * u1[index];
    } // endfor j
    q[i] = s;
  } // endfor i

  // block D
  j = DRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = DRowPtr[i+1];
    for(;j<k;j++)
    {
      index = DKCol[j];
      value1 = DEntries[j];
      s += value1 * p[index];
    } // endfor j
    q[i] += s;
  } // endfor i
  return;
}
void MatVect_CD_VMM(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  MatVectCD_VMM(A[0], B[0], B[1], A[1], x, y);
  return;
}

/** r := b - A * x */
void DefectCD_VMM(TSquareMatrix *A, TMatrix *B, TMatrix *C, TSquareMatrix *D, 
                   double *x, double *y, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, value, value1;
  double *u1, *p;
  double *v1, *q, *r1, *r2;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *CRowPtr, *DRowPtr, *CKCol, *DKCol;
  double *AEntries, *BEntries, *CEntries, *DEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B->GetRowPtr();
  BKCol = B->GetKCol();
  BEntries = B->GetEntries();

  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  DRowPtr = D->GetRowPtr();
  DKCol = D->GetKCol();
  DEntries = D->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = D->GetN_Rows();

  u1 = x;
  p  = u1+N_UDOF;

  v1 = y;
  q  = v1+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  // block A 
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u1[index];
    }
    r1[i] = s;
  } // endfor i

  // block B 
  j = BRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0; 
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value = BEntries[j];
      s += value * p[index];
    }
    r1[i] -= s;
  } // endfor i

  // block C
  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      index = CKCol[j];
      value1 = CEntries[j];
      s -= value1 * u1[index];
    } // endfor j
    r2[i] = s;
  } // endfor i

  // block D
  j = DRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = DRowPtr[i+1];
    for(;j<k;j++)
    {
      index = DKCol[j];
      value1 = DEntries[j];
      s += value1 * p[index];
    } // endfor j
    r2[i] -= s;
  } // endfor i
  return;
}

void Defect_CD_VMM(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  DefectCD_VMM(A[0], B[0], B[1], A[1], x, b, r);
  return;
}
 
/** Convection-diffusion problem with VMM [KL02] */
void MatVectCD_VMM_KL02(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, 
                        TMatrix *C1, TMatrix *C2, TSquareMatrix *D, 
                        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *p1, *p2;
  double *v1, *q1, *q2;
  const int *ARowPtr, *B1RowPtr, *AKCol, *B1KCol;
  const int *C1RowPtr, *DRowPtr, *C1KCol, *DKCol;
  double *AEntries, *B1Entries, *C1Entries, *DEntries, *B2Entries, *C2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  B1RowPtr = B1->GetRowPtr();
  B1KCol = B1->GetKCol();
  B1Entries = B1->GetEntries();

  B2Entries = B2->GetEntries();

  C1RowPtr = C1->GetRowPtr();
  C1KCol = C1->GetKCol();
  C1Entries = C1->GetEntries();

  C2Entries = C2->GetEntries();

  DRowPtr = D->GetRowPtr();
  DKCol = D->GetKCol();
  DEntries = D->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = D->GetN_Rows();

  u1 = x;
  p1  = u1+N_UDOF;
  p2  = p1+N_PDOF;

  v1 = y;
  q1  = v1+N_UDOF;
  q2  = q1+N_PDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  // block A 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
    }
    v1[i] = s;
  } // endfor i

  // block B1, B2 
  j = B1RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0; 
    t = 0;
    k = B1RowPtr[i+1];
    for(;j<k;j++)
    {
      index = B1KCol[j];
      value = B1Entries[j];
      s += value * p1[index];
      value = B2Entries[j];
      t += value * p2[index];      
    }
    v1[i] += s+t;
  } // endfor i

  // block C1, C2
  j = C1RowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    t = 0;
    k = C1RowPtr[i+1];
    for(;j<k;j++)
    {
      index = C1KCol[j];
      value1 = C1Entries[j];
      s += value1 * u1[index];
      value2 = C2Entries[j];
      t += value2 * u1[index];
   } // endfor j
    q1[i] = s;
    q2[i] = t;
  } // endfor i

  // block D
  j = DRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    t = 0;
    k = DRowPtr[i+1];
    for(;j<k;j++)
    {
      index = DKCol[j];
      value1 = DEntries[j];
      s += value1 * p1[index];
      t += value1 * p2[index];
    } // endfor j
    q1[i] += s;
    q2[i] += t;
  } // endfor i
  return;
}
void MatVect_CD_VMM_KL02(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  MatVectCD_VMM_KL02(A[0], B[0], B[1], B[2], B[3], A[1], x, y);
  return;
}
void DefectCD_VMM_KL02(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, 
                        TMatrix *C1, TMatrix *C2, TSquareMatrix *D, 
                        double *x, double *y, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *p1, *p2;
  double *v1, *q1, *q2, *r1, *r2, *r3;
  const int *ARowPtr, *B1RowPtr, *AKCol, *B1KCol;
  const int *C1RowPtr, *DRowPtr, *C1KCol, *DKCol;
  double *AEntries, *B1Entries, *C1Entries, *DEntries, *B2Entries, *C2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  B1RowPtr = B1->GetRowPtr();
  B1KCol = B1->GetKCol();
  B1Entries = B1->GetEntries();

  B2Entries = B2->GetEntries();

  C1RowPtr = C1->GetRowPtr();
  C1KCol = C1->GetKCol();
  C1Entries = C1->GetEntries();

  C2Entries = C2->GetEntries();

  DRowPtr = D->GetRowPtr();
  DKCol = D->GetKCol();
  DEntries = D->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = D->GetN_Rows();

  u1 = x;
  p1  = u1+N_UDOF;
  p2  = p1+N_PDOF;

  v1 = y;
  q1  = v1+N_UDOF;
  q2  = q1+N_PDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_PDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  // block A 
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u1[index];
    }
    r1[i] = s;
  } // endfor i

  // block B1, B2 
  j = B1RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0; 
    k = B1RowPtr[i+1];
    for(;j<k;j++)
    {
      index = B1KCol[j];
      value = B1Entries[j];
      value1 = B2Entries[j];
      s += value * p1[index] + value1 * p2[index];
    }
    r1[i] -= s;
  } // endfor i

  // block C1, C2
  j = C1RowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q1[i];
    t = q2[i];
    k = C1RowPtr[i+1];
    for(;j<k;j++)
    {
      index = C1KCol[j];
      value1 = C1Entries[j];
      s -= value1 * u1[index];
      value2 = C2Entries[j];
      t -= value2 * u1[index];
   } // endfor j
    r2[i] = s;
    r3[i] = t;
  } // endfor i

  // block D
  j = DRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    t = 0;
    k = DRowPtr[i+1];
    for(;j<k;j++)
    {
      index = DKCol[j];
      value1 = DEntries[j];
      s += value1 * p1[index];
      t += value1 * p2[index];
    } // endfor j
    r2[i] -= s;
    r3[i] -= t;
  } // endfor i
  return;
}

void Defect_CD_VMM_KL02(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  DefectCD_VMM_KL02(A[0], B[0], B[1], B[2], B[3], A[1], x, b, r);
  return;
}

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVectLV96(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *y, double delta)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

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
 
  for(i=0;i<N_UDOF;i++) // A*u
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
  for(i=0;i<N_PDOF;i++) // B^T*p and B*u
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s += value1 * u1[index] + value2 * u2[index];

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t; 
        v2[index] += value2 * t;
      }
    } // endfor j
    q[i] = s;               // B*u
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++) // B^T * B * u = B^T * q
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];

      if(index<N_Active)
      {
        t = q[i];
        v1[index] += value1 * t/delta;
        v2[index] += value2 * t/delta;
      }
    } // endfor j
  } // endfor i
  return;
}

/** r := b - A * x */
void CoupledDefectLV96(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        double *x, double *b, double *r, double delta)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q, *help;
  double *r1, *r2, *r3;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

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
  for(i=0;i<N_UDOF;i++) // -A*u
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

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
      }
    } // endfor j
    r3[i] = s;
  } // endfor i
  
  help = new double [N_PDOF];
  memset(help,0,N_PDOF*SizeOfDouble);
  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++) // B*u
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
    help[i] = s;               // help = B*u
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++) 
  {
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];

      if(index<N_Active)
      {
        t = help[i];
        r1[index] -= value1 * t/delta;
        r2[index] -= value2 * t/delta;
      }
    } // endfor j
  } // endfor i
  delete help;
}

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVectMortar(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *matrix_mortar, 
        double *x, double *y)
{
  int N_UDOF, N_PDOF,N_MORTDOF;
  int i,j,k,l,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p, *m1, *m2;
  double *v1, *v2, *q, *n1, *n2;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol, *MortarRowPtr, *MortarKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries, *MortarEntries;
  int N_Active;

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

  MortarRowPtr = matrix_mortar->GetRowPtr();
  MortarKCol =  matrix_mortar->GetKCol();
  MortarEntries = matrix_mortar->GetEntries();
  
  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();
  N_MORTDOF = matrix_mortar->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;
  m1 = p+N_PDOF;
  m2 = m1+N_MORTDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;
  n1 = q+N_PDOF;
  n2 = n1+N_MORTDOF;

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

  for (i=0;i<N_MORTDOF;i++)
    {
      l=MortarRowPtr[i];
      k=MortarRowPtr[i+1];
      s=t=0;
      for(j=l;j<k;j++)
        {
          index = MortarKCol[j];
          s+=MortarEntries[j]*u1[index];
          t+=MortarEntries[j]*u2[index];
        }
      n1[i] = s;
      n2[i] = t;
    } 

  for (i=0;i<N_MORTDOF;i++)
    {
      l=MortarRowPtr[i];
      k=MortarRowPtr[i+1];
      for(j=l;j<k;j++)
        {
          index = MortarKCol[j];
          v1[index]+=MortarEntries[j]*m1[i];
          v2[index]+=MortarEntries[j]*m2[i];
        }
    } 
 return;
}

#endif

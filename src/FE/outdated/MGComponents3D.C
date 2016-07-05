// =======================================================================
// MGComponents3D.C
//
// Purpose:     components for multigrid in 3d
//
// Author:      Gunar Matthies          27.01.1999
//              Volker John             27.10.1999  
// 
//              parallel methods  (Sashikumaar Ganesan) 19.09.2010
// =======================================================================

#ifdef _MPI
# include "mpi.h"
#endif

#include <LinAlg.h>
#include <Database.h>
#include <FEDatabase3D.h>

#include <stdlib.h>
#include <string.h>

#ifdef __3D__

#define AT(i,j) (a[(j)*LDA+(i)])
#define A(i,j) (a[(i)*LDA+(j)])

void SolveDiagonalVanka3D(double *a, double *b, int N_U, int N_P, int LDA)
// Arguments:
//    a         double array which contains the system matrix
//              A(i,j) = A[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_U       number of velocity unknowns
//    N_P       number of presure unknowns
//    LDA       leading dimension of matrix a
{
  int i,j,k, row;
  int N_Eqn;
  double pp, dp, ffp, tmp;
  double Ai[8*MaxN_BaseFunctions3D];

//  int ii, jj;

//  double S[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
//  double r[MaxN_BaseFunctions3D];
//  double q[MaxN_BaseFunctions3D];
//  int NU3;

  N_Eqn = 3*N_U+N_P;
  row = 3*N_U;

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
    for(i=0,j=N_U,k=2*N_U; i<N_U; i++,j++,k++)
    {
      tmp = Ai[i];
      dp  -= tmp * (AT(row, i) * AT(i, row) + AT(row, j) * AT(j, row)
        + AT(row, k) * AT(k, row));
      ffp -= tmp * (AT(row, i) * b[i] + AT(row, j) * b[j]+ AT(row, k) * b[k]);
    }
 
    pp = ffp / dp;
    b[row] = pp;
  
    for(i=0, j=N_U,k=2*N_U; i<N_U; i++, j++,k++)
    {
      tmp = Ai[i];
      b[i] = tmp * (b[i] - AT(i, row) * pp);
      b[j] = tmp * (b[j] - AT(j, row) * pp);
      b[k] = tmp * (b[k] - AT(k, row) * pp);
    }
  }
  else
  {
    // OutPut("SolveDiagonalVanka3D" << endl);
    //if(fabs(TDatabase::ParamDB->REACTOR_P25 - 8)<1e-10)
    {
      SolveLinearSystemLapack(a, b, N_Eqn, N_Eqn);
    }
    /* else
    {
      NU3 = 3*N_U;
      for(i=0;i<N_P;i++)
      {
        for(j=0;j<N_P;j++)
        {
          pp = 0;
          for(k=0;k<NU3;k++)
          {
            pp += A(i+NU3,k)*A(k,j+NU3)/A(k,k);
          } // endfor k
          S[i*N_P+j] = pp;
        } // endfor j
      } // endfor i

      for(i=0;i<N_P;i++)
      {
        pp = -b[NU3+i];
        for(k=0;k<NU3;k++)
        {
          pp += A(i+NU3,k)*b[k]/A(k,k);
        } // endfor k
        r[i] = pp;
      } // endfor i

      for(i=0;i<N_P;i++)
      {
        for(j=0;j<N_P;j++)
        {
          cout << setw(4) << i << setw(4) << j << setw(25) << S[i*N_P+j] << endl;
        }
      }

      SolveLinearSystem(S, r, N_P, N_P);

      for(i=0;i<N_P;i++)
        b[NU3+i] = r[i];

      for(i=0;i<NU3;i++)
      {
        pp = b[i];
        for(j=0;j<N_P;j++)
        {
          pp -= A(i,j+NU3)*r[j];
        } // endfor j
        b[i] = pp/A(i,i);
      } // endfor i

      for(i=0;i<N_Eqn;i++)
        b[i] *= TDatabase::ParamDB->REACTOR_P24;
	}*/
  }
}


double tP=0.0,tR=0.0;
/** prolongate */
void Prolongate(const TFESpace3D *CoarseSpace, 
                const TFESpace3D *FineSpace, double *CoarseFunction, 
        double *FineFunction, double *aux)

{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs; // N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
//  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
  int *DOF, Index;
  double *entry;
  double t1,t2;
#ifdef _MPI
  t1 = MPI_Wtime();
#else
  t1 = GetTime();
#endif
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

#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,j,k,l,cell,DOF,FineId,FineElement,FineBF,N_Fine, \
                                             parent, N_Children, CoarseNumber, CoarseId, CoarseElement, \
                                             CoarseBF, BaseFunctions, Ref, FineNumber,\
                                             QQ, FineDOF, CoarseDOF, Val, s, Val2, \
                                             Index, N_Coarse, entry)
{
#pragma omp for schedule(guided) 
#endif
  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
    
    DOF = FineGlobalNumbers+FineBeginIndex[i];
    FineId = FineSpace->GetFE3D(i, cell);
    FineElement = TFEDatabase3D::GetFE3D(FineId);
    FineBF = FineElement->GetBaseFunct3D_ID();
    N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();
    for(j=0;j<N_Fine;j++)
#ifdef _HYBRID      
      #pragma omp atomic 
#endif  
      aux[DOF[j]] += 1;
  }
  
#ifdef _HYBRID
#pragma omp for schedule(guided) 
#endif
  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }


#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }
  
#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
#ifdef _MPI
    if(cell->IsHaloCell())   continue;
#endif
    
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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
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
// #ifdef _HYBRID
// 	#pragma omp critical
// #endif
        {
         k = cell->GetClipBoard();
         cell->SetClipBoard(-2);
        }
        
#ifdef _HYBRID
	  if(k==-2) continue;
#endif
	FineNumber = -(k+10);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

#ifdef _HYBRID
       #pragma omp critical
#endif
	{
	  QQ = TFEDatabase3D::GetProlongationMatrix3D 
                (CoarseId, Ref, FineId, j);
	}

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];

        for(k=0;k<N_Fine;k++)
        {
          s = 0;
          entry = QQ+k*MaxN_BaseFunctions3D;
          for(l=0;l<N_Coarse;l++)
          {
            // s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
            s += entry[l] * Val[l];
            // cout << k << " " << l << " " << entry[l] << endl;
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                        ->ChangeBF(FineColl, cell, Val2);
			
        for(k=0;k<N_Fine;k++)
        {
          Index = FineDOF[k];
	 {
#ifdef _HYBRID
       #pragma omp atomic
#endif  
          FineFunction[Index] += Val2[k];
	 }
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

#ifdef _HYBRID
       #pragma omp critical
#endif
	{
	  QQ = TFEDatabase3D::GetProlongationMatrix3D 
              (CoarseId, Ref, FineId, 0);
	}

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
          s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase3D::GetBaseFunct3D(FineBF)
                      ->ChangeBF(FineColl, cell, Val2);
      
      for(k=0;k<N_Fine;k++)
      {
        Index = FineDOF[k];
	{
#ifdef _HYBRID
       #pragma omp atomic
#endif  
         FineFunction[Index] += Val2[k];
	}
      }
    } // endelse
  } // endfor i

#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif 
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] /= aux[i];
  }
#ifdef _HYBRID
} 
#endif

#ifdef _MPI
  t2 = MPI_Wtime();
#else
  t2 = GetTime();
#endif
  tP += (t2-t1);
}

void Prolongate(const TFESpace3D *CoarseSpace, const TFESpace3D *FineSpace,
        int N_Functions,
        double *CoarseFunction, double *FineFunction, double *aux)

{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
//  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
  int Index; //*DOF,
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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

        // do prolongation
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase3D::GetProlongationMatrix3D 
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
            entry = QQ+k*MaxN_BaseFunctions3D;
            for(l=0;l<N_Coarse;l++)
            {
              // s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
              s += entry[l] * Val[l];
            } // endfor l
            Val2[k] = s;
          } // endfor k

          TFEDatabase3D::GetBaseFunct3D(FineBF)
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
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      // do prolongation
      QQ = TFEDatabase3D::GetProlongationMatrix3D 
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
            s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(FineBF)
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
void DefectRestriction(const TFESpace3D *CoarseSpace,
                       const TFESpace3D *FineSpace, double *CoarseFunction,
        double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
//  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
  int *DOF, Index;
//  double *entry;
  double t1,t2;
#ifdef _MPI
  t1 = MPI_Wtime();
#else
  t1 = GetTime();
#endif
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

  memset(aux, 0, SizeOfDouble*N_FineDOFs);
  memset(CoarseFunction, 0, SizeOfDouble*N_CoarseDOFs);
  
  
  
  
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,j,k,l,cell,DOF,FineId,FineElement,FineBF,N_Fine, \
                                             parent, N_Children, CoarseNumber, CoarseId, CoarseElement, \
                                             CoarseBF, BaseFunctions, Ref, FineNumber,\
                                             QQ, FineDOF, CoarseDOF, Val, s, Val2, \
                                             Index, N_Coarse)
{
#pragma omp for schedule(guided) 
#endif
//   set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);

    DOF = FineGlobalNumbers+FineBeginIndex[i];
    FineId = FineSpace->GetFE3D(i, cell);
    FineElement = TFEDatabase3D::GetFE3D(FineId);
    FineBF = FineElement->GetBaseFunct3D_ID();
    N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();
    for(j=0;j<N_Fine;j++)
#ifdef _HYBRID      
      #pragma omp atomic 
#endif      
      aux[DOF[j]] += 1;
  }

#ifdef _HYBRID
#pragma omp for schedule(guided) nowait 
#endif
  // modify fine function values, will be repaired at end
  for(i=0;i<N_FineDOFs;i++)
  {    
    FineFunction[i] /= aux[i];
  }
   
#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
#ifdef _MPI
    if(cell->IsHaloCell())	continue;
#endif
    
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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();
  
      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
// #ifdef _HYBRID
// 	#pragma omp critical
// #endif
	{
         k = cell->GetClipBoard();
         cell->SetClipBoard(-2);
	}
	
#ifdef _HYBRID
	  if(k==-2) continue;
#endif
	  
	FineNumber = -(k+10);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();
#ifdef _HYBRID
       #pragma omp critical
#endif
	{
	  QQ = TFEDatabase3D::GetProlongationMatrix3D 
                (CoarseId, Ref, FineId, j);
	}

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions3D+k] * Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(CoarseBF)
                        ->ChangeBF(CoarseColl, parent, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          Index = CoarseDOF[k];
#ifdef _HYBRID
       #pragma omp atomic
#endif
           CoarseFunction[Index] += Val2[k];
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = NoRef;
#ifdef _HYBRID
       #pragma omp critical
#endif
	{
	  QQ = TFEDatabase3D::GetProlongationMatrix3D 
              (CoarseId, Ref, FineId, 0);
	}

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      TFEDatabase3D::GetBaseFunct3D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[l*MaxN_BaseFunctions3D+k]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase3D::GetBaseFunct3D(CoarseBF)
                      ->ChangeBF(CoarseColl, cell, Val2);

      for(k=0;k<N_Coarse;k++)
      {
#ifdef _HYBRID
       #pragma omp atomic
#endif
	CoarseFunction[CoarseDOF[k]] += Val2[k];
      }
    } // endelse
  } // endfor i

  // repair fine function values since they are modified at beginning
#ifdef _HYBRID
#pragma omp for schedule(guided)
#endif  
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] *= aux[i];
  }
#ifdef _HYBRID
} 
#endif
  
#ifdef _MPI
  t2 = MPI_Wtime();
#else
  t2 = GetTime();
#endif
  tR+=(t2-t1);
  
}

/** defect restriction from level+1 to level */
void DefectRestriction(const TFESpace3D *CoarseSpace, const TFESpace3D *FineSpace,
        int N_Functions,
        double *CoarseFunction, double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
//  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
  int *DOF, Index;
//  double *entry;
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
    FineId = FineSpace->GetFE3D(i, cell);
    FineElement = TFEDatabase3D::GetFE3D(FineId);
    FineBF = FineElement->GetBaseFunct3D_ID();
    N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();
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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase3D::GetProlongationMatrix3D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          FineOffset = IFunct*N_FineDOFs;
          CoarseOffset = IFunct*N_CoarseDOFs;

          for(l=0;l<N_Fine;l++)
            Val[l] = FineFunction[FineOffset + FineDOF[l]];

          TFEDatabase3D::GetBaseFunct3D(FineBF)
                            ->ChangeBF(FineColl, cell, Val);

          for(k=0;k<N_Coarse;k++)
          {
            s = 0;
            for(l=0;l<N_Fine;l++)
            {
              s += QQ[l*MaxN_BaseFunctions3D+k] * Val[l];
            } // endfor l
            Val2[k] = s;
          } // endfor k

          TFEDatabase3D::GetBaseFunct3D(CoarseBF)
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
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
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
      QQ = TFEDatabase3D::GetProlongationMatrix3D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        FineOffset = IFunct*N_FineDOFs;
        CoarseOffset = IFunct*N_CoarseDOFs;

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineOffset + FineDOF[l]];

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions3D+k]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(CoarseBF)
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
void RestrictFunction(const TFESpace3D *CoarseSpace, 
                      const TFESpace3D *FineSpace,
    double *CoarseFunction, double *FineFunction,
    double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_CoarseDOFs; //N_FineDOFs
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
//  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
//  int *DOF, Index;
//  double *entry;

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
  
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,cell,k)
{
#pragma omp for schedule(static) nowait 
#endif
  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }
#ifdef _HYBRID
#pragma omp for schedule(static) nowait 
#endif
  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }
#ifdef _HYBRID
#pragma omp for schedule(static) nowait 
#endif
  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }
#ifdef _HYBRID
}
#endif

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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, MaxN_BaseFunctions3D*SizeOfDouble);
#ifdef _HYBRID
// #pragma omp parallel default(shared) private(j,k,s,l,cell,FineNumber,FineId,FineElement,FineBF,N_Fine,QQ,FineDOF,Val2,Index)
// #pragma omp for schedule(guided) nowait 
#endif
      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase3D::GetRestrictionMatrix3D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions3D+l] * Val[l];
          } // endfor l
          Val2[k] += s;
        } // endfor k
      } // endfor j

      TFEDatabase3D::GetBaseFunct3D(CoarseBF)
                      ->ChangeBF(CoarseColl, parent, Val2);
		      
#ifdef _HYBRID
// #pragma omp parallel default(shared) private(k,l)
// #pragma omp for schedule(guided) nowait 
#endif
      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
        aux[l] += 1;
// 	 #pragma omp critical
        CoarseFunction[l] += Val2[k];
      } // endfor k
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
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
      QQ = TFEDatabase3D::GetRestrictionMatrix3D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

#ifdef _HYBRID
// #pragma omp parallel default(shared) private(k,s,l,Val2,Index)
{
// #pragma omp for schedule(guided) nowait 
#endif
      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      TFEDatabase3D::GetBaseFunct3D(FineBF)
                      ->ChangeBF(FineColl, cell, Val);

#ifdef _HYBRID
// #pragma omp for schedule(guided) nowait 
#endif
      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      TFEDatabase3D::GetBaseFunct3D(CoarseBF)
                      ->ChangeBF(CoarseColl, cell, Val2);

#ifdef _HYBRID
// #pragma omp for schedule(guided) nowait 
#endif
      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
#ifdef _HYBRID
	#pragma omp critical
#endif
	{
         CoarseFunction[l] += Val2[k];
         aux[l] += 1;
	}
      } // endfor k
#ifdef _HYBRID
}
#endif
    } // endelse
  } // endfor i

#ifdef _HYBRID
#pragma omp parallel default(shared) private(i)
#pragma omp for schedule(static) nowait 
#endif
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
void RestrictFunction(const TFESpace3D *CoarseSpace, const TFESpace3D *FineSpace,
    int N_Functions,
    double *CoarseFunction, double *FineFunction, double *aux)
{
  int i,j,k,l;
  TBaseCell *cell, *parent;
  TCollection *CoarseColl, *FineColl;
  FE3D CoarseId, FineId;
  TFE3D *CoarseElement, *FineElement;
  BaseFunct3D CoarseBF, FineBF;
  TBaseFunct3D *BaseFunctions;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int *CoarseBeginIndex, *FineBeginIndex;
  int *CoarseGlobalNumbers, *FineGlobalNumbers;
  int FineNumber, CoarseNumber;
  int *FineDOF, *CoarseDOF;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
//  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions3D];
  double *Val2;
//  int *DOF, Index;
//  double *entry;
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
  Val2 = new double[MaxN_BaseFunctions3D*N_Functions];

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
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, parent);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, N_Functions*MaxN_BaseFunctions3D*SizeOfDouble);

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace->GetFE3D(FineNumber, cell);
        FineElement = TFEDatabase3D::GetFE3D(FineId);
        FineBF = FineElement->GetBaseFunct3D_ID();
        N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = TFEDatabase3D::GetRestrictionMatrix3D 
                (CoarseId, Ref, FineId, j);

        FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
        CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

        for(IFunct=0;IFunct<N_Functions;IFunct++)
        {
          FineOffset = IFunct*N_FineDOFs;
          for(l=0;l<N_Fine;l++)
            Val[l] = FineFunction[FineOffset + FineDOF[l]];

          TFEDatabase3D::GetBaseFunct3D(FineBF)
                          ->ChangeBF(FineColl, cell, Val);

          Offset = IFunct*N_Coarse;
          for(k=0;k<N_Coarse;k++)
          {
            s = 0;
            for(l=0;l<N_Fine;l++)
            {
              s += QQ[k*MaxN_BaseFunctions3D+l] * Val[l];
            } // endfor l
            Val2[Offset + k] += s;
          } // endfor k
        } // endfor IFunct
      } // endfor j

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        CoarseOffset = IFunct*N_CoarseDOFs;
        Offset = IFunct*N_Coarse;

        TFEDatabase3D::GetBaseFunct3D(CoarseBF)
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
      FineId = FineSpace->GetFE3D(i, cell);
      FineElement = TFEDatabase3D::GetFE3D(FineId);
      FineBF = FineElement->GetBaseFunct3D_ID();
      N_Fine = TFEDatabase3D::GetBaseFunct3D(FineBF)->GetDimension();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace->GetFE3D(CoarseNumber, cell);

      CoarseElement = TFEDatabase3D::GetFE3D(CoarseId);
      CoarseBF = CoarseElement->GetBaseFunct3D_ID();
      BaseFunctions = TFEDatabase3D::GetBaseFunct3D(CoarseBF);
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
      QQ = TFEDatabase3D::GetRestrictionMatrix3D 
              (CoarseId, Ref, FineId, 0);

      FineDOF = FineGlobalNumbers+FineBeginIndex[FineNumber];
      CoarseDOF = CoarseGlobalNumbers+CoarseBeginIndex[CoarseNumber];

      for(IFunct=0;IFunct<N_Functions;IFunct++)
      {
        FineOffset = IFunct*N_FineDOFs;
        CoarseOffset = IFunct*N_CoarseDOFs;

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineOffset + FineDOF[l]];

        TFEDatabase3D::GetBaseFunct3D(FineBF)
                        ->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        TFEDatabase3D::GetBaseFunct3D(CoarseBF)
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
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries, *B3Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
      t += value * u2[index];
      u += value * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t;
        v2[index] += value2 * t;
        v3[index] += value3 * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
  return;
}
void MatVect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], B[2], x, y);
  return;
}
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *AEntries, *B1Entries, *B2Entries, *B3Entries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  N_Active = A->GetActiveBound();

  //OutPut("r0 " << Ddot(N_UDOF+N_PDOF,r,r) << endl);
  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      //cout << "a " << value << endl;
      s -= value * u1[index];
      t -= value * u2[index];
      u -= value * u3[index];
      //OutPut(k<< " " << value << " " << endl);
    }
    //OutPut(s << " " << t << " " << u << endl);
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
  } // endfor i
  //OutPut("r1 " << Ddot(N_UDOF+N_PDOF,r,r) << endl);

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
      value3 = B3Entries[j];
      //cout << "b " << i << " " << index << " " << value1 << " " << value2 << " " << value3 << endl;
      s -= value1 * u1[index] + value2 * u2[index] + 
        value3 * u3[index];

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
        r3[index] -= value3 * t;
      }
    } // endfor j
    r4[i] = s;
  } // endfor i
  // OutPut("r2 " << Ddot(N_UDOF+N_PDOF,r,r) << endl);
}
void Defect_NSE1(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

// #ifdef _MPI
//  N_UDOF = 0;
//  N_PDOF = 0;
// 
//  if(TDatabase::ParamDB->ActiveProcess)
// #endif
  {
   CoupledDefect(A[0], B[0], B[1], B[2], x, b, r);
   N_UDOF = A[0]->GetN_Rows();
   N_PDOF = B[0]->GetN_Rows();
  }

  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

 
/** Navier--Stokes type 2 (NSTYPE==2) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    TMatrix *B1T, TMatrix *B2T, TMatrix *B3T, 
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries, *B3Entries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  N_Active = A->GetActiveBound();
  j = ARowPtr[0];
 
  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
      t += value * u2[index];
      u += value * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    q[i] = s;
  } // endfor i
 
  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    v1[i] += s;
    v2[i] += t;
    v3[i] += u;
  } // endfor i

  return;
}
void MatVect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], B[2], B[3], B[4], B[5], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries, *B3Entries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();

  N_UDOF = A->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  N_Active = A->GetActiveBound();

  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u1[index];
      t -= value * u2[index];
      u -= value * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
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
      value3 = B3Entries[j];
      s -= value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    r4[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
    r3[i] -= u;
  } // endfor i
}
void Defect_NSE2(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], B[0], B[1], B[2], B[3], B[4], B[5], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();
  
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);

  return;
}

/** Navier--Stokes type 3 (NSTYPE==3) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                    TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                    TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                    TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  N_UDOF = A11->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  N_Active = A11->GetActiveBound();
  j = ARowPtr[0];
  // real dof
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      s += value * u1[index] + value1 * u2[index] + value2 * u3[index];
      t += value3* u1[index] + value4 * u2[index] + value5 * u3[index];
      u += value6* u1[index] + value7 * u2[index] + value8 * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
  } // endfor i
  // Dirichlet and hanging nodes
  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s += value * u1[index];
      t += value1 * u2[index];
      u += value2 * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        +  value3 * u3[index];

      if(index<N_Active)
      {
        t = p[i];
        v1[index] += value1 * t;
        v2[index] += value2 * t;
        v3[index] += value3 * t;
      }
    } // endfor j
    q[i] = s;
  } // endfor i
  return;
}
void MatVect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8],
                 B[0], B[1], B[2], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                   TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                   TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                   TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3;
  double value4, value5, value6, value7, value8;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  N_UDOF = A11->GetN_Rows();
  N_PDOF = B1->GetN_Rows();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  N_Active = A11->GetActiveBound();

  j = ARowPtr[0];
 
  for(i=0;i<N_Active;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      s -= value * u1[index] + value1 * u2[index] + value2 * u3[index];
      t -= value3* u1[index] + value4 * u2[index] + value5 * u3[index];
      u -= value6* u1[index] + value7 * u2[index] + value8 * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
  } // endfor i

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s -= value * u1[index];
      t -= value1 * u2[index];
      u -= value2 * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
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
      value3 = B3Entries[j];
      s -= value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];

      if(index<N_Active)
      {
        t = p[i];
        r1[index] -= value1 * t;
        r2[index] -= value2 * t;
        r3[index] -= value3 * t;
      }
    } // endfor j
    r4[i] = s;
  } // endfor i
}
void Defect_NSE3(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

// #ifdef _MPI
//  N_UDOF = 0;
//  N_PDOF = 0;
// 
//  if(TDatabase::ParamDB->ActiveProcess)
// #endif
  {
   CoupledDefect(A[0], A[1], A[2], A[3], A[4], A[5], A[6],A[7], A[8],
                B[0], B[1], B[2], x, b, r);
   N_UDOF = A[0]->GetN_Rows();
   N_PDOF = B[0]->GetN_Rows();
  }

  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}

/** Navier--Stokes type 4 (NSTYPE==4) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                    TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                    TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                    TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8, value9, value10, value11;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  N_UDOF = A11->GetN_Rows();
  N_Active = A11->GetActiveBound();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();
  N_PDOF = B1->GetN_Rows();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  j = ARowPtr[0];
  // real dof
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      value9 = u1[index];
      value10 = u2[index];
      value11 = u3[index];
      s += value * value9 + value1 * value10 + value2 * value11;
      t += value3* value9 + value4 * value10 + value5 * value11;
      u += value6* value9 + value7 * value10 + value8 * value11;
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
  } // endfor i
  // Dirichlet and hanging nodes

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s += value * u1[index];
      t += value1 * u2[index];
      u += value2 * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    q[i] = s;
  } // endfor i
 
  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    v1[i] += s;
    v2[i] += t;
    v3[i] += u;
  } // endfor i

  return;
}
void MatVect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8],
                 B[0], B[1], B[2], B[3], B[4], B[5], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                   TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                   TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
                   TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8, value9, value10, value11;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  N_UDOF = A11->GetN_Rows();
  N_Active = A11->GetActiveBound();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();
  N_PDOF = B1->GetN_Rows();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();
 
  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  j = ARowPtr[0];

  for(i=0;i<N_Active;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      value9 = u1[index];
      value10 = u2[index];
      value11 = u3[index];
      s -= value * value9 + value1 * value10+ value2 * value11;
      t -= value3* value9 + value4 * value10+ value5 * value11;
      u -= value6* value9 + value7 * value10 + value8 * value11;
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
  } // endfor i

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s -= value * u1[index];
      t -= value1 * u2[index];
      u -= value2 * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
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
      value3 = B3Entries[j];
      s -= value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    r4[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
    r3[i] -= u;
  } // endfor i
}
void Defect_NSE4(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

// #ifdef _MPI
//  N_UDOF = 0;
//  N_PDOF = 0;
// 
//  if(TDatabase::ParamDB->ActiveProcess)
// #endif
  {
   CoupledDefect(A[0], A[1], A[2], A[3], A[4], A[5], A[6],A[7], A[8],
                B[0], B[1], B[2], B[3], B[4], B[5], x, b, r);
   N_UDOF = A[0]->GetN_Rows();
   N_PDOF = B[0]->GetN_Rows();
  }


  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}


/** Navier--Stokes type 4 (NSTYPE==4) */
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                    TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                    TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
		    TSquareMatrix *C,
                    TMatrix *B1, TMatrix *B2, TMatrix *B3,
                    TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8, value9, value10, value11;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol, *CRowPtr, *CKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  double *CEntries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  N_UDOF = A11->GetN_Rows();
  N_Active = A11->GetActiveBound();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();
  N_PDOF = B1->GetN_Rows();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();

  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  j = ARowPtr[0];
  // real dof
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      value9 = u1[index];
      value10 = u2[index];
      value11 = u3[index];
      s += value * value9 + value1 * value10 + value2 * value11;
      t += value3* value9 + value4 * value10 + value5 * value11;
      u += value6* value9 + value7 * value10 + value8 * value11;
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
  } // endfor i
  // Dirichlet and hanging nodes

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s += value * u1[index];
      t += value1 * u2[index];
      u += value2 * u3[index];
    }
    v1[i] = s;
    v2[i] = t;
    v3[i] = u;
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
      value3 = B3Entries[j];
      s += value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    q[i] = s;
  } // endfor i
 
  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    v1[i] += s;
    v2[i] += t;
    v3[i] += u;
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
void MatVect_NSE14(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8],
                 A[9], B[0], B[1], B[2], B[3], B[4], B[5], x, y);
  return;
}
 
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A13,
                   TSquareMatrix *A21, TSquareMatrix *A22, TSquareMatrix *A23,
                   TSquareMatrix *A31, TSquareMatrix *A32, TSquareMatrix *A33,
		   TSquareMatrix *C,
                   TMatrix *B1, TMatrix *B2, TMatrix *B3,
                   TMatrix *B1T, TMatrix *B2T, TMatrix *B3T,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, u, value, value1, value2, value3, value4;
  double value5, value6, value7, value8, value9, value10, value11;
  double *u1, *u2, *u3, *p;
  double *v1, *v2, *v3, *q;
  double *r1, *r2, *r3, *r4;
  int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  int *BTRowPtr, *BTKCol, *CRowPtr, *CKCol;
  double *B1Entries, *B2Entries, *B3Entries;
  double *A11Entries, *A12Entries, *A13Entries;
  double *A21Entries, *A22Entries, *A23Entries;
  double *A31Entries, *A32Entries, *A33Entries;
  double *CEntries;
  double *B1TEntries, *B2TEntries, *B3TEntries;
  int N_Active;

  ARowPtr = A11->GetRowPtr();
  AKCol = A11->GetKCol();
  N_UDOF = A11->GetN_Rows();
  N_Active = A11->GetActiveBound();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A13Entries = A13->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();
  A23Entries = A23->GetEntries();
  A31Entries = A31->GetEntries();
  A32Entries = A32->GetEntries();
  A33Entries = A33->GetEntries();
  CRowPtr = C->GetRowPtr();
  CKCol = C->GetKCol();
  CEntries = C->GetEntries();

  BRowPtr = B1->GetRowPtr();
  BKCol = B1->GetKCol();
  N_PDOF = B1->GetN_Rows();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();
  B3Entries = B3->GetEntries();

  BTRowPtr = B1T->GetRowPtr();
  BTKCol = B1T->GetKCol();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();
  B3TEntries = B3T->GetEntries();
 
  u1 = x;
  u2 = u1+N_UDOF;
  u3 = u2+N_UDOF;
  p  = u3+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  v3 = v2+N_UDOF;
  q  = v3+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;
  r4 = r3+N_UDOF;

  j = ARowPtr[0];

  for(i=0;i<N_Active;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A13Entries[j];
      value3 = A21Entries[j];
      value4 = A22Entries[j];
      value5 = A23Entries[j];
      value6 = A31Entries[j];
      value7 = A32Entries[j];
      value8 = A33Entries[j];
      value9 = u1[index];
      value10 = u2[index];
      value11 = u3[index];
      s -= value * value9 + value1 * value10+ value2 * value11;
      t -= value3* value9 + value4 * value10+ value5 * value11;
      u -= value6* value9 + value7 * value10 + value8 * value11;
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
  } // endfor i

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    u = v3[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      value2 = A33Entries[j];
      s -= value * u1[index];
      t -= value1 * u2[index];
      u -= value2 * u3[index];
    }
    r1[i] = s;
    r2[i] = t;
    r3[i] = u;
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
      value3 = B3Entries[j];
      s -= value1 * u1[index] + value2 * u2[index]
        + value3 * u3[index];
    } // endfor j
    r4[i] = s;
    
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    u = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value3 = B3TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
      u += value3 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
    r3[i] -= u;
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
	  s += value1 * p[index];
      } // endfor j
      q[i] -= s;
  } // endfor i
}
void Defect_NSE14(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int N_UDOF,N_PDOF;

  CoupledDefect(A[0], A[1], A[2], A[3], A[4], A[5], A[6],A[7], A[8],
                A[9], B[0], B[1], B[2], B[3], B[4], B[5], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();

  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(r+3*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  return;
}




#endif 

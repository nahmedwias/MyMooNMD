// =======================================================================
// @(#)DiscreteForm2D.C        1.6 10/18/99
// 
// Class:       TDiscreteForm2D
// Purpose:     assemble a couple of matrices and right-hand side at once
//
// Author:      Gunar Matthies (02.10.98)
//
// History:     start of implementation 02.10.98 (Gunar Matthies)
//        :     added moving domains and 2phaseflows 20.09.09 (Sashikumaar Ganesan)
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <Database.h>
#include <FEDatabase2D.h>
#include <DiscreteForm2D.h>
#include <string.h>
#include <stdlib.h>


/** constructor with vector initialization */
TDiscreteForm2D::TDiscreteForm2D(char *name, char *description,
        int n_terms, MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFct2D *assemble, CoeffFct2D coeffs,
        ManipulateFct2D *manipulate)
{
  int i, j, max;
  MultiIndex2D alpha;

  Name = strdup(name);
  Description = strdup(description);

  N_Terms = n_terms;
  Derivatives = derivatives;
  FESpaceNumber = fespacenumber;

  N_Matrices = n_matrices;
  N_Rhs = n_rhs;
  RowSpace = rowspace;
  ColumnSpace = columnspace;
  RhsSpace = rhsspace;

  Coeffs = coeffs;

  Assemble = assemble;
  AssembleParam = nullptr;

  Manipulate = manipulate;

  AllOrigValues = new double** [N_Terms];
  OrigValues = new double* [N_Terms];

  // find number of spaces
  max = -1;
  for(i=0;i<N_Terms;i++)
  {
    j = FESpaceNumber[i];
    if(j > max) max = j;
  }

  N_Spaces = max+1;

  Needs2ndDerivatives = new bool[N_Spaces];
  for(i=0;i<N_Spaces;i++)
    Needs2ndDerivatives[i] = false;

  for(i=0;i<N_Terms;i++)
  {
    alpha = Derivatives[i];
    j = FESpaceNumber[i];
    if(alpha == D20 || alpha == D11 || alpha == D02)
      Needs2ndDerivatives[j] = true;
  }

  #ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0)
  #endif 
	Output::print("---------------------");
    Output::print("number of spaces: ", N_Spaces );
    for(i=0;i<N_Spaces;i++)
    	Output::print(i," ",Needs2ndDerivatives[i]);
    Output::print("---------------------");
}

/** constructor with assembling using parameters */
TDiscreteForm2D::TDiscreteForm2D(char *name, char *description,
        int n_terms, MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFctParam assembleparam, CoeffFct2D coeffs,
        ManipulateFct2D *manipulate)
{
  int i, j, max;
  MultiIndex2D alpha;

  Name = strdup(name);
  Description = strdup(description);

  N_Terms = n_terms;
  Derivatives = derivatives;
  FESpaceNumber = fespacenumber;

  N_Matrices = n_matrices;
  N_Rhs = n_rhs;
  RowSpace = rowspace;
  ColumnSpace = columnspace;
  RhsSpace = rhsspace;

  Coeffs = coeffs;

  Assemble = nullptr;
  AssembleParam = assembleparam;

  Manipulate = manipulate;

  AllOrigValues = new double** [N_Terms];
  OrigValues = new double* [N_Terms];

  // find number of spaces
  max = -1;
  for(i=0;i<N_Terms;i++)
  {
    j = FESpaceNumber[i];
    if(j > max) max = j;
  }

  N_Spaces = max+1;

  Needs2ndDerivatives = new bool[N_Spaces];
  for(i=0;i<N_Spaces;i++)
    Needs2ndDerivatives[i] = false;

  for(i=0;i<N_Terms;i++)
  {
    alpha = Derivatives[i];
    j = FESpaceNumber[i];
    if(alpha == D20 || alpha == D11 || alpha == D02)
      Needs2ndDerivatives[j] = true;
  }

#ifdef _MPI
int rank;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if(rank==0)
#endif
  Output::print("---------------------");
  Output::print("number of spaces: ", N_Spaces );
  for(i=0;i<N_Spaces;i++)
  	Output::print(i," ",Needs2ndDerivatives[i]);
  Output::print("---------------------");
}

TDiscreteForm2D::~TDiscreteForm2D()
{
  delete AllOrigValues;
  delete OrigValues;
  delete Needs2ndDerivatives;
  delete Name;
  delete Description;
}

void TDiscreteForm2D::GetLocalForms(int N_Points, double *weights, 
                                    double *AbsDetjk, double hK, 
                                    double *X, double *Y,
                                    int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                                    double **Parameters, double **AuxArray,
                                    TBaseCell *Cell, int n_matrices, int n_rhs,
                                    double ***LocMatrix, double **LocRhs,
                                    double factor)
{
  int i,j, N_Rows, N_Columns;
  double **CurrentMatrix, *MatrixRow;
  double Mult, *Coeff, *Param;

  // cout << "in TDiscreteForm2D::GetLocalForms" << endl;

  for(i=0;i<n_matrices;i++)
  {
    CurrentMatrix = LocMatrix[i];
    N_Rows = N_BaseFuncts[RowSpace[i]];
    N_Columns = N_BaseFuncts[ColumnSpace[i]];
    for(j=0;j<N_Rows;j++)
    {
      MatrixRow = CurrentMatrix[j];
      memset(MatrixRow, 0, SizeOfDouble*N_Columns);
    } // endfor j
  } // endfor i

  for(i=0;i<n_rhs;i++)
  {
    N_Rows = N_BaseFuncts[RhsSpace[i]];
    memset(LocRhs[i], 0, SizeOfDouble*N_Rows);
  }

// *****************************************************
// for 2Phase flow problems (Sashikumaar Ganesan)
  AuxArray[0][0] = Cell->GetPhase_ID();
// *****************************************************

  if(Coeffs)
    Coeffs(N_Points, X, Y, Parameters, AuxArray);

  if(Manipulate)
    Manipulate(N_Points, AuxArray, Parameters, Cell);

  for(i=0;i<N_Terms;i++)
  {
    AllOrigValues[i] = 
      TFEDatabase2D::GetOrigElementValues(BaseFuncts[FESpaceNumber[i]], 
                                        Derivatives[i]);
  }

  for(i=0;i<N_Points;i++)
  {
    Mult = weights[i] * AbsDetjk[i] * factor;
    Coeff = AuxArray[i];
    Coeff[19] = AbsDetjk[i];
    
    if(TDatabase::ParamDB->Axial3DAxis==1)
    {
     Coeff[20] = Y[i];  // r in axial3D (X: symmetric) problems (Sashikumaar Ganesan)      
    }
    else
    {
    Coeff[20] = X[i];  // r in axial3D (Y: symmetric) problems (Sashikumaar Ganesan)      
    }

    Param = Parameters[i];

    for(j=0;j<N_Terms;j++)
      OrigValues[j] = AllOrigValues[j][i];

    if(Assemble)
      Assemble(Mult, Coeff, hK, OrigValues, N_BaseFuncts, 
               LocMatrix, LocRhs);

    if(AssembleParam)
      AssembleParam(Mult, Coeff, Param, hK, OrigValues, N_BaseFuncts, 
                    LocMatrix, LocRhs);
  } // endfor i
}

void TDiscreteForm2D::GetLocalForms(int N_Points, double *weights, 
                        double *AbsDetjk, double hK, 
                        double *X, double *Y,
                        int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                        TBaseCell *Cell,
                        double ***LocMatrix, double **LocRhs)
{
  double Mult;
  double *Coefficients[N_Points];
  double *aux = new double [N_Points*20]; // do not change below 20
  for(int j=0;j<N_Points;j++)
    Coefficients[j] = aux + j*20;
  
  if(Coeffs)
    Coeffs(N_Points, X, Y, nullptr, Coefficients);

  if(Manipulate)
    Manipulate(N_Points, Coefficients, nullptr, Cell);
  for(int j=0;j<N_Terms;j++)
  {
    AllOrigValues[j] = 
      TFEDatabase2D::GetOrigElementValues(BaseFuncts[FESpaceNumber[j]], 
                                        Derivatives[j]);
  }
  
  for(int i=0;i<N_Points;i++)
  {
 
    Mult = weights[i]*AbsDetjk[i];
    Coefficients[i][19] = AbsDetjk[i];
    
    for(int j=0;j<N_Terms;j++) {
       OrigValues[j] = AllOrigValues[j][i];
    }

 
    if(Assemble)
      Assemble(Mult, Coefficients[i], hK, OrigValues, N_BaseFuncts, 
               LocMatrix, LocRhs);
    if(AssembleParam)
      AssembleParam(Mult, Coefficients[i], nullptr, hK, OrigValues, N_BaseFuncts, 
                    LocMatrix, LocRhs);
  } // endfor i
  delete [] aux;
}

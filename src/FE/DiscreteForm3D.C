// =======================================================================
// %W% %G%
// 
// Class:       TDiscreteForm3D
// Purpose:     assemble a couple of matrices and right-hand side at once
//
// Author:      Gunar Matthies (02.10.98)
//
// History:     start of implementation 02.10.98 (Gunar Matthies)
//
// =======================================================================

#include <Database.h>
#include <FEDatabase3D.h>
#include <DiscreteForm3D.h>
#include <string.h>
#include <stdlib.h>

/** constructor with vector initialization */
TDiscreteForm3D::TDiscreteForm3D(char *name, char *description, int n_terms,
        MultiIndex3D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFct3D *assemble, const CoeffFct3D& coeffs,
        ManipulateFct3D *manipulate)
{
  int i, j, max;
  MultiIndex3D alpha;

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
      if(alpha == D200 || alpha == D110 || alpha == D101 || 
         alpha == D020 || alpha == D011 || alpha == D002)
        Needs2ndDerivatives[j] = true;
    }

#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==TDatabase::ParamDB->Par_P0)
#endif
    Output::print<3>("---------------------");
    Output::print<3>("number of spaces: ", N_Spaces);
    for(i=0;i<N_Spaces;i++)
    	Output::print<3>(i ," ", Needs2ndDerivatives[i]);
    Output::print<3>("---------------------");
}

/** constructor with assembling using parameters */
TDiscreteForm3D::TDiscreteForm3D(char *name, char *description, int n_terms,
        MultiIndex3D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        const AssembleFctParam& assembleparam, const CoeffFct3D& coeffs,
        ManipulateFct3D *manipulate)
{
  int i, j, max;
  MultiIndex3D alpha;

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
    if(alpha == D200 || alpha == D110 || alpha == D101 || 
       alpha == D020 || alpha == D011 || alpha == D002) 
      Needs2ndDerivatives[j] = true;
  }


#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==TDatabase::ParamDB->Par_P0)
#endif 
  Output::print<3>("---------------------");
  Output::print<3>("number of spaces: ", N_Spaces);
  for(i=0;i<N_Spaces;i++)
    Output::print<3>(i ," ", Needs2ndDerivatives[i]);
  Output::print<3>("---------------------");
}

TDiscreteForm3D::~TDiscreteForm3D()
{
  delete AllOrigValues;
  delete OrigValues;
  delete Needs2ndDerivatives;
  delete Name;
  delete Description;
}

void TDiscreteForm3D::GetLocalForms(int N_Points, double *weights, 
                                    double *AbsDetjk, double hK, 
                                    double *X, double *Y, double *Z,
                                    int *N_BaseFuncts, BaseFunct3D *BaseFuncts, 
                                    double **Parameters, double **AuxArray,
                                    TBaseCell *Cell,
                                    int n_matrices, int n_rhs,
                                    double ***LocMatrix, double **LocRhs)
{
  int i,j,N_Rows, N_Columns;
  double **CurrentMatrix, *MatrixRow;
  double Mult, *Coeff, *Param;

  // cout << "in TDiscreteForm::GetLocalForms" << endl;

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
  AuxArray[0][1] = Cell->GetRegionID();
  AuxArray[0][2] = hK;
  
//   AuxArray[0][2] = Cell->GetVertex(0)-> GetX(); 
//   AuxArray[0][3] = Cell->GetVertex(1)-> GetX();  
//   AuxArray[0][4] = Cell->GetVertex(2)-> GetX();
//   AuxArray[0][5] = Cell->GetVertex(3)-> GetX();   
// *****************************************************

  if(Coeffs)
    Coeffs(N_Points, X, Y, Z, Parameters, AuxArray);

  if(Manipulate)
    Manipulate(N_Points, AuxArray, Parameters, Cell);

  for(i=0;i<N_Terms;i++)
  {
    AllOrigValues[i] = 
      TFEDatabase3D::GetOrigElementValues(BaseFuncts[FESpaceNumber[i]], 
                                          Derivatives[i]);
    if(!(AllOrigValues[i]))
    {
      ErrMsg("second derivatives not yet supported. Exiting");
      exit(1);
    }
  }


  for(i=0;i<N_Points;i++)
  {
    Mult = weights[i]*AbsDetjk[i];
//      if(TDatabase::ParamDB->P14<AbsDetjk[i]) TDatabase::ParamDB->P14=AbsDetjk[i];        
    Coeff = AuxArray[i];
    Coeff[19] = AbsDetjk[i];
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


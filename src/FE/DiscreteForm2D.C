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
#include <MovingNavierStokes.h>

#ifdef __2D__
#include <NSE2D_Param.h>
#include <NSE2D_EquOrd_FixPo.h>
#include <NSE2D_FixPo.h>
#include <NSE2D_Friction_FixPo.h>
#include <NSE2D_FixPoRot.h>
#include <NSE2D_FixPoSkew.h>
#include <NSE2D_Newton.h>
#include <NSE2D_AxialSymm3D_FixPo.h>
#include <TNSE2D_FixPo.h>
#include <TNSE2D_FixPoRot.h>
#include <TNSE2D_FixPo_SSMUM.h>
#include <TNSE2D_Routines.h>
#include <TCD2D.h>

#include <MainUtilities.h>
#include <ConvDiff2D.h>
#endif

/** constructor with vector initialization */
TDiscreteForm2D::TDiscreteForm2D(char *name, char *description,
        int n_terms, MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFct2D *assemble, CoeffFct2D *coeffs,
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
  AssembleParam = NULL;

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
    Needs2ndDerivatives[i] = FALSE;

  for(i=0;i<N_Terms;i++)
  {
    alpha = Derivatives[i];
    j = FESpaceNumber[i];
    if(alpha == D20 || alpha == D11 || alpha == D02)
      Needs2ndDerivatives[j] = TRUE;
  }

  #ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>0)
  #endif 
  {
  cout << "---------------------" << endl;
  cout << "number of spaces: " << N_Spaces << endl;
  for(i=0;i<N_Spaces;i++)
    cout << i << " " << Needs2ndDerivatives[i] << endl;
  cout << "---------------------" << endl;
  }

}

/** constructor with assembling using parameters */
TDiscreteForm2D::TDiscreteForm2D(char *name, char *description,
        int n_terms, MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFctParam2D *assembleparam, CoeffFct2D *coeffs,
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

  Assemble = NULL;
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
    Needs2ndDerivatives[i] = FALSE;

  for(i=0;i<N_Terms;i++)
  {
    alpha = Derivatives[i];
    j = FESpaceNumber[i];
    if(alpha == D20 || alpha == D11 || alpha == D02)
      Needs2ndDerivatives[j] = TRUE;
  }

  #ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0 && TDatabase::ParamDB->SC_VERBOSE>0)
  #endif 
  {    
  cout << "---------------------" << endl;
  cout << "number of spaces: " << N_Spaces << endl;
  for(i=0;i<N_Spaces;i++)
    cout << i << " " << Needs2ndDerivatives[i] << endl;
  cout << "---------------------" << endl;
  }
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
                        TBaseCell *Cell,
                        int n_matrices, int n_rhs,
                        double ***LocMatrix, double **LocRhs)
{
  int i,j,k,l, N_Rows, N_Columns;
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
    Mult = weights[i]*AbsDetjk[i];
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

/******************************************************************************/
//
// Routines for initializing discrete forms
//
/******************************************************************************/

// since this routine will also compiled for 3d:
#ifdef __2D__


void InitializeDiscreteFormsScalar(TDiscreteForm2D *&DiscreteFormMatrixMRhs, TDiscreteForm2D *&DiscreteFormMatrixARhs, TDiscreteForm2D *&DiscreteFormMatrixMRhs_SUPG, TDiscreteForm2D *&DiscreteFormMatrixARhs_SUPG, CoeffFct2D *LinCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char SDFEMString[] = "SDFEM";

  if(TDatabase::ParamDB->Axial3D)
   {
    // discrete form for assembling mass matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixMRhs = new TDiscreteForm2D(GalerkinString, allString, N_Terms_MatrixMRhs, Derivatives_MatrixMRhs,
      SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
      RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
      MatrixMRhsAssemble_Axial3D, LinCoeffs, NULL);

    // discrete form for assembling stiffness matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixARhs = new TDiscreteForm2D
      (GalerkinString, allString, N_Terms_MatrixARhs, Derivatives_MatrixARhs,
      SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs, N_Rhs_MatrixARhs,
      RowSpace_MatrixARhs, ColumnSpace_MatrixARhs, RhsSpace_MatrixARhs,
      MatrixARhsAssemble_Axial3D, LinCoeffs, NULL);

    // discrete form for assembling mass matrix and rhs (Galerkin FEM)
//     DiscreteFormMatrixMRhs_SUPG = NULL;
    DiscreteFormMatrixMRhs_SUPG = new TDiscreteForm2D
      (SDFEMString, allString, N_Terms_MatrixMRhs_SUPG, Derivatives_MatrixMRhs_SUPG,
      SpacesNumbers_MatrixMRhs_SUPG, N_Matrices_MatrixMRhs_SUPG, N_Rhs_MatrixMRhs_SUPG,
      RowSpace_MatrixMRhs_SUPG, ColumnSpace_MatrixMRhs_SUPG, RhsSpace_MatrixMRhs_SUPG,
      MatrixMRhsAssemble_SUPG_Axial3D, LinCoeffs, NULL);
      
      
    // discrete form for assembling stiffness matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixARhs_SUPG = new TDiscreteForm2D
      (SDFEMString, allString, N_Terms_MatricesAKRhs_SUPG, Derivatives_MatricesAKRhs_SUPG,
      SpacesNumbers_MatricesAKRhs_SUPG, N_Matrices_MatricesAKRhs_SUPG, N_Rhs_MatricesAKRhs_SUPG,
      RowSpace_MatricesAKRhs_SUPG, ColumnSpace_MatricesAKRhs_SUPG, RhsSpace_MatricesAKRhs_SUPG,
      MatricesAKRhsAssemble_SUPG_Axial3D, LinCoeffs, NULL);  
 
//     OutPut("InitializeDiscreteFormsScalar !!! Axial3D SUPG not implemented yet !!!!!!!!" << endl);
//     exit(0);
   }
  else
   {
    // discrete form for assembling mass matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixMRhs = new TDiscreteForm2D(GalerkinString, allString, N_Terms_MatrixMRhs, Derivatives_MatrixMRhs,
      SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
      RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
      MatrixMRhsAssemble, LinCoeffs, NULL);

    // discrete form for assembling stiffness matrix and rhs (Galerkin FEM)
    DiscreteFormMatrixARhs = new TDiscreteForm2D
      (GalerkinString, allString, N_Terms_MatrixARhs, Derivatives_MatrixARhs,
      SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs, N_Rhs_MatrixARhs,
      RowSpace_MatrixARhs, ColumnSpace_MatrixARhs, RhsSpace_MatrixARhs,
      MatrixARhsAssemble, LinCoeffs, NULL);

    // discrete form for assembling mass matrix and rhs ( )
    DiscreteFormMatrixMRhs_SUPG = new TDiscreteForm2D
      (SDFEMString, allString, N_Terms_MatrixMRhs_SUPG, Derivatives_MatrixMRhs_SUPG,
      SpacesNumbers_MatrixMRhs_SUPG, N_Matrices_MatrixMRhs_SUPG, N_Rhs_MatrixMRhs_SUPG,
      RowSpace_MatrixMRhs_SUPG, ColumnSpace_MatrixMRhs_SUPG, RhsSpace_MatrixMRhs_SUPG,
      MatrixMRhsAssemble_SUPG, LinCoeffs, NULL);

    // discrete form for assembling stiffness matrix and rhs ( )
    DiscreteFormMatrixARhs_SUPG = new TDiscreteForm2D
      (SDFEMString, allString, N_Terms_MatricesAKRhs_SUPG, Derivatives_MatricesAKRhs_SUPG,
      SpacesNumbers_MatricesAKRhs_SUPG, N_Matrices_MatricesAKRhs_SUPG, N_Rhs_MatricesAKRhs_SUPG,
      RowSpace_MatricesAKRhs_SUPG, ColumnSpace_MatricesAKRhs_SUPG, RhsSpace_MatricesAKRhs_SUPG,
      MatricesAKRhsAssemble_SUPG, LinCoeffs, NULL);  
   }
  
}


void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormSDFEM,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormVMSProjection,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLSDFEM,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLVMSProjection,
  TDiscreteForm2D *&DiscreteFormPressSep,
  TDiscreteForm2D *&DiscreteFormPressSepAuxProb,
  TDiscreteForm2D *&DiscreteFormNSRFBRhs,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";


  DiscreteFormVMSProjection = NULL;
  DiscreteFormNLVMSProjection = NULL;

  DiscreteFormPressSep = new TDiscreteForm2D(GalerkinString, allString,
					     NSPressSepN_Terms, NSPressSepDerivatives, 
					     NSPressSepSpaceNumbers,
					     NSPressSepN_Matrices, NSPressSepN_Rhs, 
					     NSPressSepRowSpace, NSPressSepColumnSpace,
					     NSPressSepRhsSpace, NSPressSep, LinCoeffs, NULL);

  DiscreteFormPressSepAuxProb = new TDiscreteForm2D(GalerkinString, allString,
                                                    NSPressSepAuxProbN_Terms, NSPressSepAuxProbDerivatives, 
                                                    NSPressSepAuxProbSpaceNumbers,
                                                    NSPressSepAuxProbN_Matrices, NSPressSepAuxProbN_Rhs,
                                                    NSPressSepAuxProbRowSpace, NSPressSepAuxProbColumnSpace,
                                                    NSPressSepAuxProbRhsSpace, NSPressSepAuxProb, LinCoeffs, NULL);

  DiscreteFormNSRFBRhs = new TDiscreteForm2D(GalerkinString, allString,
					   NSRFBRhsN_Terms, NSRFBRhsDerivatives, 
					   NSRFBRhsSpaceNumbers,
					   NSRFBRhsN_Matrices, NSRFBRhsN_Rhs, 
					   NSRFBRhsRowSpace, NSRFBRhsColumnSpace,
					   NSRFBRhsRhsSpace, NSRFBRhs, LinCoeffs, NULL);

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
  {
    switch(NSTYPE)
    {
      case 1:
        if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Galerkin, LinCoeffs, NULL);
                //NSType1RhsSpace, NSType1GalerkinFJMT07, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1SDFEM, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Smagorinsky, LinCoeffs, NULL);

         DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  NSType1VMSProjectionN_Terms, NSType1VMSProjectionDerivatives, 
                  NSType1VMSProjectionSpaceNumbers,
                  NSType1VMSProjectionN_Matrices, NSType1N_Rhs, 
                  NSType1VMSProjectionRowSpace, NSType1VMSProjectionColumnSpace,
                  NSType1RhsSpace, NSType1VMSProjection, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkin, LinCoeffs, NULL);
                //NSType1NLRhsSpace, NSType1_2NLGalerkinFJMT07, LinCoeffs, NULL);

          DiscreteFormNLSDFEM = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLSDFEMN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLSDFEMRhsSpace, NSType1NLSDFEM, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLSmagorinsky, LinCoeffs, NULL);

         DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  NSType1_2NLVMSProjectionN_Terms, NSType1_2NLVMSProjectionDerivatives, 
                  NSType1_2NLVMSProjectionSpaceNumbers,
                  NSType1_2NLVMSProjectionN_Matrices, NSType1NLN_Rhs,
                  NSType1_2NLVMSProjectionRowSpace, NSType1_2NLVMSProjectionColumnSpace,
                  NSType1NLRhsSpace, NSType1_2NLVMSProjection, LinCoeffs, NULL);
        }
	else
	    // skew-symmetric case
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1GalerkinSkew, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1SDFEMSkew, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1SmagorinskySkew, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkinSkew, LinCoeffs, NULL);

          DiscreteFormNLSDFEM = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLSDFEMN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLSDFEMRhsSpace, NSType1NLSDFEMSkew, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLSmagorinskySkew, LinCoeffs, NULL);
        }
        break;

      case 2:
        if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Galerkin, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEM, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkin, LinCoeffs, NULL);
      
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEM, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinsky, LinCoeffs, NULL);
        }
        else 
        { // skew symmetric case
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2GalerkinSkew, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEMSkew, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2SmagorinskySkew, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkinSkew, LinCoeffs, NULL);
      
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEMSkew, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinskySkew, LinCoeffs, NULL);
        }
        break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Galerkin, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Smagorinsky, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkin, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky, LinCoeffs, NULL);
          }
	if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
	{ // skew symmetric case
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskySkew, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew, LinCoeffs, NULL);
          }
	if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
	{ // rotation form
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinRot, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyRot, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinRot, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyRot, LinCoeffs, NULL);
	}
	}
	else
        {
          // (D(u):D(v))
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
          { // skew symmetric case
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskySkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkewDD, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
          { // rotation form
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinRotDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyRotDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinRotDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyRotDD, LinCoeffs, NULL);
	  }
	}
        break;

      case 4:
	case 14:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
	    if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Galerkin, LinCoeffs, NULL);
	    if (TDatabase::ParamDB->NSTYPE!=14)
		DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
							NSType4SDN_Terms, NSType4SDDerivatives, 
							NSType4SDSpaceNumbers,
							NSType4SDN_Matrices, NSType4SDN_Rhs, 
							NSType4SDRowSpace, NSType4SDColumnSpace,
							NSType4SDRhsSpace, NSType4SDFEM, LinCoeffs, NULL);
	    else
		DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
							NSType4EquOrdN_Terms, NSType4EquOrdDerivatives, 
							NSType4EquOrdSpaceNumbers,
							NSType4EquOrdN_Matrices, NSType4EquOrdN_Rhs, 
							NSType4EquOrdRowSpace, NSType4EquOrdColumnSpace,
							NSType4EquOrdRhsSpace, NSType4SDFEMEquOrd, LinCoeffs, NULL);
	      
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Smagorinsky, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkin, LinCoeffs, NULL);
  
	    if (TDatabase::ParamDB->NSTYPE!=14)
		DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
							  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
							  NSType4NLSDSpaceNumbers,
							  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
							  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
							  NSType4NLSDRhsSpace, NSType4NLSDFEM, LinCoeffs, NULL);
	    else
		DiscreteFormNLSDFEM = DiscreteFormSDFEM;
/*
	    new TDiscreteForm2D(SDFEMString, nonlinearString,
							  NSType4NLEquOrdN_Terms, NSType4NLEquOrdDerivatives, 
							  NSType4NLEquOrdSpaceNumbers,
							  NSType4NLEquOrdN_Matrices, NSType4NLEquOrdN_Rhs, 
							  NSType4NLEquOrdRowSpace, NSType4NLEquOrdColumnSpace,
							  NSType4NLEquOrdRhsSpace, NSType4NLSDFEMEquOrd, LinCoeffs, NULL);
*/
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky, LinCoeffs, NULL);

	    // not yet implemented 
	    if (TDatabase::ParamDB->NSTYPE==14)
	    {
		DiscreteFormGalerkin = NULL;
		DiscreteFormUpwind = NULL;
		DiscreteFormSmagorinsky = NULL;
		DiscreteFormVMSProjection = NULL;
		DiscreteFormNLGalerkin = NULL;
		DiscreteFormNLUpwind = NULL;
		DiscreteFormNLSmagorinsky = NULL;
		DiscreteFormNLVMSProjection = NULL;		
	    }

          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMSkew, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskySkew, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMSkew, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinRot, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMRot, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyRot, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinRot, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMRot, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyRot, LinCoeffs, NULL);
          }
        }
        else
        {
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDD, LinCoeffs, NULL);

            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
          {
            
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMSkewDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskySkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMSkewDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkewDD, LinCoeffs, NULL);
          }
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
          {
            
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinRotDD, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMRotDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyRotDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinRotDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMRotDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyRotDD, LinCoeffs, NULL);
          }
        }
        break;
    } // endswitch
  } 
  else //Newton iteration
  {
    switch(NSTYPE)
    {
      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
        // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyNewton, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLGalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLUpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3NLSmagorinskyNewton, LinCoeffs, NULL);
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDDNewton, LinCoeffs, NULL);
  

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLGalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLUpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3NLSmagorinskyDDNewton, LinCoeffs, NULL);
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyNewton, LinCoeffs, NULL);
 
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLGalerkinNewton, LinCoeffs, NULL);

          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDNewtonN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLUpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType4NLSmagorinskyNewton, LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDDNewton, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLGalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDNewtonN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDDNewton, LinCoeffs, NULL);
    
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLUpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType4NLSmagorinskyDDNewton, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  }
}
void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLGalerkinDuese,
  CoeffFct2D *LinCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char nonlinearString[] = "nonlinear";
  char UpwindString[] = "Upwind";

  DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                                             NSType4N_Terms, NSType4Derivatives, 
                                             NSType4SpaceNumbers,
                                             NSType4N_Matrices, NSType4N_Rhs, 
                                             NSType4RowSpace, NSType4ColumnSpace,
                                             NSType4RhsSpace, NSType4GalerkinAxialSymm3D, LinCoeffs, NULL);

  DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                                               NSType4NLN_Terms, NSType4NLDerivatives, 
                                               NSType4NLSpaceNumbers,
                                               NSType4NLN_Matrices, NSType4NLN_Rhs, 
                                               NSType4NLRowSpace, NSType4NLColumnSpace,
                                               NSType4NLRhsSpace, NSType3_4NLGalerkinAxialSymm3D, LinCoeffs, NULL);

  DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
					   NSType4N_Terms, NSType4Derivatives, 
					   NSType4SpaceNumbers,
					   NSType4N_Matrices, NSType4N_Rhs, 
					   NSType4RowSpace, NSType4ColumnSpace,
					   NSType4RhsSpace, NSType4UpwindAxialSymm3D, LinCoeffs, NULL);
  
  DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
					     NSType4NLN_Terms, NSType4NLDerivatives, 
					     NSType4NLSpaceNumbers,
					     NSType4NLN_Matrices, NSType4NLN_Rhs, 
					     NSType4NLRowSpace, NSType4NLColumnSpace,
					     NSType4NLRhsSpace, NSType3_4NLUpwindAxialSymm3D, LinCoeffs, NULL);

  DiscreteFormNLGalerkinDuese = new TDiscreteForm2D(GalerkinString, nonlinearString,
                                               NSType4N_Terms,NSType4Derivatives,
					       NSType4SpaceNumbers,
                                               NSType4NLSDNewtonN_Matrices, NSType4NLN_Rhs, 
					       NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                                               NSType4NLRhsSpace, NSType3_4NLGalerkinAxialSymm3D_Duese, LinCoeffs, NULL);
}

// friction
void InitializeDiscreteFormsFriction(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormSDFEM,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLSDFEM,
  TDiscreteForm2D *&DiscreteFormNLUpwind,
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormGalerkinFriction, 
  TDiscreteForm2D *&DiscreteFormNLGalerkinFriction,
  TDiscreteForm2D *&DiscreteFormGalerkinFrictionLocal, 
  TDiscreteForm2D *&DiscreteFormNLGalerkinFrictionLocal,
  TDiscreteForm2D *&DiscreteFormSDFEMFriction,
  TDiscreteForm2D *&DiscreteFormNLSDFEMFriction,
  TDiscreteForm2D *&DiscreteFormSDFEMDDFriction,
  TDiscreteForm2D *&DiscreteFormNLSDFEMDDFriction,
  TDiscreteForm2D *&DiscreteFormSDFEMDDFrictionRhs,
  TDiscreteForm2D *&DiscreteFormNLSDFEMDDFrictionRhs,
  TDiscreteForm2D *&DiscreteFormSDFEMFrictionRST,
  TDiscreteForm2D *&DiscreteFormNLSDFEMFrictionRST,
  TDiscreteForm2D *&DiscreteFormUpwindFriction,
  TDiscreteForm2D *&DiscreteFormNLUpwindFriction,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";


  DiscreteFormUpwindFriction = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives,
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindFriction, LinCoeffs, NULL);

  DiscreteFormNLUpwindFriction = new TDiscreteForm2D(UpwindString, 
                  nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindFriction, LinCoeffs, NULL);
  
  DiscreteFormGalerkinFrictionLocal 
           = new TDiscreteForm2D(GalerkinString, 
                  allString,
		  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDDFrictionLocal, LinCoeffs, NULL);
  
  DiscreteFormNLGalerkinFrictionLocal 
           = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDDFrictionLocal, 
                  LinCoeffs, NULL);


  DiscreteFormGalerkinFriction 
           = new TDiscreteForm2D(GalerkinString, 
                  allString,
		  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDDFriction, LinCoeffs, NULL);
  
  DiscreteFormNLGalerkinFriction 
           = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDDFriction, 
                  LinCoeffs, NULL);

  DiscreteFormSDFEMFriction = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMFriction, LinCoeffs, NULL);

  DiscreteFormNLSDFEMFriction = new TDiscreteForm2D(SDFEMString, 
		  nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMFriction, LinCoeffs, 
                  NULL);
  DiscreteFormSDFEMDDFriction = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDDFriction, LinCoeffs, NULL);

  DiscreteFormNLSDFEMDDFriction = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDDFriction, LinCoeffs, NULL);
    
  DiscreteFormSDFEMDDFrictionRhs = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDDFrictionRhs, LinCoeffs, NULL);

  DiscreteFormNLSDFEMDDFrictionRhs = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDDFrictionRhs, LinCoeffs, NULL);


  DiscreteFormSDFEMFrictionRST = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMFrictionRST, LinCoeffs, NULL);

  DiscreteFormNLSDFEMFrictionRST = new TDiscreteForm2D(SDFEMString, 
		  nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMFrictionRST, LinCoeffs, 
                  NULL);
  
  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
  {
    switch(NSTYPE)
    {
      case 1:
        if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Galerkin, LinCoeffs, NULL);

	  DiscreteFormUpwind= new TDiscreteForm2D(UpwindString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkin, LinCoeffs, NULL);

	  DiscreteFormNLUpwindFriction = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLSmagorinsky, LinCoeffs, NULL);
        }
	// skew--symmetric
        else
	{
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1GalerkinSkew, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType1N_Terms, NSType1Derivatives, 
                NSType1SpaceNumbers,
                NSType1N_Matrices, NSType1N_Rhs, 
                NSType1RowSpace, NSType1ColumnSpace,
                NSType1RhsSpace, NSType1SmagorinskySkew, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLGalerkinSkew, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType1NLN_Terms, NSType1NLDerivatives, 
                NSType1NLSpaceNumbers,
                NSType1NLN_Matrices, NSType1NLN_Rhs, 
                NSType1NLRowSpace, NSType1NLColumnSpace,
                NSType1NLRhsSpace, NSType1_2NLSmagorinskySkew, LinCoeffs, NULL);
        }
        break;

      case 2:
        if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Galerkin, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEM, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkin, LinCoeffs, NULL);
      
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEM, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinsky, LinCoeffs, NULL);
        }
        else 
        { // skew symmetric case
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2GalerkinSkew, LinCoeffs, NULL);

          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                NSType2SDN_Terms, NSType2SDDerivatives, 
                NSType2SDSpaceNumbers,
                NSType2SDN_Matrices, NSType2SDN_Rhs, 
                NSType2SDRowSpace, NSType2SDColumnSpace,
                NSType2SDRhsSpace, NSType2SDFEMSkew, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                NSType2N_Terms, NSType2Derivatives, 
                NSType2SpaceNumbers,
                NSType2N_Matrices, NSType2N_Rhs, 
                NSType2RowSpace, NSType2ColumnSpace,
                NSType2RhsSpace, NSType2SmagorinskySkew, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLGalerkinSkew, LinCoeffs, NULL);
      
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                NSType2NLSDN_Terms, NSType2NLSDDerivatives, 
                NSType2NLSDSpaceNumbers,
                NSType2NLSDN_Matrices, NSType2NLSDN_Rhs, 
                NSType2NLSDRowSpace, NSType2NLSDColumnSpace,
                NSType2NLSDRhsSpace, NSType2NLSDFEMSkew, LinCoeffs, NULL);

          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                NSType2NLN_Terms, NSType2NLDerivatives, 
                NSType2NLSpaceNumbers,
                NSType2NLN_Matrices, NSType2NLN_Rhs, 
                NSType2NLRowSpace, NSType2NLColumnSpace,
                NSType2NLRhsSpace, NSType1_2NLSmagorinskySkew, LinCoeffs, NULL);
        }
        break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Galerkin, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Smagorinsky, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkin, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky, LinCoeffs, NULL);
          }
          else
          { // skew symmetric case
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskySkew, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew, LinCoeffs, NULL);
          }
        }
        else
        {
          // (D(u):D(v))
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          }
          else
          { // skew symmetric case
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskySkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLGalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLN_Matrices, NSType3NLN_Rhs, 
                  NSType3NLRowSpace, NSType3NLColumnSpace,
                  NSType3NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkewDD, LinCoeffs, NULL);
          }
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Galerkin, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEM, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Smagorinsky, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkin, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEM, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinsky, LinCoeffs, NULL);
          }
          else
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMSkew, LinCoeffs, NULL);
  
            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4Upwind, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskySkew, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinSkew, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMSkew, LinCoeffs, NULL);
  
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwind, LinCoeffs, NULL);

            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkew, LinCoeffs, NULL);
          }
        }
        else
        {
          if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
          {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          }
          else
          {
            
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMSkewDD, LinCoeffs, NULL);

            DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDD, LinCoeffs, NULL);
  
            DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskySkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLGalerkinSkewDD, LinCoeffs, NULL);
  
            DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDRowSpace, NSType4NLSDColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMSkewDD, LinCoeffs, NULL);
    
            DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType3_4NLUpwindDD, LinCoeffs, NULL);
  
            DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLSmagorinskyN_Terms, NSType4NLSmagorinskyDerivatives, 
                  NSType4NLSmagorinskySpaceNumbers,
                  NSType4NLSmagorinskyN_Matrices, NSType4NLSmagorinskyN_Rhs, 
                  NSType4NLSmagorinskyRowSpace, NSType4NLSmagorinskyColumnSpace,
                  NSType4NLSmagorinskyRhsSpace, NSType3_4NLSmagorinskySkewDD, LinCoeffs, NULL);
          }
        }
        break;
    } // endswitch
  } 
  else //Newton iteration
  {
    switch(NSTYPE)
    {
      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
        // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyNewton, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLGalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLUpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3NLSmagorinskyNewton, LinCoeffs, NULL);
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3GalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3UpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString, 
                  NSType3N_Terms, NSType3Derivatives, 
                  NSType3SpaceNumbers,
                  NSType3N_Matrices, NSType3N_Rhs, 
                  NSType3RowSpace, NSType3ColumnSpace,
                  NSType3RhsSpace, NSType3SmagorinskyDDNewton, LinCoeffs, NULL);
  

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLGalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType3NLN_Terms, NSType3NLDerivatives, 
                  NSType3NLSpaceNumbers,
                  NSType3NLNewtonN_Matrices, NSType3NLNewtonN_Rhs, 
                  NSType3NLNewtonRowSpace, NSType3NLNewtonColumnSpace,
                  NSType3NLNewtonRhsSpace, NSType3NLUpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString, 
                  NSType3NLSmagorinskyN_Terms, NSType3NLSmagorinskyDerivatives, 
                  NSType3NLSmagorinskySpaceNumbers,
                  NSType3NLSmagorinskyN_Matrices, NSType3NLSmagorinskyN_Rhs, 
                  NSType3NLSmagorinskyRowSpace, NSType3NLSmagorinskyColumnSpace,
                  NSType3NLSmagorinskyRhsSpace, NSType3NLSmagorinskyDDNewton, LinCoeffs, NULL);
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinNewton, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMNewton, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyNewton, LinCoeffs, NULL);
 
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLGalerkinNewton, LinCoeffs, NULL);

          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDNewtonN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMNewton, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLUpwindNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType4NLSmagorinskyNewton, LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4GalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSDFEM = new TDiscreteForm2D(SDFEMString, allString,
                  NSType4SDN_Terms, NSType4SDDerivatives, 
                  NSType4SDSpaceNumbers,
                  NSType4SDN_Matrices, NSType4SDN_Rhs, 
                  NSType4SDRowSpace, NSType4SDColumnSpace,
                  NSType4SDRhsSpace, NSType4SDFEMDDNewton, LinCoeffs, NULL);

          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4UpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(UpwindString, allString, 
                  NSType4N_Terms, NSType4Derivatives, 
                  NSType4SpaceNumbers,
                  NSType4N_Matrices, NSType4N_Rhs, 
                  NSType4RowSpace, NSType4ColumnSpace,
                  NSType4RhsSpace, NSType4SmagorinskyDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLGalerkinDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSDFEM = new TDiscreteForm2D(SDFEMString, nonlinearString,
                  NSType4NLSDN_Terms, NSType4NLSDDerivatives, 
                  NSType4NLSDSpaceNumbers,
                  NSType4NLSDNewtonN_Matrices, NSType4NLSDN_Rhs, 
                  NSType4NLSDNewtonRowSpace, NSType4NLSDNewtonColumnSpace,
                  NSType4NLSDRhsSpace, NSType4NLSDFEMDDNewton, LinCoeffs, NULL);
    
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLNewtonN_Matrices, NSType4NLNewtonN_Rhs, 
                  NSType4NLNewtonRowSpace, NSType4NLNewtonColumnSpace,
                  NSType4NLNewtonRhsSpace, NSType4NLUpwindDDNewton, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  NSType4NLN_Terms, NSType4NLDerivatives, 
                  NSType4NLSpaceNumbers,
                  NSType4NLN_Matrices, NSType4NLN_Rhs, 
                  NSType4NLRowSpace, NSType4NLColumnSpace,
                  NSType4NLRhsSpace, NSType4NLSmagorinskyDDNewton, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  }
}

void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormColetti,
  TDiscreteForm2D *&DiscreteFormGL00Convolution, 
  TDiscreteForm2D *&DiscreteFormGL00AuxProblem, 
  TDiscreteForm2D *&DiscreteFormVMSProjection,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind, 
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLColetti, 
  TDiscreteForm2D *&DiscreteFormNLGL00Convolution,
  TDiscreteForm2D *&DiscreteFormNLGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormNLVMSProjection,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormRHSColetti,
  TDiscreteForm2D *&DiscreteFormRHSLESModel,
  TDiscreteForm2D *&DiscreteFormMatrixGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormGL00AuxProblemRHS,
  TDiscreteForm2D *&DiscreteFormRHSSmagorinskyExpl,
  TDiscreteForm2D *&DiscreteFormMatrixAuxProblemU,
  TDiscreteForm2D *&DiscreteFormRHSAuxProblemU,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";
  char rhsString[] = "rhs";
  char ColettiString[] = "Coletti";
  char GaldiLaytonString[] = "Galdi-Layton";
  char auxprobString[] = "aux prob";
  char Layton96String[] = "Layton96";

  DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHS, LinCoeffs, NULL);

  DiscreteFormRHSAuxProblemU = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHSAuxProblemU, LinCoeffs, NULL);

  DiscreteFormMatrixAuxProblemU = new TDiscreteForm2D(auxprobString, allString,
            MatrixAuxiliaryProblemN_Terms, MatrixAuxiliaryProblemDerivatives, 
            MatrixAuxiliaryProblemSpaceNumbers,
            MatrixAuxiliaryProblemN_Matrices, MatrixAuxiliaryProblemN_Rhs, 
            MatrixAuxiliaryProblemRowSpace, MatrixAuxiliaryProblemColumnSpace,
            MatrixAuxiliaryProblemRhsSpace, MatrixAuxiliaryProblem, LinCoeffs, NULL);

  DiscreteFormRHSColetti = new TDiscreteForm2D(ColettiString, rhsString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSColetti, LinCoeffs, NULL);


   DiscreteFormRHSLESModel =  new TDiscreteForm2D(GaldiLaytonString, rhsString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSLESModel, LinCoeffs, NULL);

  DiscreteFormGL00AuxProblemRHS = new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSGL00AuxProblemRHSN_Terms, TimeNSGL00AuxProblemRHSDerivatives, 
            TimeNSGL00AuxProblemRHSSpaceNumbers,
            TimeNSGL00AuxProblemRHSN_Matrices, TimeNSGL00AuxProblemRHSN_Rhs, 
            TimeNSGL00AuxProblemRHSRowSpace, TimeNSGL00AuxProblemRHSColumnSpace,
            TimeNSGL00AuxProblemRHSRhsSpace, TimeNSGL00AuxProblemRHS,
            LinCoeffs, NULL);

  // this is temporarily changed to the discrete form for the defect correction 
	    if ( TDatabase::ParamDB->DEFECT_CORRECTION_TYPE==0)
		DiscreteFormRHSSmagorinskyExpl= new TDiscreteForm2D(GaldiLaytonString, auxprobString,
								    TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, 
								    TimeNSRHSColSpaceNumbers,
								    TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
								    TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
								    TimeNSRHSColRhsSpace, TimeNSRHSSmagorinskyExplicit,
								    LinCoeffs, NULL);
	    if ( TDatabase::ParamDB->DEFECT_CORRECTION_TYPE==1)
		DiscreteFormRHSSmagorinskyExpl= new TDiscreteForm2D(GaldiLaytonString, auxprobString,
								    TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, 
								    TimeNSRHSColSpaceNumbers,
								    TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
								    TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
								    TimeNSRHSColRhsSpace, TimeNSRHSDefectCorrectionU2,
								    LinCoeffs, NULL);
	    
	    if ( TDatabase::ParamDB->DEFECT_CORRECTION_TYPE==2)
		DiscreteFormRHSSmagorinskyExpl= new TDiscreteForm2D(GaldiLaytonString, auxprobString,
								    TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, 
								    TimeNSRHSColSpaceNumbers,
								    TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
								    TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
								    TimeNSRHSColRhsSpace, TimeNSRHSDefectCorrectionU2_1,
								    LinCoeffs, NULL);
		

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Upwind, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);
     
       DiscreteFormSmagorinsky = new TDiscreteForm2D(Layton96String, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Smagorinsky, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, 
                TimeNSType1SpaceNumbers,
                TimeNSType1GL00AuxProblemN_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1GL00AuxProblemRowSpace, TimeNSType1GL00AuxProblemColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GL00AuxProblem, LinCoeffs, NULL);

        // same as Galerkin
        DiscreteFormColetti = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        // same as Smagorinsky
        DiscreteFormNLColetti =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType2N_Terms, TimeNSType2Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Upwind, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Smagorinsky, LinCoeffs, NULL);


        // same as Smagorinsky
        DiscreteFormColetti = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2GL00AuxProblemN_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2GL00AuxProblemRowSpace, TimeNSType2GL00AuxProblemColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2GL00AuxProblem, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        // same as Galerkin
        DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
       break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind, LinCoeffs, NULL);

          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Smagorinsky, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

	  if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
	      DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
		  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
 
          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblem, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
          
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;

	  if (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2)
	      DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyRot, LinCoeffs, NULL);

        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3UpwindDD, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD, LinCoeffs, NULL);
 
          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblemDD, LinCoeffs, NULL);

         DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3VMSProjectionN_Terms, TimeNSType3VMSProjectionDerivatives, 
                  TimeNSType3VMSProjectionSpaceNumbers,
                  TimeNSType3VMSProjectionN_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3VMSProjectionRowSpace, TimeNSType3VMSProjectionColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3VMSProjectionDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          
          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
 
         DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLVMSProjectionN_Terms, TimeNSType3NLVMSProjectionDerivatives, 
                  TimeNSType3NLVMSProjectionSpaceNumbers,
                  TimeNSType3NLVMSProjectionN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLVMSProjectionRowSpace, TimeNSType3NLVMSProjectionColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLVMSProjectionDD, LinCoeffs, NULL);
        
       }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Smagorinsky, LinCoeffs, NULL);


          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblem, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        }
        else
        {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4UpwindDD, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4SmagorinskyDD, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;
        
          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblemDD, LinCoeffs, NULL);

          DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4VMSProjectionN_Terms, TimeNSType4VMSProjectionDerivatives, 
                  TimeNSType4VMSProjectionSpaceNumbers,
                  TimeNSType4VMSProjectionN_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4VMSProjectionRowSpace, TimeNSType4VMSProjectionColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4VMSProjectionDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
   
          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;

          DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLVMSProjectionN_Terms, TimeNSType3NLVMSProjectionDerivatives, 
                  TimeNSType3NLVMSProjectionSpaceNumbers,
                  TimeNSType3NLVMSProjectionN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLVMSProjectionRowSpace, TimeNSType3NLVMSProjectionColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLVMSProjectionDD, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  else
  {
    OutPut("Newton method not implemented yet !!!!!!!!" << endl);
    exit(4711);
  }
}
void InitializeDiscreteForms(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormColetti,
  TDiscreteForm2D *&DiscreteFormGL00Convolution, 
  TDiscreteForm2D *&DiscreteFormGL00AuxProblem, 
  TDiscreteForm2D *&DiscreteFormVMSProjection,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind, 
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLColetti, 
  TDiscreteForm2D *&DiscreteFormNLGL00Convolution,
  TDiscreteForm2D *&DiscreteFormNLGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormNLVMSProjection,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormRHSColetti,
  TDiscreteForm2D *&DiscreteFormRHSLESModel,
  TDiscreteForm2D *&DiscreteFormMatrixGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormGL00AuxProblemRHS,
  TDiscreteForm2D *&DiscreteFormRHSSmagorinskyExpl,
  TDiscreteForm2D *&DiscreteFormMatrixAuxProblemU,
  TDiscreteForm2D *&DiscreteFormRHSAuxProblemU,
  TDiscreteForm2D *&DiscreteFormC,
  TDiscreteForm2D *&DiscreteFormJ,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";
  char rhsString[] = "rhs";
  char ColettiString[] = "Coletti";
  char GaldiLaytonString[] = "Galdi-Layton";
  char auxprobString[] = "aux prob";
  char Layton96String[] = "Layton96";

  DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHS, LinCoeffs, NULL);

  DiscreteFormRHSAuxProblemU = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHSAuxProblemU, LinCoeffs, NULL);

  DiscreteFormMatrixAuxProblemU = new TDiscreteForm2D(auxprobString, allString,
            MatrixAuxiliaryProblemN_Terms, MatrixAuxiliaryProblemDerivatives, 
            MatrixAuxiliaryProblemSpaceNumbers,
            MatrixAuxiliaryProblemN_Matrices, MatrixAuxiliaryProblemN_Rhs, 
            MatrixAuxiliaryProblemRowSpace, MatrixAuxiliaryProblemColumnSpace,
            MatrixAuxiliaryProblemRhsSpace, MatrixAuxiliaryProblem, LinCoeffs, NULL);

  DiscreteFormRHSColetti = new TDiscreteForm2D(ColettiString, rhsString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSColetti, LinCoeffs, NULL);


   DiscreteFormRHSLESModel =  new TDiscreteForm2D(GaldiLaytonString, rhsString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSLESModel, LinCoeffs, NULL);

  DiscreteFormGL00AuxProblemRHS = new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSGL00AuxProblemRHSN_Terms, TimeNSGL00AuxProblemRHSDerivatives, 
            TimeNSGL00AuxProblemRHSSpaceNumbers,
            TimeNSGL00AuxProblemRHSN_Matrices, TimeNSGL00AuxProblemRHSN_Rhs, 
            TimeNSGL00AuxProblemRHSRowSpace, TimeNSGL00AuxProblemRHSColumnSpace,
            TimeNSGL00AuxProblemRHSRhsSpace, TimeNSGL00AuxProblemRHS,
            LinCoeffs, NULL);

  DiscreteFormRHSSmagorinskyExpl= new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, 
            TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSSmagorinskyExplicit,
            LinCoeffs, NULL);

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        cout << "Super!!!" << endl;
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Upwind, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);
     
       DiscreteFormSmagorinsky = new TDiscreteForm2D(Layton96String, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Smagorinsky, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, 
                TimeNSType1SpaceNumbers,
                TimeNSType1GL00AuxProblemN_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1GL00AuxProblemRowSpace, TimeNSType1GL00AuxProblemColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GL00AuxProblem, LinCoeffs, NULL);

        // same as Galerkin
        DiscreteFormColetti = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        // same as Smagorinsky
        DiscreteFormNLColetti =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution =  DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        
        // Discrete Forms for Rosenbrock Methods
        DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                0, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GalerkinRHS, LinCoeffs, NULL);
                
        DiscreteFormC = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                0, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GalerkinC, LinCoeffs, NULL);
        
        DiscreteFormJ = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1GalerkinJ, LinCoeffs, NULL);
        
        break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType2N_Terms, TimeNSType2Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Upwind, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Smagorinsky, LinCoeffs, NULL);


        // same as Smagorinsky
        DiscreteFormColetti = DiscreteFormSmagorinsky;
        DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

        DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2GL00AuxProblemN_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2GL00AuxProblemRowSpace, TimeNSType2GL00AuxProblemColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2GL00AuxProblem, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        // same as Galerkin
        DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
        DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        
        // Discrete Forms for Rosenbrock Methods
        DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                0, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType1GalerkinRHS, LinCoeffs, NULL);
                
        DiscreteFormC = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                0, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType1GalerkinC, LinCoeffs, NULL);
        
        DiscreteFormJ = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType1GalerkinJ, LinCoeffs, NULL);
        
       break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Smagorinsky, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblem, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
          
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        }
        else
        {
          // (D(u):D(v))
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3UpwindDD, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD, LinCoeffs, NULL);
 
          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, 
                TimeNSType3SpaceNumbers,
                TimeNSType3GL00AuxProblemN_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3GL00AuxProblemRowSpace, TimeNSType3GL00AuxProblemColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GL00AuxProblemDD, LinCoeffs, NULL);

         DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3VMSProjectionN_Terms, TimeNSType3VMSProjectionDerivatives, 
                  TimeNSType3VMSProjectionSpaceNumbers,
                  TimeNSType3VMSProjectionN_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3VMSProjectionRowSpace, TimeNSType3VMSProjectionColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3VMSProjectionDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);
  
          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
          
          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
 
         DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLVMSProjectionN_Terms, TimeNSType3NLVMSProjectionDerivatives, 
                  TimeNSType3NLVMSProjectionSpaceNumbers,
                  TimeNSType3NLVMSProjectionN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLVMSProjectionRowSpace, TimeNSType3NLVMSProjectionColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLVMSProjectionDD, LinCoeffs, NULL);
        
       }
       
         // Discrete Forms for Rosenbrock Methods
        DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                0, TimeNSType3N_Rhs, 
                TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                TimeNSType3RhsSpace, TimeNSType1GalerkinRHS, LinCoeffs, NULL);
                
        DiscreteFormC = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                0, TimeNSType3N_Rhs, 
                TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                TimeNSType3RhsSpace, TimeNSType1GalerkinC, LinCoeffs, NULL);
        
        DiscreteFormJ = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                TimeNSType3RhsSpace, TimeNSType3GalerkinJ, LinCoeffs, NULL);
       
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Smagorinsky, LinCoeffs, NULL);


          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;

          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblem, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;
        }
        else
        {
            DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4UpwindDD, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4SmagorinskyDD, LinCoeffs, NULL);

          // same as Smagorinsky
          DiscreteFormColetti = DiscreteFormSmagorinsky;
          DiscreteFormGL00Convolution = DiscreteFormSmagorinsky;
        
          DiscreteFormGL00AuxProblem =  new TDiscreteForm2D(GaldiLaytonString, allString,
                TimeNSType4N_Terms, TimeNSType4Derivatives, 
                TimeNSType4SpaceNumbers,
                TimeNSType4GL00AuxProblemN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4GL00AuxProblemRowSpace, TimeNSType4GL00AuxProblemColumnSpace,
                TimeNSType4RhsSpace, TimeNSType4GL00AuxProblemDD, LinCoeffs, NULL);

          DiscreteFormVMSProjection = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4VMSProjectionN_Terms, TimeNSType4VMSProjectionDerivatives, 
                  TimeNSType4VMSProjectionSpaceNumbers,
                  TimeNSType4VMSProjectionN_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4VMSProjectionRowSpace, TimeNSType4VMSProjectionColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4VMSProjectionDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
   
          // same as Smagorinsky
          DiscreteFormNLColetti = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00Convolution = DiscreteFormNLSmagorinsky;
          DiscreteFormNLGL00AuxProblem =  DiscreteFormNLSmagorinsky;

          DiscreteFormNLVMSProjection = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLVMSProjectionN_Terms, TimeNSType3NLVMSProjectionDerivatives, 
                  TimeNSType3NLVMSProjectionSpaceNumbers,
                  TimeNSType3NLVMSProjectionN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLVMSProjectionRowSpace, TimeNSType3NLVMSProjectionColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLVMSProjectionDD, LinCoeffs, NULL);
        }
        
        // Discrete Forms for Rosenbrock Methods
        DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                0, TimeNSType4N_Rhs, 
                TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                TimeNSType4RhsSpace, TimeNSType1GalerkinRHS, LinCoeffs, NULL);
                
        DiscreteFormC = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                0, TimeNSType4N_Rhs, 
                TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                TimeNSType4RhsSpace, TimeNSType1GalerkinC, LinCoeffs, NULL);
        
        DiscreteFormJ = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                TimeNSType4NLN_Matrices, TimeNSType4N_Rhs, 
                TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                TimeNSType4RhsSpace, TimeNSType3GalerkinJ, LinCoeffs, NULL); 
        
        break;
    } // endswitch
  else
  {
    OutPut("Newton method not implemented yet !!!!!!!!" << endl);
    exit(4711);
  }
}
void InitializeDiscreteFormsPaper2(  
  TDiscreteForm2D *&DiscreteFormRHSGL00AuxProblem,
  TDiscreteForm2D *&DiscreteFormGL00AuxProblemRHSPaper2,
  CoeffFct2D *LinCoeffs, CoeffFct2D *Coeffs)
{
  char GaldiLaytonString[] = "GaldiLayton";
  char auxprobString[] = "aux prob";

  DiscreteFormRHSGL00AuxProblem = new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSRHSColN_Terms, TimeNSRHSColDerivatives, TimeNSRHSColSpaceNumbers,
            TimeNSRHSColN_Matrices, TimeNSRHSColN_Rhs, 
            TimeNSRHSColRowSpace, TimeNSRHSColColumnSpace,
            TimeNSRHSColRhsSpace, TimeNSRHSGL00AuxProblemPaper2, LinCoeffs, NULL);

  DiscreteFormGL00AuxProblemRHSPaper2 = new TDiscreteForm2D(GaldiLaytonString, auxprobString,
            TimeNSGL00AuxProblemRHSN_Terms, TimeNSGL00AuxProblemRHSDerivatives, 
            TimeNSGL00AuxProblemRHSSpaceNumbers,
            TimeNSGL00AuxProblemRHSN_Matrices, TimeNSGL00AuxProblemRHSN_Rhs, 
            TimeNSGL00AuxProblemRHSRowSpace, TimeNSGL00AuxProblemRHSColumnSpace,
            TimeNSGL00AuxProblemRHSRhsSpace,  TimeNSGL00AuxProblemRHSPaper2 ,
            Coeffs, NULL);
}
void InitializeDiscreteFormsVMS(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormUpwind,
  TDiscreteForm2D *&DiscreteFormSmagorinsky,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormNLUpwind, 
  TDiscreteForm2D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteForm_ho_RHS,
  TDiscreteForm2D *&DiscreteForm_ls_RHS,
  CoeffFct2D *LinCoeffs, int NSTYPE)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char SmagorinskyString[] = "Smagorinsky";
  char nonlinearString[] = "nonlinear";
  char SDFEMString[] = "SDFEM";
  char rhsString[] = "rhs";
  char ColettiString[] = "Coletti";
  char GaldiLaytonString[] = "Galdi-Layton";
  char auxprobString[] = "aux prob";
  char Layton96String[] = "Layton96";

  DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHS, LinCoeffs, NULL);

  DiscreteForm_ho_RHS = new TDiscreteForm2D(
    GalerkinString, rhsString,
    TimeNS_ho_RHSN_Terms, TimeNS_ho_RHSDerivatives, TimeNS_ho_RHSSpaceNumbers,
    TimeNS_ho_RHSN_Matrices, TimeNS_ho_RHSN_Rhs, 
    TimeNS_ho_RHSRowSpace, TimeNS_ho_RHSColumnSpace,
    TimeNS_ho_RHSRhsSpace, TimeNS_VMS_SmallRhs2D, LinCoeffs, NULL);

  // THIS SHOULD BE DELETED IN THE MAIN PROGRAM
  DiscreteForm_ls_RHS = NULL;

  //fixed point iteration
  if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Upwind, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);
     
       DiscreteFormSmagorinsky = new TDiscreteForm2D(Layton96String, allString,
                TimeNSType1N_Terms, TimeNSType1Derivatives, TimeNSType1SpaceNumbers,
                TimeNSType1N_Matrices, TimeNSType1N_Rhs, 
                TimeNSType1RowSpace, TimeNSType1ColumnSpace,
                TimeNSType1RhsSpace, TimeNSType1Smagorinsky, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType1NLN_Terms, TimeNSType1NLDerivatives, TimeNSType1NLSpaceNumbers,
                TimeNSType1NLN_Matrices, TimeNSType1NLN_Rhs, 
                TimeNSType1NLRowSpace, TimeNSType1NLColumnSpace,
                TimeNSType1NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);
        break;

      case 2:
        DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Galerkin, LinCoeffs, NULL);

        DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                TimeNSType2N_Terms, TimeNSType2Derivatives, TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Upwind, LinCoeffs, NULL);

        DiscreteFormSmagorinsky = new TDiscreteForm2D(GalerkinString, allString,
                TimeNSType2N_Terms, TimeNSType2Derivatives, 
                TimeNSType2SpaceNumbers,
                TimeNSType2N_Matrices, TimeNSType2N_Rhs, 
                TimeNSType2RowSpace, TimeNSType2ColumnSpace,
                TimeNSType2RhsSpace, TimeNSType2Smagorinsky, LinCoeffs, NULL);

        DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLGalerkin, LinCoeffs, NULL);

        DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLUpwind, LinCoeffs, NULL);

        DiscreteFormNLSmagorinsky = new TDiscreteForm2D(GalerkinString, nonlinearString,
                TimeNSType2NLN_Terms, TimeNSType2NLDerivatives, TimeNSType2NLSpaceNumbers,
                TimeNSType2NLN_Matrices, TimeNSType2NLN_Rhs, 
                TimeNSType2NLRowSpace, TimeNSType2NLColumnSpace,
                TimeNSType2NLRhsSpace, TimeNSType1_2NLSmagorinsky, LinCoeffs, NULL);

        break;

      case 3:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          // (grad, grad)
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Upwind, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                  TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                  TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                  TimeNSType3RhsSpace, TimeNSType3Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                  TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                  TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                  TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
        }
        else
        {
	    //switch(TDatabase::ParamDB->VMS_TYPE)
	    //{
	    //case 0:
                 // (D(u):D(v))
                 DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                                                            TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                                                            TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                                                            TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                                                            TimeNSType3RhsSpace, TimeNSType3GalerkinDD, LinCoeffs, NULL);
                 
                 DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                                                          TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                                                          TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                                                          TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                                                          TimeNSType3RhsSpace, TimeNSType3UpwindDD, LinCoeffs, NULL);
                 
                 DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                                                               TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                                                               TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                                                               TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                                                               TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                                                              TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                              TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                                                              TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                                                              TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                                                            TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                            TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                                                            TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                                                            TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                                                                 TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                                 TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                                                                 TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                                                                 TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
		 /*     break;
             case 1:
                 // (D(u):D(v))
                 DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                         TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                         TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                         TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                         TimeNSType3RhsSpace, TimeNSType3GalerkinDD, LinCoeffs, NULL);
                 
                 DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                         TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                         TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                         TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                         TimeNSType3RhsSpace, TimeNSType3UpwindDD, LinCoeffs, NULL);
                 
                 DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                         TimeNSType3N_Terms, TimeNSType3Derivatives, TimeNSType3SpaceNumbers,
                         TimeNSType3N_Matrices, TimeNSType3N_Rhs, 
                         TimeNSType3RowSpace, TimeNSType3ColumnSpace,
                         TimeNSType3RhsSpace, TimeNSType3SmagorinskyDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                         TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                         TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                         TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                         TimeNSType3NLRhsSpace, TimeNSType3_4NLGalerkin_VMS_1_DD, LinCoeffs, NULL);
                 
                 DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                                                            TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                            TimeNSType3NLN_Matrices, TimeNSType3NLN_Rhs, 
                                                            TimeNSType3NLRowSpace, TimeNSType3NLColumnSpace,
                                                            TimeNSType3NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);
                 
                 DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                                                                 TimeNSType3NLN_Terms, TimeNSType3NLDerivatives, TimeNSType3NLSpaceNumbers,
                                                                 TimeNSType3NLSmagorinskyN_Matrices, TimeNSType3NLN_Rhs, 
                                                                 TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                                                                 TimeNSType3NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
                 break;
		 }*/
        }
        break;

      case 4:
        if(TDatabase::ParamDB->LAPLACETYPE == 0)
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Galerkin, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Upwind, LinCoeffs, NULL);
  
         DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4Smagorinsky, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkin, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwind, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinsky, LinCoeffs, NULL);
        }
        else
        {
          DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4UpwindDD, LinCoeffs, NULL);
  
          DiscreteFormSmagorinsky = new TDiscreteForm2D(SmagorinskyString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4SmagorinskyDD, LinCoeffs, NULL);

          DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);
  
          DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLUpwindDD, LinCoeffs, NULL);

          DiscreteFormNLSmagorinsky = new TDiscreteForm2D(SmagorinskyString, nonlinearString,
                  TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType3NLSmagorinskyN_Matrices, TimeNSType4NLN_Rhs, 
                  TimeNSType3NLSmagorinskyRowSpace, TimeNSType3NLSmagorinskyColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLSmagorinskyDD, LinCoeffs, NULL);
        }
        break;
    } // endswitch
  else
  {
    OutPut("Newton method not implemented yet !!!!!!!!" << endl);
    exit(4711);
  }
}

void InitializeDiscreteForms_SSMUM(  
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormRHS1,
  CoeffFct2D *LinCoeffs)
{
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char nonlinearString[] = "nonlinear";
  char rhsString[] = "rhs";

  DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNSRHS_SSMUM_ALE, LinCoeffs, NULL);

  DiscreteFormRHS1 = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs, 
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, TimeNS_REL_VELO_SSMUM_WITH_ROTFRAME, LinCoeffs, NULL);

  DiscreteFormGalerkin = 
      new TDiscreteForm2D(GalerkinString, allString,
			  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
			  TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
			  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
			  TimeNSType4RhsSpace, TimeNSType4Galerkin_SSMUM_ALE, LinCoeffs, NULL);
  
  DiscreteFormNLGalerkin = 
      new TDiscreteForm2D(GalerkinString, nonlinearString,
			  TimeNSType4_SSMUM_NLN_Terms, TimeNSType4_SSMUM_NLDerivatives, TimeNSType4_SSMUM_NLSpaceNumbers,
			  TimeNSType4_SSMUM_NLN_Matrices, TimeNSType4_SSMUM_NLN_Rhs, 
			  TimeNSType4_SSMUM_NLRowSpace, TimeNSType4_SSMUM_NLColumnSpace,
			  TimeNSType4_SSMUM_NLRhsSpace, TimeNSType3_4NLGalerkin_SSMUM_ALE, LinCoeffs, NULL);
}

void  InitializeDiscreteForms_Moving(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormNLGalerkin,
                                  TDiscreteForm2D *&DiscreteFormGrid, CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs)
{

  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char nonlinearString[] = "nonlinear";
  char GridString[] = "Grid";
  
  DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
                         GridN_Terms, GridDerivatives, GridSpaceNumbers,
                         GridN_Matrices, GridN_Rhs,
                         GridRowSpace, GridColumnSpace,
                         GridRhsSpace, GridAssemble4,
                         GridCoeffs, NULL);


 if(TDatabase::ParamDB->Axial3D)
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                             TimeNSType4N_Matrices, TimeNSType4N_Rhs,
                             TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                             TimeNSType4RhsSpace, TimeNSType4GalerkinDD_Axial3D, LinCoeffs, NULL);

   DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                               TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                               TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs,
                               TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                               TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD_Axial3D, LinCoeffs, NULL);
  } 
 else
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                             TimeNSType4N_Matrices, TimeNSType4N_Rhs, 
                             TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                             TimeNSType4RhsSpace, TimeNSType4GalerkinDD, LinCoeffs, NULL);

   DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                               TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                               TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs, 
                               TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                               TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD, LinCoeffs, NULL);

 }

}


void InitializeDiscreteForms_ScalarMoving(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormGrid,
                                    TDiscreteForm2D *&DiscreteFormSUPG, CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char SUPGString[] = "SUPG";
  char GridString[] = "Grid";
  
    DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
                         GridN_Terms, GridDerivatives, GridSpaceNumbers,
                         GridN_Matrices, GridN_Rhs,
                         GridRowSpace, GridColumnSpace,
                         GridRhsSpace, GridAssemble4,
                         GridCoeffs, NULL);
  
  
  
 if(TDatabase::ParamDB->Axial3D)
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixMARhs, Derivatives_MatrixMARhs, SpacesNumbers_MatrixMARhs,
                             N_Matrices_MatrixMARhs, N_Rhs_MatrixMARhs,
                             RowSpace_MatrixMARhs, ColumnSpace_MatrixMARhs,
                             RhsSpace_MatrixMARhs, MatrixMARhsAssemble_Axial3D, LinCoeffs, NULL);  
   
   
    cout<< "Only DiscreteFormGalerkin for Axial3D is implemented !!!!!!"<<endl;  
    cout<< "Only DiscreteFormGalerkin for Axial3D is implemented InitializeDiscreteForms !!!!!!"<<endl;     
    DiscreteFormSUPG = NULL; 
  } 
 else
 {   
  DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixMARhs, Derivatives_MatrixMARhs, SpacesNumbers_MatrixMARhs,
                             N_Matrices_MatrixMARhs, N_Rhs_MatrixMARhs,
                             RowSpace_MatrixMARhs, ColumnSpace_MatrixMARhs,
                             RhsSpace_MatrixMARhs, MatrixMARhsAssemble, LinCoeffs, NULL);
  
  
  DiscreteFormSUPG = new TDiscreteForm2D(SUPGString, allString,
                             N_Terms_MatrixMARhs_SUPG, Derivatives_MatrixMARhs_SUPG, SpacesNumbers_MatrixMARhs_SUPG,
                             N_Matrices_MatrixMARhs_SUPG, N_Rhs_MatrixMARhs_SUPG,
                             RowSpace_MatrixMARhs_SUPG, ColumnSpace_MatrixMARhs_SUPG,
                             RhsSpace_MatrixMARhs_SUPG, MatrixMARhsAssemble_SUPG, LinCoeffs, NULL);
 }
  
  
} // InitializeDiscreteForms
   
   
   
   


void InitializeDiscreteForms_Moving(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormSUPG, CoeffFct2D *LinCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char SUPGString[] = "SUPG";
  
 if(TDatabase::ParamDB->Axial3D)
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixMARhs, Derivatives_MatrixMARhs, SpacesNumbers_MatrixMARhs,
                             N_Matrices_MatrixMARhs, N_Rhs_MatrixMARhs,
                             RowSpace_MatrixMARhs, ColumnSpace_MatrixMARhs,
                             RhsSpace_MatrixMARhs, MatrixMARhsAssemble_Axial3D, LinCoeffs, NULL);  
   
   
    cout<< "Only DiscreteFormGalerkin for Axial3D is implemented !!!!!!"<<endl;  
    cout<< "Only DiscreteFormGalerkin for Axial3D is implemented InitializeDiscreteForms !!!!!!"<<endl;     
    DiscreteFormSUPG = NULL; 
  } 
 else
 {   
  DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                             N_Terms_MatrixMARhs, Derivatives_MatrixMARhs, SpacesNumbers_MatrixMARhs,
                             N_Matrices_MatrixMARhs, N_Rhs_MatrixMARhs,
                             RowSpace_MatrixMARhs, ColumnSpace_MatrixMARhs,
                             RhsSpace_MatrixMARhs, MatrixMARhsAssemble, LinCoeffs, NULL);
  
  
  DiscreteFormSUPG = new TDiscreteForm2D(SUPGString, allString,
                             N_Terms_MatrixMARhs_SUPG, Derivatives_MatrixMARhs_SUPG, SpacesNumbers_MatrixMARhs_SUPG,
                             N_Matrices_MatrixMARhs_SUPG, N_Rhs_MatrixMARhs_SUPG,
                             RowSpace_MatrixMARhs_SUPG, ColumnSpace_MatrixMARhs_SUPG,
                             RhsSpace_MatrixMARhs_SUPG, MatrixMARhsAssemble_SUPG, LinCoeffs, NULL);
 }
  
  
} // InitializeDiscreteForms
   
   
void InitializeDiscreteForms_Stationary( TDiscreteForm2D *&DiscreteFormUpwind,  TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormSDFEM,
                                TDiscreteForm2D *&DiscreteFormGLS, CoeffFct2D *LinCoeffs)
{
  
  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char SUPGString[] = "SUPG";
  
    ManipulateFct2D *manipulate;
  if (TDatabase::ParamDB->SDFEM_NORM_B==0)
    manipulate = linfb;
  else
    manipulate = ave_l2b_quad_points;
  
 if(TDatabase::ParamDB->Axial3D)
  {
   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                                  N_Terms, Derivatives, SpacesNumbers,
                                  CD_N_Matrices, CD_N_Rhs, CD_RowSpace,
                                  CD_ColumnSpace, CD_RhsSpace, BilinearAssemble_Axial3D,
                                  LinCoeffs, NULL); 
   
  
   cout<< "Only DiscreteFormGalerkin for Axial3D is implemented !!!!!!"<<endl;
   cout<< "Only DiscreteFormGalerkin for Axial3D is implemented !!!!!!"<<endl;   
   cout<< "InitializeDiscreteForms_Stationary"<<endl;  
   DiscreteFormUpwind = NULL;
   DiscreteFormSDFEM = NULL;
   DiscreteFormGLS = NULL;
  } 
 else
 {
   
     DiscreteFormUpwind = new TDiscreteForm2D(GalerkinString, allString,
                                  N_Terms, Derivatives, SpacesNumbers,
                                  CD_N_Matrices, CD_N_Rhs, CD_RowSpace,
                                  CD_ColumnSpace, CD_RhsSpace, BilinearAssemble_UPW2,
                                  LinCoeffs, NULL);  
   
     DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                                  N_Terms, Derivatives, SpacesNumbers,
                                  CD_N_Matrices, CD_N_Rhs, CD_RowSpace,
                                  CD_ColumnSpace, CD_RhsSpace, BilinearAssemble,
                                  LinCoeffs, NULL);
     
     DiscreteFormSDFEM = new TDiscreteForm2D(SUPGString, SUPGString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD, 
                                             CD_N_Matrices, CD_N_Rhs, CD_RowSpace, CD_ColumnSpace, CD_RhsSpace,
                                             BilinearAssemble_SD, LinCoeffs,  manipulate);
     
     DiscreteFormGLS = new TDiscreteForm2D(SUPGString, SUPGString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD, 
                                             CD_N_Matrices, CD_N_Rhs, CD_RowSpace, CD_ColumnSpace, CD_RhsSpace,
                                             BilinearAssemble_GLS, LinCoeffs,  manipulate);     
     
     
 }  
} // InitializeDiscreteForms
     

void  InitializeDiscreteForms(TDiscreteForm2D *&DiscreteFormGalerkin, TDiscreteForm2D *&DiscreteFormNLGalerkin,
                              CoeffFct2D *LinCoeffs)
{

  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char nonlinearString[] = "nonlinear";




//   DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
//                                              NSType4N_Terms, NSType4Derivatives, 
//                                              NSType4SpaceNumbers,
//                                              NSType4N_Matrices, NSType4N_Rhs, 
//                                              NSType4RowSpace, NSType4ColumnSpace,
//                                              NSType4RhsSpace, NSType4GalerkinAxialSymm3D, LinCoeffs, NULL);
// 
//   DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
//                                                NSType4NLN_Terms, NSType4NLDerivatives, 
//                                                NSType4NLSpaceNumbers,
//                                                NSType4NLN_Matrices, NSType4NLN_Rhs, 
//                                                NSType4NLRowSpace, NSType4NLColumnSpace,
//                                                NSType4NLRhsSpace, NSType3_4NLGalerkinAxialSymm3D, LinCoeffs, NULL);
// 
//   DiscreteFormUpwind = new TDiscreteForm2D(UpwindString, allString, 
// 					   NSType4N_Terms, NSType4Derivatives, 
// 					   NSType4SpaceNumbers,
// 					   NSType4N_Matrices, NSType4N_Rhs, 
// 					   NSType4RowSpace, NSType4ColumnSpace,
// 					   NSType4RhsSpace, NSType4UpwindAxialSymm3D, LinCoeffs, NULL);
//   
//   DiscreteFormNLUpwind = new TDiscreteForm2D(UpwindString, nonlinearString, 
// 					     NSType4NLN_Terms, NSType4NLDerivatives, 
// 					     NSType4NLSpaceNumbers,
// 					     NSType4NLN_Matrices, NSType4NLN_Rhs, 
// 					     NSType4NLRowSpace, NSType4NLColumnSpace,
// 					     NSType4NLRhsSpace, NSType3_4NLUpwindAxialSymm3D, LinCoeffs, NULL);
}


void InitializeDiscreteFormGrid(TDiscreteForm2D *&DiscreteFormGrid, CoeffFct2D *GridCoeffs)
{
  char GridString[] = "Grid";
  char allString[] = "all";
  
  DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
            GridN_Terms, GridDerivatives, GridSpaceNumbers,
            GridN_Matrices, GridN_Rhs,
            GridRowSpace, GridColumnSpace,
            GridRhsSpace, GridAssemble4,
             GridCoeffs, NULL);    
    
    
}
    
void InitializeDiscreteForms_2PhaseAxial3D(
  TDiscreteForm2D *&DiscreteFormGalerkin,
  TDiscreteForm2D *&DiscreteFormNLGalerkin,
  TDiscreteForm2D *&DiscreteFormRHS,
  TDiscreteForm2D *&DiscreteFormGrid,
  CoeffFct2D *LinCoeffs, CoeffFct2D *GridCoeffs, int NSTYPE)
{

  char GalerkinString[] = "Galerkin";
  char allString[] = "all";
  char UpwindString[] = "Upwind";
  char nonlinearString[] = "nonlinear";
  char rhsString[] = "rhs";
  char GridString[] = "Grid";

DiscreteFormRHS = new TDiscreteForm2D(GalerkinString, rhsString,
            TimeNSRHSN_Terms, TimeNSRHSDerivatives, TimeNSRHSSpaceNumbers,
            TimeNSRHSN_Matrices, TimeNSRHSN_Rhs,
            TimeNSRHSRowSpace, TimeNSRHSColumnSpace,
            TimeNSRHSRhsSpace, MovingNSRHS, LinCoeffs, NULL);

 DiscreteFormGrid = new TDiscreteForm2D(GridString, allString,
            GridN_Terms, GridDerivatives, GridSpaceNumbers,
            GridN_Matrices, GridN_Rhs,
            GridRowSpace, GridColumnSpace,
            GridRhsSpace, GridAssemble4,
             GridCoeffs, NULL);


  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
    case 2:
    case 3:
    cout<< " DiscreteFormGalerkin_Psep not implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
    cout<< " DiscreteFormUpwind_Psep not implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;


    cout<< " DiscreteForm noting is implemented for 2Phase NSTYPE 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
    exit(-1);
  
    break;

    case 4:
      if(TDatabase::ParamDB->LAPLACETYPE == 0)
      {

         cout<< " DiscreteForm noting is implemented for 2Phase NSTYPE 4 grad:grad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
        exit(-1);
  
      }
      else
      {

       DiscreteFormGalerkin = new TDiscreteForm2D(GalerkinString, allString,
                  TimeNSType4N_Terms, TimeNSType4Derivatives, TimeNSType4SpaceNumbers,
                   TimeNSType4N_Matrices, TimeNSType4N_Rhs,
                  TimeNSType4RowSpace, TimeNSType4ColumnSpace,
                  TimeNSType4RhsSpace, TimeNSType4GalerkinDD_2PhaseAxial3D, LinCoeffs, NULL);

       DiscreteFormNLGalerkin = new TDiscreteForm2D(GalerkinString, nonlinearString,
                   TimeNSType4NLN_Terms, TimeNSType4NLDerivatives, TimeNSType4NLSpaceNumbers,
                  TimeNSType4NLN_Matrices, TimeNSType4NLN_Rhs,
                  TimeNSType4NLRowSpace, TimeNSType4NLColumnSpace,
                  TimeNSType4NLRhsSpace, TimeNSType3_4NLGalerkinDD_2PhaseAxial3D, LinCoeffs, NULL);

      }
    break;
  } // endswitch


}







#endif

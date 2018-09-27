// =======================================================================
// @(#)MultiGridIte.C        1.24 06/27/00
//
// Class:       TMultiGridIte
// Purpose:     iteration methods
//
// Author:      Volker John 24.10.2000
//
// History:     24.10.2000 start of implementation
//
// =======================================================================
#include <ItMethod.h>
#include <MultiGridIte.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** constructor with initialization */
TMultiGridIte::TMultiGridIte(MatVecProc *MatVec, 
                             DefectProc *Defect, 
                             TItMethod *Prec,
                             int n_aux, int n_dof,
                             TNSE_MultiGrid *MG, int zero_start)
  : TItMethod(MatVec, Defect, Prec, n_aux, n_dof)
                                    
{
  // initialize only mg
  // every other thing is initialized in the definition of the
  // multigrid
  mg = MG;
  N_DOF = n_dof;
  N_Aux = 0;
  Zero_Start=zero_start;
}

TMultiGridIte::~TMultiGridIte()
{
}

int TMultiGridIte::Iterate (TSquareMatrix **sqmat,
                            TMatrix **mat, double *sol, 
                            double *rhs)
{
  double res, *mgsol, *mgrhs;
//  int umfpack_flag;
  
  // umfpack_flag = 0 ==> allocieren und LU-Zerlegung
  // umfpack_flag = 1 ==> nur vorw./rueckw.
  // umfpack_flag = 2 ==> vorw./rueckw. und speicher wieder freigeben
  // umfpack_flag = 3 ==> allocieren, LU-Zerl, Freigabe
  //OutPut("substep and Substeps: " << substep << " - " << Substeps << endl);
/*  if (substep==0)
  {
     umfpack_flag = 0;
  }
  else if (substep<=Substeps)
  {
     umfpack_flag = 1;
  }
  else if (Substeps==substep-1)
  {
     umfpack_flag = 2;
  }
  else
  {
     OutPut("This error should never be happen!!!" << endl);
     exit(4711);
  }
  //OutPut("umfpack_flag: " << umfpack_flag << endl);
  //OutPut("TMultiGridIte::IterateNe ... Here we are!!!" << endl);
  */
  // set data for multigrid cycle
  mg->GetLevel(mg->GetN_Levels()-1)->GetSolutionVector(mgsol);
  mg->GetLevel(mg->GetN_Levels()-1)->GetRhsVector(mgrhs);
  // initialize start solution
  if (Zero_Start)
    memset(mgsol, 0, N_DOF*SizeOfDouble);
  else
    memcpy(mgsol, sol, N_DOF*SizeOfDouble);
  memcpy(mgrhs, rhs, N_DOF*SizeOfDouble);
  mg->SetDirichletNodes(mg->GetN_Levels()-1);
  mg->SetRecursion(mg->GetN_Levels()-1);
  // one multigrid cycle
  mg->Cycle(mg->GetN_Levels()-1, res);
  // store solution on rhs
  memcpy(sol, mgsol, N_DOF*SizeOfDouble);
  return(0);
}                        

// =======================================================================
// @(#)ParMGComponents3D.C
//

// Purpose:    components for parallel multigrid in 3d
//
// Author:     Sashikumaar Ganesan (03.06.12)
//
// History:    Start of implementation 03.06.12 (Sashikumaar Ganesan)
//
// ======================================================================= 
#ifdef _MPI

#include "mpi.h"
#include <ParVector3D.h>
#include <ParFECommunicator3D.h>
#include <ParVectorNSE3D.h>
#include <LinAlg.h>
#include <Database.h>

void ParDefect_NSE2(TSquareMatrix **A, TMatrix **B, TParVectorNSE3D *ParSol, TParVectorNSE3D *ParRhs, TParVectorNSE3D *ParDef)
{
  int N_UDOF,N_PDOF;
  double *x, *b, *r;
  
  x = ParSol->GetValues();
  b = ParRhs->GetValues();
  r = ParDef->GetValues();  

  CoupledDefect(A[0], B[0], B[1], B[2], B[3], B[4], B[5], x, b, r);
  N_UDOF = A[0]->GetN_Rows();
  N_PDOF = B[0]->GetN_Rows();

  ParDef->AssembleByADD();
  
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
      ParDef->IntoL20Pressure();
}


#endif
// =======================================================================
// @(#)ParDirectSolver.h
//
// Class:      TParDirectSolver
// Purpose:    Class for interfacing ParMooN with MumpsSolver
//
// Author:     Sashikumaar Ganesan & Abdus Shamim (27.04.15)
//
// History:    Start of implementation 27.04.15 (Sashikumaar Ganesan & Abdus Shamim)
//
// ======================================================================= 
#ifdef _MPI
#  include "mpi.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <LinAlg.h>

#include <ParDirectSolver.h>
#include <Database.h>
#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>
#include <MumpsSolver.h>

TParDirectSolver::TParDirectSolver(TParFECommunicator3D *parcomm,TSquareMatrix3D *mat)
{
  ParComm = parcomm;
  NDof    = ParComm->GetNDof();
  Comm    = TDatabase::ParamDB->Comm;
  N_rhs   = 1;
  Mat     = mat;
  DSType  = TDatabase::ParamDB->DSType;
  
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  if(DSType == 1)
  {
    if(rank == 0)
      printf("ParDirectSolver Type ---(Select _OMPONLY in makefile)---> ParDiso\n");
    MPI_Finalize();
    exit(0);
  }
  else if(DSType == 2)
  {
    if(rank == 0)
      printf("ParDirectSolver Type ------> Mumps\n");
    
    InitMumps(); 
    Mumps = new TMumpsSolver(N_Eqns, N_Nz, I_rn, J_cn, N_rhs);
  }
  else
  {
    printf("Select ParDirectSolver Type\n");
    MPI_Finalize();
    exit(0);
  }
  
}

TParDirectSolver::~TParDirectSolver()
{
  if(DSType == 2)
    Mumps->Clean();

  delete [] MatLoc;
  delete [] OwnRhs;
  delete [] I_rn;
  delete [] J_cn;
  delete [] GlobalRhs;
//   delete Mumps;
}

void TParDirectSolver::InitMumps()
{
  int i,j,k,t;
  int *Master       = ParComm->GetMaster();
  int *local2global = ParComm->Get_Local2Global();
   
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  RowPtr = Mat->GetRowPtr();
  KCol   = Mat->GetKCol();
  
  N_Nz = 0;
  for(i=0;i<NDof;i++)
    if(Master[i] == rank)
      N_Nz += RowPtr[i+1] - RowPtr[i];
  
  N_Master = ParComm->GetN_Master();
  
  I_rn   = new int[N_Nz];
  J_cn   = new int[N_Nz];
  MatLoc = new double[N_Nz];
  
  k = 0;
  for(i=0;i<NDof;i++)
  {
    if(Master[i] == rank)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
      {
	I_rn[k] = local2global[i] + 1;              //fortran format
	J_cn[k] = local2global[KCol[j]] + 1;        //fortran format
	k++;
      }
    }
  }

  OwnRhs = new double[N_Master];
  
  MPI_Allreduce(&N_Master, &N_Eqns, 1, MPI_INT, MPI_SUM, Comm);
  
  GlobalRhsSize = N_Eqns;
  
  if(rank == 0)
    GlobalRhs = new double[GlobalRhsSize];
 
}

void TParDirectSolver::AssembleMatrix(TSquareMatrix3D *matrix)
{ 
  int i,j,k,t;
  int *Master      = ParComm->GetMaster();
  double *EntriesA = matrix->GetEntries();
 
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  k = 0;
  for(i=0;i<NDof;i++)
  {
    if(Master[i] == rank)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
      {
	MatLoc[k] = EntriesA[j];
	k++;
      }
    }
  }

}

void TParDirectSolver::GetRhs(double *Rhs)
{
  int i,j,k,t;
  int *Master = ParComm->GetMaster();
  
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  t = 0;
  for(i=0;i<NDof;i++)
  {
    if(Master[i] == rank)
    {
      OwnRhs[t++] = Rhs[i];
    }
  }
  
  ParComm->GatherToRoot(GlobalRhs, GlobalRhsSize, OwnRhs, N_Master, 0);
    
}

void TParDirectSolver::UpdateSol(double *Sol)
{
  int i,j,k,t;
  int *Master = ParComm->GetMaster();
  
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  ParComm->ScatterFromRoot(GlobalRhs, GlobalRhsSize, OwnRhs, N_Master, 0);
  
  t = 0;
  for(i=0;i<NDof;i++)
  {
    if(Master[i] == rank)
    {
      Sol[i] = OwnRhs[t++];
    }
  }
  ParComm->CommUpdate(Sol);
}

void TParDirectSolver::Solve(double *Sol, double *Rhs, bool Factorize)
{
  double t = MPI_Wtime();
  int i,j;
 
  if(DSType == 2)
  {
    GetRhs(Rhs);
    
    if(Factorize)
    {
      Mumps->FactorizeAndSolve(MatLoc,GlobalRhs);
    }
    else
    {
      Mumps->Solve(MatLoc,GlobalRhs);
    }
    
    UpdateSol(Sol);
  }
  else
  {
    printf("Select ParDirectSolver Type\n");
    MPI_Finalize();
    exit(0);
  }
  
  printf("time taken for solving::%lf\n",MPI_Wtime()-t);

}

#else

//--------------------------------------------------------------------------------------------------------//
//						OMP ONLY
//--------------------------------------------------------------------------------------------------------//


#ifdef _OMPONLY

#include "omp.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <ParDirectSolver.h>
#include <Database.h>
#include <SquareMatrix3D.h>
#include <ParDiso.h>

TParDirectSolver::TParDirectSolver(TSquareMatrix3D *mat)
{
  N_rhs   = 1;
  Mat     = mat;
  DSType  = TDatabase::ParamDB->DSType;

  if(DSType == 1)
  {
      printf("ParDirectSolver Type ------> ParDiso\n");
    
    InitPardiso();
    ParDiso = new TParDiso(N_Eqns, N_Nz, RowPtr, KCol);
  }
  else
  {
    printf("Select ParDirectSolver Type as 1 for ParDiso\n");
    exit(0);
  }
  
}

TParDirectSolver::~TParDirectSolver()
{
//   ParDiso->Clean(Mat_Values, Rhs, Sol);
  delete [] Mat_Values; Mat_Values = NULL;
}

void TParDirectSolver::InitPardiso()
{
  int i,j,k,t;
   
  RowPtr = Mat->GetRowPtr();
  KCol   = Mat->GetKCol();

  N_Nz   = Mat->GetN_Entries();
  
  N_Eqns = Mat->GetN_Rows();
  
  Mat_Values = new double[N_Nz];
}

void TParDirectSolver::AssembleMatrix(TSquareMatrix3D *matrix)
{
  int i;
  double *EntriesA = matrix->GetEntries();

  for(i=0;i<N_Nz;i++)
    Mat_Values[i] = EntriesA[i];
}

void TParDirectSolver::Solve(double *Sol, double *Rhs, bool Factorize)
{
  int i,j;
  double t = GetTime();
  
  if(DSType == 1)
  {
    ParDiso->FactorizeAndSolve(Mat_Values, Rhs, Sol, Factorize);
  }
  else
  {
    printf("Select ParDirectSolver Type\n");
    exit(0);
  }
  
  printf("time taken for solving::%lf\n",GetTime()-t);

   ParDiso->Clean(Mat_Values, Rhs, Sol);
}

#endif

#endif















   
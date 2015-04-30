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
//   Mapper = mapper;
  ParComm = parcomm;
  NDof    = ParComm->GetNDof();
  Comm    = TDatabase::ParamDB->Comm;
  N_rhs   = 1;
  Mat     = mat;
  
  InitMumps();

  //check if the parameters are global or local????
  Mumps = new TMumpsSolver(N_Eqns, N_Nz, I_rn, J_cn, N_rhs);
  
}

TParDirectSolver::~TParDirectSolver()
{
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
  {
    if(Master[i] == rank)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
      {
	N_Nz++;
      }
    }
  }
  
  I_rn   = new int[N_Nz];
  J_cn   = new int[N_Nz];
  
  k = 0;t = 0;
  N_Master = 0;
  for(i=0;i<NDof;i++)
  {
    if(Master[i] == rank)
    {
      N_Master++;
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
      {
	I_rn[k] = local2global[i] + 1;              //fortran format
	J_cn[k] = local2global[KCol[j]] + 1;        //fortran format
// 	cout<<"local2global[i]::"<<local2global[i]<<endl;
// 	cout<<"local2global[KCol[j]]::"<<local2global[KCol[j]]<<endl;
	k++;
      }
    }
  }

  OwnRhs = new double[N_Master];
  
  MPI_Allreduce(&N_Master, &N_Eqns,    1, MPI_INT, MPI_SUM, Comm);
  
  GlobalRhsSize = N_Eqns;
//   MPI_Allreduce(&N_Master, &GlobalRhsSize, 1, MPI_INT, MPI_SUM, Comm);
  
  if(rank == 0)
    GlobalRhs = new double[GlobalRhsSize];
 
//         printf("N_Nz :: %d\t k :: %d\n",N_Nz,k);
}

void TParDirectSolver::AssembleLocMatrix(TSquareMatrix3D *matrix)
{
  int i,j,k;
  int *Master = ParComm->GetMaster();
  
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  MatLoc = new double[N_Nz];

  double *EntriesA = matrix->GetEntries();  
  k =0;
  for(i=0;i<NDof;i++)
  {
    if(Master[i] == rank)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
      {
	MatLoc[k++] = EntriesA[j];
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
  
  GetRhs(Rhs);
  
  if(Factorize)
  {
    Mumps->FactorizeAndSolve(MatLoc,GlobalRhs);
//     int t;
  }
  else
  {
    Mumps->Solve(MatLoc,GlobalRhs);
//     int t;
  }
  
  UpdateSol(Sol);
  
  
  printf("time taken for solving::%lf\n",MPI_Wtime()-t);
}

#endif



















   
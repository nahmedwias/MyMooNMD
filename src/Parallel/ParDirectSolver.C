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

TParDirectSolver::TParDirectSolver(TParFECommunicator3D *parcomm,TParFECommunicator3D *parcomm_p,TSquareMatrix3D **mat,TMatrix3D **matB)
{
  ParComm   = parcomm;
  NDof      = ParComm->GetNDof();
  Comm      = TDatabase::ParamDB->Comm;
  N_rhs     = 1;
  Mat       = mat[0];
  DSType    = TDatabase::ParamDB->DSType;
  
  //for NSTYPE 4
  if(matB != NULL)
  {
    ParComm_P   = parcomm_p;
    NDof_P      = ParComm_P->GetNDof();
    NDof_U      = ParComm->GetNDof();
    NDof        = 3*NDof_U + NDof_P;
    MatB        = matB[0];   
    MatBT       = matB[3];
    
    MatA11 = mat[0];    MatA12 = mat[1];    MatA13 = mat[2];	MatBT1 = matB[3];
    MatA21 = mat[3];    MatA22 = mat[4];    MatA23 = mat[5];	MatBT2 = matB[4];
    MatA31 = mat[6];    MatA32 = mat[7];    MatA33 = mat[8];	MatBT3 = matB[5];
    MatB1  = matB[0];   MatB2  = matB[1];   MatB3  = matB[2];      //[0]
    
  }
  else
  {
    MatB = NULL;
  }
    
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
//     exit(0);
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
  int i,j,k,l,m,t;
  int *Master_P,*local2global_P;
  int *Master         = ParComm->GetMaster();
  int *local2global   = ParComm->Get_Local2Global();

  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  if(MatB == NULL)
  {
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
  else
  {
    Master_P       = ParComm_P->GetMaster();
    local2global_P = ParComm_P->Get_Local2Global();
    
    //compute N_Nz in sqmatrices i.e. A blocks
    RowPtr = Mat->GetRowPtr();
    KCol   = Mat->GetKCol();
    
    N_Nz_U = 0;
    for(i=0;i<NDof_U;i++)
      if(Master[i] == rank)
	N_Nz_U += RowPtr[i+1] - RowPtr[i];

    //compute N_Nz in B blocks
    RowPtr_P = MatB->GetRowPtr();
    KCol_P   = MatB->GetKCol();

    RowPtr_PT = MatBT->GetRowPtr();
    KCol_PT   = MatBT->GetKCol();
    
    N_Nz_P = 0;
    for(i=0;i<NDof_P;i++)
      if(Master_P[i] == rank)
	N_Nz_P += RowPtr_P[i+1] - RowPtr_P[i];  

    N_Nz = 6*N_Nz_P + 9*N_Nz_U; 	//total non zeros in system matrix
      
    N_Master_U  = ParComm->GetN_Master();
    N_Master_P  = ParComm_P->GetN_Master();
    N_Master    = N_Master_P + 3*N_Master_U;
    
    OwnRhs = new double[N_Master];
    
    MPI_Allreduce(&N_Master, &N_Eqns, 1, MPI_INT, MPI_SUM, Comm);
    
    GlobalRhsSize = N_Eqns;
    
    if(rank == 0)
      GlobalRhs = new double[GlobalRhsSize];
    
    I_rn   = new int[N_Nz];
    J_cn   = new int[N_Nz];
    MatLoc = new double[N_Nz];
    
    k = 0;
    for(i=0;i<NDof_U;i++)
    {
      
      if(Master[i] == rank)
      {
	for(m=0;m<3;m++)
	{
	  //Mat A(m,l)
	  for(l=0;l<3;l++)
	  {
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      I_rn[k] = NDof_U*m + local2global[i] + 1;              //fortran format
	      J_cn[k] = NDof_U*l + local2global[KCol[j]] + 1;        //fortran format
	      k++;
	    }//for(j=
	  }//for(l=
	
	  //Mat BT(m)
	  for(j=RowPtr_PT[i];j<RowPtr_PT[i+1];j++)
	  {
	    I_rn[k] = NDof_U*m + local2global_P[i] + 1;                //fortran format
	    J_cn[k] = NDof_U*3 + local2global_P[KCol_PT[j]] + 1;        //fortran format
	    k++;
	  }//for(j=
	  
        }//for(m=
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
    for(i=0;i<NDof_P;i++)
    {
      if(Master_P[i] == rank)
      {
	//Mat B(m)
	for(l=0;l<3;l++)
	{
	  for(j=RowPtr_P[i];j<RowPtr_P[i+1];j++)
	  {
	    I_rn[k] = NDof_U*3 + local2global_P[i] + 1;                //fortran format
	    J_cn[k] = NDof_U*l + local2global_P[KCol_P[j]] + 1;        //fortran format
	    k++;
	  }//for(j=
	}//for(l=
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
  }
}

void TParDirectSolver::AssembleMatrix()
{ 
  int i,j,k,l,m,t;
  int *Master = ParComm->GetMaster();;
  double *EntriesA;
 
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  if(MatB == NULL)
  {
    EntriesA = Mat->GetEntries();
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
  else
  {
    //test
    for(i=0;i<N_Nz;i++)
    {
      if(I_rn[i]==J_cn[i])
	MatLoc[i]=1;
      else
	MatLoc[i]=0; 
    }
    return;
    
    double *Entries[15];
    Entries[0]  = MatA11->GetEntries();
    Entries[1]  = MatA12->GetEntries();
    Entries[2]  = MatA13->GetEntries();
    
    Entries[3]  = MatA21->GetEntries();
    Entries[4]  = MatA22->GetEntries();
    Entries[5]  = MatA23->GetEntries();
    
    Entries[6]  = MatA31->GetEntries();
    Entries[7]  = MatA32->GetEntries();
    Entries[8]  = MatA33->GetEntries();
    
    Entries[9]  = MatB1->GetEntries();
    Entries[10] = MatB2->GetEntries();
    Entries[11] = MatB3->GetEntries();
    
    Entries[12] = MatBT1->GetEntries();
    Entries[13] = MatBT2->GetEntries();
    Entries[14] = MatBT3->GetEntries();
    
    int *Master_P = ParComm_P->GetMaster();
    
    k = 0;
    for(i=0;i<NDof_U;i++)
    {
      if(Master[i] == rank)
      {
	for(m=0;m<3;m++)
	{
	  //Mat A(m,l)
	  for(l=0;l<3;l++)
	  {
	    for(j=RowPtr[i];j<RowPtr[i+1];j++)
	    {
	      MatLoc[k] = Entries[m*3+l][j]; 
	      k++;
	    }//for(j=
	  }//for(l=
	
	  //Mat BT(m)
	  for(j=RowPtr_PT[i];j<RowPtr_PT[i+1];j++)
	  {
	    MatLoc[k] = Entries[12+m][j];
	    k++;
	  }//for(j=
	  
        }//for(m=
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
    for(i=0;i<NDof_P;i++)
    {
      if(Master_P[i] == rank)
      {
	//Mat B(m)
	for(l=0;l<3;l++)
	{
	  for(j=RowPtr_P[i];j<RowPtr_P[i+1];j++)
	  {
	    MatLoc[k] = Entries[9+l][j];
	    k++;
	  }//for(j=
	}//for(l=
      }//if(Master_P[i] == rank)
    }//for(i=0;i<
    
  }//else

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
      AssembleMatrix();
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
  ParDiso->Clean(Mat_Values);
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
 
}

#endif

#endif















   
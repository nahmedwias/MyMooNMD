// =======================================================================
// @(#)Scalar_ParSolver.C
//
// Class:      TScalar_ParSolver 
// Purpose:    Super class for all Scalar_ParSolver
//
// Author:     Sashikumaar Ganesan (22.12.10)
//
// History:    Start of implementation 22.12.10 (Sashikumaar Ganesan)
//        :    Iterative solved added 05.07.2012 (Sashikumaar Ganesan)
// =======================================================================
#ifdef _MPI
// #  include "mpi.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <LinAlg.h>
 
#include <Scalar_ParSolver.h>
#include <Database.h>
#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>
 
 
 /** constructor */
TScalar_ParSolver::TScalar_ParSolver(TParFECommunicator3D *fEComm, 
                                     TSquareStructure3D *squareStructure, int n_rhs,
                                     TParVector3D *parSol, TParVector3D *parRhs)
{
   Comm = TDatabase::ParamDB->Comm;
   FEComm = fEComm;
   Global_DOF = FEComm->GetN_GlobalDegreesOfFreedom();
   SquareStructure = squareStructure;
   N_Rhs = n_rhs;
   MUMPS_Solver = NULL;
   ParSol = parSol;
   ParRhs = parRhs;
   
   //set up the parallel solver
   switch(int(TDatabase::ParamDB->SOLVER_TYPE))
    {
      case 101:
        // MUMPS Parallel solver
       this->SetUpDistMat();
       this->InitMumps();
      break;

      case 201:
        // Parallel Iterativer solver
        this->SetUpIterDistMat(); 
 
      break;
      
      default: 
        cout << "wrong  Parallel solver type !!!!!!!!!!!!!" << endl;
        exit(0);
      break;

     }   
   
  //allocate memory for system matrix
  DistMat_Entries = new double[N_DistMatEntries];
  
//   printf("Mumps called!  SOLVER  \n");
//  MPI_Finalize();
//  exit(0); 
} //TScalar_ParSolver


void TScalar_ParSolver::Solve(TSquareMatrix3D *Mat, bool FACTORIZE)
{
 int i, rank, distbegin, locbegin, len;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
      
 MPI_Comm_rank(Comm, &rank);  
 
   //set up the parallel solver
   switch(int(TDatabase::ParamDB->SOLVER_TYPE))
    {
      case 101:      
        //centralized rhs at host for MUMPS
        
       this->AssembleDistMatrix(Mat);       
        
        if(rank==0)
         {   
          RHS = new double[N_Rhs*N_Eqns];
          // memset(RHS, 0, N_Rhs*N_Eqns*SizeOfDouble); 
         }      

       ParRhs->AssembleAtRootByADD(RHS);
     
         /** call the solver to solve the system */
       if(FACTORIZE)
        { MUMPS_Solver->FactorizeAndSolve(DistMat_Entries, RHS); 
 
       
/*          if(rank==0)
         {   
            OutPut("MUMPs Norm of sol " <<  sqrt(Ddot(N_Eqns,RHS,RHS))  << endl); 
 
          for(i=0; i<N_Eqns; i++)
              if(RHS[i]<0.)
                OutPut(i<<" RHS " <<   RHS[i]  << endl); 

          memset(RHS, 0, N_Rhs*N_Eqns*SizeOfDouble); 
         }     */   
       
	}
       else
        { MUMPS_Solver->Solve(DistMat_Entries, RHS); }

       ParSol->ScatterFromRoot(RHS);
       
        if(rank==0)
          delete [] RHS;

      break;

      
      case 201:
        // Parallel Iterativer solver
        this->AssembleIterDistMatrix(Mat);
      break;
      
      
      
     default: 
      cout << "wrong  Parallel solver type !!!!!!!!!!!!!" << endl;
      MPI_Finalize();
      exit(0);
     break;     
      
    }//  switch(int

//   printf("TScalar_ParSolver \n");  
//   MPI_Finalize();
//   exit(0);
  
}//Solve



void TScalar_ParSolver::AssembleDistMatrix(TSquareMatrix3D *Mat)
{
 int i, j, k, l, m, M, N, S, Master_ID, Neib_ID, dof, distdof, rank, size, begin, end;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *DeptDofs, *N_DeptDofNeibs;
 int MaxN_DeptDofs_All, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocDof;
 int MaxN_MasterDofs_All, MaxN_SlaveDofs_All, *DeptDofMasterID; //, *N_DofRankIndex;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;    
 int nu, *RowPtrA, *KColA, N_Active;
 int **LocalDofOfRecbuf, **LocalDofOfRecbuf_send;
 int *N_SendDofs, *N_RecevDofs, index_dep;
 int **sendbuf_Kcol, **recevbuf_KcolA, *KColBegin, N_UCol, len_A, rowbegin;

 double **sendbuf, **recevbuf, *EntriesA, *row_Entries;
 
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  FEComm->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs,
                     DeptDofs, N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks,
                     IndexOfNeibRank, LocDOFofNeib, NeibLocDof);

  FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  FEComm->GetDeptDofMasterInfo(MaxN_MasterDofs_All, MaxN_SlaveDofs_All, DeptDofMasterID);
  
  nu = Mat->GetN_Rows();
  RowPtrA = Mat->GetRowPtr();
  KColA = Mat->GetKCol();
  N_Active = Mat->GetActiveBound();
  EntriesA =  Mat->GetEntries();  
  
  LocalDofOfRecbuf = new int*[2];
  LocalDofOfRecbuf_send = new int*[2];
  
  if(N_Neibs)
   {
    N_SendDofs =  new int[N_Neibs];
    N_RecevDofs =  new int[N_Neibs];
    memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    sendbuf = new double*[1]; 
    recevbuf = new double*[1];
    sendbuf_Kcol = new int*[1];
    recevbuf_KcolA = new int*[1];

    sendbuf[0] = new double[N_Neibs*MaxN_SlaveDofs_All*MaxN_NzInDepRow_All_A];
    recevbuf[0]  = new double[N_Neibs*MaxN_SlaveDofs_All*MaxN_NzInDepRow_All_A];
    sendbuf_Kcol[0] = new int[N_Neibs*MaxN_SlaveDofs_All*MaxN_NzInDepRow_All_A];
    recevbuf_KcolA[0] = new int[N_Neibs*MaxN_SlaveDofs_All*MaxN_NzInDepRow_All_A];
    
    LocalDofOfRecbuf[0] = new int[N_Neibs*MaxN_SlaveDofs_All];
    LocalDofOfRecbuf[1] = new int[N_Neibs*MaxN_SlaveDofs_All];
    LocalDofOfRecbuf_send[0]  = new int[N_Neibs*MaxN_SlaveDofs_All];  
    LocalDofOfRecbuf_send[1]  = new int[N_Neibs*MaxN_SlaveDofs_All];
            
     for(i=0;i<N_DeptDofs;i++)
     {
      N = N_DeptDofNeibs[i];
      Master_ID = DeptDofMasterID[i];

      if(Master_ID==rank)
       {
        for(j=0;j<N;j++)
         {
          S = i*MaxSubDomainPerDof + j;
          Neib_ID = DeptDofNeibRanks[S];
          k = IndexOfNeibRank[Neib_ID];
          N_RecevDofs[k]++;

         } // for(j=0;j<N;j++)
       }
      else 
       {
        for(j=0;j<N;j++)
         {
          S = i*MaxSubDomainPerDof + j;
          Neib_ID = DeptDofNeibRanks[S];
          if(Neib_ID==Master_ID)
           break;
         }

        dof = DeptDofs[i];
        k = IndexOfNeibRank[Master_ID];
        M = k*MaxN_SlaveDofs_All + N_SendDofs[k];

        begin = RowPtrA[dof];
        end=RowPtrA[dof+1];

        LocalDofOfRecbuf_send[0][M] = NeibLocDof[S];
        LocalDofOfRecbuf_send[1][M] = end - begin;

        M *=MaxN_NzInDepRow_All_A;
        for(l=begin;l<end;l++)
         {
          sendbuf[0][M] =  EntriesA[l];
          sendbuf_Kcol[0][M++] =  GlobalDofOfLocalDof[KColA[l]];
         }

        N_SendDofs[k]++;
       } //else / if(Master_ID==rank)
      } //  for(i=0;i<N_DeptDofs;i++)     
   } //  if(N_Neibs)
// =====================================================================================
//   communicate 
  FEComm->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_SlaveDofs_All, N_SendDofs,
                            LocalDofOfRecbuf,  N_RecevDofs, 2);

  for(i=0; i<N_Neibs; i++)
   {
    N_SendDofs[i] *= MaxN_NzInDepRow_All_A;
    N_RecevDofs[i] *= MaxN_NzInDepRow_All_A;
   }
  FEComm->FECommunicateNeib(sendbuf, MaxN_SlaveDofs_All*MaxN_NzInDepRow_All_A,
                            N_SendDofs, recevbuf, N_RecevDofs, 1);
  FEComm->FECommunicateNeib(sendbuf_Kcol, MaxN_SlaveDofs_All*MaxN_NzInDepRow_All_A,
                            N_SendDofs, recevbuf_KcolA, N_RecevDofs, 1);
  for(i=0; i<N_Neibs; i++)
    N_RecevDofs[i] /= MaxN_NzInDepRow_All_A;
// =====================================================================================

  if(N_Neibs)
   {
    delete [] N_SendDofs;     
    delete [] sendbuf[0];
    delete [] sendbuf_Kcol[0];
    delete [] LocalDofOfRecbuf_send[0];
    delete [] LocalDofOfRecbuf_send[1];

    delete [] sendbuf;
    delete [] sendbuf_Kcol;    
    delete [] LocalDofOfRecbuf_send;
   }//if(N_Neibs)


  // assemble system sidt matrix directly
  memset(DistMat_Entries, 0, N_DistMatEntries*SizeOfDouble);

  // first add entries at Master dof from Neibs
  for(i=0;i<N_Neibs;i++)
   {
    N = N_RecevDofs[i];
    Neib_ID = NeibsRank[i];

    for(j=0;j<N;j++)
     {
      M = i*MaxN_SlaveDofs_All + j;
      dof = LocalDofOfRecbuf[0][M];
      N_UCol = LocalDofOfRecbuf[1][M]; 
      distdof = DistDofofLocDof[dof];      
      
      // A
      begin = DistA_RowPtr[distdof];
      len_A = DistA_RowPtr[distdof+1] - begin;
      KColBegin = DistA_KCol+begin;
      
      //system mat
      rowbegin = DistMat_RowPtr[distdof];
      row_Entries = DistMat_Entries+rowbegin;  
      
      M *= MaxN_NzInDepRow_All_A;     
      for(l=0; l<N_UCol; l++)
       {
        // find the col pos in the dist Kcol
        index_dep = GetIndex(KColBegin, len_A, recevbuf_KcolA[0][M]); 
        row_Entries[index_dep] += recevbuf[0][M++]; 
       }//    for(l=0; l<N_UCol; l++)   
     }//  for(j=0;j<N;j++)     
   } // for(i=0;i<N_Neibs;i++)
  
  if(N_Neibs)
   {
   delete [] recevbuf[0]; 
   delete [] recevbuf_KcolA[0]; 
   delete [] LocalDofOfRecbuf[0];
   delete [] LocalDofOfRecbuf[1];

   delete [] recevbuf; 
   delete [] N_RecevDofs;
   delete [] recevbuf_KcolA;
   delete [] LocalDofOfRecbuf;
  }

  // now copy own values
  for(i=0; i<nu; i++)
   {
    distdof=DistDofofLocDof[i];

    if(distdof==-1) continue; // this row not belongs to the dist mat

    row_Entries = DistMat_Entries+DistMat_RowPtr[distdof]; //u1 row

    begin = DistA_RowPtr[distdof];
    len_A = DistA_RowPtr[distdof+1] - begin;
    KColBegin = DistA_KCol+begin;

    end =  RowPtrA[i+1];
    for(j=RowPtrA[i]; j<end; j++)
     {
      // store the index for non-linear iteration
      index_dep = GetIndex(KColBegin, len_A,  GlobalDofOfLocalDof[KColA[j]]); 
      row_Entries[index_dep] +=  EntriesA[j];
     } // for(j=begin; j<end; j++)
   }

//   printf("AssembleDistMatrix Done \n");  
//   MPI_Finalize();
//   exit(0);
  
}// AssembleDistMatrix



/** generate all data needed for system dist matrix with local numbering*/
void TScalar_ParSolver::SetUpIterDistMat()
{
 int rank, i, j, len;
  
 this->SetUpDistMat();
  
  if(N_Rhs>1)
   {
    printf("Rank %d More than one rhs not yet implemented in iterative solver %d\n", rank,   N_Rhs);    
    MPI_Finalize();
    exit(0);   
   }
 
  N_Eqns = Global_DOF; 
  N_DistMatEntries = DistA_N_Entries;  
  N_DistMatEqns = N_distDof;
  DistMat_RowPtr = new int[N_DistMatEqns+1];  
  
  DistMat_RowPtr[0] = 0;
  for(i=0;i<N_distDof;i++)
   { 
    len = DistA_RowPtr[i+1] - DistA_RowPtr[i];
    DistMat_RowPtr[i+1] = DistMat_RowPtr[i] +  len;     
   }  
   
  RHS = new double[N_distDof];
  TCLS_SOL = new double[N_Eqns];
 
  memset(RHS, 0, N_distDof*SizeOfDouble); 
  memset(TCLS_SOL, 0, N_Eqns*SizeOfDouble); 
   
//   MPI_Comm_rank(Comm, &rank);  
//   printf("Rank %d SetUpIterDistMat N_Active %d\n", rank,   N_distActiveDof);    
//   MPI_Finalize();
//   exit(0);
      
} //  TScalar_ParSolver::SetUpIterDistMat()



/** generate all data needed for system dist matrix with global numbering*/
void TScalar_ParSolver::SetUpDistMat()
{
 int i, ii, j, k, l, M, N, S, Neib_ID, rank, size, tmp;
 int n_dof, *RowPtrA, *KColA, N_Active;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *DeptDofs, *N_DeptDofNeibs;
 int MaxN_DeptDofs_All, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocDof;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
 int MaxN_MasterDofs_All, MaxN_SlaveDofs_All, *DeptDofMasterID, *N_DofRankIndex;
 int MaxN_NzInRow, dof; 
 int **LocalDofOfRecbuf, **LocalDofOfRecbuf_send;
 int *N_SendDofs, *N_RecevDofs, **sendbuf, **recevbuf, disp, Master_ID;
 int *tmp_RowPtrA, MaxN_NzInRow_All_A, *aux, *aux_loc, N_Col;
 int *tmp_KColA, last, *DepDofIndexOfLocDOF;
 int N_MastDofToNeib_max, *aux_ColDof_SlaveToMaster, coldof, *neibaux_ColDof_SlaveToMaster;
 int *N_MastDofToNeib, N_MastDofToNeib_All, *MastDofToNeib;
 bool UPDATE;
 
  MPI_Comm_rank(Comm, &rank); 
  MPI_Comm_size(Comm, &size);

  n_dof = SquareStructure->GetN_Rows();
  RowPtrA = SquareStructure->GetRowPtr();
  KColA = SquareStructure->GetKCol();
  N_Active = SquareStructure->GetActiveBound();

  FEComm->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs,
                     DeptDofs, N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks,
                     IndexOfNeibRank, LocDOFofNeib, NeibLocDof);

  FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  FEComm->GetDeptDofMasterInfo(MaxN_MasterDofs_All, MaxN_SlaveDofs_All, DeptDofMasterID);

  N_DofRankIndex = FEComm->GetN_DofRankIndex();
  DepDofIndexOfLocDOF = FEComm->GetDepDofIndexOfLocDof();
  

  Global_beginDof = 0;
  for(i=0;i<rank;i++)
   Global_beginDof += Dist_Dof_All[i];
MPI_Finalize();
exit(0);
  MaxN_NzInRow = -1;
  for(i=0;i<N_DeptDofs;i++)
   {
    dof = DeptDofs[i];
    if(MaxN_NzInRow < (RowPtrA[dof+1] - RowPtrA[dof] )  )
     MaxN_NzInRow = RowPtrA[dof+1] - RowPtrA[dof];
   }  

  MPI_Allreduce(&MaxN_NzInRow, &MaxN_NzInDepRow_All_A, 1, MPI_INT, MPI_MAX, Comm);
  // printf("Rank %d MaxN_NzInDepRow_All_A %d  \n", rank, MaxN_NzInDepRow_All_A);

  LocalDofOfRecbuf = new int*[1];
  LocalDofOfRecbuf_send = new int*[1];
  
  /** generate data for A */  
  if(N_Neibs)
   {
    N_SendDofs =  new int[N_Neibs];
    N_RecevDofs =  new int[N_Neibs];
    N_MastDofToNeib = new int[N_Neibs];
 
    
    
    memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);
    memset(N_MastDofToNeib, 0, N_Neibs*SizeOfInt);
    
    sendbuf = new int*[1];
    recevbuf = new int*[1];

    // 1st column stores N_columns of U 
    disp = MaxN_NzInDepRow_All_A + 1;
    
    sendbuf[0] = new int[N_Neibs*MaxN_SlaveDofs_All*disp];
    recevbuf[0]  = new int[N_Neibs*MaxN_SlaveDofs_All*disp];

    LocalDofOfRecbuf[0] = new int[N_Neibs*MaxN_SlaveDofs_All];
    LocalDofOfRecbuf_send[0]  = new int[N_Neibs*MaxN_SlaveDofs_All];

    // max possible size  
    aux_ColDof_SlaveToMaster =  new int[N_Neibs*MaxN_SlaveDofs_All*disp];    
    
    for(i=0;i<N_DeptDofs;i++)
     {
      N = N_DeptDofNeibs[i];
      Master_ID = DeptDofMasterID[i];

      if(Master_ID==rank) // master of this row, info will be recieved from all neibs
       {
        for(j=0;j<N;j++)
         {
          S = i*MaxSubDomainPerDof + j;
          Neib_ID = DeptDofNeibRanks[S];
          k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
          N_RecevDofs[k]++;

         } // for(j=0;j<N;j++)
       }
      else  //slave of this dof, info of this row have to be send to the master
       {
        // fill the neib's (i.e., Master's) local dof
        for(j=0;j<N;j++)
         {
          S = i*MaxSubDomainPerDof + j;
          Neib_ID = DeptDofNeibRanks[S];
          if(Neib_ID==Master_ID)
           break;
         }
         
        k = IndexOfNeibRank[Master_ID]; // LocindexOfNeib_ID
        M = k*MaxN_SlaveDofs_All + N_SendDofs[k];
        LocalDofOfRecbuf_send[0][M] = NeibLocDof[S];

        //1st column to store N_columns of U  
        //further are global colum index of U  
        dof = DeptDofs[i];
        M *=disp;

        last = RowPtrA[dof+1];
        sendbuf[0][M++] = last - RowPtrA[dof]; // number of columns

        for(l=RowPtrA[dof];l<last;l++)
         {
          coldof = KColA[l];
          sendbuf[0][M++] = GlobalDofOfLocalDof[coldof];
  
/**generate info for updating sol from master to slave in iterative solver */
/**send all master dofs assiciated with this slave row to the master of the slave row */
          UPDATE = FALSE;
          if(N_DofRankIndex[coldof]==1) // self dof, i.e., Master coldof
           { 
            UPDATE = TRUE; 
           }
          else
           {
            ii = DepDofIndexOfLocDOF[coldof];
            if(DeptDofMasterID[ii]==rank)
             {
              UPDATE = TRUE; 
             } 
            } //  else
            
          if(UPDATE)  // Master coldof, this coldof has to be send to master of the slave row
          {
           tmp = N_MastDofToNeib[k];  
           neibaux_ColDof_SlaveToMaster = aux_ColDof_SlaveToMaster+(k*MaxN_SlaveDofs_All*disp);   
           for(ii=0; ii<tmp; ii++)
            if(coldof==neibaux_ColDof_SlaveToMaster[ii]) // already added
             {
              UPDATE = FALSE;
              break;
             }     
  
           if(UPDATE)   
            { 
             neibaux_ColDof_SlaveToMaster[N_MastDofToNeib[k]] = coldof;
             N_MastDofToNeib[k]++;
            }
           } // if(UPDATE)     
         } // for(l=RowPtrA[dof];l<last;l++)
  
        N_SendDofs[k]++;
       } //  else, if(Master_ID==rank)       
     }//  for(i=0;i<N_DeptDofs;i++)   
     
     
    N_MastDofToNeib_max = 0;
    for(i=0;i<N_Neibs;i++)
     if(N_MastDofToNeib_max < N_MastDofToNeib[i]) 
       N_MastDofToNeib_max = N_MastDofToNeib[i];

    MPI_Allreduce(&N_MastDofToNeib_max, &N_MastDofToNeib_All, 1, MPI_INT, MPI_MAX, Comm);  

    MastDofToNeib = new int[N_Neibs*N_MastDofToNeib_All];
    memset(MastDofToNeib, -1, N_Neibs*N_MastDofToNeib_All*SizeOfInt);
    
    for(i=0;i<N_Neibs;i++)
     {
      aux_loc = MastDofToNeib+(i*N_MastDofToNeib_All);
      neibaux_ColDof_SlaveToMaster  = aux_ColDof_SlaveToMaster+(i*MaxN_SlaveDofs_All*disp);
      last = N_MastDofToNeib[i];
      
      for(l=0;l<last;l++)     
       aux_loc[l] = neibaux_ColDof_SlaveToMaster[l];
     } 
     
     delete [] aux_ColDof_SlaveToMaster;    
    
    // set all these info in ParSol array
    ParSol->SetMastDofToHalloNeib(N_Neibs, N_MastDofToNeib, N_MastDofToNeib_All, MastDofToNeib);  
     
   }//   if(N_Neibs)
   
// ==============================================================
//   communicate
  FEComm->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_SlaveDofs_All, N_SendDofs,
                            LocalDofOfRecbuf,  N_RecevDofs, 1);

  for(i=0;i<N_Neibs;i++)
   {
    N_SendDofs[i] *= disp;
    N_RecevDofs[i] *= disp;
   }
  FEComm->FECommunicateNeib(sendbuf, MaxN_SlaveDofs_All*disp, N_SendDofs, recevbuf, N_RecevDofs, 1);
  for(i=0; i<N_Neibs; i++)
    N_RecevDofs[i] /= disp;

  if(N_Neibs)
   {
    delete [] N_SendDofs;
    delete [] sendbuf[0];
    delete [] sendbuf;
   }

// ==============================================================
// start to fill all necessary info 
// =================================================================
  tmp_RowPtrA = new int[n_dof];
  memset(tmp_RowPtrA, 0, n_dof*SizeOfInt);

  //find total number of entries in the dist system matrix
  DistA_N_Entries = 0;

  // first fill self dofs 
  for(i=0;i<n_dof;i++)  
   if(N_DofRankIndex[i]==1)
    {
     tmp = RowPtrA[i+1] - RowPtrA[i];
     DistA_N_Entries +=  tmp;  
     tmp_RowPtrA[i] += tmp;    
    } // if(N_DofRankIndex[i]==1)
  
   // fill master dofs
   for(i=0;i<N_DeptDofs;i++)
    if(DeptDofMasterID[i]==rank)
     {
      dof = DeptDofs[i];
      tmp = RowPtrA[dof+1] - RowPtrA[dof];

      DistA_N_Entries += tmp;
      tmp_RowPtrA[dof] += tmp;
     }  

  // now fill the master dof info from neibs
  for(i=0;i<N_Neibs;i++)
   {
    N = N_RecevDofs[i];

    for(j=0;j<N;j++)
     {
      M = i*MaxN_SlaveDofs_All + j;
      dof = LocalDofOfRecbuf[0][M];
      M *=disp;

      DistA_N_Entries +=  recevbuf[0][M];
      tmp_RowPtrA[dof] += recevbuf[0][M];
     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)

  // fill col with Global dofs
  MaxN_NzInRow_All_A = -1;
  for(i=0;i<n_dof;i++)
   if(MaxN_NzInRow_All_A <  tmp_RowPtrA[i]   )
     MaxN_NzInRow_All_A = tmp_RowPtrA[i];

  aux= new int[n_dof*MaxN_NzInRow_All_A];
  memset(aux, 0, n_dof*MaxN_NzInRow_All_A*SizeOfInt);
  memset(tmp_RowPtrA, 0, n_dof*SizeOfInt);

  // now put column info of self dof
  for(i=0;i<n_dof;i++)
   {
    aux_loc = aux + i*MaxN_NzInRow_All_A;

    if(N_DofRankIndex[i]==1)
     {
      tmp = RowPtrA[i+1];
      for(j=RowPtrA[i]; j<tmp; j++)
       aux_loc[tmp_RowPtrA[i]++] = GlobalDofOfLocalDof[KColA[j]];
     } //if(N_DofRankIndex[i]==1)
    }
  
   // fill master dofs own columns  
   for(i=0;i<N_DeptDofs;i++)
    if(DeptDofMasterID[i]==rank)
     {
      dof = DeptDofs[i];
      aux_loc = aux + dof*MaxN_NzInRow_All_A; 
      tmp = RowPtrA[dof+1];
      
      for(j=RowPtrA[dof]; j<tmp; j++)
       aux_loc[tmp_RowPtrA[dof]++] = GlobalDofOfLocalDof[KColA[j]];
     }  // f(DeptDofMasterID[
  
  // now fill the master dof info from neibs
  for(i=0;i<N_Neibs;i++)
   {
    N = N_RecevDofs[i];

    for(j=0;j<N;j++)
     {
      M = i*MaxN_SlaveDofs_All + j;
      dof = LocalDofOfRecbuf[0][M];
      M *=disp;

      aux_loc = aux + dof*MaxN_NzInRow_All_A;
      N_Col = recevbuf[0][M++];

      for(l=0; l<N_Col; l++)
        aux_loc[tmp_RowPtrA[dof]++] = recevbuf[0][M++];

     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)

  // neib info nolonger needed free the memory
  if(N_Neibs)
   {
    delete [] N_RecevDofs;
    delete [] recevbuf[0];
    delete [] recevbuf;
    delete [] LocalDofOfRecbuf[0];
    delete [] LocalDofOfRecbuf_send[0];
   }
  
   // process system matrix
   tmp_KColA = new int[DistA_N_Entries];

   N_distDof = Dist_Dof_All[rank];        /** init class variable */ 
   DistA_RowPtr = new int[N_distDof+1];   /** init class variable */ 
   DistDofofLocDof = new int[n_dof];      /** init class variable */ 
   memset(DistDofofLocDof, -1, n_dof*SizeOfInt);
  
   //find total number of entries in the dist system matrix
   DistA_N_Entries = 0;                    /** init class variable */
   DistA_RowPtr[0] = 0;                    /** init class variable */
   N_distActiveDof = 0;                    /** init class variable */
   N=0;

   for(i=0;i<n_dof;i++)
    {
     aux_loc = aux + i*MaxN_NzInRow_All_A;
     M = tmp_RowPtrA[i];

     if(M>1)
       Sort(aux_loc, M);

     last=-1;
     for(j=0;j<M;j++)
      if(aux_loc[j] != last)
       {
        last = aux_loc[j];
        tmp_KColA[DistA_N_Entries++] = last;
       }

     //set the row pointer
     if(last!=-1)
      {
       DistDofofLocDof[i] = N;
       DistA_RowPtr[++N] = DistA_N_Entries;

       if(i<N_Active)
        {
         N_distActiveDof++;
// 	 if(rank==0)
// 	   printf("Rank %d N_distActiveDof %d i %d M %d \n", rank, N_distActiveDof, i, GlobalDofOfLocalDof[i]);  
	}
      }
    } //  for(i=0;i<n_dof;i++)

  delete [] aux;
  delete [] tmp_RowPtrA;
  
  DistA_KCol = new int[DistA_N_Entries]; 
  memcpy(DistA_KCol, tmp_KColA, DistA_N_Entries*SizeOfInt);

  delete [] tmp_KColA;  
  
//   int count, N_OwnDof, *OwnDofs;
//   FEComm->GetOwnDofs(N_OwnDof, OwnDofs);    
  
//   count=1;
  
//   if(rank==1)
//   for(i=0;i<N_OwnDof;i++)
//      printf("Rank %d index %d OwnDofs  %d GlobalDofOfLocalDof   %d \n", rank, count++, OwnDofs[i], GlobalDofOfLocalDof[OwnDofs[i]]);  
//       if(OwnDofs[i]<N_Active)
  
//   printf("Rank %d SetUpDistMat n_dof %d DistA_N_Entries %d N_Active %d\n", rank, n_dof , N_OwnDof, N_distActiveDof);  
//   MPI_Finalize();
//   exit(0);
    
} //TScalar_ParSolver::SetUpDistMat()



void TScalar_ParSolver::InitMumps()
{
 int i, I, j, end, pos, rank;

  N_Eqns = Global_DOF; 
  N_DistMatEntries = DistA_N_Entries;
  DistMat_Irn = new int[N_DistMatEntries];
  DistMat_Jcn = new int[N_DistMatEntries];
  
  
  N_DistMatEqns = N_distDof;
  DistMat_RowPtr = new int[N_DistMatEqns+1];
    
  DistMat_RowPtr[0] = 0;
  pos = 0;
  for(i=0;i<N_distDof;i++)
   { 
    I = Global_beginDof + i + 1; // fortran style
    end = DistA_RowPtr[i+1];
    for(j=DistA_RowPtr[i]; j<end; j++)
     {
      DistMat_Irn[pos] = I;    
      DistMat_Jcn[pos] = DistA_KCol[pos] + 1; // fortran style
      pos++;
     }
     
    DistMat_RowPtr[i+1] = pos;     
   }

/*  MPI_Comm_rank(Comm, &rank);   
  if(rank==TDatabase::ParamDB->Par_P5)
   {
    pos=0;
    for(i=0;i<N_distDof;i++)
    {
      printf("Rank %d  u1row %d    ",  rank, DistMat_Irn[pos]  );

      for(j=DistA_RowPtr[i]; j<DistA_RowPtr[i+1]; j++)
       printf("%d, ", DistMat_Jcn[pos++]);

     printf(" \n");
    }  
  }*/
   
/*  printf("%d Mumps pos %d  N_DistMatEntries %d \n", rank, Global_beginDof, DistA_N_Entries);
   MPI_Finalize();
 exit(0); */ 


 MUMPS_Solver = new TMumpsSolver(N_Eqns, N_DistMatEntries, DistMat_Irn, DistMat_Jcn, N_Rhs);
} //TScalar_ParSolver::InitMumps()



void TScalar_ParSolver::AssembleIterDistMatrix(TSquareMatrix3D *Mat)
{
 int rank, i, N_LocDof, distdof;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
 
  double *rhs_loc, *sol_loc;
  
  
  FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  
  
  this->AssembleDistMatrix(Mat); 
  
  N_LocDof = Mat->GetN_Rows();
  
  rhs_loc = new double[N_LocDof];    
  memcpy(rhs_loc, ParRhs->GetValues(), N_LocDof*SizeOfDouble);    
  ParRhs->AssembleByADD(rhs_loc); 
  sol_loc = ParSol->GetValues();
  
  
  // now copy self + master values
  for(i=0; i<N_LocDof; i++)
   {
    distdof=DistDofofLocDof[i];
    if(distdof==-1) continue; // this row not belongs to the dist mat
  
    RHS[distdof] = rhs_loc[i];
    TCLS_SOL[GlobalDofOfLocalDof[i]] = sol_loc[i];
   }//for(i=0; i<N_LocDof; i++)
   
   
   // now copy the slave dep. dofs and Hallo dofs value from neib masters
   ParSol->UpdateDepAndHalloDof(TCLS_SOL);
   
  MPI_Comm_rank(Comm, &rank);  
  printf("Rank %d AssembleIterDistMatrix \n", rank);    
  MPI_Finalize();
  exit(0); 
  
}


TScalar_ParSolver::~TScalar_ParSolver()
{

}

#endif

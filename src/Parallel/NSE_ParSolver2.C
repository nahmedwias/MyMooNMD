// =======================================================================
// @(#)NSE_ParSolver2.C
//
// Class:      TNSE_ParSolver2
// Purpose:    Class containing all info needed for communication between subdomains
//             for square matrices
//
// Author:     Sashikumaar Ganesan (19.10.10)
//
// History:    Start of implementation 19.10.10 (Sashikumaar Ganesan)
//
// =======================================================================

#ifdef _MPI
#  include "mpi.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <NSE_ParSolver2.h>
#include <Database.h>
#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>
#include <ParVectorNSE3D.h>
#include <LinAlg.h>

    /** constructor */
TNSE_ParSolver2::TNSE_ParSolver2(TParFECommunicator3D *velo_FEComm, 
        TSquareStructure3D *squareStructure, TParFECommunicator3D *p_FECom, 
        TStructure3D * structureBT, TStructure3D *structureB)
            : TNSE_ParSolver(velo_FEComm, squareStructure, p_FECom,
                            structureBT, structureB, 2)
 {

   //set up the parallel solver
   switch(int(TDatabase::ParamDB->SOLVER_TYPE))
    {
      case 101:
        // MUMPS Parallel solver
        InitMumps();

       break;

       default: 
        cout << "wrong  Parallel solver type !!!!!!!!!!!!!" << endl;
        exit(0);
       break;

     }

 }


void TNSE_ParSolver2::Solve(TSquareMatrix3D *MatA, TMatrix3D *MatB1T, TMatrix3D *MatB2T,
                       TMatrix3D *MatB3T, TMatrix3D  *MatB1, TMatrix3D *MatB2, TMatrix3D * MatB3,
                       TParVectorNSE3D  *rhs, TParVectorNSE3D  *sol)
{
 int i, rank, size, locbegin_P;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All, distbegin, locbegin, len, len_P;
 int N_GlobDof_P, *GlobalDofOfLocalDof_P, *Dist_Dof_All_P;
 int nu, np, N_Active, N_NonActive, distdof;

 double *tmp_rhs, *SOL, *OLDSOL, *tmp_sol;

 this->AssembleDistMatrix(MatA, MatB1T, MatB2T, MatB3T, MatB1, MatB2, MatB3, rhs->GetValues());  
   
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

   //set up the parallel solver
   switch(int(TDatabase::ParamDB->SOLVER_TYPE))
    {
      case 101:

        //centralized rhs at host for MUMPS
        if(rank==0)
         {   
          tmp_rhs = new double[N_Eqns];
           memset(tmp_rhs, 0, N_Eqns*SizeOfDouble);
         }

       rhs->AssembleVeloAtRootByADD(tmp_rhs);

       if(rank==0)
        {
         RHS = new double[N_Eqns];
         memset(RHS, 0, N_Eqns*SizeOfDouble);     
       
         Velo_FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
         P_FEComm->GetGlobalDistDofInfo(N_GlobDof_P, GlobalDofOfLocalDof_P, Dist_Dof_All_P);

         distbegin = 0;
         locbegin = 0;
         for(i=0;i<size;i++)
          {
           len = Dist_Dof_All[i];
           memcpy(RHS+distbegin, tmp_rhs+locbegin, len*SizeOfDouble);
           memcpy(RHS+(len+distbegin), tmp_rhs+(Global_N_U+locbegin), len*SizeOfDouble);
           memcpy(RHS+(2*len+distbegin), tmp_rhs+(2*Global_N_U+locbegin), len*SizeOfDouble);

           distbegin += 3*len + Dist_Dof_All_P[i];
           locbegin += len;
          }
          
//            printf("RHs Dot  %e\n", Ddot(N_Eqns, RHS, RHS));     
          
        }
 
//         cout << "wrong  Parallel solver type !!!!!!!!!!!!!" << endl;

         /** call the solver to solve the system */
         MUMPS_Solver->FactorizeAndSolve(DistMat_Entries, RHS);

         // put the sol in MooNMD format
         if(rank==0)
          {   
//            printf("Sol Dot  %e\n", Ddot(3*Global_N_U, RHS, RHS));     		    
//            printf("Sol Dot  %e\n", Ddot(Global_N_P, RHS+3*Global_N_U, RHS+3*Global_N_U));     	    

           distbegin = 0;
           locbegin = 0;
           locbegin_P = 0;
           for(i=0;i<size;i++)
            {
             len = Dist_Dof_All[i];
             len_P = Dist_Dof_All_P[i];
             memcpy(tmp_rhs+locbegin, RHS+distbegin, len*SizeOfDouble);
             memcpy(tmp_rhs+(Global_N_U+locbegin), RHS+(len+distbegin), len*SizeOfDouble);
             memcpy(tmp_rhs+(2*Global_N_U+locbegin), RHS+(2*len+distbegin), len*SizeOfDouble);
             memcpy(tmp_rhs+(3*Global_N_U+locbegin_P), RHS+(3*len+distbegin), len_P*SizeOfDouble);

             distbegin += 3*len + len_P;
             locbegin += len;
             locbegin_P += len_P;
            }
           }
/*          MPI_Finalize();
        exit(0);   */    
        // scatter the sol vector from root to all sub domains
        sol->ScatterFromRoot(tmp_rhs);

        if(rank==0)
         {
          delete [] RHS;
          delete [] tmp_rhs;
         }
         
       break;

       default: 
        cout << "wrong  Parallel solver type !!!!!!!!!!!!!" << endl;
        MPI_Finalize();
        exit(0);
       break;

     }

}


void TNSE_ParSolver2::Solve(TSquareMatrix3D *MatA, TParVectorNSE3D  *rhs, TParVectorNSE3D  *sol)
{
 int i, rank, size, locbegin_P;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All, distbegin, locbegin, len, len_P;
 int N_GlobDof_P, *GlobalDofOfLocalDof_P, *Dist_Dof_All_P;
 int nu, np, N_Active, N_NonActive, distdof;

 double *tmp_rhs, *SOL;

 this->AssembleDistMatrix_A(MatA);

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

   //set up the parallel solver
   switch(int(TDatabase::ParamDB->SOLVER_TYPE))
    {
      case 101:

        //centralized rhs at host for MUMPS
        if(rank==0)
          tmp_rhs = new double[N_Eqns];

        rhs->AssembleVeloAtRootByADD(tmp_rhs);

        if(rank==0)
         {
          RHS = new double[N_Eqns];
          memset(RHS, 0, N_Eqns*SizeOfDouble);

          Velo_FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
          P_FEComm->GetGlobalDistDofInfo(N_GlobDof_P, GlobalDofOfLocalDof_P, Dist_Dof_All_P);

          distbegin = 0;
          locbegin = 0;
          for(i=0;i<size;i++)
           {
            len = Dist_Dof_All[i];
            memcpy(RHS+distbegin, tmp_rhs+locbegin, len*SizeOfDouble);
            memcpy(RHS+(len+distbegin), tmp_rhs+(Global_N_U+locbegin), len*SizeOfDouble);
            memcpy(RHS+(2*len+distbegin), tmp_rhs+(2*Global_N_U+locbegin), len*SizeOfDouble);

            distbegin += 3*len + Dist_Dof_All_P[i];
            locbegin += len;
           }
         }

         /** call the solver to solve the system */
 


         MUMPS_Solver->FactorizeAndSolve(DistMat_Entries, RHS);
//          MUMPS_Solver->Solve(DistMat_Entries, RHS);


         // put the sol in MooNMD format
         if(rank==0)
          {
           distbegin = 0;
           locbegin = 0;
           locbegin_P = 0;
           for(i=0;i<size;i++)
            {
             len = Dist_Dof_All[i];
             len_P = Dist_Dof_All_P[i];
             memcpy(tmp_rhs+locbegin, RHS+distbegin, len*SizeOfDouble);
             memcpy(tmp_rhs+(Global_N_U+locbegin), RHS+(len+distbegin), len*SizeOfDouble);
             memcpy(tmp_rhs+(2*Global_N_U+locbegin), RHS+(2*len+distbegin), len*SizeOfDouble);
             memcpy(tmp_rhs+(3*Global_N_U+locbegin_P), RHS+(3*len+distbegin), len_P*SizeOfDouble);

             distbegin += 3*len + len_P;
             locbegin += len;
             locbegin_P += len_P;
            }
           }

         // scatter the sol vector from root to all sub domains
         sol->ScatterFromRoot(tmp_rhs);

         if(rank==0)
          {
           delete [] RHS;
           delete [] tmp_rhs;
          }

       break;

       default: 
        cout << "wrong  Parallel solver type !!!!!!!!!!!!!" << endl;
        exit(0);
       break;

     }



}


/** Assemble dist non-linear dioganal matrix */
void TNSE_ParSolver2::AssembleDistMatrix_A(TSquareMatrix3D *MatA)
{
 int i, j, k, l, M, N, S, begin,  dispA, Master_ID, Neib_ID, dof, distdof, rank, size, end;
 int pos, N_UCol, len_A, col, index, index_dep;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *DeptDofs, *N_DeptDofNeibs;
 int MaxN_DeptDofs_All, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocDof;
 int MaxN_MasterDofs_All, MaxN_SlaveDofs_All, *DeptDofMasterID;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
 int **sendbuf_Kcol, **recevbuf_KcolA,  **recevbuf_KcolBT, *KColBegin;
 int nu, *RowPtrA, *KColA, N_Active;
 int **LocalDofOfRecbuf, **LocalDofOfRecbuf_send;
 int *N_SendDofs, *N_RecevDofs, Dist3U, Dist2U;
 
 double **sendbuf, **recevbuf, *EntriesA, *u1row_Entries, *u2row_Entries, *u3row_Entries;
 
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Velo_FEComm->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs,
                     DeptDofs, N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks,
                     IndexOfNeibRank, LocDOFofNeib, NeibLocDof);

  Velo_FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  Velo_FEComm->GetDeptDofMasterInfo(MaxN_MasterDofs_All, MaxN_SlaveDofs_All, DeptDofMasterID);

  nu = MatA->GetN_Rows();
  RowPtrA = MatA->GetRowPtr();
  KColA = MatA->GetKCol();
  N_Active = MatA->GetActiveBound();
  EntriesA =  MatA->GetEntries();  
  
 
  LocalDofOfRecbuf = new int*[2];
  LocalDofOfRecbuf_send = new int*[2];
  
   /** first generate data for A  matrix then BT and after B */
  if(N_Neibs)
   {
    N_SendDofs =  new int[N_Neibs];
    N_RecevDofs =  new int[N_Neibs];
    memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    sendbuf = new double*[1]; // u1, u2, u3 part will be send together
    recevbuf = new double*[1];
    sendbuf_Kcol = new int*[1];
    recevbuf_KcolA = new int*[1];   

    // N_columns of U 
    dispA = MaxN_NzInDepRow_All_A;
    sendbuf[0] = new double[N_Neibs*MaxN_SlaveDofs_All*dispA];
    recevbuf[0]  = new double[N_Neibs*MaxN_SlaveDofs_All*dispA];
    sendbuf_Kcol[0] = new int[N_Neibs*MaxN_SlaveDofs_All*dispA];
    recevbuf_KcolA[0] = new int[N_Neibs*MaxN_SlaveDofs_All*dispA];

     for(i=0;i<2;i++)
     {
      LocalDofOfRecbuf[i] = new int[N_Neibs*MaxN_SlaveDofs_All];
      LocalDofOfRecbuf_send[i]  = new int[N_Neibs*MaxN_SlaveDofs_All];
     }   

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

        M *=dispA;
        for(l=begin;l<end;l++)
         {
          sendbuf[0][M] =  EntriesA[l];
          sendbuf_Kcol[0][M++] =  GlobalDofOfLocalDof[KColA[l]];
         }

        N_SendDofs[k]++;
       }
     } //  for(i=0;i<N_DeptDofs;i++)
   } //   if(N_Neibs)

// =====================================================================================
//   communicate 
  Velo_FEComm->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_SlaveDofs_All, N_SendDofs,
                                 LocalDofOfRecbuf,  N_RecevDofs, 2);

  for(i=0; i<N_Neibs; i++)
   {
    N_SendDofs[i] *= dispA;
    N_RecevDofs[i] *= dispA;
   }
  Velo_FEComm->FECommunicateNeib(sendbuf, MaxN_SlaveDofs_All*dispA,
                                 N_SendDofs, recevbuf, N_RecevDofs, 1);
  Velo_FEComm->FECommunicateNeib(sendbuf_Kcol, MaxN_SlaveDofs_All*dispA,
                                 N_SendDofs, recevbuf_KcolA, N_RecevDofs, 1);
  for(i=0; i<N_Neibs; i++)
    N_RecevDofs[i] /= dispA;
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
   } 
    

  // assemble system dist matrix directly
  //set diognal mat entries as zeros
  Dist3U = 3*N_distU;
  Dist2U = 2*N_distU;
  index = 0;  
  for(i=0; i<nu; i++)
   {
    distdof=DistDofofLocDof_U[i];
    if(distdof==-1) continue; // this row not belongs to the dist mat

    u1row_Entries = DistMat_Entries+DistMat_RowPtr[distdof]; //u1 row
    u2row_Entries = DistMat_Entries+DistMat_RowPtr[N_distU + distdof]; //u2 row
    u3row_Entries = DistMat_Entries+DistMat_RowPtr[Dist2U + distdof]; //u3 row

    end =  RowPtrA[i+1];
    for(j=RowPtrA[i]; j<end; j++)
     {
      u1row_Entries[Aii_pos[index]] =  0.; // A11
      u2row_Entries[Aii_pos[index]] =  0.; // A22
      u3row_Entries[Aii_pos[index]] =  0.; // A33
      index++;
     } // for(j= RowPtrA[i]; j<end; j++)
   } //for(i=0; i<nu; i++  

  index_dep = 0; 
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
      distdof = DistDofofLocDof_U[dof];

      //system mat
      u1row_Entries = DistMat_Entries+DistMat_RowPtr[distdof];
      u2row_Entries = DistMat_Entries+DistMat_RowPtr[N_distU + distdof];
      u3row_Entries = DistMat_Entries+DistMat_RowPtr[Dist2U + distdof];

      pos = M*dispA;
      for(l=0; l<N_UCol; l++)
       {
        u1row_Entries[Aii_Deptpos[index_dep]] =  0.; // A11
        u2row_Entries[Aii_Deptpos[index_dep]] =  0.; // A22
        u3row_Entries[Aii_Deptpos[index_dep]] =  0.; // A33
        index_dep++;
        pos++;
       }
     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)   
 
 
  // first add entries at Master dof from Neibs
  index_dep = 0; 
  for(i=0;i<N_Neibs;i++)
   {
    N = N_RecevDofs[i];
    Neib_ID = NeibsRank[i];
    for(j=0;j<N;j++)
     {
      M = i*MaxN_SlaveDofs_All + j;
      dof = LocalDofOfRecbuf[0][M];
      N_UCol = LocalDofOfRecbuf[1][M];
      distdof = DistDofofLocDof_U[dof];

      // A
      begin = DistA_RowPtr[distdof];
      len_A = DistA_RowPtr[distdof+1] - begin;
      KColBegin = DistA_KCol+begin;

      //system mat
      u1row_Entries = DistMat_Entries+DistMat_RowPtr[distdof];
      u2row_Entries = DistMat_Entries+DistMat_RowPtr[N_distU + distdof];
      u3row_Entries = DistMat_Entries+DistMat_RowPtr[Dist2U + distdof];

      pos = M*dispA;
      for(l=0; l<N_UCol; l++)
       {
        col = recevbuf_KcolA[0][pos];
        index = GetIndex(KColBegin, len_A, col); // find the col pos in the dist Kcol
        // printf("AssembleDist  Rank %d  Neib %d N_SendDofs %d \n", rank, index, Aii_Deptpos[index_dep]);
        u1row_Entries[Aii_Deptpos[index_dep]] += recevbuf[0][pos]; // A11
        u2row_Entries[Aii_Deptpos[index_dep]] += recevbuf[0][pos]; // A22
        u3row_Entries[Aii_Deptpos[index_dep]] += recevbuf[0][pos]; // A33
        index_dep++;     
        pos++;
       }
     } // for(j=0;j<N;j++)
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
  index = 0;
  for(i=0; i<nu; i++)
   {
    distdof=DistDofofLocDof_U[i];
    if(distdof==-1) continue; // this row not belongs to the dist mat
    u1row_Entries = DistMat_Entries+DistMat_RowPtr[distdof]; //u1 row
    u2row_Entries = DistMat_Entries+DistMat_RowPtr[N_distU + distdof]; //u2 row
    u3row_Entries = DistMat_Entries+DistMat_RowPtr[Dist2U + distdof]; //u3 row

    end =  RowPtrA[i+1];
    for(j=RowPtrA[i]; j<end; j++)
     {
      u1row_Entries[Aii_pos[index]] +=  EntriesA[j]; // A11
      u2row_Entries[Aii_pos[index]] +=  EntriesA[j]; // A22
      u3row_Entries[Aii_pos[index]] +=  EntriesA[j]; // A33
      index++;
     } // for(j=begin; j<end; j++)
   } //for(i=0; i<nu; i++  
  
} // AssembleDistMatrix_A

void TNSE_ParSolver2::AssembleDistMatrix(TSquareMatrix3D *MatA, TMatrix3D *MatB1T, 
                       TMatrix3D *MatB2T, TMatrix3D *MatB3T, TMatrix3D  *MatB1, 
                       TMatrix3D *MatB2, TMatrix3D * MatB3, double *Rhs)
{
 int i, j, k, l, m, M, N, S, dispA, dispBT, dispB, Master_ID, Neib_ID, dof, distdof, rank, size, begin, end;
 int pos, N_UCol, N_PCol, len_A, len_BT, len_B, col, index, index_pos, index_dep;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *DeptDofs, *N_DeptDofNeibs;
 int MaxN_DeptDofs_All, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocDof;
 int MaxN_MasterDofs_All, MaxN_SlaveDofs_All, *DeptDofMasterID; //, *N_DofRankIndex;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
 int u1_rowbegin, u2_rowbegin, u3_rowbegin, J1, J2, J3;

 int MaxSubDomainPerDof_P, N_Neibs_P, *NeibsRank_P, N_DeptDofs_P, *DeptDofs_P, *N_DeptDofNeibs_P;  
 int MaxN_DeptDofs_All_P, *DeptDofNeibRanks_P, *IndexOfNeibRank_P, *LocDOFofNeib_P, *NeibLocDof_P;
 int N_GlobDof_P, *GlobalDofOfLocalDof_P, *Dist_Dof_All_P;
 int MaxN_MasterDofs_All_P, MaxN_SlaveDofs_All_P, *DeptDofMasterID_P;

 int **LocalDofOfRecbuf, **LocalDofOfRecbuf_send;
 int *N_SendDofs, *N_RecevDofs, Dist3U, Dist2U;
 int nu, *RowPtrA, *KColA, N_Active, np, *RowPtrBT, *KColBT, *RowPtrB, *KColB;
 int **sendbuf_Kcol, **recevbuf_KcolA,  **recevbuf_KcolBT, *KColBegin;

 double **sendbuf, **recevbuf, *EntriesA, *EntriesB1T, *EntriesB2T, *EntriesB3T;
 double *EntriesB1, *EntriesB2, *EntriesB3;
 double **recebufBT, *u1row_Entries, *u2row_Entries, *u3row_Entries, *prow_Entries;


  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Velo_FEComm->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs,
                     DeptDofs, N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks,
                     IndexOfNeibRank, LocDOFofNeib, NeibLocDof);

  Velo_FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  Velo_FEComm->GetDeptDofMasterInfo(MaxN_MasterDofs_All, MaxN_SlaveDofs_All, DeptDofMasterID);


  P_FEComm->GetDofNeibInfo(MaxSubDomainPerDof_P, N_Neibs_P, NeibsRank_P, N_DeptDofs_P, DeptDofs_P, 
                     N_DeptDofNeibs_P, MaxN_DeptDofs_All_P, DeptDofNeibRanks_P, IndexOfNeibRank_P,
                     LocDOFofNeib_P, NeibLocDof_P);
  P_FEComm->GetGlobalDistDofInfo(N_GlobDof_P, GlobalDofOfLocalDof_P, Dist_Dof_All_P);
  P_FEComm->GetDeptDofMasterInfo(MaxN_MasterDofs_All_P, MaxN_SlaveDofs_All_P, DeptDofMasterID_P);

  nu = MatA->GetN_Rows();
  RowPtrA = MatA->GetRowPtr();
  KColA = MatA->GetKCol();
  N_Active = MatA->GetActiveBound();
  EntriesA =  MatA->GetEntries();

  RowPtrBT = MatB1T->GetRowPtr();
  KColBT = MatB1T->GetKCol();
  EntriesB1T =  MatB1T->GetEntries();
  EntriesB2T =  MatB2T->GetEntries();
  EntriesB3T =  MatB3T->GetEntries();

  np = MatB1->GetN_Rows();
  RowPtrB = MatB1->GetRowPtr();
  KColB = MatB1->GetKCol();
  EntriesB1 = MatB1->GetEntries();
  EntriesB2 = MatB2->GetEntries();
  EntriesB3 = MatB3->GetEntries();

  LocalDofOfRecbuf = new int*[3];
  LocalDofOfRecbuf_send = new int*[3];

  /** first generate data for A  matrix then BT and after B */
  if(N_Neibs)
   {
    N_SendDofs =  new int[N_Neibs];
    N_RecevDofs =  new int[N_Neibs];
    memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    sendbuf = new double*[1]; // u1, u2, u3 part will be send together
    recevbuf = new double*[1];
    sendbuf_Kcol = new int*[1];
    recevbuf_KcolA = new int*[1];

    // N_columns of U 
    dispA = MaxN_NzInDepRow_All_A;
    sendbuf[0] = new double[N_Neibs*MaxN_SlaveDofs_All*dispA];
    recevbuf[0]  = new double[N_Neibs*MaxN_SlaveDofs_All*dispA];
    sendbuf_Kcol[0] = new int[N_Neibs*MaxN_SlaveDofs_All*dispA];
    recevbuf_KcolA[0] = new int[N_Neibs*MaxN_SlaveDofs_All*dispA];

    for(i=0;i<3;i++)
     {
      LocalDofOfRecbuf[i] = new int[N_Neibs*MaxN_SlaveDofs_All];
      LocalDofOfRecbuf_send[i]  = new int[N_Neibs*MaxN_SlaveDofs_All];
     }

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

//         if(dof<N_Active)
//          {
          begin = RowPtrA[dof];
          end=RowPtrA[dof+1];

          LocalDofOfRecbuf_send[0][M] = NeibLocDof[S];
          LocalDofOfRecbuf_send[1][M] = end - begin;
          LocalDofOfRecbuf_send[2][M] = RowPtrBT[dof+1] - RowPtrBT[dof];

          M *=dispA;
          for(l=begin;l<end;l++)
           {
            sendbuf[0][M] =  EntriesA[l];
            sendbuf_Kcol[0][M++] =  GlobalDofOfLocalDof[KColA[l]];

//            if(Master_ID==1 &&  NeibLocDof[S]==1072)
//              printf("Rank %d slave dof %d  KColA %d GlobKColA %d M %d total %d\n", rank,   dof,    KColA[l], GlobalDofOfLocalDof[KColA[l]], M-1, end - begin); 
           }

//           }
//          else // Dirichlet row values wil not be added to the master row
//           {
//            LocalDofOfRecbuf_send[0][M] = -1;
//           }

        N_SendDofs[k]++;
       }
     } //  for(i=0;i<N_DeptDofs;i++)
   } //   if(N_Neibs)

// =====================================================================================
//   communicate 
  Velo_FEComm->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_SlaveDofs_All, N_SendDofs,
                                 LocalDofOfRecbuf,  N_RecevDofs, 3);

  for(i=0; i<N_Neibs; i++)
   {
    N_SendDofs[i] *= dispA;
    N_RecevDofs[i] *= dispA;
   }
  Velo_FEComm->FECommunicateNeib(sendbuf, MaxN_SlaveDofs_All*dispA,
                                 N_SendDofs, recevbuf, N_RecevDofs, 1);
  Velo_FEComm->FECommunicateNeib(sendbuf_Kcol, MaxN_SlaveDofs_All*dispA,
                                 N_SendDofs, recevbuf_KcolA, N_RecevDofs, 1);
  for(i=0; i<N_Neibs; i++)
    N_RecevDofs[i] /= dispA;
// =====================================================================================

  if(N_Neibs)
   {
    delete [] sendbuf[0];
    delete [] sendbuf_Kcol[0];
    delete [] LocalDofOfRecbuf_send[0];
    delete [] LocalDofOfRecbuf_send[1];
    delete [] LocalDofOfRecbuf_send[2];

    delete [] sendbuf;
    delete [] LocalDofOfRecbuf_send;

    sendbuf = new double*[3];
    recebufBT = new double*[3];
    recevbuf_KcolBT = new int*[1];

    // N_columns of B1T, B2T, B3T
    dispBT = MaxN_NzInDepRow_All_BT;
    for(i=0;i<3;i++)
     {
      sendbuf[i] = new double[N_Neibs*MaxN_SlaveDofs_All*dispBT];
      recebufBT[i] = new double[N_Neibs*MaxN_SlaveDofs_All*dispBT];
     }

    sendbuf_Kcol[0] = new int[N_Neibs*MaxN_SlaveDofs_All*dispBT];
    recevbuf_KcolBT[0] = new int[N_Neibs*MaxN_SlaveDofs_All*dispBT];

    // npw send dep dof values of BT to the master row
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    for(i=0;i<N_DeptDofs;i++)
     {
      Master_ID = DeptDofMasterID[i];

      if(Master_ID!=rank)
       {
        dof = DeptDofs[i];
        k = IndexOfNeibRank[Master_ID];

//         if(dof<N_Active)
//          {
          M = (k*MaxN_SlaveDofs_All + N_SendDofs[k])*dispBT;

          end=RowPtrBT[dof+1];
          for(l=RowPtrBT[dof];l<end;l++)
           {
            sendbuf[0][M] =  EntriesB1T[l];
            sendbuf[1][M] =  EntriesB2T[l];
            sendbuf[2][M] =  EntriesB3T[l];
            sendbuf_Kcol[0][M] = GlobalDofOfLocalDof_P[KColBT[l]];
            M++;
           }
//           }
//          else
//           {
//            Rhs[dof] = 0.; //Dep only Master Dirichlet value will be be taken
//           }


        N_SendDofs[k]++;
       }
     } //  for(i=0;i<N_DeptDofs;i++)
   }

// =====================================================================================
//   communicate 
  for(i=0; i<N_Neibs; i++)
   {
    N_SendDofs[i] *= dispBT;
    N_RecevDofs[i] *= dispBT;
   }

//   for(i=0; i<N_Neibs; i++)
// printf("AssembleDistMatrix Rank %d  Neib %d N_SendDofs %d \n", rank, NeibsRank[i], N_SendDofs[i]);

//  TDatabase::ParamDB->Par_P7 = 1;

  Velo_FEComm->FECommunicateNeib(sendbuf, MaxN_SlaveDofs_All*dispBT,
                                 N_SendDofs, recebufBT, N_RecevDofs, 3);

  Velo_FEComm->FECommunicateNeib(sendbuf_Kcol, MaxN_SlaveDofs_All*dispBT,
                                 N_SendDofs, recevbuf_KcolBT, N_RecevDofs, 1);
  for(i=0; i<N_Neibs; i++)
    N_RecevDofs[i] /= dispBT;
// =====================================================================================
// now assemble the dist mat A and BT
  if(N_Neibs)
   {
    delete [] N_SendDofs;
    delete [] sendbuf[0];
    delete [] sendbuf_Kcol[0];

    delete [] sendbuf;
    delete [] sendbuf_Kcol;
   }


  // assemble system sidt matrix directly
  DistMat_Entries = new double[N_DistMatEntries];
  memset(DistMat_Entries, 0, N_DistMatEntries*SizeOfDouble);

  Dist3U = 3*N_distU;
  Dist2U = 2*N_distU;

  index_dep = 0;
  // first add entries at Master dof from Neibs
  for(i=0;i<N_Neibs;i++)
   {
    N = N_RecevDofs[i];
    Neib_ID = NeibsRank[i];
    for(j=0;j<N;j++)
     {
      M = i*MaxN_SlaveDofs_All + j;
      dof = LocalDofOfRecbuf[0][M];

//       if(dof==-1) // Dirichlet row
//        {
//         continue;
// //  if(rank==TDatabase::ParamDB->Par_P0)
//  printf("Rank  %d main  \n",rank );
//  MPI_Finalize();
//  exit(0);
//        }

      N_UCol = LocalDofOfRecbuf[1][M];
      N_PCol = LocalDofOfRecbuf[2][M];

      distdof = DistDofofLocDof_U[dof];

      // A
      begin = DistA_RowPtr[distdof];
      len_A = DistA_RowPtr[distdof+1] - begin;
      KColBegin = DistA_KCol+begin;

      //system mat
      u1_rowbegin = DistMat_RowPtr[distdof];
      u2_rowbegin = DistMat_RowPtr[N_distU + distdof];
      u3_rowbegin = DistMat_RowPtr[Dist2U + distdof];

      u1row_Entries = DistMat_Entries+u1_rowbegin;
      u2row_Entries = DistMat_Entries+u2_rowbegin;
      u3row_Entries = DistMat_Entries+u3_rowbegin;

      pos = M*dispA;
      for(l=0; l<N_UCol; l++)
       {
        col = recevbuf_KcolA[0][pos];
        Aii_Deptpos[index_dep] = GetIndex(KColBegin, len_A, col); // find the col pos in the dist Kcol

        u1row_Entries[Aii_Deptpos[index_dep]] += recevbuf[0][pos]; // A11
        u2row_Entries[Aii_Deptpos[index_dep]] += recevbuf[0][pos]; // A22
        u3row_Entries[Aii_Deptpos[index_dep]] += recevbuf[0][pos]; // A33
        index_dep++;
        pos++;
       }

      if(distdof<N_distActiveU)
       {
        // B1T, B2T, B3T
        begin = DistBT_RowPtr[distdof];
        len_BT = DistBT_RowPtr[distdof+1] - begin;
        KColBegin = DistBT_KCol+begin;

        //system mat
        u1_rowbegin = DistMat_RowPtr[distdof] + len_A;
        u2_rowbegin = DistMat_RowPtr[N_distU + distdof] + len_A;
        u3_rowbegin = DistMat_RowPtr[Dist2U + distdof] + len_A;

        u1row_Entries = DistMat_Entries+u1_rowbegin;
        u2row_Entries = DistMat_Entries+u2_rowbegin;
        u3row_Entries = DistMat_Entries+u3_rowbegin;

        pos = M*dispBT;
        for(l=0; l<N_PCol; l++)
         {
          col = recevbuf_KcolBT[0][pos];
          index = GetIndex(KColBegin, len_BT, col); // find the col pos in the dist Kcol

//   if(rank==TDatabase::ParamDB->Par_P5)
//     printf("Mumps Rank %d  index  %d val  %f \n", rank,index, recebufBT[0][pos]);
          u1row_Entries[index] += recebufBT[0][pos]; // B1T
          u2row_Entries[index] += recebufBT[1][pos]; // B2T
          u3row_Entries[index] += recebufBT[2][pos]; // B2T
          pos++;
         }
        } // if(distdof<N_distActiveU)
     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)


  if(N_Neibs)
   {
   delete [] recevbuf[0];
   delete [] recebufBT[0];
   delete [] recebufBT[1];
   delete [] recebufBT[2];
   delete [] recevbuf_KcolA[0];
   delete [] recevbuf_KcolBT[0];
   delete [] LocalDofOfRecbuf[0];
   delete [] LocalDofOfRecbuf[1];
   delete [] LocalDofOfRecbuf[2];

   delete [] recevbuf;
   delete [] recebufBT;
   delete [] N_RecevDofs;
   delete [] recevbuf_KcolA;
   delete [] recevbuf_KcolBT;
   delete [] LocalDofOfRecbuf;
  }

  // add master values of B1, B2, B3 matrices
  if(N_Neibs_P)
   {
    N_SendDofs =  new int[N_Neibs_P];
    N_RecevDofs =  new int[N_Neibs_P];

    memset(N_RecevDofs, 0, N_Neibs_P*SizeOfInt);
    memset(N_SendDofs, 0, N_Neibs_P*SizeOfInt);

    sendbuf = new double*[3];
    recevbuf = new double*[3];

    LocalDofOfRecbuf_send = new int*[2];
    LocalDofOfRecbuf = new int*[2];

    sendbuf_Kcol = new int*[1];
    recevbuf_KcolA = new int*[1];

    // N_columns of B1, B2, B3
    dispB = MaxN_NzInDepRow_All_B;
    for(i=0;i<3;i++)
     {
      sendbuf[i] = new double[N_Neibs_P*MaxN_SlaveDofs_All_P*dispB];
      recevbuf[i] = new double[N_Neibs_P*MaxN_SlaveDofs_All_P*dispB];
     }

   LocalDofOfRecbuf_send[0]  = new int[N_Neibs_P*MaxN_SlaveDofs_All_P];
   LocalDofOfRecbuf_send[1]  = new int[N_Neibs_P*MaxN_SlaveDofs_All_P];

   LocalDofOfRecbuf[0]  = new int[N_Neibs_P*MaxN_SlaveDofs_All_P];
   LocalDofOfRecbuf[1]  = new int[N_Neibs_P*MaxN_SlaveDofs_All_P];

    sendbuf_Kcol[0] = new int[N_Neibs_P*MaxN_SlaveDofs_All_P*dispB];
    recevbuf_KcolA[0] = new int[N_Neibs_P*MaxN_SlaveDofs_All_P*dispB];

    for(i=0;i<N_DeptDofs_P;i++)
     {
      Master_ID = DeptDofMasterID_P[i];
      N = N_DeptDofNeibs_P[i];

      if(Master_ID==rank)
       {
        for(j=0;j<N;j++)
         {
          S = i*MaxSubDomainPerDof_P + j;
          Neib_ID = DeptDofNeibRanks_P[S];
          k = IndexOfNeibRank_P[Neib_ID];
          N_RecevDofs[k]++;
         } // for(j=0;j<N;j++)
       }
      else 
       {
        for(j=0;j<N;j++)
         {
          S = i*MaxSubDomainPerDof_P + j;
          Neib_ID = DeptDofNeibRanks_P[S];
          if(Neib_ID==Master_ID)
           break;
         }

        dof = DeptDofs_P[i];
        k = IndexOfNeibRank_P[Master_ID];
        M = k*MaxN_SlaveDofs_All_P + N_SendDofs[k];

        begin = RowPtrB[dof];
        end = RowPtrB[dof+1];

        LocalDofOfRecbuf_send[0][M] = NeibLocDof_P[S];
        LocalDofOfRecbuf_send[1][M] = end - begin;

        M *= dispB;

        for(l=begin;l<end;l++)
         {
          sendbuf[0][M] =  EntriesB1[l];
          sendbuf[1][M] =  EntriesB2[l];
          sendbuf[2][M] =  EntriesB3[l];
          sendbuf_Kcol[0][M] = GlobalDofOfLocalDof[KColB[l]];
//         printf("Rank %d  dof %d  Neibdof %d col %d\n",  rank, dof, NeibLocDof_P[S], sendbuf_Kcol[0][M]);
          M++;
         }

        N_SendDofs[k]++;
       }
     } //  for(i=0;i<N_DeptDofs;i++)
   } // if(N_Neibs_P)

// =====================================================================================
//   communicate 
  P_FEComm->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_SlaveDofs_All_P, N_SendDofs,
                              LocalDofOfRecbuf,  N_RecevDofs, 2);

  for(i=0; i<N_Neibs_P; i++)
   {
    N_SendDofs[i] *= dispB;
    N_RecevDofs[i] *= dispB;
   }
  P_FEComm->FECommunicateNeib(sendbuf, MaxN_SlaveDofs_All_P*dispB,
                                 N_SendDofs, recevbuf, N_RecevDofs, 3);
  P_FEComm->FECommunicateNeib(sendbuf_Kcol, MaxN_SlaveDofs_All_P*dispB,
                                 N_SendDofs, recevbuf_KcolA, N_RecevDofs, 1);
  for(i=0; i<N_Neibs_P; i++)
   {
    N_RecevDofs[i] /= dispB;
   }
// =====================================================================================


 if(N_Neibs_P)
  {
   for(i=0;i<3;i++)
    delete [] sendbuf[i];
   delete [] LocalDofOfRecbuf_send[0];
   delete [] LocalDofOfRecbuf_send[1];
   delete [] sendbuf_Kcol[0];

   delete [] N_SendDofs;
   delete [] sendbuf;
   delete [] LocalDofOfRecbuf_send;
   delete [] sendbuf_Kcol;
  } // if(N_Neibs_P)

  //  add entries to Master dof rows from Neibs
  for(i=0;i<N_Neibs_P;i++)
   {
    N = N_RecevDofs[i];

    for(j=0;j<N;j++)
     {
      M = i*MaxN_SlaveDofs_All_P + j;
      dof = LocalDofOfRecbuf[0][M];
      N_PCol = LocalDofOfRecbuf[1][M];

      distdof = DistDofofLocDof_P[dof];
      prow_Entries = DistMat_Entries + DistMat_RowPtr[Dist3U + distdof];

      if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE && distdof==Pressure_Const_DistDof)
        continue;

      // B1, B2, B3
      begin = DistB_RowPtr[distdof];
      len_B = DistB_RowPtr[distdof+1] - begin;
      KColBegin = DistB_KCol+begin;

      M *= dispB;
      for(k=0; k<N_PCol; k++)
       {
        col = recevbuf_KcolA[0][M];
        index = GetIndex(KColBegin, len_B, col);

        prow_Entries[index] += recevbuf[0][M];
        prow_Entries[len_B + index] += recevbuf[1][M];
        prow_Entries[2*len_B + index] += recevbuf[2][M];
        M++;
       }

     } // for(j=0;j<N;j++)
   } //

 if(N_Neibs_P)
  {

   for(i=0;i<3;i++)
    delete [] recevbuf[i];

   delete [] LocalDofOfRecbuf[0];
   delete [] LocalDofOfRecbuf[1];
   delete [] recevbuf_KcolA[0];

   delete [] N_RecevDofs;
   delete [] recevbuf;
   delete [] LocalDofOfRecbuf;
   delete [] recevbuf_KcolA;
  } // if(N_Neibs_P)

  index_pos = 0;
  // now copy own values
  for(i=0; i<nu; i++)
   {
    distdof=DistDofofLocDof_U[i];

    if(distdof==-1) continue; // this row not belongs to the dist mat

    u1row_Entries = DistMat_Entries+DistMat_RowPtr[distdof]; //u1 row
    u2row_Entries = DistMat_Entries+DistMat_RowPtr[N_distU + distdof]; //u2 row
    u3row_Entries = DistMat_Entries+DistMat_RowPtr[Dist2U + distdof]; //u3 row

    begin = DistA_RowPtr[distdof];
    len_A = DistA_RowPtr[distdof+1] - begin;
    KColBegin = DistA_KCol+begin;

    end =  RowPtrA[i+1];
//   if(rank==TDatabase::ParamDB->Par_P5)
//       printf("Rank %d  urow %d GlobRow %d  distdof %d %f COL:",  rank, i,  GlobalDofOfLocalDof[i], DistDofofLocDof_U[i], u1row_Entries[0] );


    for(j=RowPtrA[i]; j<end; j++)
     {
      col = GlobalDofOfLocalDof[KColA[j]];
      Aii_pos[index_pos] = GetIndex(KColBegin, len_A, col); // store the index for non-linear iteration

      u1row_Entries[Aii_pos[index_pos]] +=  EntriesA[j]; // A11
      u2row_Entries[Aii_pos[index_pos]] +=  EntriesA[j]; // A22
      u3row_Entries[Aii_pos[index_pos]] +=  EntriesA[j]; // A33
      index_pos++;      
//   if(rank==TDatabase::ParamDB->Par_P5)
//         printf("%d: %d:  %f:%f:%f  ", KColA[j], DistMat_RowPtr[distdof]+index,
//             EntriesA[j], u1row_Entries[index], DistMat_Entries[DistMat_RowPtr[distdof]+index]);


     } // for(j=begin; j<end; j++)
//   if(rank==TDatabase::ParamDB->Par_P5)
//      printf(" \n");

    if(distdof>=N_distActiveU)
      continue;



    //system mat
    u1row_Entries = u1row_Entries+len_A; // push A11 col
    u2row_Entries = u2row_Entries+len_A; // push A22 col
    u3row_Entries = u3row_Entries+len_A; // push A33 col

//   if(rank==TDatabase::ParamDB->Par_P5)
//       printf("Rank %d  i %d   distdof %d  N_distActiveU %d  \n",  rank, i,  distdof, N_distActiveU);

    // B1T, B2T, B3T
    begin = DistBT_RowPtr[distdof];
    len_BT = DistBT_RowPtr[distdof+1] - begin;
    KColBegin = DistBT_KCol+begin;

    end =  RowPtrBT[i+1];

//   if(rank==TDatabase::ParamDB->Par_P5)
//       printf("Rank %d  urow %d GlobRow %d  distdof %d %f COL:",  rank, i,  GlobalDofOfLocalDof[i], DistDofofLocDof_U[i], u1row_Entries[0] );



    for(j=RowPtrBT[i]; j<end; j++)
     {
      col = GlobalDofOfLocalDof_P[KColBT[j]];
      index = GetIndex(KColBegin, len_BT, col);

      u1row_Entries[index] +=  EntriesB1T[j];
      u2row_Entries[index] +=  EntriesB2T[j];
      u3row_Entries[index] +=  EntriesB3T[j];
     } // for(j=begin; j<end; j++)
   } //for(i=0; i<nu; i++

  for(i=0; i<np; i++)
   {
    distdof=DistDofofLocDof_P[i];

    if(distdof==-1) continue; // this row not belongs to the dist mat

    prow_Entries = DistMat_Entries + DistMat_RowPtr[Dist3U + distdof];

    begin = DistB_RowPtr[distdof];
    len_B = DistB_RowPtr[distdof+1] - begin;
    KColBegin = DistB_KCol+begin;

    end =  RowPtrB[i+1];

    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE && distdof==Pressure_Const_DistDof)
     {
      prow_Entries[0] = 1.;
      Rhs[3*nu] = 0.;
      continue;
     }

    for(j=RowPtrB[i]; j<end; j++)
     {
      col = GlobalDofOfLocalDof[KColB[j]];
      index = GetIndex(KColBegin, len_B, col);

      prow_Entries[index] +=  EntriesB1[j];
      prow_Entries[len_B + index] +=  EntriesB2[j];
      prow_Entries[2*len_B + index] +=  EntriesB3[j];

// // //   if(rank==TDatabase::ParamDB->Par_P5)
// // //         printf("%d: %d:  %f:%f:%f  ", KColB[j], DistMat_RowPtr[distdof]+index,
// // //             EntriesB1[j], prow_Entries[index], DistMat_Entries[DistMat_RowPtr[Dist3U + distdof]+index]);
     } // for(j=begin; j<end; j++)

//   if(rank==TDatabase::ParamDB->Par_P5)
//      printf(" \n");

   } //for(i=0; i<np; i++

//   if(rank==TDatabase::ParamDB->Par_P5)
//    {
//      // A mat orig
//     for(i=0;i<nu;i++)
//     {
//      if(DistDofofLocDof_U[i]==-1)
//       continue;
// 
//       printf("Rank %d  urow %d GlobRow %d  distdof %d ",  rank, i,  GlobalDofOfLocalDof[i], DistDofofLocDof_U[i] );
// 
//       for(j=RowPtrA[i]; j<RowPtrA[i+1]; j++)
//        {
//         printf("%d:%f  ", KColA[j], EntriesA[j]);
//         pos++;
//        }
//      printf(" \n");
//     }


//     for(i=0;i<N_distU;i++)
//     {
// 
//     if(i<N_distActiveU)
//       continue;
// 
// 
//       pos=DistMat_RowPtr[i];
//       printf("Rank %d  u1row %d    ",  rank, DistMat_Irn[pos]  );
// 
//       for(j=DistMat_RowPtr[i]; j<DistMat_RowPtr[i+1]; j++)
//        {
//         printf("%d:%f  ", DistMat_Jcn[pos], DistMat_Entries[pos]);
//         pos++;
//        }
//      printf(" \n");
//     }
// 
//      printf(" \n");
//     for(i=0;i<N_distU;i++)
//     {
// 
//      if(i<N_distActiveU)
//       continue;
// 
// 
//       pos=DistMat_RowPtr[N_distU+i];
//       printf("Rank %d  u2row %d     ",  rank, DistMat_Irn[pos]  );
// 
//       for(j=DistMat_RowPtr[N_distU+i]; j<DistMat_RowPtr[N_distU+i+1]; j++)
//        {
//         printf("%d:%f  ", DistMat_Jcn[pos], DistMat_Entries[pos]);
//         pos++;
//        }
// 
//      printf(" \n");
//     }
// 
//      printf(" \n");
//    for(i=0;i<N_distU;i++)
//     {
//      if(i<N_distActiveU)
//       continue;
// 
//       pos=DistMat_RowPtr[2*N_distU+i];
//       printf("Rank %d  u3row %d    ",  rank,  DistMat_Irn[pos]   );
// 
//       for(j=DistMat_RowPtr[2*N_distU+i]; j<DistMat_RowPtr[2*N_distU+i+1]; j++)
//        {
//         printf("%d:%f  ", DistMat_Jcn[pos], DistMat_Entries[pos]);
//         pos++;
//        }
// 
//      printf(" \n");
//     }

//      printf(" \n");
//    for(i=0;i<N_distP;i++)
//     {
// //      if(GlobalDofOfLocalDof_P[i] != 0)
// //       continue;
// 
//       pos=DistMat_RowPtr[3*N_distU+i];
//       printf("Rank %d DistDof %d pos %d ",  rank,  DistDofofLocDof_P[i], pos );
// 
//       for(j=DistMat_RowPtr[3*N_distU+i]; j<DistMat_RowPtr[3*N_distU+i+1]; j++)
//        {
//         printf("%d:%f  ", DistMat_Jcn[pos], DistMat_Entries[pos]);
//         pos++;
//        }
// 
//      printf(" \n");
//     }


//    } //  if(rank==TDatabase::ParamDB->Par_P5)

// MPI_Finalize();
// exit(0);

}

/** init in such a way that Process 0: u1_dist0, u2_dist0, p_dist0; 
  Process 1: u1_dist1, u2_dist1, p_dist1 and so on */ 
void TNSE_ParSolver2::InitMumps()
{
 int i, j, k, l, I1, I2, I3, J, J1, J2, J3, step, end;
 int pos1, pos2, pos3, index, rank, size, U_BlockN_Entries;
 int u1_rowbegin, u2_rowbegin, u3_rowbegin, p_rowbegin;
 int *distu_coldisp, *distp_coldisp;
 int *distu_range, *distp_range;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
 int N_GlobDof_P, *GlobalDofOfLocalDof_P, *Dist_Dof_All_P;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Velo_FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  P_FEComm->GetGlobalDistDofInfo(N_GlobDof_P, GlobalDofOfLocalDof_P, Dist_Dof_All_P);

  N_Eqns = 3*Global_N_U + Global_N_P;
  u1_rowbegin = 3*Global_beginU + Global_beginP + 1; // + 1 => fortran style
  u2_rowbegin = u1_rowbegin + N_distU;
  u3_rowbegin = u2_rowbegin + N_distU;
  p_rowbegin = u3_rowbegin + N_distU;

  distu_coldisp = new int[size+1];
  distp_coldisp = new int[size+1];
  distu_coldisp[0] = 0;
  distp_coldisp[0] = 0;

  for(i=0;i<size;i++)
   {
    distu_coldisp[i+1] = distu_coldisp[i] + 2*Dist_Dof_All[i] + Dist_Dof_All_P[i];
    distp_coldisp[i+1] = distp_coldisp[i] + 3*Dist_Dof_All[i];
   }

  distu_range = new int[size+1];
  distp_range = new int[size+1];

  distu_range[0] = 0;
  distp_range[0] = 0;

  for(i=0;i<size;i++)
   {
    distu_range[i+1] =  distu_range[i] + Dist_Dof_All[i];
    distp_range[i+1] =  distp_range[i] + Dist_Dof_All_P[i];
   }

  N_DistMatEqns = 3*N_distU + N_distP;  
  DistMat_RowPtr = new int[N_DistMatEqns+1];

  N_DistMatEntries = 3*DistA_RowPtr[N_distU] + 3*DistB_RowPtr[N_distP]
                     + 3*DistBT_RowPtr[N_distActiveU];
  DistMat_Irn = new int[N_DistMatEntries];
  DistMat_Jcn = new int[N_DistMatEntries];
  Aii_pos = new int[DistA_RowPtr[N_distU]];
  Aii_Deptpos = new int[DistA_RowPtr[N_distU]];
  
  U_BlockN_Entries = DistA_RowPtr[N_distU] + DistBT_RowPtr[N_distActiveU];

  index = 0;
  pos1 = 0;
  pos2 = U_BlockN_Entries;
  pos3 = 2*U_BlockN_Entries;
  DistMat_RowPtr[0] = 0;

  for(i=0;i<N_distU;i++)
   {
    I1 = u1_rowbegin+i;
    I2 = u2_rowbegin+i;
    I3 = u3_rowbegin+i;

    end = DistA_RowPtr[i+1];
    for(j=DistA_RowPtr[i]; j<end; j++)
     {
      J = DistA_KCol[j];

      // find col belongs to which rank, and make disp accordingly
      for(k=0;k<size;k++)
       if(J>=distu_range[k]  && J<distu_range[k+1] )
        {
         J1 = J + distu_coldisp[k];
         J2 = J + Dist_Dof_All[k] + distu_coldisp[k];
         J3 = J + 2*Dist_Dof_All[k] + distu_coldisp[k];
         break;
        }

      DistMat_Irn[pos1] = I1;
      DistMat_Irn[pos2] = I2;
      DistMat_Irn[pos3] = I3;

      DistMat_Jcn[pos1] = J1+1; // fortran style
      DistMat_Jcn[pos2] = J2+1; // fortran style
      DistMat_Jcn[pos3] = J3+1; // fortran style

//       Aii_pos[index++] = pos1;

      pos1++;
      pos2++;
      pos3++;
     }


    if(i<N_distActiveU)
     {
      end = DistBT_RowPtr[i+1];
      for(j=DistBT_RowPtr[i]; j<end; j++)
       {
        J = DistBT_KCol[j];

        // find col belongs to which rank, and make disp accordingly
        for(k=0;k<size;k++)
         if(J>=distp_range[k]  && J<distp_range[k+1] )
          {
           J += (3*Dist_Dof_All[k] + distp_coldisp[k]);
           break;
          }

        DistMat_Irn[pos1] = I1;
        DistMat_Irn[pos2] = I2;
        DistMat_Irn[pos3] = I3;

        DistMat_Jcn[pos1] = J + 1; // fortran style
        DistMat_Jcn[pos2] = J + 1; // fortran style
        DistMat_Jcn[pos3] = J + 1; // fortran style

        pos1++;
        pos2++;
        pos3++;
       }
      }

      DistMat_RowPtr[i+1] = pos1;
      DistMat_RowPtr[N_distU + i+1] = pos2;
      DistMat_RowPtr[2*N_distU + i+1] = pos3;
     } // for(i=0;i<N_distU;i++)


  pos3 = 3*U_BlockN_Entries;
  for(i=0;i<N_distP;i++)
   {
     I3  =  p_rowbegin +  i;

    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE && i==Pressure_Const_DistDof)
     {
     // printf("%d MUMPS Rank %d DistDof %d pos %d \n", TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE, rank,  i, pos3 );
      N_DistMatEntries -=3*(DistB_RowPtr[i+1] - DistB_RowPtr[i]) - 1;
      DistMat_Irn[pos3] = I3;
      DistMat_Jcn[pos3] = I3;
      pos3++;
     }
    else
     {
      end = DistB_RowPtr[i+1];
      step = end - DistB_RowPtr[i];
      l = pos3;
      for(j=DistB_RowPtr[i];j<end;j++)
       {
        J = DistB_KCol[j];

        // find col belongs to which rank, and make disp accordingly
        for(k=0;k<size;k++)
         if(J>=distu_range[k]  && J<distu_range[k+1] )
          {
           J1 = J + distu_coldisp[k];
           J2 = J + Dist_Dof_All[k] + distu_coldisp[k];
           J3 = J + 2*Dist_Dof_All[k] + distu_coldisp[k];
           break;
          }

        DistMat_Irn[l] = I3; //B1
        DistMat_Irn[step + l] = I3; //B2
        DistMat_Irn[2*step + l] = I3; //B3

        DistMat_Jcn[l] = J1+1;
        DistMat_Jcn[step + l] = J2+1;
        DistMat_Jcn[2*step + l] = J3+1;
        l++;
       } //  for(j=DistB_RowPtr[i];j<

      pos3 +=3*step;
     }

    DistMat_RowPtr[3*N_distU + i+1] = pos3;
   } // for(i=0;i<N_distP;i++)

//   if(rank==TDatabase::ParamDB->Par_P5)
//    {
// 
// //     printf("Mumps Rank %d  u1_rowbegin %d u2_rowbegin %d u3_rowbegin %d \n", rank, u1_rowbegin,
// //             u2_rowbegin, u3_rowbegin);
// 
// 
//     pos1=0;
//     for(i=0;i<N_distU;i++)
//     {
//       printf("Rank %d  u1row %d    ",  rank, DistMat_Irn[pos1]  );
// 
//       for(j=DistMat_RowPtr[i]; j<DistMat_RowPtr[i+1]; j++)
//        printf("%d, ", DistMat_Jcn[pos1++]);
// 
//      printf(" \n");
//     }
//      printf(" \n");
//     for(i=0;i<N_distU;i++)
//     {
//       printf("Rank %d  u2row %d     ",  rank, DistMat_Irn[pos1]  );
// 
//       for(j=DistMat_RowPtr[N_distU+i]; j<DistMat_RowPtr[N_distU+i+1]; j++)
//        printf("%d, ", DistMat_Jcn[pos1++]);
// 
//      printf(" \n");
//     }
// 
//      printf(" \n");
//    for(i=0;i<N_distU;i++)
//     {
//       printf("Rank %d  u3row %d    ",  rank,  DistMat_Irn[pos1]   );
// 
//       for(j=DistMat_RowPtr[2*N_distU+i]; j<DistMat_RowPtr[2*N_distU+i+1]; j++)
//        printf("%d, ", DistMat_Jcn[pos1++]);
// 
//      printf(" \n");
//     }
// 
//      printf(" \n");
//    for(i=0;i<N_distP;i++)
//     {
//       printf("Rank %d  prow %d    ",  rank,   DistMat_Irn[pos1]  );
// 
//       for(j=DistMat_RowPtr[3*N_distU+i]; j<DistMat_RowPtr[3*N_distU+i+1]; j++)
//        printf("%d, ", DistMat_Jcn[pos1++]);
// 
//      printf(" \n");
//     }
//    }

//  MPI_Allreduce(&N_DistMatEntries, &pos3, 1, MPI_INT, MPI_SUM, Comm);
// 
//  if(rank==0)
//   printf("Rank %d Mumps Inti  %d\n", rank, pos3);
 
//    printf("Rank %d Mumps Inti  %d\n", rank, N_DistMatEntries);
 
// MPI_Finalize();
// exit(0);


 delete [] distu_range;
 delete [] distu_coldisp;

 delete [] distp_range;
 delete [] distp_coldisp;

 MUMPS_Solver = new TMumpsSolver(N_Eqns, N_DistMatEntries, DistMat_Irn, DistMat_Jcn, 1);
 
} // InitMumps


/** init in such a way that Process 0: u1_dist0, u2_dist0, p_dist0; 
  Process 1: u1_dist1, u2_dist1, p_dist1 and so on */ 
void TNSE_ParSolver2::InitSuperLU()
{
 int i, j, k, l, I1, I2, I3, J, J1, J2, J3, step, end;
 int pos1, pos2, pos3, index, rank, size, U_BlockN_Entries;
 int u1_rowbegin, u2_rowbegin, u3_rowbegin, p_rowbegin;
 int *distu_coldisp, *distp_coldisp;
 int *distu_range, *distp_range;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
 int N_GlobDof_P, *GlobalDofOfLocalDof_P, *Dist_Dof_All_P;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Velo_FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  P_FEComm->GetGlobalDistDofInfo(N_GlobDof_P, GlobalDofOfLocalDof_P, Dist_Dof_All_P);

  N_Eqns = 3*Global_N_U + Global_N_P;
  u1_rowbegin = 3*Global_beginU + Global_beginP;
  u2_rowbegin = u1_rowbegin + N_distU;
  u3_rowbegin = u2_rowbegin + N_distU;
  p_rowbegin = u3_rowbegin + N_distU;

  distu_coldisp = new int[size+1];
  distp_coldisp = new int[size+1];
  distu_coldisp[0] = 0;
  distp_coldisp[0] = 0;

  for(i=0;i<size;i++)
   {
    distu_coldisp[i+1] = distu_coldisp[i] + 2*Dist_Dof_All[i] + Dist_Dof_All_P[i];
    distp_coldisp[i+1] = distp_coldisp[i] + 3*Dist_Dof_All[i];
   }

  distu_range = new int[size+1];
  distp_range = new int[size+1];

  distu_range[0] = 0;
  distp_range[0] = 0;

  for(i=0;i<size;i++)
   {
    distu_range[i+1] =  distu_range[i] + Dist_Dof_All[i];
    distp_range[i+1] =  distp_range[i] + Dist_Dof_All_P[i];
   }

  N_DistMatEqns = 3*N_distU + N_distP;  
  DistMat_RowPtr = new int[N_DistMatEqns+1];

  N_DistMatEntries = 3*DistA_RowPtr[N_distU] + 3*DistB_RowPtr[N_distP]
                     + 3*DistBT_RowPtr[N_distActiveU];
  DistMat_Irn = new int[N_DistMatEntries];
  DistMat_Jcn = new int[N_DistMatEntries];
  Aii_pos = new int[DistA_RowPtr[N_distU]];
  Aii_Deptpos = new int[DistA_RowPtr[N_distU]];
  
  U_BlockN_Entries = DistA_RowPtr[N_distU] + DistBT_RowPtr[N_distActiveU];

  index = 0;
  pos1 = 0;
  pos2 = U_BlockN_Entries;
  pos3 = 2*U_BlockN_Entries;
  DistMat_RowPtr[0] = 0;

  for(i=0;i<N_distU;i++)
   {
    I1 = u1_rowbegin+i;
    I2 = u2_rowbegin+i;
    I3 = u3_rowbegin+i;

    end = DistA_RowPtr[i+1];
    for(j=DistA_RowPtr[i]; j<end; j++)
     {
      J = DistA_KCol[j];

      // find col belongs to which rank, and make disp accordingly
      for(k=0;k<size;k++)
       if(J>=distu_range[k]  && J<distu_range[k+1] )
        {
         J1 = J + distu_coldisp[k];
         J2 = J + Dist_Dof_All[k] + distu_coldisp[k];
         J3 = J + 2*Dist_Dof_All[k] + distu_coldisp[k];
         break;
        }

      DistMat_Irn[pos1] = I1;
      DistMat_Irn[pos2] = I2;
      DistMat_Irn[pos3] = I3;

      DistMat_Jcn[pos1] = J1;
      DistMat_Jcn[pos2] = J2;
      DistMat_Jcn[pos3] = J3;

      pos1++;
      pos2++;
      pos3++;
     }


    if(i<N_distActiveU)
     {
      end = DistBT_RowPtr[i+1];
      for(j=DistBT_RowPtr[i]; j<end; j++)
       {
        J = DistBT_KCol[j];

        // find col belongs to which rank, and make disp accordingly
        for(k=0;k<size;k++)
         if(J>=distp_range[k]  && J<distp_range[k+1] )
          {
           J += (3*Dist_Dof_All[k] + distp_coldisp[k]);
           break;
          }

        DistMat_Irn[pos1] = I1;
        DistMat_Irn[pos2] = I2;
        DistMat_Irn[pos3] = I3;

        DistMat_Jcn[pos1] = J;
        DistMat_Jcn[pos2] = J;
        DistMat_Jcn[pos3] = J;

        pos1++;
        pos2++;
        pos3++;
       }
      }

      DistMat_RowPtr[i+1] = pos1;
      DistMat_RowPtr[N_distU + i+1] = pos2;
      DistMat_RowPtr[2*N_distU + i+1] = pos3;
     } // for(i=0;i<N_distU;i++)


  pos3 = 3*U_BlockN_Entries;
  for(i=0;i<N_distP;i++)
   {
     I3  =  p_rowbegin +  i;

    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE && i==Pressure_Const_DistDof)
     {
      // printf("MUMPS Rank %d DistDof %d pos %d ",  rank,  i, pos3 );
      N_DistMatEntries -=3*(DistB_RowPtr[i+1] - DistB_RowPtr[i]) - 1;
      DistMat_Irn[pos3] = I3;
      DistMat_Jcn[pos3] = I3;
      pos3++;
     }
    else
     {
      end = DistB_RowPtr[i+1];
      step = end - DistB_RowPtr[i];
      l = pos3;
      for(j=DistB_RowPtr[i];j<end;j++)
       {
        J = DistB_KCol[j];

        // find col belongs to which rank, and make disp accordingly
        for(k=0;k<size;k++)
         if(J>=distu_range[k]  && J<distu_range[k+1] )
          {
           J1 = J + distu_coldisp[k];
           J2 = J + Dist_Dof_All[k] + distu_coldisp[k];
           J3 = J + 2*Dist_Dof_All[k] + distu_coldisp[k];
           break;
          }

        DistMat_Irn[l] = I3; //B1
        DistMat_Irn[step + l] = I3; //B2
        DistMat_Irn[2*step + l] = I3; //B3

        DistMat_Jcn[l] = J1;
        DistMat_Jcn[step + l] = J2;
        DistMat_Jcn[2*step + l] = J3;
        l++;
       } //  for(j=DistB_RowPtr[i];j<

      pos3 +=3*step;
     }

    DistMat_RowPtr[3*N_distU + i+1] = pos3;
   } // for(i=0;i<N_distP;i++)

//   if(rank==TDatabase::ParamDB->Par_P5)
//    {
// 
// //     printf("Mumps Rank %d  u1_rowbegin %d u2_rowbegin %d u3_rowbegin %d \n", rank, u1_rowbegin,
// //             u2_rowbegin, u3_rowbegin);
// 
// 
//     pos1=0;
//     for(i=0;i<N_distU;i++)
//     {
//       printf("Rank %d  u1row %d    ",  rank, DistMat_Irn[pos1]  );
// 
//       for(j=DistMat_RowPtr[i]; j<DistMat_RowPtr[i+1]; j++)
//        printf("%d, ", DistMat_Jcn[pos1++]);
// 
//      printf(" \n");
//     }
//      printf(" \n");
//     for(i=0;i<N_distU;i++)
//     {
//       printf("Rank %d  u2row %d     ",  rank, DistMat_Irn[pos1]  );
// 
//       for(j=DistMat_RowPtr[N_distU+i]; j<DistMat_RowPtr[N_distU+i+1]; j++)
//        printf("%d, ", DistMat_Jcn[pos1++]);
// 
//      printf(" \n");
//     }
// 
//      printf(" \n");
//    for(i=0;i<N_distU;i++)
//     {
//       printf("Rank %d  u3row %d    ",  rank,  DistMat_Irn[pos1]   );
// 
//       for(j=DistMat_RowPtr[2*N_distU+i]; j<DistMat_RowPtr[2*N_distU+i+1]; j++)
//        printf("%d, ", DistMat_Jcn[pos1++]);
// 
//      printf(" \n");
//     }
// 
//      printf(" \n");
//    for(i=0;i<N_distP;i++)
//     {
//       printf("Rank %d  prow %d    ",  rank,   DistMat_Irn[pos1]  );
// 
//       for(j=DistMat_RowPtr[3*N_distU+i]; j<DistMat_RowPtr[3*N_distU+i+1]; j++)
//        printf("%d, ", DistMat_Jcn[pos1++]);
// 
//      printf(" \n");
//     }
//    }
// printf("Mumps Inti  \n");
// MPI_Finalize();
// exit(0);


 delete [] distu_range;
 delete [] distu_coldisp;

 delete [] distp_range;
 delete [] distp_coldisp;

#ifdef _SUPERLU_
 SuperLU_Solver = new TSuperLUSolver(N_Eqns, N_DistMatEqns, N_DistMatEntries, 
                                     DistMat_RowPtr, DistMat_Jcn, u1_rowbegin);
#else
printf("SuperLU_Solver cannot be called! change the SOLVER flag in makefile  \n");
MPI_Finalize();
exit(0);
#endif
} // InitSuperLU




/** init in such a way that Process 0: u1_dist0, u2_dist0, p_dist0; 
  Process 1: u1_dist1, u2_dist1, p_dist1 and so on */ 

void TNSE_ParSolver2::InitHips()
{
 int i, j, k, l, I1, I2, I3, J, J1, J2, J3, step, end;
 int pos1, pos2, pos3, index, rank, size, U_BlockN_Entries;
 int u1_rowbegin, u2_rowbegin, u3_rowbegin, p_rowbegin;
 int *distu_coldisp, *distp_coldisp;
 int *distu_range, *distp_range;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
 int N_GlobDof_P, *GlobalDofOfLocalDof_P, *Dist_Dof_All_P, *GlobalEqIndex;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Velo_FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  P_FEComm->GetGlobalDistDofInfo(N_GlobDof_P, GlobalDofOfLocalDof_P, Dist_Dof_All_P);

  N_Eqns = 3*Global_N_U + Global_N_P;
  u1_rowbegin = 3*Global_beginU + Global_beginP;
  u2_rowbegin = u1_rowbegin + N_distU;
  u3_rowbegin = u2_rowbegin + N_distU;
  p_rowbegin = u3_rowbegin + N_distU;

  distu_coldisp = new int[size+1];
  distp_coldisp = new int[size+1];
  distu_coldisp[0] = 0;
  distp_coldisp[0] = 0;

  for(i=0;i<size;i++)
   {
    distu_coldisp[i+1] = distu_coldisp[i] + 2*Dist_Dof_All[i] + Dist_Dof_All_P[i];
    distp_coldisp[i+1] = distp_coldisp[i] + 3*Dist_Dof_All[i];
   }

  distu_range = new int[size+1];
  distp_range = new int[size+1];

  distu_range[0] = 0;
  distp_range[0] = 0;

  for(i=0;i<size;i++)
   {
    distu_range[i+1] =  distu_range[i] + Dist_Dof_All[i];
    distp_range[i+1] =  distp_range[i] + Dist_Dof_All_P[i];
   }

  N_DistMatEqns = 3*N_distU + N_distP;
  DistMat_RowPtr = new int[N_DistMatEqns+1];

  N_DistMatEntries = 3*DistA_RowPtr[N_distU] + 3*DistB_RowPtr[N_distP]
                     + 3*DistBT_RowPtr[N_distActiveU];
//   DistMat_Irn = new int[N_DistMatEqns+1];
  DistMat_Jcn = new int[N_DistMatEntries];
  Aii_pos = new int[DistA_RowPtr[N_distU]];
  Aii_Deptpos = new int[DistA_RowPtr[N_distU]];
  
  U_BlockN_Entries = DistA_RowPtr[N_distU] + DistBT_RowPtr[N_distActiveU];
  GlobalEqIndex = new int[N_DistMatEqns];  
  
  index = 0;
  pos1 = 0;
  pos2 = U_BlockN_Entries;
  pos3 = 2*U_BlockN_Entries;
  DistMat_RowPtr[0] = 0;
//   DistMat_Irn[0] = 0;

  for(i=0;i<N_distU;i++)
   {
    I1 = u1_rowbegin+i;
    I2 = u2_rowbegin+i;
    I3 = u3_rowbegin+i; 

    end = DistA_RowPtr[i+1];
    for(j=DistA_RowPtr[i]; j<end; j++)
     {
      J = DistA_KCol[j];

      // find col belongs to which rank, and make disp accordingly
      for(k=0;k<size;k++)
       if(J>=distu_range[k]  && J<distu_range[k+1] )
        {
         J1 = J + distu_coldisp[k];
         J2 = J + Dist_Dof_All[k] + distu_coldisp[k];
         J3 = J + 2*Dist_Dof_All[k] + distu_coldisp[k];
         break;
        }

      DistMat_Jcn[pos1] = J1;
      DistMat_Jcn[pos2] = J2;
      DistMat_Jcn[pos3] = J3;

      pos1++;
      pos2++;
      pos3++;
     }


    if(i<N_distActiveU)
     {
      end = DistBT_RowPtr[i+1];
      for(j=DistBT_RowPtr[i]; j<end; j++)
       {
        J = DistBT_KCol[j];

        // find col belongs to which rank, and make disp accordingly
        for(k=0;k<size;k++)
         if(J>=distp_range[k]  && J<distp_range[k+1] )
          {
           J += (3*Dist_Dof_All[k] + distp_coldisp[k]);
           break;
          }

        DistMat_Jcn[pos1] = J;
        DistMat_Jcn[pos2] = J;
        DistMat_Jcn[pos3] = J;

        pos1++;
        pos2++;
        pos3++;
       }
      }

      DistMat_RowPtr[i+1] = pos1;
      DistMat_RowPtr[N_distU + i+1] = pos2;
      DistMat_RowPtr[2*N_distU + i+1] = pos3;

//       DistMat_Irn[i+1] = pos1;
//       DistMat_Irn[N_distU + i+1] = pos2;
//       DistMat_Irn[2*N_distU + i+1] = pos3;

      GlobalEqIndex[i] = I1;
      GlobalEqIndex[N_distU + i] = I2;
      GlobalEqIndex[2*N_distU + i] = I3;
     } // for(i=0;i<N_distU;i++)


  pos3 = 3*U_BlockN_Entries;
  for(i=0;i<N_distP;i++)
   {
     I3  =  p_rowbegin +  i;

    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE && i==Pressure_Const_DistDof)
     {
      // printf("MUMPS Rank %d DistDof %d pos %d ",  rank,  i, pos3 );
      N_DistMatEntries -=3*(DistB_RowPtr[i+1] - DistB_RowPtr[i]) - 1;
      DistMat_Jcn[pos3] = I3;
      pos3++;
     }
    else
     {
      end = DistB_RowPtr[i+1];
      step = end - DistB_RowPtr[i];
      l = pos3;
      for(j=DistB_RowPtr[i];j<end;j++)
       {
        J = DistB_KCol[j];

        // find col belongs to which rank, and make disp accordingly
        for(k=0;k<size;k++)
         if(J>=distu_range[k]  && J<distu_range[k+1] )
          {
           J1 = J + distu_coldisp[k];
           J2 = J + Dist_Dof_All[k] + distu_coldisp[k];
           J3 = J + 2*Dist_Dof_All[k] + distu_coldisp[k];
           break;
          }

        DistMat_Jcn[l] = J1;
        DistMat_Jcn[step + l] = J2;
        DistMat_Jcn[2*step + l] = J3;
        l++;
       } //  for(j=DistB_RowPtr[i];j<

      pos3 +=3*step;
     }

    DistMat_RowPtr[3*N_distU + i+1] = pos3;
//     DistMat_Irn[3*N_distU + i+1] = pos3;
    GlobalEqIndex[3*N_distU + i] = I3;
   } // for(i=0;i<N_distP;i++)

 delete [] distu_range;
 delete [] distu_coldisp;

 delete [] distp_range;
 delete [] distp_coldisp;
 
#ifdef _HIPS_
int *mapptr, *mapp;
   if(rank==0)
    {
     mapptr = new int[size+1];
     mapp = new int[N_Eqns];
     mapptr[0] = 0;
     for(i=0;i<size;i++)
       mapptr[i+1] = 3*Dist_Dof_All[i] + Dist_Dof_All_P[i];

    for(i=0;i<N_Eqns;i++)
      mapp[i] = i;
    }
    
  Hips_Solver = new THipsSolver(N_Eqns, N_DistMatEqns, N_DistMatEntries, DistMat_RowPtr,
                      DistMat_Jcn, mapptr, mapp, GlobalEqIndex);
     
  if(rank==0)
   {  
    delete [] mapptr;
    delete [] mapp;
   }     
   
#else
printf("Hips cannot be called! change the SOLVER flag in makefile  \n");
MPI_Finalize();
exit(0);
#endif
  delete [] GlobalEqIndex;
 

} // InitHips()



void TNSE_ParSolver2::InitpARMS()
{

  
 int i, j, k, l, I1, I2, I3, J, J1, J2, J3, step, end;
 int pos1, pos2, pos3, index, rank, size, U_BlockN_Entries;
 int u1_rowbegin, u2_rowbegin, u3_rowbegin, p_rowbegin;
 int *distu_coldisp, *distp_coldisp;
 int *distu_range, *distp_range;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
 int N_GlobDof_P, *GlobalDofOfLocalDof_P, *Dist_Dof_All_P, *GlobalEqIndex;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Velo_FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  P_FEComm->GetGlobalDistDofInfo(N_GlobDof_P, GlobalDofOfLocalDof_P, Dist_Dof_All_P);

  N_Eqns = 3*Global_N_U + Global_N_P;
  u1_rowbegin = 3*Global_beginU + Global_beginP;
  u2_rowbegin = u1_rowbegin + N_distU;
  u3_rowbegin = u2_rowbegin + N_distU;
  p_rowbegin = u3_rowbegin + N_distU;

  distu_coldisp = new int[size+1];
  distp_coldisp = new int[size+1];
  distu_coldisp[0] = 0;
  distp_coldisp[0] = 0;

  for(i=0;i<size;i++)
   {
    distu_coldisp[i+1] = distu_coldisp[i] + 2*Dist_Dof_All[i] + Dist_Dof_All_P[i];
    distp_coldisp[i+1] = distp_coldisp[i] + 3*Dist_Dof_All[i];
   }

  distu_range = new int[size+1];
  distp_range = new int[size+1];

  distu_range[0] = 0;
  distp_range[0] = 0;

  for(i=0;i<size;i++)
   {
    distu_range[i+1] =  distu_range[i] + Dist_Dof_All[i];
    distp_range[i+1] =  distp_range[i] + Dist_Dof_All_P[i];
   }

  N_DistMatEqns = 3*N_distU + N_distP;  
  DistMat_RowPtr = new int[N_DistMatEqns+1];

  N_DistMatEntries = 3*DistA_RowPtr[N_distU] + 3*DistB_RowPtr[N_distP]
                     + 3*DistBT_RowPtr[N_distActiveU];
//   DistMat_Irn = new int[N_DistMatEqns+1];
  DistMat_Jcn = new int[N_DistMatEntries];
  Aii_pos = new int[DistA_RowPtr[N_distU]];
  Aii_Deptpos = new int[DistA_RowPtr[N_distU]];
  
  U_BlockN_Entries = DistA_RowPtr[N_distU] + DistBT_RowPtr[N_distActiveU];
  GlobalEqIndex = new int[N_DistMatEqns];  
  
  index = 0;
  pos1 = 0;
  pos2 = U_BlockN_Entries;
  pos3 = 2*U_BlockN_Entries;
  DistMat_RowPtr[0] = 0;
//   DistMat_Irn[0] = 0;
      
  for(i=0;i<N_distU;i++)
   {
    I1 = u1_rowbegin+i;
    I2 = u2_rowbegin+i;
    I3 = u3_rowbegin+i; 
    
    end = DistA_RowPtr[i+1];
    for(j=DistA_RowPtr[i]; j<end; j++)
     {
      J = DistA_KCol[j];

      // find col belongs to which rank, and make disp accordingly
      for(k=0;k<size;k++)
       if(J>=distu_range[k]  && J<distu_range[k+1] )
        {
         J1 = J + distu_coldisp[k];
         J2 = J + Dist_Dof_All[k] + distu_coldisp[k];
         J3 = J + 2*Dist_Dof_All[k] + distu_coldisp[k];
         break;
        }

      DistMat_Jcn[pos1] = J1;
      DistMat_Jcn[pos2] = J2;
      DistMat_Jcn[pos3] = J3;

      pos1++;
      pos2++;
      pos3++;
     }


    if(i<N_distActiveU)
     {
      end = DistBT_RowPtr[i+1];
      for(j=DistBT_RowPtr[i]; j<end; j++)
       {
        J = DistBT_KCol[j];

        // find col belongs to which rank, and make disp accordingly
        for(k=0;k<size;k++)
         if(J>=distp_range[k]  && J<distp_range[k+1] )
          {
           J += (3*Dist_Dof_All[k] + distp_coldisp[k]);
           break;
          }

        DistMat_Jcn[pos1] = J;
        DistMat_Jcn[pos2] = J;
        DistMat_Jcn[pos3] = J;

        pos1++;
        pos2++;
        pos3++;
       }
      }

      DistMat_RowPtr[i+1] = pos1;
      DistMat_RowPtr[N_distU + i+1] = pos2;
      DistMat_RowPtr[2*N_distU + i+1] = pos3;

//       DistMat_Irn[i+1] = pos1;
//       DistMat_Irn[N_distU + i+1] = pos2;
//       DistMat_Irn[2*N_distU + i+1] = pos3;

      GlobalEqIndex[i] = I1;
      GlobalEqIndex[N_distU + i] = I2;
      GlobalEqIndex[2*N_distU + i] = I3;
     } // for(i=0;i<N_distU;i++)


  pos3 = 3*U_BlockN_Entries;
  for(i=0;i<N_distP;i++)
   {
     I3  =  p_rowbegin +  i;

    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE && i==Pressure_Const_DistDof)
     {
      // printf("MUMPS Rank %d DistDof %d pos %d ",  rank,  i, pos3 );
      N_DistMatEntries -=3*(DistB_RowPtr[i+1] - DistB_RowPtr[i]) - 1;
      DistMat_Jcn[pos3] = I3;
      pos3++;
     }
    else
     {
      end = DistB_RowPtr[i+1];
      step = end - DistB_RowPtr[i];
      l = pos3;
      for(j=DistB_RowPtr[i];j<end;j++)
       {
        J = DistB_KCol[j];

        // find col belongs to which rank, and make disp accordingly
        for(k=0;k<size;k++)
         if(J>=distu_range[k]  && J<distu_range[k+1] )
          {
           J1 = J + distu_coldisp[k];
           J2 = J + Dist_Dof_All[k] + distu_coldisp[k];
           J3 = J + 2*Dist_Dof_All[k] + distu_coldisp[k];
           break;
          }

        DistMat_Jcn[l] = J1;
        DistMat_Jcn[step + l] = J2;
        DistMat_Jcn[2*step + l] = J3;
        l++;
       } //  for(j=DistB_RowPtr[i];j<

      pos3 +=3*step;
     }

    DistMat_RowPtr[3*N_distU + i+1] = pos3;
//     DistMat_Irn[3*N_distU + i+1] = pos3;
    GlobalEqIndex[3*N_distU + i] = I3;
   } // for(i=0;i<N_distP;i++)

 delete [] distu_range;
 delete [] distu_coldisp;

 delete [] distp_range;
 delete [] distp_coldisp;

#ifdef _PARMS_
  Parms_Solver = new TParmsSolver(N_Eqns, N_DistMatEqns, N_DistMatEntries, DistMat_RowPtr,
                      DistMat_Jcn, N_distU, N_distP, Dist_Dof_All, Dist_Dof_All_P, GlobalEqIndex); 
#else
printf("Parms cannot be called! change the SOLVER flag in makefile  \n");
MPI_Finalize();
exit(0);
#endif
  delete [] GlobalEqIndex;
} // InitpARMS()




















#endif
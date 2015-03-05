// =======================================================================
// @(#)NSE_ParSolver.h
//
// Class:      TNSE_ParSolver
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

#include <NSE_ParSolver.h>
#include <Database.h>
#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>
#include <ParVectorNSE3D.h>

    /** constructor */
TNSE_ParSolver::TNSE_ParSolver(TParFECommunicator3D *velo_FEComm, 
        TSquareStructure3D *squareStructure, TParFECommunicator3D *p_FECom, 
        TStructure3D * structureBT, TStructure3D *structureB, int type)
 {
   Comm = TDatabase::ParamDB->Comm;

   Velo_FEComm = velo_FEComm;
   P_FEComm = p_FECom;

   Global_N_U = Velo_FEComm->GetN_GlobalDegreesOfFreedom();
   Global_N_P = P_FEComm->GetN_GlobalDegreesOfFreedom();

   SquareStructure = squareStructure;
   Structure = structureB;
   StructureT = structureBT;

   Pressure_Const_DistDof = -1;

   // NSType
   if(type==2 || type==4)
    SetUpDistMat();

   MUMPS_Solver = NULL;
//    Parms_Solver = NULL;
 }


/** generate all data needed for Square Matrix */
void TNSE_ParSolver::SetUpDistMat()
{
 int i, j, k, l,  M, N, P, S, T, dof, Neib_ID, rank, size;
 int nu, *RowPtrA, *KColA, N_Active, *KColBT, *RowPtrBT, MaxN_NzInRow;
 int MaxN_NzInRow_All_A, MaxN_NzInRow_All_BT;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *DeptDofs, *N_DeptDofNeibs;
 int MaxN_DeptDofs_All, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocDof;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All;
 int MaxN_MasterDofs_All, MaxN_SlaveDofs_All, *DeptDofMasterID, *N_DofRankIndex;

 int MaxSubDomainPerDof_P, N_Neibs_P, *NeibsRank_P, N_DeptDofs_P, *DeptDofs_P, *N_DeptDofNeibs_P;  
 int MaxN_DeptDofs_All_P, *DeptDofNeibRanks_P, *IndexOfNeibRank_P, *LocDOFofNeib_P, *NeibLocDof_P;
 int N_GlobDof_P, *GlobalDofOfLocalDof_P, *Dist_Dof_All_P;
 int MaxN_MasterDofs_All_P, MaxN_SlaveDofs_All_P, *DeptDofMasterID_P, *N_DofRankIndex_P;

 int *N_SendDofs, *N_RecevDofs, **sendbuf, **recevbuf, disp, Master_ID;
 int *tmp_RowPtrA, *tmp_RowPtrBT, tmp, tmp1;
 int **LocalDofOfRecbuf, **LocalDofOfRecbuf_send;
 int *aux1, *aux2, *aux1_loc, *aux2_loc, N_UCol, N_PCol, last;
 int *tmp_KColA, *tmp_KColBT, *tmp_RowPtrB;

 int np, *RowPtrB, *KColB, MaxN_NzInRow_All_B, *tmp_KColB;

 double x, y, z;

  MPI_Comm_rank(Comm, &rank); 
  MPI_Comm_size(Comm, &size);

  nu = SquareStructure->GetN_Rows();
  RowPtrA = SquareStructure->GetRowPtr();
  KColA = SquareStructure->GetKCol();
  N_Active = SquareStructure->GetActiveBound();

  KColBT = StructureT->GetKCol();
  RowPtrBT = StructureT->GetRowPtr();

  Velo_FEComm->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs,
                     DeptDofs, N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks,
                     IndexOfNeibRank, LocDOFofNeib, NeibLocDof);

  Velo_FEComm->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  Velo_FEComm->GetDeptDofMasterInfo(MaxN_MasterDofs_All, MaxN_SlaveDofs_All, DeptDofMasterID);

  N_DofRankIndex = Velo_FEComm->GetN_DofRankIndex();

  Global_beginU = 0;
  for(i=0;i<rank;i++)
   Global_beginU += Dist_Dof_All[i];

  MaxN_NzInRow = -1;
  for(i=0;i<N_DeptDofs;i++)
   {
    dof = DeptDofs[i];
    if(MaxN_NzInRow < (RowPtrA[dof+1] - RowPtrA[dof] )  )
     MaxN_NzInRow = RowPtrA[dof+1] - RowPtrA[dof];
   }

  MPI_Allreduce(&MaxN_NzInRow, &MaxN_NzInDepRow_All_A, 1, MPI_INT, MPI_MAX, Comm);
  // printf("Rank %d MaxN_NzInDepRow_All_A %d  \n", rank, MaxN_NzInDepRow_All_A);


  P_FEComm->GetDofNeibInfo(MaxSubDomainPerDof_P, N_Neibs_P, NeibsRank_P, N_DeptDofs_P, DeptDofs_P, 
                     N_DeptDofNeibs_P, MaxN_DeptDofs_All_P, DeptDofNeibRanks_P, IndexOfNeibRank_P,
                     LocDOFofNeib_P, NeibLocDof_P);

  P_FEComm->GetGlobalDistDofInfo(N_GlobDof_P, GlobalDofOfLocalDof_P, Dist_Dof_All_P);
  P_FEComm->GetDeptDofMasterInfo(MaxN_MasterDofs_All_P, MaxN_SlaveDofs_All_P, DeptDofMasterID_P);

  N_DofRankIndex_P = P_FEComm->GetN_DofRankIndex();

  Global_beginP = 0;
  for(i=0;i<rank;i++)
   Global_beginP += Dist_Dof_All_P[i];

  MaxN_NzInRow = -1;
  for(i=0;i<N_DeptDofs;i++)
   {
    dof = DeptDofs[i];
    if(MaxN_NzInRow < (RowPtrBT[dof+1] - RowPtrBT[dof] )  )
     MaxN_NzInRow = RowPtrBT[dof+1] - RowPtrBT[dof];
   }

  MPI_Allreduce(&MaxN_NzInRow, &MaxN_NzInDepRow_All_BT, 1, MPI_INT, MPI_MAX, Comm);

  LocalDofOfRecbuf = new int*[1];
  LocalDofOfRecbuf_send = new int*[1];

  /** first generate data for A & BT matrices, then for B */
  if(N_Neibs)
   {
    N_SendDofs =  new int[N_Neibs];
    N_RecevDofs =  new int[N_Neibs];
    memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    sendbuf = new int*[1];
    recevbuf = new int*[1];

    // 1st & 2nd columns to store N_columns of U & P
    disp = MaxN_NzInDepRow_All_A*MaxN_NzInDepRow_All_BT+2;

    sendbuf[0] = new int[N_Neibs*MaxN_SlaveDofs_All*disp];
    recevbuf[0]  = new int[N_Neibs*MaxN_SlaveDofs_All*disp];

    LocalDofOfRecbuf[0] = new int[N_Neibs*MaxN_SlaveDofs_All];
    LocalDofOfRecbuf_send[0]  = new int[N_Neibs*MaxN_SlaveDofs_All];

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

        //1st & 2nd columns to store N_columns of U & P
        //further are global colum index of U and P
        dof = DeptDofs[i];
        M *=disp;

        sendbuf[0][M++] = RowPtrA[dof+1] - RowPtrA[dof]; // number of columns
        sendbuf[0][M++] = RowPtrBT[dof+1] - RowPtrBT[dof]; // number of columns

        for(l=RowPtrA[dof];l<RowPtrA[dof+1];l++)
          sendbuf[0][M++] = GlobalDofOfLocalDof[KColA[l]];

        for(l= RowPtrBT[dof];l<RowPtrBT[dof+1];l++)
         sendbuf[0][M++] = GlobalDofOfLocalDof_P[KColBT[l]];

        N_SendDofs[k]++;
       }
     } //  for(i=0;i<N_DeptDofs;i++)
   } // if(N_Neibs)
// ==============================================================
//   communicate
  Velo_FEComm->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_SlaveDofs_All, N_SendDofs,
                                 LocalDofOfRecbuf,  N_RecevDofs, 1);

  for(i=0; i<N_Neibs; i++)
   {
    N_SendDofs[i] *= disp;
    N_RecevDofs[i] *= disp;
   }
  Velo_FEComm->FECommunicateNeib(sendbuf, MaxN_SlaveDofs_All*disp,
                                 N_SendDofs, recevbuf, N_RecevDofs, 1);
  for(i=0; i<N_Neibs; i++)
    N_RecevDofs[i] /= disp;

  if(N_Neibs)
   {
    delete [] N_SendDofs;
    delete [] sendbuf[0];
    delete [] sendbuf;
   }

// ==============================================================
// start to fill all necessary info, first velo
// =================================================================
  tmp_RowPtrA = new int[nu];
  tmp_RowPtrBT = new int[nu];

  memset(tmp_RowPtrA, 0, nu*SizeOfInt);
  memset(tmp_RowPtrBT, 0, nu*SizeOfInt);

  //find total number of entries in the dist system matrix
  DistA_N_Entries = 0;
  DistBT_N_Entries = 0;
 
  // first fill self dofs for first veo component
  // total No. velo Comp. Entries = 3*DistA_N_Entries;
  for(i=0;i<nu;i++)
   if(N_DofRankIndex[i]==1)
    {
     tmp = DistA_N_Entries;
     tmp1 = DistBT_N_Entries;

     DistA_N_Entries +=  (RowPtrA[i+1] - RowPtrA[i]);
     DistBT_N_Entries +=  (RowPtrBT[i+1] - RowPtrBT[i]);

     tmp_RowPtrA[i] += DistA_N_Entries - tmp;
     tmp_RowPtrBT[i] += DistBT_N_Entries - tmp1;
    }

   // fill master dofs
   for(i=0;i<N_DeptDofs;i++)
    if(DeptDofMasterID[i]==rank)
     {
      dof = DeptDofs[i];
      tmp = DistA_N_Entries;
      tmp1 = DistBT_N_Entries;

      DistA_N_Entries += (RowPtrA[dof+1] - RowPtrA[dof]);
      DistBT_N_Entries += (RowPtrBT[dof+1] - RowPtrBT[dof]);

      tmp_RowPtrA[dof] += DistA_N_Entries - tmp;
      tmp_RowPtrBT[dof] += DistBT_N_Entries - tmp1;
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
      DistBT_N_Entries += recevbuf[0][M+1];

      tmp_RowPtrA[dof] += recevbuf[0][M];
      tmp_RowPtrBT[dof] += recevbuf[0][M+1];
     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)

  // fill col with Global dofs
  MaxN_NzInRow_All_A = -1;
  for(i=0;i<nu;i++)
   if(MaxN_NzInRow_All_A <  tmp_RowPtrA[i]   )
     MaxN_NzInRow_All_A = tmp_RowPtrA[i];

  MaxN_NzInRow_All_BT = -1;
  for(i=0;i<nu;i++)
   if(MaxN_NzInRow_All_BT <  tmp_RowPtrBT[i]   )
     MaxN_NzInRow_All_BT = tmp_RowPtrBT[i];

   aux1= new int[nu*MaxN_NzInRow_All_A];
   aux2= new int[nu*MaxN_NzInRow_All_BT];

   memset(aux1, 0, nu*MaxN_NzInRow_All_A*SizeOfInt);
   memset(aux2, 0, nu*MaxN_NzInRow_All_BT*SizeOfInt);

   memset(tmp_RowPtrA, 0, nu*SizeOfInt);
   memset(tmp_RowPtrBT, 0, nu*SizeOfInt);

   // now put column info
   for(i=0;i<nu;i++)
    {
     aux1_loc = aux1 + i*MaxN_NzInRow_All_A;
     aux2_loc = aux2 + i*MaxN_NzInRow_All_BT;

     if(N_DofRankIndex[i]==1)
      {
       for(j=RowPtrA[i]; j<RowPtrA[i+1]; j++)
        aux1_loc[tmp_RowPtrA[i]++] = GlobalDofOfLocalDof[KColA[j]];

       for(j=RowPtrBT[i]; j<RowPtrBT[i+1]; j++)
        aux2_loc[tmp_RowPtrBT[i]++] = GlobalDofOfLocalDof_P[KColBT[j]];
    } //if(N_DofRankIndex[i]==1)
   }


   // fill master dofs own columns
   for(i=0;i<N_DeptDofs;i++)
    if(DeptDofMasterID[i]==rank)
     {
      dof = DeptDofs[i];
      aux1_loc = aux1 + dof*MaxN_NzInRow_All_A;
      aux2_loc = aux2 + dof*MaxN_NzInRow_All_BT;

      for(j=RowPtrA[dof]; j<RowPtrA[dof+1]; j++)
       aux1_loc[tmp_RowPtrA[dof]++] = GlobalDofOfLocalDof[KColA[j]];

      for(j=RowPtrBT[dof]; j<RowPtrBT[dof+1]; j++)
       aux2_loc[tmp_RowPtrBT[dof]++] = GlobalDofOfLocalDof_P[KColBT[j]];
      }

// now fill the master dof info from neibs
  for(i=0;i<N_Neibs;i++)
   {
    N = N_RecevDofs[i];
//     Neib_ID = NeibsRank[i];
    for(j=0;j<N;j++)
     {
      M = i*MaxN_SlaveDofs_All + j;
      dof = LocalDofOfRecbuf[0][M];
      M *=disp;

      aux1_loc = aux1 + dof*MaxN_NzInRow_All_A;
      aux2_loc = aux2 + dof*MaxN_NzInRow_All_BT;

      N_UCol = recevbuf[0][M++];
      N_PCol = recevbuf[0][M++];

      for(l=0; l<N_UCol; l++)
        aux1_loc[tmp_RowPtrA[dof]++] = recevbuf[0][M++];
//            if(rank==1 &&  dof==1072 && Neib_ID==0)
//             printf("Rank %d Master slave dof %d  KColA %d N_SendDofs[k] %d\n", rank,   dof,    recevbuf[0][M-1] ,  N); 
//        }

      for(l=0; l<N_PCol; l++)
       aux2_loc[tmp_RowPtrBT[dof]++] = recevbuf[0][M++];

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

   // process A mat
   tmp_KColA = new int[DistA_N_Entries];

   N_distU = Dist_Dof_All[rank];        /** init class variable */ 
   DistA_RowPtr = new int[N_distU+1];
   DistDofofLocDof_U = new int[nu];
   memset(DistDofofLocDof_U, -1, nu*SizeOfInt);

   //find total number of entries in the dist system matrix
   DistA_N_Entries = 0;                    /** init class variable */
   DistA_RowPtr[0] = 0;                    /** init class variable */
   N_distActiveU = 0;
   N=0;

   for(i=0;i<nu;i++)
    {
     aux1_loc = aux1 + i*MaxN_NzInRow_All_A;
     M = tmp_RowPtrA[i];

     if(M>1)
      {
       Sort(aux1_loc, M);
      }

     last=-1;
     for(j=0;j<M;j++)
       if(aux1_loc[j] != last)
        {
         last = aux1_loc[j];
         tmp_KColA[DistA_N_Entries] = aux1_loc[j];
         DistA_N_Entries++;
        }

     if(last!=-1)
      {
       DistDofofLocDof_U[i] = N;
       DistA_RowPtr[++N] = DistA_N_Entries;

//       if( (i== 1072) &&  rank==1)
//        {
//         printf("Rank %d dof %d N_Active %d distdof %d A  distlen  %d len %d \n",
//                rank, i, N_Active, N-1,  DistA_RowPtr[N] - DistA_RowPtr[N-1], RowPtrA[i+1]-RowPtrA[i] );
//  
//      for(j=DistA_RowPtr[N-1];j<DistA_RowPtr[N];j++)
//            printf("Rank %d index %d dof %d tmp_KColA %d \n",
//                rank,  j-DistA_RowPtr[N-1], i,    tmp_KColA[j] );
// 
// 
//        }

       if(i<N_Active)
        N_distActiveU++;
      }

    } //  for(i=0;i<nu;i++)

  delete [] aux1;
  delete [] tmp_RowPtrA;

  DistA_KCol = new int[DistA_N_Entries]; 
  memcpy(DistA_KCol, tmp_KColA, DistA_N_Entries*SizeOfInt);

  delete [] tmp_KColA;

// MPI_Finalize();
// exit(0);


// //   print matrix
//  if(rank==TDatabase::ParamDB->Par_P5)
//   for(j=0;j<N_distU;j++)
//    {
//     printf("Row %d: ", Global_beginU + j);
// 
//     for(k=DistA_RowPtr[j];k<DistA_RowPtr[j+1];k++)
//       printf(" %d ",   DistA_KCol[k] );
//       printf(" \n");
// //      printf(" N_Columns %d \n", DistA_RowPtr[j+1] - DistA_RowPtr[j]);
//    }


   /** now process dist BT mat */
   tmp_KColBT = new int[DistBT_N_Entries];
   DistBT_RowPtr = new int[N_distU+1]; 

   DistBT_N_Entries = 0; 
   DistBT_RowPtr[0] = 0; 
   N=0;
   for(i=0;i<nu;i++)
    {
     aux2_loc = aux2 + i*MaxN_NzInRow_All_BT;
     M = tmp_RowPtrBT[i];

     if(M>1)
      Sort(aux2_loc, M);

     last=-1;
     for(j=0;j<M;j++)
       if(aux2_loc[j] != last)
        {
         last = aux2_loc[j];
         tmp_KColBT[DistBT_N_Entries] = aux2_loc[j];
         DistBT_N_Entries++;
        }

     if(last!=-1)
       DistBT_RowPtr[++N] = DistBT_N_Entries;
    }

  delete [] aux2;
  delete [] tmp_RowPtrBT;

  DistBT_KCol = new int[DistBT_N_Entries];
  memcpy(DistBT_KCol, tmp_KColBT, DistBT_N_Entries*SizeOfInt);

  delete [] tmp_KColBT;

// //   print matrix
//  if(rank==TDatabase::ParamDB->Par_P5)
// 
//   for(j=0;j<N_distU;j++)
//    {
//     printf("Row %d: ", Global_beginU + j);
// 
//     for(k=DistBT_RowPtr[j];k<DistBT_RowPtr[j+1];k++)
//       printf(" %d ",   DistBT_KCol[k] );
//       printf(" \n");
// //      printf(" N_Columns %d \n", DistA_RowPtr[j+1] - DistA_RowPtr[j]);
//    }


 /** Process B Mat start */

  np = Structure->GetN_Rows();
  RowPtrB = Structure->GetRowPtr();
  KColB = Structure->GetKCol();

  MaxN_NzInRow = -1;
  for(i=0;i<N_DeptDofs_P;i++)
   {
    dof = DeptDofs_P[i];
    if(MaxN_NzInRow < (RowPtrB[dof+1] - RowPtrB[dof] )  )
     MaxN_NzInRow = RowPtrB[dof+1] - RowPtrB[dof];
   }

  MPI_Allreduce(&MaxN_NzInRow, &MaxN_NzInDepRow_All_B, 1, MPI_INT, MPI_MAX, Comm);

  if(N_Neibs_P)
   {
    N_SendDofs =  new int[N_Neibs_P];
    N_RecevDofs =  new int[N_Neibs_P];
    memset(N_RecevDofs, 0, N_Neibs_P*SizeOfInt);
    memset(N_SendDofs, 0, N_Neibs_P*SizeOfInt);

    sendbuf = new int*[1];
    recevbuf = new int*[1];

    // 1st  column to store N_columns of U 
    disp = MaxN_NzInDepRow_All_B+1;

    sendbuf[0] = new int[N_Neibs_P*MaxN_SlaveDofs_All_P*disp];
    recevbuf[0]  = new int[N_Neibs_P*MaxN_SlaveDofs_All_P*disp];

    LocalDofOfRecbuf[0] = new int[N_Neibs_P*MaxN_SlaveDofs_All_P];
    LocalDofOfRecbuf_send[0] = new int[N_Neibs_P*MaxN_SlaveDofs_All_P];

    for(i=0;i<N_DeptDofs_P;i++)
     {
      N = N_DeptDofNeibs_P[i];
      Master_ID = DeptDofMasterID_P[i];

      if(Master_ID==rank) // master of this row, info will be recieved from all neibs
       {
        for(j=0;j<N;j++)
         {
          S = i*MaxSubDomainPerDof_P + j;
          Neib_ID = DeptDofNeibRanks_P[S];
          k = IndexOfNeibRank_P[Neib_ID]; // LocindexOfNeib_ID
          N_RecevDofs[k]++;
         } // for(j=0;j<N;j++)
       }
      else  //info from all neibs of this row have to be send to the master
       {
        // fill the neib's (i.e., Master's) local dof
        for(j=0;j<N;j++)
         {
          S = i*MaxSubDomainPerDof_P + j;
          Neib_ID = DeptDofNeibRanks_P[S];
          if(Neib_ID==Master_ID)
           break;
         }
        k = IndexOfNeibRank_P[Master_ID]; // LocindexOfNeib_ID
        M = k*MaxN_SlaveDofs_All_P + N_SendDofs[k];
        LocalDofOfRecbuf_send[0][M] = NeibLocDof_P[S];

        dof = DeptDofs_P[i];
        M *= disp;
        sendbuf[0][M++] = RowPtrB[dof+1] - RowPtrB[dof];

        for(l=RowPtrB[dof];l<RowPtrB[dof+1];l++)
         sendbuf[0][M++] = GlobalDofOfLocalDof[KColB[l]];

        N_SendDofs[k]++;
       }

     }
   } // if(N_Neibs_P)

// ==============================================================
//   communicate 

  P_FEComm->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_SlaveDofs_All_P, N_SendDofs,
                                 LocalDofOfRecbuf,  N_RecevDofs, 1);

  for(i=0; i<N_Neibs_P; i++)
   {
    N_SendDofs[i] *= disp;
    N_RecevDofs[i] *= disp;
   }
   P_FEComm->FECommunicateNeib(sendbuf, MaxN_SlaveDofs_All_P*disp,
                             N_SendDofs, recevbuf, N_RecevDofs, 1);

 if(N_Neibs_P)
  {
   for(i=0; i<N_Neibs_P; i++)
    N_RecevDofs[i] /= disp;

   delete [] N_SendDofs;
   delete [] sendbuf[0];
   delete [] sendbuf;
  }

  // ==============================================================
  // start to fill all necessary info, first pressure
  // =================================================================
  tmp_RowPtrB = new int[np];
  memset(tmp_RowPtrB, 0, np*SizeOfInt);

  DistB_N_Entries = 0;
  // fill own dofs
  for(i=0;i<np;i++)
   if(N_DofRankIndex_P[i]==1)
    {
     tmp = RowPtrB[i+1] - RowPtrB[i];
     DistB_N_Entries += tmp;
     tmp_RowPtrB[i] += tmp;
    }

   // fill master dofs
   for(i=0;i<N_DeptDofs_P;i++)
    if(DeptDofMasterID_P[i]==rank)
     {
      dof = DeptDofs_P[i];
      tmp = RowPtrB[dof+1] - RowPtrB[dof];
      DistB_N_Entries += tmp;
      tmp_RowPtrB[dof] += tmp;
     }

// now fill the master dof info from neibs
  for(i=0;i<N_Neibs_P;i++)
   {
    N = N_RecevDofs[i];
    Neib_ID = NeibsRank_P[i];
    for(j=0;j<N;j++)
     {
      M = i*MaxN_SlaveDofs_All_P + j;
      dof = LocalDofOfRecbuf[0][M];
      M *=disp;

      DistB_N_Entries +=  recevbuf[0][M];
      tmp_RowPtrB[dof] += recevbuf[0][M];
     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)

  MaxN_NzInRow_All_B = -1;
  for(i=0;i<np;i++)
   if(MaxN_NzInRow_All_B <  tmp_RowPtrB[i]   )
     MaxN_NzInRow_All_B = tmp_RowPtrB[i];

   aux1= new int[np*MaxN_NzInRow_All_B];

   memset(tmp_RowPtrB, 0, np*SizeOfInt);

  // fill own dofs
  for(i=0;i<np;i++)
   if(N_DofRankIndex_P[i]==1)
    {
     aux1_loc = aux1 + i*MaxN_NzInRow_All_B;

     for(j=RowPtrB[i]; j<RowPtrB[i+1]; j++)
      aux1_loc[tmp_RowPtrB[i]++] = GlobalDofOfLocalDof[KColB[j]];
    }

   // fill master dofs
   for(i=0;i<N_DeptDofs_P;i++)
    if(DeptDofMasterID_P[i]==rank)
     {
      dof = DeptDofs_P[i];
      aux1_loc = aux1 + dof*MaxN_NzInRow_All_B;

     for(j=RowPtrB[dof]; j<RowPtrB[dof+1]; j++)
      aux1_loc[tmp_RowPtrB[dof]++] = GlobalDofOfLocalDof[KColB[j]];
     }

// now fill the master dof info from neibs
  for(i=0;i<N_Neibs_P;i++)
   {
    N = N_RecevDofs[i];
    Neib_ID = NeibsRank_P[i];

    for(j=0;j<N;j++)
     {
      M = i*MaxN_SlaveDofs_All_P + j;
      dof = LocalDofOfRecbuf[0][M];
      M *=disp;

      aux1_loc = aux1 + dof*MaxN_NzInRow_All_B;
      N_UCol = recevbuf[0][M++];

      for(l=0; l<N_UCol; l++)
       aux1_loc[tmp_RowPtrB[dof]++] = recevbuf[0][M++];
     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)


  // neib infi nolonger needed free the memory
  if(N_Neibs_P)
   {
    delete [] N_RecevDofs;
    delete [] recevbuf[0];
    delete [] recevbuf;
    delete [] LocalDofOfRecbuf[0];
    delete [] LocalDofOfRecbuf_send[0];
   }

   // process B mat
   tmp_KColB = new int[DistB_N_Entries];

   N_distP = Dist_Dof_All_P[rank]; 
   DistB_RowPtr = new int[N_distP+1];
   DistDofofLocDof_P = new int[np];
   memset(DistDofofLocDof_P, -1, np*SizeOfInt);

   //find total number of entries in the dist system matrix
   DistB_N_Entries = 0;
   DistB_RowPtr[0] = 0;
   N=0;
   for(i=0;i<np;i++)
    {
     aux1_loc = aux1 + i*MaxN_NzInRow_All_B;
     M = tmp_RowPtrB[i];

     if(M>1)
      Sort(aux1_loc, M);

     last=-1;
     for(j=0;j<M;j++)
       if(aux1_loc[j] != last)
        {
         last = aux1_loc[j];
         tmp_KColB[DistB_N_Entries++] = aux1_loc[j];
        }

     if(last!=-1)
      {
       DistDofofLocDof_P[i] = N;
       DistB_RowPtr[++N] = DistB_N_Entries;

       if(GlobalDofOfLocalDof_P[i]==0 &&  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        Pressure_Const_DistDof = i;
      }
    }

  delete [] aux1;
  delete [] tmp_RowPtrB;

  DistB_KCol = new int[DistB_N_Entries];
  memcpy(DistB_KCol, tmp_KColB, DistB_N_Entries*SizeOfInt);

  delete [] tmp_KColB;
  delete [] LocalDofOfRecbuf;
  delete [] LocalDofOfRecbuf_send;

// //   print matrix
//  if(rank==TDatabase::ParamDB->Par_P5)
// 
//   for(j=0;j<N_distP;j++)
//    {
//     printf("Row %d: ", Global_beginP + j);
// 
//     for(k=DistB_RowPtr[j];k<DistB_RowPtr[j+1];k++)
//       printf(" %d ",   DistB_KCol[k] );
//       printf(" \n");
// //      printf(" N_Columns %d \n", DistA_RowPtr[j+1] - DistA_RowPtr[j]);
//    }
// 


// // //  test
//  MPI_Allreduce(&DistA_N_Entries, &MaxN_NzInRow_All_A, 1, MPI_INT, MPI_SUM, Comm);  
//  MPI_Allreduce(&DistBT_N_Entries, &MaxN_NzInRow_All_BT, 1, MPI_INT, MPI_SUM, Comm);
//  MPI_Allreduce(&DistB_N_Entries, &MaxN_NzInRow_All_B, 1, MPI_INT, MPI_SUM, Comm);
// 
// 
// if(rank==0)
//  {
//   printf("Rank %d A N_Entries  %d \n",rank, MaxN_NzInRow_All_A);
//   printf("Rank %d BT N_Entries %d \n",rank, MaxN_NzInRow_All_BT);
//   printf("Rank %d B N_Entries %d \n",rank,  MaxN_NzInRow_All_B);
//  }
// 
// MPI_Finalize();
// exit(0);
}



void TNSE_ParSolver::GetDistArray(double *RHS, int nu, int np, double *tmp_rhs)
{
  int i, distdof;
  double *rhsu1, *rhsu2, *rhsu3, *tmp1;
  double *Globrhsu1, *Globrhsu2, *Globrhsu3, *tmp2;
  
    rhsu1 = tmp_rhs;
    rhsu2 = tmp_rhs+N_distU;
    rhsu3 = tmp_rhs+2*N_distU; 

    Globrhsu1 = RHS;
    Globrhsu2 = RHS+nu;
    Globrhsu3 = RHS+2*nu;
                    
    for(i=0; i<nu; i++)
     {
      distdof=DistDofofLocDof_U[i];
      if(distdof==-1) continue; // this row not belongs to the dist mat

      rhsu1[distdof] = Globrhsu1[i];
      rhsu2[distdof] = Globrhsu2[i];
      rhsu3[distdof] = Globrhsu3[i];
     }

    tmp1 = tmp_rhs+(3*N_distU);
    tmp2 = RHS+(3*nu);
                
    for(i=0; i<np; i++)
     {
      distdof=DistDofofLocDof_P[i];
      if(distdof==-1) continue; // this row not belongs to the dist mat

      tmp1[distdof] = tmp2[i];
     }
        
        
} // GetDistArray



void TNSE_ParSolver::GetLocArray(double *tmp_rhs, int nu, int np, double *RHS)
{
  int i, distdof;
  double *rhsu1, *rhsu2, *rhsu3, *tmp1;
  double *Globrhsu1, *Globrhsu2, *Globrhsu3, *tmp2;
  
    rhsu1 = tmp_rhs;
    rhsu2 = tmp_rhs+N_distU;
    rhsu3 = tmp_rhs+2*N_distU; 

    Globrhsu1 = RHS;
    Globrhsu2 = RHS+nu;
    Globrhsu3 = RHS+2*nu;
                    
    for(i=0; i<nu; i++)
     {
      distdof=DistDofofLocDof_U[i];
      if(distdof==-1) continue; // this row not belongs to the dist mat

      Globrhsu1[i] = rhsu1[distdof];
      Globrhsu2[i] = rhsu2[distdof];
      Globrhsu3[i] = rhsu3[distdof];
     }

    tmp1 = tmp_rhs+(3*N_distU);
    tmp2 = RHS+(3*nu);
                
    for(i=0; i<np; i++)
     {
      distdof=DistDofofLocDof_P[i];
      if(distdof==-1) continue; // this row not belongs to the dist mat

      tmp2[i] = tmp1[distdof];
     }
                
} // GetDistArray



TNSE_ParSolver::~TNSE_ParSolver()
{

}

#endif
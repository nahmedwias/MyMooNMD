// =======================================================================
// @(#)ParVectorNSE.C
//
// Class:      TParVectorNSE
// Purpose:    Class containing all info needed for communication between subdomains
//             for solution/rhs vectors of NSE 2D/3D
//
// Author:     Sashikumaar Ganesan (12.10.09)
//
// History:    Start of implementation 12.10.09 (Sashikumaar Ganesan)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Database.h>
#include <MumpsSolver.h>
#include <ParFECommunicator3D.h>
#include <ParVectorNSE3D.h>

/** constructor with given values of velocity and pressure*/
TParVectorNSE3D::TParVectorNSE3D(MPI_Comm comm, double *U, int nu, int np, int N_Dim,  TParFECommunicator3D *velocityCommunicator,
                       TParFECommunicator3D *pressureCommunicator) :
       TParVector(comm, U, N_Dim*nu+np, N_Dim)

{
 NU = nu;
 NP = np;

 if(velocityCommunicator)
  VelocityFESpace = velocityCommunicator->Getfespace();
 
 if(pressureCommunicator)
  PressureFESpace = pressureCommunicator->Getfespace();

 VelocityCommunicator = velocityCommunicator;
 PressureCommunicator = pressureCommunicator;
}

void TParVectorNSE3D::ParDdot(int assembletype, double &residual, double &impuls_residual)
{
 int i, N_OwnDof, *OwnDofs, m;
 double *u1, *u2, *u3, *p, tmp1, tmp2;
 
//  BYADD 0
//  BYMASTER 1
// BYOWN 2 
 
  if(assembletype==BYADD)
   { this->AssembleByADD(); }
  else if(assembletype==BYMASTER)
   {  this->AssembleWithMaster(); }
  else if(assembletype==BYOWN) 
    {
     // do nothing, since the vector is already in consistant form 
    }
  else
   {
    cout << "Wrong Assemble type !!!!!!!!!!!!! TParVectorNSE3D::ParDdot" << endl;
    MPI_Finalize();
    exit(0);       
   }
    
    
  VelocityCommunicator->GetOwnDofs(N_OwnDof, OwnDofs);

  m = N_OwnDof%3;

  u1 = Values;
  u2 = Values+NU;
  u3 = Values+2*NU;

  tmp1 = 0;
  for(i=0;i<m;i++)
   tmp1 += u1[OwnDofs[i]]*u1[OwnDofs[i]] + u2[OwnDofs[i]]*u2[OwnDofs[i]] +
           u3[OwnDofs[i]]*u3[OwnDofs[i]];

  for(i=m;i<N_OwnDof;i+=3)
   tmp1 += u1[OwnDofs[i]]*u1[OwnDofs[i]] + u2[OwnDofs[i]]*u2[OwnDofs[i]] +
           u3[OwnDofs[i]]*u3[OwnDofs[i]] +
           u1[OwnDofs[i+1]]*u1[OwnDofs[i+1]] + u2[OwnDofs[i+1]]*u2[OwnDofs[i+1]] +
           u3[OwnDofs[i+1]]*u3[OwnDofs[i+1]] +
           u1[OwnDofs[i+2]]*u1[OwnDofs[i+2]] + u2[OwnDofs[i+2]]*u2[OwnDofs[i+2]] +
           u3[OwnDofs[i+2]]*u3[OwnDofs[i+2]];

  PressureCommunicator->GetOwnDofs(N_OwnDof, OwnDofs);

  m = N_OwnDof%5;

  p = Values+3*NU;

  tmp2 = tmp1;
  for(i=0;i<m;i++)
   tmp2 += p[OwnDofs[i]]*p[OwnDofs[i]];

  for(i=m;i<N_OwnDof;i+=5)
   tmp2 += p[OwnDofs[i]]*p[OwnDofs[i]] + p[OwnDofs[i+1]]*p[OwnDofs[i+1]] +
           p[OwnDofs[i+2]]*p[OwnDofs[i+2]] + p[OwnDofs[i+3]]*p[OwnDofs[i+3]] +
           p[OwnDofs[i+4]]*p[OwnDofs[i+4]];

  impuls_residual = 0;
  residual = 0;
  MPI_Allreduce(&tmp1, &impuls_residual, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
  MPI_Allreduce(&tmp2, &residual, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);

}// ParDefectVect

void TParVectorNSE3D::AssembleByADD()
{
 int i, j, k, l, M, N, P, S, T, dof, Neib_ID, *N_SendDofs;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *N_DeptDofNeibs, MaxN_DeptDofs_All;
 int *DeptDofs, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocalDOF;
 
 double *sendbuf, *recevbuf, *pressure;


 // first velocity
 VelocityCommunicator->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs, DeptDofs, 
                     N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks, IndexOfNeibRank,
                     LocDOFofNeib, NeibLocalDOF);

  if(N_Neibs)
   {
    N_SendDofs = new int[N_Neibs];
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    sendbuf = new double[N_Neibs*MaxN_DeptDofs_All*N_DIM];
    recevbuf = new double[N_Neibs*MaxN_DeptDofs_All*N_DIM];

    for(i=0;i<N_DeptDofs;i++)
     {
      N = N_DeptDofNeibs[i];
      for(j=0;j<N;j++)
       {
        S = i*MaxSubDomainPerDof + j;
        Neib_ID = DeptDofNeibRanks[S];
        k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
        T = k*MaxN_DeptDofs_All;

        M = N_DIM*T + N_SendDofs[k];
        for(l=0;l<N_DIM;l++)
         {
          P = M + l*MaxN_DeptDofs_All;
          sendbuf[P]= Values[l*NU + DeptDofs[i]];
         }
        N_SendDofs[k]++;
       } // for(j=0;j<N;j++)
     }// for(i=0;i<N_DeptDofs;i++)
   } // if(N_Neibs)

  VelocityCommunicator->FECommunicateNeib(sendbuf, MaxN_DeptDofs_All, N_SendDofs, recevbuf, N_DIM);

  if(N_Neibs)
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

  for(i=0;i<N_DeptDofs;i++)
   {
    N = N_DeptDofNeibs[i];
    for(j=0;j<N;j++)
     {
      S = i*MaxSubDomainPerDof + j;
      Neib_ID = DeptDofNeibRanks[S];
      k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
      T = k*MaxN_DeptDofs_All;
      dof =  LocDOFofNeib[S];

      M = N_DIM*T + N_SendDofs[k];

      for(l=0;l<N_DIM;l++)
       {
        P = M + l*MaxN_DeptDofs_All;
        Values[l*NU + dof] += recevbuf[P];
       }
      N_SendDofs[k]++;
     } // for(j=0;j<N;j++)
   } //  for(i=0;i<N_DeptDofs;i++)

  if(N_Neibs)
   {
    delete [] N_SendDofs;
    delete [] sendbuf;
    delete [] recevbuf;
   }

 // pressure 
 PressureCommunicator->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs, DeptDofs, 
                     N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks, IndexOfNeibRank,
                     LocDOFofNeib, NeibLocalDOF);

  if(N_Neibs)
   {
    pressure = Values+N_DIM*NU;

    N_SendDofs = new int[N_Neibs];
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    sendbuf = new double[N_Neibs*MaxN_DeptDofs_All];
    recevbuf = new double[N_Neibs*MaxN_DeptDofs_All];

    for(i=0;i<N_DeptDofs;i++)
     {
      N = N_DeptDofNeibs[i];
      for(j=0;j<N;j++)
       {
        S = i*MaxSubDomainPerDof + j;
        Neib_ID = DeptDofNeibRanks[S];
        k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
        M = k*MaxN_DeptDofs_All + N_SendDofs[k];
        sendbuf[M]= pressure[DeptDofs[i]];
        N_SendDofs[k]++;
       } // for(j=0;j<N;j++)
     }// for(i=0;i<N_DeptDofs;i++)
   } // if(N_Neibs)

  PressureCommunicator->FECommunicateNeib(sendbuf, MaxN_DeptDofs_All, N_SendDofs, recevbuf, 1);

  if(N_Neibs)
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

  for(i=0;i<N_DeptDofs;i++)
   {
    N = N_DeptDofNeibs[i];
    for(j=0;j<N;j++)
     {
      S = i*MaxSubDomainPerDof + j;
      Neib_ID = DeptDofNeibRanks[S];
      k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
      dof =  LocDOFofNeib[S];
      M = k*MaxN_DeptDofs_All + N_SendDofs[k];
      pressure[dof] += recevbuf[M];
      N_SendDofs[k]++;
     } // for(j=0;j<N;j++)
   } //  for(i=0;i<N_DeptDofs;i++)

  if(N_Neibs)
   {
    delete [] N_SendDofs;
    delete [] sendbuf;
    delete [] recevbuf;
   }

//  MPI_Finalize();
//  exit(0);
}


void TParVectorNSE3D::AssembleVeloAtRootByADD(double *GlobalVal)
{
 int rank;

  MPI_Comm_rank(Comm, &rank);

  if(rank!=0)
   {
    MPI_Isend(Values, 3*NU, MPI_DOUBLE, 0, 100, Comm, &requestAssembleAtRoot);
   }
  else
   {
    int M, *N_LocalDofAllRank, MaxN_LocalDofAllRank, *GlobalDofOFLocalDofAllRank;
    int i, j, k, l, size, Global_NU, Global_NP;
    double *val;

    MPI_Comm_size(Comm, &size);

    Global_NU = VelocityCommunicator->GetN_GlobalDegreesOfFreedom();
    Global_NP = PressureCommunicator->GetN_GlobalDegreesOfFreedom();

    memset(GlobalVal, 0, (3*Global_NU+Global_NP)*SizeOfDouble);

    VelocityCommunicator->GetLocalDofAllRankInfo(N_LocalDofAllRank, MaxN_LocalDofAllRank,
                                                 GlobalDofOFLocalDofAllRank);
    // copy own values first
    for(i=0; i<NU; i++)
     {
      l = GlobalDofOFLocalDofAllRank[i];
      GlobalVal[l] += Values[i];
      GlobalVal[l+Global_NU] += Values[i+NU];
      GlobalVal[l+2*Global_NU] += Values[i+2*NU];
     }

    // add only 3 velo components
    val = new double[3*MaxN_LocalDofAllRank];

    for(i=1; i<size; i++)
     {
      M = N_LocalDofAllRank[i];
      MPI_Recv(val, 3*M, MPI_DOUBLE, i, 100, Comm, MPI_STATUS_IGNORE);

      k = i*MaxN_LocalDofAllRank;
      for(j=0; j<M; j++)
       {
        l = GlobalDofOFLocalDofAllRank[k + j];
        GlobalVal[l] += val[j];
        GlobalVal[l+Global_NU] += val[j+M];
        GlobalVal[l+2*Global_NU] += val[j+2*M];
       } // for(j=0; j<M; j++)
     } //for(i=1; i<N; i++)

    delete [] val;
   }
}


/** scatter the global vector from root to all sub domains */
/** must be called by all ranks in MPI_Comm */
void TParVectorNSE3D::ScatterFromRoot(double *GlobalVal)
{
 int rank;

  MPI_Comm_rank(Comm, &rank);

  if(rank!=0)
   {
    MPI_Recv(Values, 3*NU+NP, MPI_DOUBLE, 0, 100, Comm, MPI_STATUS_IGNORE);
   }
  else
   {
    int n, n1, M, *N_LocalDofAllRank, MaxN_LocalDofAllRank, *GlobalDofOFLocalDofAllRank;
    int M_P, *N_LocalDofAllRank_P, MaxN_LocalDofAllRank_P, *GlobalDofOFLocalDofAllRank_P;   
    int i, j, k, l, size, Global_NU, Global_NP, *GlobalDof;
    double *val, *val_P, *val_loc;

    MPI_Request request;
	
    MPI_Comm_size(Comm, &size);

    Global_NU = VelocityCommunicator->GetN_GlobalDegreesOfFreedom();
    Global_NP = PressureCommunicator->GetN_GlobalDegreesOfFreedom();

    VelocityCommunicator->GetLocalDofAllRankInfo(N_LocalDofAllRank, MaxN_LocalDofAllRank,
                                                 GlobalDofOFLocalDofAllRank);

    PressureCommunicator->GetLocalDofAllRankInfo(N_LocalDofAllRank_P, MaxN_LocalDofAllRank_P,
                                                 GlobalDofOFLocalDofAllRank_P);


    val = new double[size*(3*MaxN_LocalDofAllRank+ MaxN_LocalDofAllRank_P)];

    for(i=1; i<size; i++)
     {
      M = N_LocalDofAllRank[i];
      n = 2*Global_NU;
      n1 = 2*M;

      val_loc = val + i*(3*MaxN_LocalDofAllRank+ MaxN_LocalDofAllRank_P);
      GlobalDof = GlobalDofOFLocalDofAllRank+(i*MaxN_LocalDofAllRank);

      for(j=0; j<M; j++)
       {
        l = GlobalDof[j];
        val_loc[j]     = GlobalVal[l];
        val_loc[j+M]   = GlobalVal[l+Global_NU];
        val_loc[j+n1] = GlobalVal[l+n];
       } // for(j=0; j<M; j++)

       M_P = N_LocalDofAllRank_P[i];
       n = 3*Global_NU;

       val_P = val_loc + 3*M;
       GlobalDof = GlobalDofOFLocalDofAllRank_P+(i*MaxN_LocalDofAllRank_P);

       for(j=0; j<M_P; j++)
       {
        l = GlobalDof[j];
        val_P[j] = GlobalVal[n+l];
       } // for(j=0; j<M; j++)

      if(i>1)
       MPI_Wait(&request, MPI_STATUS_IGNORE);

      MPI_Isend(val_loc, 3*M+M_P, MPI_DOUBLE, i, 100, Comm, &request);

     } //for(i=1; i<size; i++)

 
    // copy own values first
    n = 2*Global_NU;
    n1 = 2*NU;
    for(j=0; j<NU; j++)
    {
     l = GlobalDofOFLocalDofAllRank[j];
     Values[j] = GlobalVal[l];
     Values[j+NU] = GlobalVal[l+Global_NU];
     Values[j+n1] = GlobalVal[l+n];
    } // for(j=0; j<M; j++)
 
    n = 3*Global_NU;
    n1 = 3*NU;
    
    for(j=0; j<NP; j++)
     {
      l = GlobalDofOFLocalDofAllRank_P[j];
      Values[n1 + j] = GlobalVal[n+l];
     } // for(j=0; j<M; j++)

   if(size>1)
    MPI_Wait(&request, MPI_STATUS_IGNORE);
   
    delete [] val;
   }
   
   
} // ScatterFromRoot(double *GlobalVal)


/** Copy the master dof value to all slave dofs */
void TParVectorNSE3D::AssembleWithMaster()
{
 int i, j, k, l, M, N, P, dof, Neib_ID, *N_SendDofs, *N_RecevDofs, Master_ID, rank, disp;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *N_DeptDofNeibs, MaxN_DeptDofs_All;
 int *DeptDofs, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocalDOF;
 int MaxN_MasterDofs_All, MaxN_SlaveDofs_All, *DeptDofMasterID;
 int **LocalDofOfRecbuf_send, **LocalDofOfRecbuf;

 double **sendbuf, **recevbuf, *u1, *u2, *u3, *p;

 MPI_Comm_rank(Comm, &rank);

 // first velocity
 VelocityCommunicator->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs, DeptDofs, 
                     N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks, IndexOfNeibRank,
                     LocDOFofNeib, NeibLocalDOF);

 VelocityCommunicator->GetDeptDofMasterInfo(MaxN_MasterDofs_All, MaxN_SlaveDofs_All,
                                             DeptDofMasterID);

  u1 = Values;
  u2 = Values+NU;
  u3 = Values+2*NU;

  if(N_Neibs)
   {
    N_SendDofs =  new int[N_Neibs];
    N_RecevDofs =  new int[N_Neibs];
    memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    LocalDofOfRecbuf = new int*[1];
    LocalDofOfRecbuf_send  = new int*[1];

    sendbuf = new double*[3];
    recevbuf = new double*[3];

    for(i=0;i<3;i++) // u1, u2, u3
     {
      sendbuf[i] = new double[N_Neibs*MaxN_MasterDofs_All];
      recevbuf[i]  = new double[N_Neibs*MaxN_MasterDofs_All];
     }

    LocalDofOfRecbuf[0] = new int[N_Neibs*MaxN_MasterDofs_All];
    LocalDofOfRecbuf_send[0]  = new int[N_Neibs*MaxN_MasterDofs_All];


    for(i=0;i<N_DeptDofs;i++)
     {
      Master_ID = DeptDofMasterID[i];

      if(Master_ID==rank)
       {
        N = N_DeptDofNeibs[i];
        dof = DeptDofs[i];

        for(j=0;j<N;j++)
         {
          P = i*MaxSubDomainPerDof + j;
          Neib_ID = DeptDofNeibRanks[P];

          k = IndexOfNeibRank[Neib_ID];
          M = k*MaxN_MasterDofs_All + N_SendDofs[k];
          N_SendDofs[k]++;


          LocalDofOfRecbuf_send[0][M] = NeibLocalDOF[P];

          sendbuf[0][M] = u1[dof];
          sendbuf[1][M] = u2[dof];
          sendbuf[2][M] = u3[dof];
         } // for(j=0;j<N;j++)
       }
      else 
       {
        k = IndexOfNeibRank[Master_ID];
        N_RecevDofs[k]++;
       } // else
     } //  for(i=0;i<N_DeptDofs;i++)
   }

  //   communicate 
  VelocityCommunicator->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_MasterDofs_All, N_SendDofs,
                                 LocalDofOfRecbuf,  N_RecevDofs, 1);

  VelocityCommunicator->FECommunicateNeib(sendbuf, MaxN_MasterDofs_All, N_SendDofs, recevbuf,
                                          N_RecevDofs, 3);


// now fill the master dof info from neibs
  for(i=0;i<N_Neibs;i++)
   {
    N = N_RecevDofs[i];

    for(j=0;j<N;j++)
     {
      M = i*MaxN_MasterDofs_All + j;
      dof = LocalDofOfRecbuf[0][M];

      u1[dof] = recevbuf[0][M];
      u2[dof] = recevbuf[1][M];
      u3[dof] = recevbuf[2][M];
     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)


  if(N_Neibs)
   {
    delete [] N_SendDofs;
    delete [] N_RecevDofs;

    for(i=0;i<3;i++) // u1, u2, u3
     {
      delete [] sendbuf[i];
      delete [] recevbuf[i];
     }

    delete [] LocalDofOfRecbuf_send[0];
    delete [] LocalDofOfRecbuf[0];

    delete [] sendbuf;
    delete [] recevbuf;
    delete [] LocalDofOfRecbuf_send;
    delete [] LocalDofOfRecbuf;
   }

  //presure
  PressureCommunicator->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs, DeptDofs, 
                     N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks, IndexOfNeibRank,
                     LocDOFofNeib, NeibLocalDOF);

  PressureCommunicator->GetDeptDofMasterInfo(MaxN_MasterDofs_All, MaxN_SlaveDofs_All,
                                             DeptDofMasterID);

  p = Values+3*NU;

  if(N_Neibs)
   {
    N_SendDofs =  new int[N_Neibs];
    N_RecevDofs =  new int[N_Neibs];
    memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    LocalDofOfRecbuf = new int*[1];
    LocalDofOfRecbuf_send  = new int*[1];

    sendbuf = new double*[1];
    recevbuf = new double*[1];

    sendbuf[0] = new double[N_Neibs*MaxN_MasterDofs_All];
    recevbuf[0]  = new double[N_Neibs*MaxN_MasterDofs_All];
    LocalDofOfRecbuf[0] = new int[N_Neibs*MaxN_MasterDofs_All];
    LocalDofOfRecbuf_send[0]  = new int[N_Neibs*MaxN_MasterDofs_All];

    for(i=0;i<N_DeptDofs;i++)
     {
      Master_ID = DeptDofMasterID[i];

      if(Master_ID==rank)
       {
        N = N_DeptDofNeibs[i];
        dof = DeptDofs[i];

        for(j=0;j<N;j++)
         {
          P = i*MaxSubDomainPerDof + j;
          Neib_ID = DeptDofNeibRanks[P];

          k = IndexOfNeibRank[Neib_ID];
          M = k*MaxN_MasterDofs_All + N_SendDofs[k];
          N_SendDofs[k]++;


          LocalDofOfRecbuf_send[0][M] = NeibLocalDOF[P];
          sendbuf[0][M] = p[dof];
         } // for(j=0;j<N;j++)
       }
      else 
       {
        k = IndexOfNeibRank[Master_ID];
        N_RecevDofs[k]++;
       } // else
     } //  for(i=0;i<N_DeptDofs;i++)
   } // if(N_Neibs)

  //   communicate 
  PressureCommunicator->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_MasterDofs_All, N_SendDofs,
                                          LocalDofOfRecbuf,  N_RecevDofs, 1);

  PressureCommunicator->FECommunicateNeib(sendbuf, MaxN_MasterDofs_All, N_SendDofs, recevbuf,
                                          N_RecevDofs, 1);


// now fill the master dof info from neibs
  for(i=0;i<N_Neibs;i++)
   {
    N = N_RecevDofs[i];

    for(j=0;j<N;j++)
     {
      M = i*MaxN_MasterDofs_All + j;
      dof = LocalDofOfRecbuf[0][M];

      p[dof] = recevbuf[0][M];
     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)

  if(N_Neibs)
   {
    delete [] N_SendDofs;
    delete [] N_RecevDofs;

    delete [] sendbuf[0];
    delete [] recevbuf[0];
    delete [] LocalDofOfRecbuf_send[0];
    delete [] LocalDofOfRecbuf[0];

    delete [] sendbuf;
    delete [] recevbuf;
    delete [] LocalDofOfRecbuf_send;
    delete [] LocalDofOfRecbuf;
   }





// printf("AssembleWithMaster  \n");
// 
// MPI_Finalize();
// exit(0);
} // AssembleWithMaster





void TParVectorNSE3D::IntoL20Pressure()
{
 int order = TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE;
 int i, totallength, N_OwnDof, *OwnDofs, Global_NP;
  
 double s, temp, *v;

   
  v = Values+3*NU;
     
  switch(order)
  {
    case -11:
      s=0;
      for(i=0;i<NP;i+=4) 
        s += v[i];

      MPI_Allreduce(&NP, &totallength, 1, MPI_INT, MPI_SUM, Comm);
      MPI_Allreduce(&s, &temp, 1, MPI_DOUBLE, MPI_SUM, Comm);

      s = temp/(double)totallength;
      s *= 4;

      for(i=0;i<NP;i+=4)
        v[i] -= s;
      break;

    case -12:
      s=0;
      for(i=0;i<NP;i+=10)
        s += v[i];

      MPI_Allreduce(&NP, &totallength, 1, MPI_INT, MPI_SUM, Comm);
      MPI_Allreduce(&s, &temp, 1, MPI_DOUBLE, MPI_SUM, Comm);
      
      s = temp/(double)totallength;
      s *= 10;

      for(i=0;i<NP;i+=10)
        v[i] -= s;

      break;

    default :
      
     PressureCommunicator->GetOwnDofs(N_OwnDof, OwnDofs);
     Global_NP = PressureCommunicator->GetN_GlobalDegreesOfFreedom();

     s=0;
     for(i=0;i<N_OwnDof;i++)
       s += v[OwnDofs[i]];

      MPI_Allreduce(&s, &temp, 1, MPI_DOUBLE, MPI_SUM, Comm);
      s = temp/(double)Global_NP;

      for(i=0;i<N_OwnDof;i++)
        v[OwnDofs[i]] -= s;
      break;
  }
}



// void TParVectorNSE::GetOwnValues(double *values_own)
// {
//  int i, k, *OwnDofIndexU, *OwnDofIndexP;
// 
//   N_own = N_DIM*NU_own + NP_own;
//   OwnDofIndexU=  VelocityCommunicator->GetOwnDofIndex();
//   OwnDofIndexP=  PressureCommunicator->GetOwnDofIndex();
// 
//   for(i=0; i<NU_own; i++)
//    {
//     k = OwnDofIndexU[i];
// 
//     values_own[i] = Values[k];
//     values_own[NU_own + i] = Values[NU + k];
//    }
// 
//   for(i=0; i<NP_own; i++)
//    {
//     k = OwnDofIndexP[i];
// 
//     values_own[2*NU_own + i] = Values[2*NU + k];
//    }
// 
// }


TParVectorNSE3D::~TParVectorNSE3D()
{
 if(Values)  delete [] Values;
}
#endif




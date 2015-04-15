// =======================================================================
// @(#)ParVector3D.C
//
// Class:      TParVector3D
// Purpose:    Class containing all info needed for communication between subdomains
//             for solution/rhs vectors
//
// Author:     Sashikumaar Ganesan (14.12.10)
//
// History:    Start of implementation 14.12.10 (Sashikumaar Ganesan)
//
// =======================================================================

#ifdef _MPI
#include "mpi.h"
#include <ParVector3D.h>

#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <LinAlg.h>

TParVector3D::TParVector3D(MPI_Comm comm, double *U, int ndof, int N_Dim,  TParFECommunicator3D *communicator)
                   :  TParVector(comm, U, N_Dim*ndof, N_Dim)
{
   
 if(communicator)
  FESpace = communicator->Getfespace();
 
 Communicator = communicator;
}

/** Set the No. of Master Dofs which is needed in other processors in iterative solver */
/** in each neib */
void TParVector3D::SetMastDofToHalloNeib(int N_Neibs, int *N_MastDofToHalloNeib, 
                                         int N_MastDofToHalloNeib_All, int *MastDofToHalloNeib)
{ 
  TCLS_N_MastDofToHalloNeib = new int[N_Neibs];
  TCLS_MastDofToHalloNeib = new int[N_Neibs*N_MastDofToHalloNeib_All];
 
  memcpy(TCLS_N_MastDofToHalloNeib, N_MastDofToHalloNeib, N_Neibs*SizeOfInt);
  memcpy(TCLS_MastDofToHalloNeib, MastDofToHalloNeib, N_Neibs*N_MastDofToHalloNeib_All*SizeOfInt);    
  TCLS_N_MastDofToHalloNeib_All = N_MastDofToHalloNeib_All;
}


void TParVector3D::UpdateDepAndHalloDof(double *GlobalVal)
{
 int i, j, N, M, P, Master_ID;  
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *DeptDofs, *N_DeptDofNeibs;
 int MaxN_DeptDofs_All, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocDof;  
 int MaxN_MasterDofs_All, MaxN_SlaveDofs_All, *DeptDofMasterID, *N_DofRankIndex;
 int N_GlobDof, *GlobalDofOfLocalDof, *Dist_Dof_All; 
 int disp, dof;
 double *sendbuf, *recevbuf, *buf;
 
 if(TCLS_N_MastDofToHalloNeib_All==-1) 
  {
   cout << "Set SetMastDofToHalloNeib before calling TParVector3D::UpdateHalloDof" << endl;
   cout << "that is after initialize the parallel solver " << endl;
   MPI_Finalize();
   exit(0);      
  }   
    
  Communicator->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs,
                     DeptDofs, N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks,
                     IndexOfNeibRank, LocDOFofNeib, NeibLocDof);   
   
  Communicator->GetDeptDofMasterInfo(MaxN_MasterDofs_All, MaxN_SlaveDofs_All, DeptDofMasterID);
  Communicator->GetGlobalDistDofInfo(N_GlobDof, GlobalDofOfLocalDof, Dist_Dof_All);
  
  /** first update the dep dof with master and then update hallo here */
  this->AssembleWithMaster();

  for(i=0;i<N_DeptDofs;i++)
   {
    dof = DeptDofs[i];
    GlobalVal[GlobalDofOfLocalDof[dof]] = Values[dof];  
   }

  /** now update the hallo dof */
   // 1st column stores N_columns of U 
   disp = TCLS_N_MastDofToHalloNeib_All + 1;    
   sendbuf = new double[N_Neibs*disp];
   recevbuf = new double[N_Neibs*disp];
  
   for(i=0;i<N_Neibs;i++)
    {
     N = TCLS_N_MastDofToHalloNeib[i];    
     
     if(N<1)
      continue;
      
     M = i*disp;
     P = i*TCLS_N_MastDofToHalloNeib_All;
     
     //first pos contains the number of calues    
     sendbuf[M++] = (double)N;
 
     for(j=0;j<N;j++)
      sendbuf[M+j] = Values[TCLS_MastDofToHalloNeib[P+j]]; 
    } //   for(i=0;i<N_DeptDofs;i++)
 
 
   // type =1, i.e., higher to lower ranks
   Communicator->FECommunicateOneWay(sendbuf, disp, recevbuf, 1);
   
   for(i=0;i<N_Neibs;i++)
    {  
     buf =  recevbuf+(i*disp);
     N = (int)buf[0];  //first pos contains the number of column    
     
//      for(j=0;j<N;j++) 
   
   
    } //   for(i=0;i<N_DeptDofs;i++)  
   cout << "TParVector3D::UpdateHalloDof no yet implemented" << endl;
   MPI_Finalize();
   exit(0);    
 
}



// Methods
void TParVector3D::ParDdot(int assembletype, double *residual)
{
 int ii, i, N_OwnDof, *OwnDofs;
 double locres=0., *tempval, *locval, res;
 
  if(assembletype==BYADD)
   {
    tempval = new double[Len];    
    memcpy(tempval, Values, Len*SizeOfDouble);        
    this->AssembleByADD(tempval);
   }
  else if(assembletype==BYMASTER)
   { 
     this->AssembleWithMaster();
     tempval=Values;
   }
  else if(assembletype==BYOWN) 
   {
     // do nothing, since the vector is already in consistant form 
     tempval=Values;
   }
  else
   {
    cout << "Wrong Assemble type !!!!!!!!!!!!! TParVector3D::ParDdot" << endl;
    MPI_Finalize();
    exit(0);    
   }
    
  Communicator->GetOwnDofs(N_OwnDof, OwnDofs);    
  
  for(ii=0;ii<N_DIM;ii++)  // more than one components
   {
    locres = 0.;
    locval = tempval + (ii*LDIM);
    
    
    for(i=0;i<N_OwnDof;i++)   
     locres += locval[OwnDofs[i]]*locval[OwnDofs[i]];    
    
    
     MPI_Allreduce(&locres, &res, 1, MPI_DOUBLE, MPI_SUM, Comm);
    
    residual[ii] = res;
   }
   
  if(assembletype==BYADD)
   delete [] tempval;  
}

void TParVector3D::AssembleByADD(double *tempval)
{
 int i, j, k, l, M, N, P, S, T, dof, Neib_ID, *N_SendDofs, LocDofofDepDof;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *N_DeptDofNeibs, MaxN_DeptDofs_All;
 int *DeptDofs, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocalDOF;
 
 double *sendbuf, *recevbuf;  
  
  
 Communicator->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs, DeptDofs, 
                              N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks, IndexOfNeibRank,
                              LocDOFofNeib, NeibLocalDOF);  
 
 // add subdomain interface dof tempval
 if(N_Neibs)
  {
    N_SendDofs = new int[N_Neibs];
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    sendbuf = new double[N_Neibs*MaxN_DeptDofs_All*N_DIM];
    recevbuf = new double[N_Neibs*MaxN_DeptDofs_All*N_DIM];
    
    for(i=0;i<N_DeptDofs;i++)
     {
      N = N_DeptDofNeibs[i];
      LocDofofDepDof = DeptDofs[i];
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
          sendbuf[P]= tempval[l*LDIM + LocDofofDepDof];
         }
        N_SendDofs[k]++;
       } // for(j=0;j<N;j++)
     }// for(i=0;i<N_DeptDofs;i++)       
  }// if(N_Neibs)
  
  Communicator->FECommunicateNeib(sendbuf, MaxN_DeptDofs_All, N_SendDofs, recevbuf, N_DIM); 
  
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
        tempval[l*LDIM + dof] += recevbuf[P];
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
   
}

/** Copy the master dof value to all slave dofs */
void TParVector3D::AssembleWithMaster()
{
 int rank, i, j, k, l, M, N, P, S, T, dof, Neib_ID, *N_SendDofs, *N_RecevDofs;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *N_DeptDofNeibs, MaxN_DeptDofs_All;
 int *DeptDofs, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocalDOF;
 int MaxN_MasterDofs_All, MaxN_SlaveDofs_All, *DeptDofMasterID, Master_ID;
 int **LocalDofOfRecbuf_send, **LocalDofOfRecbuf;
 
 double *sendbuf, *recevbuf;  
  
 MPI_Comm_rank(Comm, &rank);
 
 Communicator->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs, DeptDofs, 
                              N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks, IndexOfNeibRank,
                              LocDOFofNeib, NeibLocalDOF);    
 
 Communicator->GetDeptDofMasterInfo(MaxN_MasterDofs_All, MaxN_SlaveDofs_All, DeptDofMasterID);
 // add subdomain interface dof tempval
 if(N_Neibs)
  {
   N_SendDofs =  new int[N_Neibs];
   N_RecevDofs =  new int[N_Neibs];
   memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);
   memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

   sendbuf = new double[N_Neibs*MaxN_MasterDofs_All*N_DIM];
   recevbuf = new double[N_Neibs*MaxN_MasterDofs_All*N_DIM];
   
   LocalDofOfRecbuf = new int*[1];
   LocalDofOfRecbuf_send  = new int*[1];
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
          T = k*MaxN_MasterDofs_All;  

          LocalDofOfRecbuf_send[0][T + N_SendDofs[k]] = NeibLocalDOF[P];
          M = N_DIM*T + N_SendDofs[k]; 
          N_SendDofs[k]++;

          for(l=0;l<N_DIM;l++)
           {
            S = M + l*MaxN_MasterDofs_All;
            sendbuf[S]= Values[l*LDIM + dof];
           }

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
  Communicator->FECommunicateNeib(LocalDofOfRecbuf_send, MaxN_MasterDofs_All, N_SendDofs,
                                  LocalDofOfRecbuf,  N_RecevDofs, 1);

  Communicator->FECommunicateNeib(sendbuf, MaxN_MasterDofs_All, N_SendDofs, recevbuf, N_RecevDofs, N_DIM);
 
  // now fill the master dof info from neibs
  for(i=0;i<N_Neibs;i++)
   {
    N = N_RecevDofs[i];

    for(j=0;j<N;j++)
     {
      T = i*MaxN_MasterDofs_All;
      M = N_DIM*T + j;       
      dof = LocalDofOfRecbuf[0][T+j];

      for(l=0;l<N_DIM;l++)
       {
        S = M + l*MaxN_MasterDofs_All;
        Values[l*LDIM + dof] = recevbuf[S];     
       }
     } // for(j=0;j<N;j++)
   } // for(i=0;i<N_Neibs;i++)
   
   
//    printf("AssembleWithMaster not yet implemented \n"); 
//    MPI_Finalize();
//    exit(0);
}

/** assemble the vector by averaging dof values on subdomain interface */   
void TParVector3D::AssembleByAverage()
{
 int i, j, k, l, M, N, P, S, T, dof, Neib_ID, *N_SendDofs, LocDofofDepDof;
 int MaxSubDomainPerDof, N_Neibs, *NeibsRank, N_DeptDofs, *N_DeptDofNeibs, MaxN_DeptDofs_All;
 int *DeptDofs, *DeptDofNeibRanks, *IndexOfNeibRank, *LocDOFofNeib, *NeibLocalDOF;
 
 double *sendbuf, *recevbuf;  
  
  
 Communicator->GetDofNeibInfo(MaxSubDomainPerDof, N_Neibs, NeibsRank, N_DeptDofs, DeptDofs, 
                              N_DeptDofNeibs, MaxN_DeptDofs_All, DeptDofNeibRanks, IndexOfNeibRank,
                              LocDOFofNeib, NeibLocalDOF);  
 
 // first add subdomain interface dof values and then average it
 if(N_Neibs)
  {
    N_SendDofs = new int[N_Neibs];
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    sendbuf = new double[N_Neibs*MaxN_DeptDofs_All*N_DIM];
    recevbuf = new double[N_Neibs*MaxN_DeptDofs_All*N_DIM];
    
    for(i=0;i<N_DeptDofs;i++)
     {
      N = N_DeptDofNeibs[i];
      LocDofofDepDof = DeptDofs[i];
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
          sendbuf[P]= Values[l*LDIM + LocDofofDepDof];
         }
        N_SendDofs[k]++;
       } // for(j=0;j<N;j++)
     }// for(i=0;i<N_DeptDofs;i++)       
  }// if(N_Neibs)
  
  Communicator->FECommunicateNeib(sendbuf, MaxN_DeptDofs_All, N_SendDofs, recevbuf, N_DIM); 
  
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
        Values[l*LDIM + dof] += recevbuf[P];
       }
      N_SendDofs[k]++;
     } // for(j=0;j<N;j++)
   } //  for(i=0;i<N_DeptDofs;i++)

 // averaging
  for(i=0;i<N_DeptDofs;i++)
   {
    N = N_DeptDofNeibs[i];    
    dof = DeptDofs[i];  
    
    for(l=0;l<N_DIM;l++)
     Values[l*LDIM + dof] /= (double)(N+1.);    
   }//  for(i=0;i<N_DeptDofs;i++)

  if(N_Neibs)
   {
    delete [] N_SendDofs;
    delete [] sendbuf;
    delete [] recevbuf;
   }
//   printf("AssembleByAverage \n");
}// AssembleByAverage


void TParVector3D::AssembleAtRootByADD(double *GlobalVal)
{
 int rank;

  MPI_Comm_rank(Comm, &rank);
  
  if(rank!=0)
   {
    // Len=N_Dim*ndof, see the above call to TParVector constructor 
    MPI_Send(Values, Len, MPI_DOUBLE, 0, 100, Comm);
   }
  else
   {
    int M, *N_LocalDofAllRank, MaxN_LocalDofAllRank, *GlobalDofOFLocalDofAllRank;
    int i, j, k, l, size, Global_N, *GDOF;
    double *val, *GlobalVal_Loc, *Values_Loc;

    MPI_Comm_size(Comm, &size);

    Global_N = Communicator->GetN_GlobalDegreesOfFreedom();     
    memset(GlobalVal, 0, N_DIM*Global_N*SizeOfDouble); 
    
    Communicator->GetLocalDofAllRankInfo(N_LocalDofAllRank, MaxN_LocalDofAllRank,
                                         GlobalDofOFLocalDofAllRank);

    // copy own values first
    for(j=0; j<N_DIM; j++)   
     {
      GlobalVal_Loc = GlobalVal +(j*Global_N);
      Values_Loc = Values +(j*LDIM);
      
      for(i=0; i<LDIM; i++)
       {
        l = GlobalDofOFLocalDofAllRank[i];
// 	if(l==542)
// 	   printf("l %d GlobalVal_Loc[l] %e  \n", l,   Values_Loc[i]);	 
        GlobalVal_Loc[l] += Values_Loc[i];
       } // for(i=0; i<LDIM; i+                                      
     } // for(j=0; j<N_DIM; 
     
     
    val = new double[N_DIM*MaxN_LocalDofAllRank];                                        
                                         
    for(i=1; i<size; i++)
     {
      M = N_LocalDofAllRank[i];
      MPI_Recv(val, M*N_DIM, MPI_DOUBLE, i, 100, Comm, MPI_STATUS_IGNORE);

      GDOF = GlobalDofOFLocalDofAllRank + (i*MaxN_LocalDofAllRank);
      
      for(k=0; k<N_DIM; k++)    
       {
        GlobalVal_Loc = GlobalVal +(k*Global_N);
        Values_Loc = val +(k*M); 
 
        for(j=0; j<M; j++)
         {
          l = GDOF[j];
 
          GlobalVal_Loc[l] += Values_Loc[j];
         } // for(j=0; j<M; j++)
       }// for(k=0; k<N_DIM; k++)
       
       
     } //for(i=1; i< size; i++)

    delete [] val;                                
   }  
}// AssembleAtRootByADD


/** scatter the global vector from root to all sub domains */
void TParVector3D::ScatterFromRoot(double *GlobalVal)
{
 int rank;

  MPI_Comm_rank(Comm, &rank);

  if(rank!=0)
   {
    // Len=N_Dim*ndof, see the above call to TParVector constructor      
    MPI_Recv(Values, Len, MPI_DOUBLE, 0, 100, Comm, MPI_STATUS_IGNORE);
    
    
    
//  int k;  
//      if(rank==5)
//       {
//        if(N_DIM>1)
//         { 
//          for(k=0; k<N_DIM; k++)
//           OutPut(k << " PSD ParDdot "<< setprecision(15) << Ddot(LDIM, Values+k*LDIM, Values+k*LDIM)<< endl);  }
//        else
//         { OutPut("ParDdot "<< setprecision(15) << Ddot(LDIM, Values, Values)<< endl);  }       
//       }
      
   }
  else
   {
    int M, *N_LocalDofAllRank, MaxN_LocalDofAllRank, *GlobalDofOFLocalDofAllRank; 
    int i, j, k, l, size, Global_N, *GlobalDof;
    double *val, *val_loc, *GlobalVal_Loc;

    MPI_Request request;

    MPI_Comm_size(Comm, &size);

    Global_N = Communicator->GetN_GlobalDegreesOfFreedom();
    Communicator->GetLocalDofAllRankInfo(N_LocalDofAllRank, MaxN_LocalDofAllRank,
                                         GlobalDofOFLocalDofAllRank);  
     
    val = new double[N_DIM*MaxN_LocalDofAllRank];
    memset(val, 0, N_DIM*MaxN_LocalDofAllRank*SizeOfDouble); 
    
    for(i=1; i<size; i++)
     {
      M = N_LocalDofAllRank[i];   
      GlobalDof = GlobalDofOFLocalDofAllRank+(i*MaxN_LocalDofAllRank);

      for(k=0; k<N_DIM; k++)    
       {
        val_loc = val + (k*M);
        GlobalVal_Loc = GlobalVal +(k*Global_N);

        for(j=0; j<M; j++)
         {
          l = GlobalDof[j];
          val_loc[j] = GlobalVal_Loc[l]; 
         } // for(j=0; j<M; j++)
       } // for(k=0; k<N_DIM; k++)    
       
//       if(i>1)
//        MPI_Wait(&request, MPI_STATUS_IGNORE);

//       MPI_Isend(val, M*N_DIM, MPI_DOUBLE, i, 100, Comm, &request);
       // val should not be accessed untill MPI_Wait, so changed to MPI_Send

      MPI_Send(val, M*N_DIM, MPI_DOUBLE, i, 100, Comm);

     } //for(i=1; i<size; i++)
   
    // copy own values
    for(k=0; k<N_DIM; k++)  
     {
      GlobalVal_Loc = GlobalVal +(k*Global_N);
      val_loc = Values +(k*LDIM);
      
      for(j=0; j<LDIM; j++)
       {
        l = GlobalDofOFLocalDofAllRank[j];
//         printf(" locDof %d GlobalDofOFLocalDofAllRank %d \n", j, l); 	
        val_loc[j] = GlobalVal_Loc[l];
       } // for(j=0; j<M; j++)    
             
       
     } // for(k=0; k<N_DIM; k++) 
         
//     if(size>1)
//      MPI_Wait(&request, MPI_STATUS_IGNORE);
   
    delete [] val;
   }
  
}// void TParVector3D::ScatterFromRoot

TParVector3D::~TParVector3D()
{
 if(Values)  delete [] Values;

}

#endif
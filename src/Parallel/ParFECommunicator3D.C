// =======================================================================
// @(#)TParFECommunicator3D.C
//
// Class:      TParFECommunicator3D
// Purpose:    Class containing all info needed for communication between subdomains
//
// Author:     Sashikumaar Ganesan / Abdus Shamim (24.04.15)
//
// History:    Start of implementation 24.04.15 (Sashikumaar Ganesan/Abdus Shamim)
//
// ======================================================================= 
#ifdef _MPI

#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#include <FEDatabase3D.h>
#include <Database.h>
#include <SubDomainJoint.h>
#include <Edge.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#define GLOBAL_NO 0
#define DOF_NO 1

#define HALOCELL 0
#define NONHALO  1

extern double timeC;

TParFECommunicator3D::TParFECommunicator3D(TParFEMapper3D *mapper)
{
  Mapper = mapper;
  
  Comm = TDatabase::ParamDB->Comm;
  
  Master = Mapper->GetMaster();
  
  Mapper->GetCommInfo(N_Dim, N_Dof,
		      N_SendDof, N_SendDofMS, N_SendDofH1, N_SendDofH2,
		      Send_Info, Send_InfoMS, Send_InfoH1, Send_InfoH2,
		      Recv_Info, Recv_InfoMS, Recv_InfoH1, Recv_InfoH2,
		      N_DofSend, N_DofSendMS, N_DofSendH1, N_DofSendH2,
		      N_DofRecv, N_DofRecvMS, N_DofRecvH1, N_DofRecvH2,
		      sdispl,    sdisplMS,    sdisplH1,    sdisplH2,
		      rdispl,    rdisplMS,    rdisplH1,    rdisplH2,
		      DofSend,   DofSendMS,   DofSendH1,   DofSendH2,
		      DofRecv,   DofRecvMS,   DofRecvH1,   DofRecvH2,
		      N_Slave,   N_InterfaceS,N_Halo1,     N_Halo2);

}

TParFECommunicator3D::TParFECommunicator3D()
{
  Mapper = NULL;
  
  Comm = TDatabase::ParamDB->Comm;

}

void TParFECommunicator3D::CommUpdateMS(double *sol)
{
  if(!Mapper)
  {
    printf("Set The Mapper for Communicator routines\n");
    MPI_Finalize();
    exit(0);
  }

  int i,j,k;
  double t1,t2;
  t1=MPI_Wtime();
  
  for(i=0;i<N_SendDofMS;i++)
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      Send_InfoMS[k]=sol[DofSendMS[i]+j*N_Dof];
    }
  }

  MPI_Alltoallv(Send_InfoMS,N_DofSendMS,sdisplMS,MPI_DOUBLE,Recv_InfoMS,N_DofRecvMS,rdisplMS,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_InterfaceS;i++)  
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      sol[DofRecvMS[i]+j*N_Dof] = Recv_InfoMS[k];
    }
  }
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

void TParFECommunicator3D::CommUpdateH1(double *sol)
{
  if(!Mapper)
  {
    printf("Set The Mapper for Communicator routines\n");
    MPI_Finalize();
    exit(0);
  }
  
  int i,j,k;
  double t1,t2;
  t1=MPI_Wtime();
  
  for(i=0;i<N_SendDofH1;i++)
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      Send_InfoH1[k]=sol[DofSendH1[i]+j*N_Dof];
    }
  }

  MPI_Alltoallv(Send_InfoH1,N_DofSendH1,sdisplH1,MPI_DOUBLE,Recv_InfoH1,N_DofRecvH1,rdisplH1,MPI_DOUBLE,Comm);
  
//   for(i=0;i<N_Halo1;i++)  
//     sol[DofRecvH1[i]] = Recv_InfoH1[i];
  
  for(i=0;i<N_Halo1;i++)  
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      sol[DofRecvH1[i]+j*N_Dof] = Recv_InfoH1[k];
    }
  }
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

void TParFECommunicator3D::CommUpdateH2(double *sol)
{
  if(!Mapper)
  {
    printf("Set The Mapper for Communicator routines\n");
    MPI_Finalize();
    exit(0);
  }
  
  int i,j,k;
  double t1,t2;
  t1=MPI_Wtime();
  
  for(i=0;i<N_SendDofH2;i++)
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      Send_InfoH2[k]=sol[DofSendH2[i]+j*N_Dof];
    }
  }
  
  MPI_Alltoallv(Send_InfoH2,N_DofSendH2,sdisplH2,MPI_DOUBLE,Recv_InfoH2,N_DofRecvH2,rdisplH2,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_Halo2;i++)  
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      sol[DofRecvH2[i]+j*N_Dof] = Recv_InfoH2[k];
    }
  }
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

void TParFECommunicator3D::CommUpdate_M_H1(double *sol)
{
  if(!Mapper)
  {
    printf("Set The Mapper for Communicator routines\n");
    MPI_Finalize();
    exit(0);
  }
  
  
  CommUpdateMS(sol);
  CommUpdateH1(sol);
  
}

void TParFECommunicator3D::CommUpdate(double *sol)
{
  
  if(!Mapper)
  {
    printf("Set The Mapper for Communicator routines\n");
    MPI_Finalize();
    exit(0);
  }
  
  int i,j,k;
  double t1,t2;
  
  if(TDatabase::ParamDB->MapperType != 2)
  {
    CommUpdateMS(sol);
    CommUpdateH1(sol);
    CommUpdateH2(sol);
  }
  else
  {
    t1=MPI_Wtime();
    
    for(i=0;i<N_SendDof;i++)
    {
      for(j=0;j<N_Dim;j++)
      {
	k = i*N_Dim + j;
	Send_Info[k]=sol[DofSend[i]+j*N_Dof];
      }
    }

    MPI_Alltoallv(Send_Info,N_DofSend,sdispl,MPI_DOUBLE,Recv_Info,N_DofRecv,rdispl,MPI_DOUBLE,Comm);
  
    for(i=0;i<N_Slave;i++)  
    {
      for(j=0;j<N_Dim;j++)
      {
        k = i*N_Dim + j;
        sol[DofRecv[i]+j*N_Dof] = Recv_Info[k];
      }
    }
    
    t2=MPI_Wtime(); 
    timeC+=(t2-t1);
  } 

}

void TParFECommunicator3D::CommUpdateReduce(double *rhs)
{
  if(!Mapper)
  {
    printf("Set The Mapper for Communicator routines\n");
    MPI_Finalize();
    exit(0);
  }
  
  double t1,t2;
  t1=MPI_Wtime();
  
  int i,j,k,rank,size;
  MPI_Status status;
  MPI_Comm_rank(Comm,&rank);		
  MPI_Comm_size(Comm, &size);
  
  int N_S, N_SendDof_temp;
  double *Recv_Info_temp, *Send_Info_temp;
  int *DofSend_temp, *DofRecv_temp, *sdisp_temp, *rdisp_temp;
  int *N_DofRecv_temp, *N_DofSend_temp;
  
  if(TDatabase::ParamDB->MapperType !=2 )
  {
    N_S             = N_InterfaceS;
    N_SendDof_temp  = N_SendDofMS;
    Recv_Info_temp  = Recv_InfoMS;
    Send_Info_temp  = Send_InfoMS;
    DofSend_temp    = DofSendMS;
    DofRecv_temp    = DofRecvMS;
    sdisp_temp      = sdisplMS;
    rdisp_temp      = rdisplMS;
    N_DofRecv_temp  = N_DofRecvMS;
    N_DofSend_temp  = N_DofSendMS;
  }
  else
  {
    N_S             = N_Slave;
    N_SendDof_temp  = N_SendDof;
    Recv_Info_temp  = Recv_Info;
    Send_Info_temp  = Send_Info;
    DofSend_temp    = DofSend;
    DofRecv_temp    = DofRecv;
    sdisp_temp      = sdispl;
    rdisp_temp      = rdispl;
    N_DofRecv_temp  = N_DofRecv;
    N_DofSend_temp  = N_DofSend;
  }

  for(i=0;i<N_S;i++)
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      Recv_Info_temp[k]=rhs[DofRecv_temp[i]+j*N_Dof];
    }
  }
 
  MPI_Alltoallv(Recv_Info_temp,N_DofRecv_temp,rdisp_temp,MPI_DOUBLE,Send_Info_temp,N_DofSend_temp,sdisp_temp,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_SendDof_temp;i++)
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      rhs[DofSend_temp[i]+j*N_Dof] += Send_Info_temp[k];
    }
  }

  /** MASTER HAS FINAL VALUE **/  
  for(i=0;i<N_SendDof_temp;i++)
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      Send_Info_temp[k]=rhs[DofSend_temp[i]+j*N_Dof];
    }
  }
  
  MPI_Alltoallv(Send_Info_temp,N_DofSend_temp,sdisp_temp,MPI_DOUBLE,Recv_Info_temp,N_DofRecv_temp,rdisp_temp,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_S;i++)  
  {
    for(j=0;j<N_Dim;j++)
    {
      k = i*N_Dim + j;
      rhs[DofRecv_temp[i]+j*N_Dof] = Recv_Info_temp[k];
    }
  }
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}



void TParFECommunicator3D::CommUpdate(double *sol, double *rhs)
{
  printf("GMRES not yet verified. Check if rhs update is reqd????\n");
  MPI_Finalize();
  exit(0);
}

#endif














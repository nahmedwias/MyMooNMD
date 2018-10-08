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
  Mapper = nullptr;
  
  Comm = TDatabase::ParamDB->Comm;

}

size_t TParFECommunicator3D::get_n_global_dof() const
{
  int n_m_total = 0;
  int sendbuf = Mapper->GetN_Master();
  MPI_Allreduce(&sendbuf, &n_m_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(n_m_total >= 0)
    return size_t(n_m_total);
  ErrThrow("negative number of global dofs computed ??? ", n_m_total);
}


void TParFECommunicator3D::print_info() const
{
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int root = 0;
  //TODO this method is open for extension and enhancement in every direction

  //1) get number of master dofs per process
  std::vector<int> ns_masters(size, 0);
  {
   int sendbuf = Mapper->GetN_Master();
   MPI_Gather(&sendbuf, 1, MPI_INT, 		 //send
              &ns_masters.at(0), 1, MPI_INT, //receive
              root, MPI_COMM_WORLD);         //control
  }

  //only root does the printing
  if(my_rank == 0)
  {
    Output::stat("ParFECommunicator3D of FESpace ",
                  Mapper->get_fe_space()->GetName());
    int n_m_total = 0;
    for( int i=0 ; i<size ; ++i )
    {
      Output::dash("Rank ", i, "\t n_masters: ", ns_masters.at(i)); //TODO further fille this output line
      n_m_total += ns_masters.at(i);
    }
    Output::dash("Total number of dofs in the communicator: ", n_m_total);
  }
}

void TParFECommunicator3D::consistency_update(double* vector, size_t level) const
{
   switch (level)
   {
     case 1:
       CommUpdateMS(vector); //restore level 1 consistency by updating all interface slaves
       break;
     case 2:
       CommUpdateMS(vector); //restore level 2 consistency by updating
       CommUpdateH1(vector); // all interface slaves and Halo1 dofs
       break;
     case 3:
       CommUpdateMS(vector); // restore level 2 consistency by updating
       CommUpdateH1(vector); // all interface slaves, Halo1 and Halo2 dofs
       CommUpdateH2(vector);
       break;
     default:
       ErrThrow("Unknown consistency level ", level, ". Choose between 1, 2 "
           "and 3 (full consistency).")
   }
}

void TParFECommunicator3D::average_at_interface(double* vector) const
{
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(N_Dim != 1)
    ErrThrow("ParFECommunicators of dimension greater than 1 are a "
             " not-yet-reintroduced feature!");

  /* ****************************************************** */
  /* *** Report from interface slaves to their masters **** */
  /* ****************************************************** */
  {
    //send - all values of interface slaves
    double* sbuf = new double[N_InterfaceS];
    for(int i=0;i<N_InterfaceS;i++)
    {
      sbuf[i]=vector[DofRecvMS[i]]; //fill all values of interface slaves into sendbuffer
    }
    int* scounts = N_DofRecvMS;
    int* sdispls = rdisplMS;

    //receive - all values of interface slaves belonging to one of my_ranks masters
    double* rbuf = new double[N_SendDofMS];
    int* rcounts = N_DofSendMS;
    int* rdispls = sdisplMS;
    //control

    // this is a topsy-turvy call of MPI_Alltoallv
    MPI_Alltoallv(
        sbuf, scounts, sdispls, MPI_DOUBLE, //send
        rbuf, rcounts, rdispls, MPI_DOUBLE, //receive
        MPI_COMM_WORLD); //control

    // take the received values, interpret them via DofSendMS, and add them to the
    // correct interface master dof value
    std::vector<int> n_updates(N_Dof, 0); //TODO there must be a more clever way to do this!
    for( int i = 0 ; i < N_SendDofMS ; i++)
    {
      vector[DofSendMS[i]] += rbuf[i];
      n_updates[DofSendMS[i]] += 1;
    }

    // now this is the averaging
    for( int j = 0 ; j < N_Dof ; ++j)
    {
      if( n_updates.at(j) > 0)
        vector[j] /=  n_updates[j] + 1; // "+1": this is the master dof itself!
    }

    delete[] sbuf;
    delete[] rbuf;
  }
}

void TParFECommunicator3D::update_from_additive_to_consistent_storage(
    double* vector, size_t consist_level) const
{
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(N_Dim != 1)
    ErrThrow("ParFECommunicators of dimension greater than 1 are an outdated feature!");

  /* ****************************************************** */
  /* *** Report from interface slaves to their masters **** */
  /* ****************************************************** */
  {
    //send - all values of interface slaves
    double* sbuf = new double[N_InterfaceS];
    for(int i=0;i<N_InterfaceS;i++)
    {
      sbuf[i]=vector[DofRecvMS[i]]; //fill all values of interface slaves into sendbuffer
    }
    int* scounts = N_DofRecvMS;
    int* sdispls = rdisplMS;

    //receive - all values of interface slaves belonging to one of my_ranks masters
    double* rbuf = new double[N_SendDofMS];
    int* rcounts = N_DofSendMS;
    int* rdispls = sdisplMS;
    //control

    // this is a topsy-turvy call of MPI_Alltoallv
    MPI_Alltoallv(
        sbuf, scounts, sdispls, MPI_DOUBLE, //send
        rbuf, rcounts, rdispls, MPI_DOUBLE, //receive
        MPI_COMM_WORLD); //control

    // take the received values, interpret them via DofSendMS, and add them to the
    // correct interface master dof value
    for( int i = 0 ; i < N_SendDofMS ; i++)
    {
      vector[DofSendMS[i]] += rbuf[i];
    }
    delete[] sbuf;
    delete[] rbuf;
  }

  /* ****************************************************** */
  /* **** Report from HALO 1  D.O.F to their masters ****** */
  /* ****************************************************** */
  {
    //send - all values of halo 1 d.o.f.
    double* sbuf = new double[N_Halo1];
    for( int i = 0 ; i < N_Halo1 ; i++ )
    {
      sbuf[i]=vector[DofRecvH1[i]]; //fill all values of halo1 d.o.f. into send buffer
    }
    int* scounts = N_DofRecvH1;
    int* sdispls = rdisplH1;

    //receive
    double* rbuf = new double[N_SendDofH1];
    int* rcounts = N_DofSendH1;
    int* rdispls = sdisplH1;
    //control

    // this is a topsy-turvy call of MPI_Alltoallv
    MPI_Alltoallv(
        sbuf, scounts, sdispls, MPI_DOUBLE, //send
        rbuf, rcounts, rdispls, MPI_DOUBLE, //receive
        MPI_COMM_WORLD); //control

    // take the received values, interpret them via DofSendH1
    for( int i = 0 ; i <  N_SendDofH1; i++)
    {
      vector[DofSendH1[i]] += rbuf[i];
    }
    delete[] sbuf;
    delete[] rbuf;
  }

  /* ****************************************************** */
  /* **** Report from HALO 2  D.O.F to their masters ****** */
  /* ****************************************************** */
  {
    //send - all values of halo2 d.o.f.
    double* sbuf = new double[N_Halo2];
    for( int i = 0 ; i < N_Halo2 ; i++ )
    {
      sbuf[i]=vector[DofRecvH2[i]]; //fill all values of halo2 d.o.f. into send buffer
    }
    int* scounts = N_DofRecvH2;
    int* sdispls = rdisplH2;

    //receive
    double* rbuf = new double[N_SendDofH2];
    int* rcounts = N_DofSendH2;
    int* rdispls = sdisplH2;
    //control

    // this is a topsy-turvy call of MPI_Alltoallv
    MPI_Alltoallv(
        sbuf, scounts, sdispls, MPI_DOUBLE, //send
        rbuf, rcounts, rdispls, MPI_DOUBLE, //receive
        MPI_COMM_WORLD); //control

    // take the received values, interpret them via DofSendH2
    for( int i = 0 ; i <  N_SendDofH2; i++)
    {
      vector[DofSendH2[i]] += rbuf[i];
    }
    delete[] sbuf;
    delete[] rbuf;
  }
  /* ****************************************************** */
  /* ****    Restore a certain consistency level     ****** */
  /* ****************************************************** */

  // Now the vector is stored level-0 consistent. If another
  // consistency level is required, perform the update here.
  if( consist_level > 0 )
  {
    consistency_update(vector, consist_level);
  }

}




void TParFECommunicator3D::CommUpdateReduce(double *rhs) const
{
  
  int i,j,k,rank,size;
//  MPI_Status status;
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
  
}

int TParFECommunicator3D::dof_ping(size_t process, size_t dof) const
{
  int size_gl, rank_gl;
  MPI_Comm_size(MPI_COMM_WORLD, &size_gl);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_gl);

  //Check input (on 'process' only)
  bool invalid_input= false;
  if(rank_gl == (int) process)
  {
    if((int) dof >= N_Dof)
      ErrThrow("There is no dof ", dof, " on process ", process, ", which knows ",N_Dof , " dof.");
    if(Master[dof] != (int) process)
    {
      Output::warn("DOF_PING", "Dof ", dof, " is not master on process ",
                   process, ", therefore no ping is sent");
      invalid_input = true;
    }
  }
  // Let everyone know if the input was valid (i.e. dof a master on 'process')
  // or not - if not so, everyone returns -1.
  MPI_Bcast(&invalid_input, 1, MPI_C_BOOL, process, MPI_COMM_WORLD);
  if(invalid_input)
    return(-1);

  // Input is valid - let's start working!
  std::vector<double> ping(N_Dof, 0);
  if(rank_gl == (int) process)
    ping.at(dof) = 1; //the actual ping is just a 1.0

  // ... update ping vector to full consistency...
  consistency_update(&ping[0], 3);

  // ...and start talking!
  const char* markers = get_dof_markers();

  int local_dof = -1;
  bool pinged = false;
  //look through the vector, find out if and where the ping was received
  for( size_t i = 0 ; i < ping.size() ; ++i )
  {
    if(ping[i] == 1.0)
    {
      if(pinged) //ping found, but already found another one!
        ErrThrow("Dof found multiple times through ping!");
      local_dof = i;
      pinged = true;
    }
  }

  /* ****** */
  // Set up a new communicator consisting of all those ranks, which contain
  // the pinged dof. With this communicator, gather information in that rank
  // which sent the ping and print the information.
  MPI_Comm mpi_comm_ping;
  int color = pinged ? 1 : MPI_UNDEFINED;
  int key = (rank_gl == (int) process) ? -1 : rank_gl ; //rank 0 of new comm will be 'process', the others maintain their order
  MPI_Comm_split( MPI_COMM_WORLD, color, key, &mpi_comm_ping);

  // This code is to be executed only by ranks which received the ping.
  // - gather information in root and print it (TODO and check its validity)
  if(color != MPI_UNDEFINED)
  {
    int rank_ping, size_ping;
    MPI_Comm_rank(mpi_comm_ping, &rank_ping);
    MPI_Comm_size(mpi_comm_ping, &size_ping);

    int* recv_proc = new int[size_ping];
    int* recv_dof = new int[size_ping];
    double* recv_pos = new double[3*size_ping]; //3D positions
    char* recv_mark = new char[size_ping];

    // Gather global process numbers.
    {
      int sbuf[1] = {rank_gl};
      MPI_Gather(sbuf, 1, MPI_INT,      //send
                 recv_proc, 1, MPI_INT, //receive
                 0, mpi_comm_ping);     //control
    }
    // Gather the local dof numbers where the ping was received
    {
      int sbuf[1] = {local_dof};
      MPI_Gather(sbuf, 1, MPI_INT,      //send
                 recv_dof, 1, MPI_INT, //receive
                 0, mpi_comm_ping);     //control
    }
    //Gather dof positions.
    {
      double x,y,z;
      get_fe_space()->GetDOFPosition(local_dof, x,y,z);
      double sbuf[3] = {x,y,z};
      MPI_Gather(sbuf, 3, MPI_DOUBLE,      //send
                 recv_pos, 3, MPI_DOUBLE, //receive
                 0, mpi_comm_ping);     //control
    }
    // Gather markers.
    {
      char sbuf[1] = {markers[local_dof]};
      MPI_Gather(sbuf, 1, MPI_CHAR,    //send
               recv_mark, 1, MPI_CHAR, //receive
               0, mpi_comm_ping);      //control

    }


    if(rank_ping == 0) //only one process does the printing
    {
      Output::info("PING","Ping sent by process ", recv_proc[0] ," \t ",
                   recv_dof[0], " \t ", recv_mark[0], " \t ",
                   recv_pos[0]," ", recv_pos[0+1], " ", recv_pos[0 + 2] );
      for(int i = 1; i < size_ping; ++i)
      {
        Output::dash("Ping received by process ", recv_proc[i]," \t ",
                     recv_dof[i], " \t " ,recv_mark[i], " \t ",
                     recv_pos[3*i]," ", recv_pos[3*i+1], " ", recv_pos[3*i + 2] );
      }
    }

    //free memory
    delete[] recv_proc;
    delete[] recv_dof;
    delete[] recv_pos;
    delete[] recv_mark;
    MPI_Comm_free(&mpi_comm_ping);
  }

  /* ****** */

  return local_dof;
  //TODO Gather ping's Output::print in one process which then prints the whole chunk of info.

}

void TParFECommunicator3D::CommUpdateMS(double *sol) const
{
  int i,j,k;

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

}

void TParFECommunicator3D::CommUpdateH1(double *sol) const
{
  int i,j,k;

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
}

void TParFECommunicator3D::CommUpdateH2(double *sol) const
{
  int i,j,k;

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
}


#endif














// =======================================================================
// @(#)ParFECommunicator3D.C
//
// Class:      TParFECommunicator3D
// Purpose:    Class containing all info needed for communication between subdomains
//
// Author:     Sashikumaar Ganesan (01.10.09)
//
// History:    Start of implementation 01.10.09 (Sashikumaar Ganesan)
//
// ======================================================================= 
#ifdef _MPI

#include "mpi.h"
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

TParFECommunicator3D::TParFECommunicator3D(MPI_Comm comm, TFESpace3D *fespace, TSquareStructure3D* Sqstruct)
{
 int N_U;

 Comm    = comm;
 FESpace = fespace;
 sqstruct = Sqstruct;
 MaxSubDomainPerDof = fespace->GetMaxSubDomainPerDof();

 if(MaxSubDomainPerDof<0)
 {
   printf("Error: SetMaxSubDomainPerDof in FeSpace before calling ParFECommunicator2D \n");
   MPI_Finalize();
   exit(0);
 }

 N_U = FESpace->GetN_DegreesOfFreedom();

 MPI_Allreduce(&N_U, &MaxN_LocalDofAllRank, 1, MPI_INT, MPI_MAX, Comm);

 if((TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR != 6 && TDatabase::ParamDB->SC_SMOOTHER_SCALAR==6) ||
     (TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR == 6 && TDatabase::ParamDB->SC_SMOOTHER_SCALAR!=6)){
   printf("use consistent smoother types\n (set 6 smoother for both smoother and coarese smoother if using 6 for any of them)\n");
   MPI_Finalize();
   exit(0);
 }
   
 if(TDatabase::ParamDB->SC_SMOOTHER_SCALAR==6){
   ConstructDofMap_light();
//    Color(N_CInt,ptrCInt,'i');
//    Color(N_CMaster,ptrCMaster,'m');
//    Color(N_CDept1,ptrCDept1,'D');
//    Color(N_CDept2,ptrCDept2,'d');
//    Color(N_CDept3,ptrCDept3,'x');
 }
 else{
   ConstructDofMap();
   SetFENeibCommunicationSteps();
 }
}

static int GetLocalIndex(int N, int *array, int val)
{
 int m=0;
 while(array[m] != val)
    m++;
 return m;
}

int TParFECommunicator3D::find_min(int *arr, int N, char *temp_arr){
  int i,j;
  int min = 1000000; 
  j = 0;
  for(i=0;i<N;i++){
    if(temp_arr[i] == 'x')   continue;
    if(arr[i]<min){
      min = arr[i];
      j = i;
    }
  }
  if(min == 1000000){
    printf("\n................this shudnt happen...................\n");
    MPI_Finalize();
    exit(0);
  }
  return j;
}

void TParFECommunicator3D::ConstructDofMap_light(){
//#################################################################  Variable Declarations  ##################################################################################//  
  int rank,size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  int i,j,k,aa,bb,N;
  int N_Cells, N_OwnCells, N_Active, N_Dof, N_LocDof, ID;
  int *GlobalNumbers, *BeginIndex, *LocalIndex, *DOF;
  int **MappingData;
  int temp,temp_globalno,temp_dofno;
  double global_start_time, start_time, end_time;
  
  TCollection *Coll;
  TBaseCell *cell;
//#################################################################  Variable Declarations  ##################################################################################//  

  
  Coll       = FESpace->GetCollection();
  N_Cells    = Coll->GetN_Cells();
  N_OwnCells = Coll->GetN_OwnCells();
  N_Dof      = FESpace->GetN_DegreesOfFreedom();
  N_Active   = FESpace->GetN_ActiveDegrees();
  
  start_time = MPI_Wtime();
  global_start_time = start_time;
  
  //check if the first N_OwnCells are OwnCells
  for(i=0;i<N_OwnCells;i++){
    cell = Coll->GetCell(i);
    if(cell->IsHaloCell()){ 
      printf("This shudnt happen---- Why am i(rank=%d) halo?? ----\n",rank);
      printf("rank=%d OwnCells infected\n",rank);break;
    }
  }
  
  //check if the remaining are halo cells
  for(i=N_OwnCells;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(!cell->IsHaloCell()){
      printf("This shudnt happen--- Why am i(rank=%d, cell_id=%d) own cell?? -----\n",rank,cell->GetSubDomainNo());
//       if(cell->IsOwnCell())
// 	printf("Recheck::Yes I am own cell\n");
//       else
// 	printf("Recheck::No I am not own cell\n");
      printf("rank=%d HaloCells infected\n",rank);break;
    }
  }
  
  /** *************************************************************************************/
  /** ARRAY CONTAINING GLOBAL DOF numbers (GLOBAL means LOCALLY GLOBAL for the subdomain) */ 
  /** BeginIndex: START ID IN GlobalNumbers ARRAY Cellwise                                */
  /** *************************************************************************************/
  BeginIndex = FESpace->GetBeginIndex();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  
  /** *************************************************************************************/
  /** LocalIndex :: Array containing global number of all Local cells                     */ 
  /** MappingData for mapping the dofs across processors                                  */
  /** *************************************************************************************/
  LocalIndex = new int[N_Cells];
  MappingData = new int*[2];
  for(i=0;i<2;i++){
    MappingData[i] = new int[N_Dof];
    memset(MappingData[i], 0, N_Dof*SizeOfInt);
  }
  
  /** *************************************************************************************/
  /** Master:: Array containing the rank to which the Dof belongs                         */ 
  /** Verify:: Array to help finding Master, Slave, Dependent, Halo                       */
  /** Verify:: ALL the Dofs are marked 'i'                                                */
  /** *************************************************************************************/
  Master = new int[N_Dof];
  for(i=0;i<N_Dof;i++) Master[i]=rank;
  char *Verify = new char[N_Dof];
  memset(Verify, 'i', N_Dof*sizeof(char)); 
  
  /** *************************************************************************************/
  /** Local Index Array is filled with GlobalNumbers of each cell                         */ 
  /** Verify::Dofs of Dependent cells are marked 'd'                                      */
  /** *************************************************************************************/
  for(i=0; i<N_Cells; i++){
     cell = Coll->GetCell(i);
     LocalIndex[i] = cell->GetGlobalCellNo();
     
     if(cell->IsDependentCell() && !cell->IsHaloCell()){
       DOF = GlobalNumbers + BeginIndex[i];
       N_LocDof = BeginIndex[i+1] - BeginIndex[i];
       
       for(j=0; j<N_LocDof; j++){
	 N = DOF[j];
	 Verify[N] = 'd';
       }
     }
  }

  /** *************************************************************************************/
  /** MappingData for dofs to be updated across domains are stored                        */ 
  /** Verify::Dofs of Halo cells not a part of dependent cells are marked 'h'             */
  /** *************************************************************************************/  
  for(i=N_OwnCells;i<N_Cells;i++){
    cell = Coll->GetCell(i);
    ID = cell->GetSubDomainNo();
    DOF = GlobalNumbers + BeginIndex[i];
    N_LocDof = BeginIndex[i+1] - BeginIndex[i];
    
    for(j=0; j<N_LocDof; j++){
      N = DOF[j];
      if(Verify[N] == 'i'){
	Verify[N] = 'h';
	Master[N] = ID;
      }
      
      if(ID < Master[N]) Master[N] = ID;
      MappingData[GLOBAL_NO][N]    = cell->GetGlobalCellNo();
      MappingData[DOF_NO][N]       = j; 
   }   
 }

//#################################################################  Master verification  ##################################################################################//  
 //------------------------------------------//
 /** Master DOF verification by other ranks **/
 //------------------------------------------//
 
 int *N_DOFtobeverified_otherRank;		//Array containing how many DOF's need to be verified from other processors(ranks) 
 int Total_DOFtobeverified_otherRank=0;		//Array containing how many DOF's need to be verified by this rank for other processors(ranks)
 int *N_DOFtobeverified_thisRank;		//Total no of DOF's that need to be verified from other ranks
 int Total_DOFtobeverified_thisRank=0;		//Total no of DOF verified by this rank
 int **sendbuf;
 int **recvbuf;
 int *verrecvbuf;
 int *versendbuf;
 int *temp_arr;
 
 sendbuf  = new int*[2];
 recvbuf  = new int*[2];
 sdispl   = new int[size];
 rdispl   = new int[size];
 temp_arr = new int[size];

 N_DOFtobeverified_otherRank = new int[size];
 memset(N_DOFtobeverified_otherRank,0,size*SizeOfInt);
 N_DOFtobeverified_thisRank = new int[size];
 memset(N_DOFtobeverified_thisRank,0,size*SizeOfInt); 
 
 //only the dofs in HALO cells (not dependent cells) needs to be verified
 for(i=0;i<N_Dof;i++){
    if(Verify[i] == 'h'){
      N_DOFtobeverified_otherRank[Master[i]]++;
      Total_DOFtobeverified_otherRank++;
    }
 }
 
 MPI_Alltoall(N_DOFtobeverified_otherRank,1, MPI_INT,N_DOFtobeverified_thisRank ,1, MPI_INT, Comm);
 for(i=0;i<size;i++)
    Total_DOFtobeverified_thisRank +=  N_DOFtobeverified_thisRank[i];
 
 for(i=0;i<2;i++){
   sendbuf[i] = new int[Total_DOFtobeverified_otherRank];
   recvbuf[i] = new int[Total_DOFtobeverified_thisRank];
 }
 
  //Send dof info to other processors and verify for master
  sdispl[0] = 0;
  for(i=1;i<size;i++){
    sdispl[i] = N_DOFtobeverified_otherRank[i-1] + sdispl[i-1];
  }
      
  rdispl[0] = 0;
  for(i=1;i<size;i++)
    rdispl[i] = N_DOFtobeverified_thisRank[i-1] + rdispl[i-1];
  
  memcpy (temp_arr,sdispl, size*SizeOfInt );
  
  for(i=0;i<N_Dof;i++){
     if(Verify[i] == 'h'){
       //if(rank==0) printf("%d ",temp_arr[Master[i]]);
       sendbuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
       sendbuf[DOF_NO][temp_arr[Master[i]]] = MappingData[DOF_NO][i];
       temp_arr[Master[i]]++;
     }
  }
  
  MPI_Alltoallv(sendbuf[GLOBAL_NO],N_DOFtobeverified_otherRank, sdispl, MPI_INT,recvbuf[GLOBAL_NO], N_DOFtobeverified_thisRank, rdispl, MPI_INT, Comm);
  MPI_Alltoallv(sendbuf[DOF_NO],   N_DOFtobeverified_otherRank, sdispl, MPI_INT,recvbuf[DOF_NO],    N_DOFtobeverified_thisRank, rdispl, MPI_INT, Comm);
  
  verrecvbuf  	= new int[Total_DOFtobeverified_thisRank];
  versendbuf  	= new int[Total_DOFtobeverified_otherRank];

  for(i=0;i<Total_DOFtobeverified_thisRank;i++){
    temp = GetLocalIndex(N_Cells,LocalIndex,recvbuf[GLOBAL_NO][i]);
    verrecvbuf[i] = Master[GlobalNumbers[temp * N_LocDof + recvbuf[DOF_NO][i]]];
  }
  
  MPI_Alltoallv(verrecvbuf, N_DOFtobeverified_thisRank, rdispl, MPI_INT,versendbuf,N_DOFtobeverified_otherRank, sdispl, MPI_INT, Comm);

  for(i=0;i<Total_DOFtobeverified_otherRank;i++){
    temp_globalno = sendbuf[GLOBAL_NO][i];
    temp_dofno    = sendbuf[DOF_NO][i];
    temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
    temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
    if(Verify[temp] != 'h')
      printf("Error : This degree of Freedom (%d) didn't require verification\n",temp);
    Master[temp] = versendbuf[i];
  }
  
 //-----------------------------------------------------//
 /** Master DOF verification by other ranks completed  **/
 //-----------------------------------------------------//
  
  delete [] verrecvbuf;                       verrecvbuf = NULL;
  delete [] versendbuf;                       versendbuf = NULL; 
  delete [] N_DOFtobeverified_otherRank;      N_DOFtobeverified_otherRank = NULL;
  delete [] N_DOFtobeverified_thisRank;       N_DOFtobeverified_thisRank  = NULL;
  
  for(i=0;i<2;i++){
    delete [] sendbuf[i];                     sendbuf[i] = NULL;
    delete [] recvbuf[i];                     recvbuf[i] = NULL;
  }
  delete [] sendbuf;                          sendbuf = NULL;
  delete [] recvbuf;                          recvbuf = NULL;  
  
  end_time = MPI_Wtime();
  if(rank == 0)
    printf("total time taken for master verification = %lf\n",end_time-start_time);
//#################################################################  Master verification  ##################################################################################//   

    {  
//################################################################# Redistribution of interface dofs #######################################################################//
//*/
  start_time = MPI_Wtime();
//---------------------------------------------------------------------------------------------------------------------------------------------------------------//  
  int total_interface_dofs        = 0;                    //these are the dofs which lie on the interface of sub domains (total over all sub domains)
  int N_own_interface_dofs        = 0;                    //these are the own interface dofs
  int N_interface_dofs            = 0;                    //these are the interface dofs including shared ones
  int T_interface_dofs            = 0;                    //these are the total of N_interface_dofs over all sub domains
  int total_own_dofs              = 0;                    //these are the dofs which belongs to my rank (total implies total over all sub domains)
  int start = 0,           end    = 0;                    //temporary variables used for numbering the dofs
  int max_n_ranks_interface_dofs  = 0;                    //maximum number of ranks sharing any interface dof
  
  int *GlobalDofNo                = new int[N_Dof];       //this array is used to assign unique global dof no. over all sub domains 
      memset(GlobalDofNo, -1, N_Dof*sizeof(int));         //all is set to -1 for a default value
  int *N_Dof_Slave                = new int[size];        //array of number of slave interface dofs to be verified by other ranks
      memset(N_Dof_Slave, 0, size*sizeof(int));
  int *N_Dof_Master               = new int[size];        //array of interface dofs to be verified by this rank for other ranks 
  int **MasterPos                 = new int*[2];          //a double array for storing the position info for master interface dofs
  int **SlavePos                  = new int*[2];          //a double array for storing the position info for slave interface dofs 
  int *masterInfo, *slaveInfo;                            //these are used to update the global dof no for interface master-slave dofs
  int *GlobalDofNo_interface;                             //this array contains only the globaldof no. of interface dofs
  int *all_own_dofs_info          = new int[size];        //array of N_own dofs by each rank
  int *all_interface_dofs_info    = new int[size];        //array of N_own interface dofs by each rank
  int *all_T_interface_dofs_info  = new int[size];        //array of interface dofs(including shared ones) by each rank
  int *all_GlobalDofNo;                                   //array of global interface dof no from all ranks
  int *N_ranks_per_interface_dofs;                        //array of count of number of ranks sharing an interface dof
  int *N_allocated_masters        = new int[size];        //array containing number of masters(interface dofs) allocated to each processors
  char *tempc                     = new char[size];       //temporary array
  char **Master_Table             = new char*[size];      //a Table to mark the shared interface dofs

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//   
  
  //mark all (interface) dofs as 'z'
  for(i=N_OwnCells;i<N_Cells;i++){
    cell = Coll->GetCell(i);
    ID = cell->GetSubDomainNo();
    DOF = GlobalNumbers + BeginIndex[i];
    N_LocDof = BeginIndex[i+1] - BeginIndex[i];
    
    //z is both interface and dependent dof, while d is only dependent dof
    for(j=0; j<N_LocDof; j++){
      N = DOF[j];
      if(Verify[N] == 'd'){
	Verify[N] = 'z';
      }
    }
  }
  
  //count total number of OWN interface_dofs and own dofs
  N_OwnDof = 0;
  for(i=0;i<N_Dof;i++){
    if(Verify[i] == 'z')
      if(Master[i] == rank)
        N_own_interface_dofs++;
    
    if(Master[i] == rank)
      N_OwnDof++;
  }
  
  //Send the count info to all ranks
  MPI_Allgather(&N_own_interface_dofs, 1, MPI_INT, all_interface_dofs_info, 1, MPI_INT, Comm);
  MPI_Allgather(&N_OwnDof,             1, MPI_INT, all_own_dofs_info,       1, MPI_INT, Comm);
    
  for(aa=0;aa<size;aa++){
    total_interface_dofs += all_interface_dofs_info[aa];
    total_own_dofs       += all_own_dofs_info[aa];
  }
  
  //start- numbering position for interface dofs
  aa=0;
  while(aa<rank){
    start += all_interface_dofs_info[aa];
    aa++;
  } 
  //printf("###########   rank = %d       start=%d\n",rank,start);
  
  Total_DOFtobeverified_otherRank = 0;                                               //this counts the total number of interface dofs(slave) that needs to be verified from other ranks(master) 
  Total_DOFtobeverified_thisRank  = 0;                                               //this counts the total number of interface dofs(masters) that needs to be verified for other ranks(slaves)
  
  //number the own interface dofs
  //count the total dofs to be verified by this and other rank
  for(i=0;i<N_Dof;i++){
    if(Verify[i]=='z'){
      if(Master[i]==rank){
	GlobalDofNo[i] = start;
	start++;
      }
      else{
	N_Dof_Slave[Master[i]]++;
	Total_DOFtobeverified_otherRank++;
      }
    }
  }
 
 MPI_Alltoall(N_Dof_Slave, 1, MPI_INT, N_Dof_Master, 1, MPI_INT, Comm);
 for(i=0;i<size;i++)
    Total_DOFtobeverified_thisRank +=  N_Dof_Master[i];
 
  for(i=0;i<2;i++){
   MasterPos[i] = new int[Total_DOFtobeverified_thisRank];
   SlavePos[i]  = new int[Total_DOFtobeverified_otherRank];
  }
  
  masterInfo  = new int[Total_DOFtobeverified_thisRank];
  slaveInfo   = new int[Total_DOFtobeverified_otherRank];
  
  rdispl[0] = 0;
  sdispl[0] = 0;
  for(i=1;i<size;i++){
    rdispl[i] = rdispl[i-1] + N_Dof_Master[i-1];
    sdispl[i] = sdispl[i-1] + N_Dof_Slave[i-1];                                 //slaves will be sent to gather info about their global dof no.
  }
 
  memcpy (temp_arr, sdispl, size*SizeOfInt );
 
  //store the position info of the slave dofs
  for(i=0;i<N_Dof;i++){
    if(Verify[i]=='z'){
      if(Master[i]!=rank){
       SlavePos[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
       SlavePos[DOF_NO][temp_arr[Master[i]]]    = MappingData[DOF_NO][i];
       temp_arr[Master[i]]++;
      }
    }
  }
  
//   for(aa=0;aa<size;aa++){
//     if(rank==aa){
//       printf("rank = %d\n",rank);	
//       for(i=0;i<size;i++){
// 	printf("N_Dof_Master[%d] = %d\t N_Dof_Slave[%d]=%d\n",i,N_Dof_Master[i],i,N_Dof_Slave[i]);
//       }
//     }
//     MPI_Barrier(Comm);
//   }
//   exit(0);

  //send the slave info to masters
  MPI_Alltoallv(SlavePos[GLOBAL_NO], N_Dof_Slave, sdispl, MPI_INT, MasterPos[GLOBAL_NO], N_Dof_Master, rdispl, MPI_INT, Comm);
  MPI_Alltoallv(SlavePos[DOF_NO],    N_Dof_Slave, sdispl, MPI_INT, MasterPos[DOF_NO],    N_Dof_Master, rdispl, MPI_INT, Comm);
 
  start = 0;
  for(aa=0;aa<size;aa++){
    end = start + N_Dof_Master[aa];
    for(i=start;i<end;i++){
      //get the cell with global cell no = MasterPos[GLOBAL_NO][i]
      temp = GetLocalIndex(N_Cells,LocalIndex,MasterPos[GLOBAL_NO][i]);
      
      //the master of this dof, i should be rank
      if(Master[GlobalNumbers[temp * N_LocDof + MasterPos[DOF_NO][i]]] != rank){
	printf("....................Wrong Master.....................\n");
	MPI_Finalize();
	exit(0);
      }
      
      //this dof should be a category z dof
      if(Verify[GlobalNumbers[temp * N_LocDof + MasterPos[DOF_NO][i]]] != 'z'){
	printf("....................Wrong Dof.........................\n");
	MPI_Finalize();
	exit(0);
      }
      
      //this dof should have a GlobalDofNo
      if(GlobalDofNo[GlobalNumbers[temp * N_LocDof + MasterPos[DOF_NO][i]]] == -1){
	printf("....................Wrong GlobalDofNo.........................\n");
	MPI_Finalize();
	exit(0);
      }
      
      masterInfo[i] = GlobalDofNo[GlobalNumbers[temp * N_LocDof + MasterPos[DOF_NO][i]]];
    }
    start = end;
  }

  MPI_Alltoallv(masterInfo, N_Dof_Master, rdispl, MPI_INT, slaveInfo, N_Dof_Slave, sdispl, MPI_INT, Comm);

  start = 0;
  for(aa=0;aa<size;aa++){
    end = start + N_Dof_Slave[aa];
    for(i=start;i<end;i++){
      temp_globalno = SlavePos[GLOBAL_NO][i];
      temp_dofno    = SlavePos[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp * N_LocDof + temp_dofno];
      
      //the master of this dof, i should be aa
      if(Master[temp] != aa){
	printf("...................2.Wrong Master.....master[%d]=%d....rank=%d............\n",temp,Master[temp],rank);
	MPI_Finalize();
	exit(0);
      }
      
      //this dof should be a category z dof
      if(Verify[temp] != 'z'){
	printf("...................2.Wrong Dof.........................\n");
	MPI_Finalize();
	exit(0);
      }
      
      //this dof should not have a GlobalDofNo
      if(GlobalDofNo[temp] != -1){
	printf("...................2.Wrong GlobalDofNo.........................\n");
	MPI_Finalize();
	exit(0);
      }
      //assign the global dof no to the slave obtained from master
      GlobalDofNo[temp] = slaveInfo[i];
    }
    start = end;
  }
  
  //now we mark the own dofs apart from the interface ones
  aa=0;
  start = total_interface_dofs;
  while(aa<rank){
    start += all_own_dofs_info[aa];
    aa++;
  }

  for(i=0;i<N_Dof;i++){
    if(Verify[i] != 'z'){
      if(Master[i]==rank){
	GlobalDofNo[i] = start;
	start++;
      }
    }
    
    if(Master[i]==rank){
      if(GlobalDofNo[i] == -1){
	printf("......................This shudnt happen.................\n");
	MPI_Finalize();
	exit(0);
      }
    }
    
    if(Verify[i]=='z'){
      if(GlobalDofNo[i]>=total_interface_dofs){
	printf("......................This shudnt happen.................\n");
	MPI_Finalize();
	exit(0);
      }
      
      N_interface_dofs++;
    }
  }
  
  GlobalDofNo_interface = new int[N_interface_dofs];
  aa = 0;
  for(i=0;i<N_Dof;i++){
     if(Verify[i]=='z'){
      GlobalDofNo_interface[aa] = GlobalDofNo[i];
      aa++;
     }
  }
  
  MPI_Allreduce(&N_interface_dofs, &T_interface_dofs, 1, MPI_INT, MPI_SUM, Comm);
  
  delete [] N_Dof_Master;
  delete [] N_Dof_Slave;
  
  for(i=0;i<2;i++){
    delete [] MasterPos[i];                    MasterPos[i] = NULL;
    delete [] SlavePos[i];                     SlavePos[i]  = NULL;
  }
  delete [] MasterPos;                         MasterPos = NULL;
  delete [] SlavePos;                          SlavePos  = NULL;  
  
  delete [] masterInfo;                        masterInfo = NULL;
  delete [] slaveInfo;                         slaveInfo  = NULL;
  
  MPI_Allgather(&N_interface_dofs, 1, MPI_INT, all_T_interface_dofs_info, 1, MPI_INT, Comm);
  all_GlobalDofNo = new int[T_interface_dofs];
  
//   for(aa=0;aa<size;aa++){
//     if(rank==aa){
//       printf("rank = %d      #total interface(own) dofs = %d       #interface dofs(own) = %d\n",rank,total_interface_dofs,N_own_interface_dofs);
//       printf("              #total interface dofs      = %d      #interface dofs      = %d\n",T_interface_dofs,N_interface_dofs);
//     }
//     MPI_Barrier(MPI_COMM_WORLD);
//   }
  
  rdispl[0] = 0;
  for(aa=1;aa<size;aa++)
    rdispl[aa] = rdispl[aa-1] + all_T_interface_dofs_info[aa-1];
  MPI_Allgatherv(GlobalDofNo_interface, N_interface_dofs, MPI_INT, all_GlobalDofNo, all_T_interface_dofs_info, rdispl, MPI_INT, Comm);
  
  //table is created only over the total own interface dofs over the sub domains
  for(i=0;i<size;i++){
    Master_Table[i] = new char[total_interface_dofs];
    memset(Master_Table[i], 'x', total_interface_dofs*sizeof(char));
  }
  
  start = 0;end = 0;
  for(i=0;i<size;i++){
    end += all_T_interface_dofs_info[i];
    for(aa=start;aa<end;aa++){
      Master_Table[i][all_GlobalDofNo[aa]] = 'y';
    }
    start = end; 
  }
  
  N_ranks_per_interface_dofs = new int[total_interface_dofs];
  memset(N_ranks_per_interface_dofs, 0, total_interface_dofs*sizeof(int));
  
  for(j=0;j<total_interface_dofs;j++){
    //check how many ranks share the interface dof
    for(i=0;i<size;i++){
      if(Master_Table[i][j] == 'y'){
	N_ranks_per_interface_dofs[j]++;
      }
    }
    //compute the max munber of rank sharing any interface dof
    if(max_n_ranks_interface_dofs<N_ranks_per_interface_dofs[j])
      max_n_ranks_interface_dofs = N_ranks_per_interface_dofs[j];
  }
  
//   if(rank==1)
//   for(i=0;i<size;i++){
//     for(aa=0;aa<20;aa++){
//       printf("%c\t",Master_Table[i][aa]);
//     }
//     printf("\n");
//   }
  
  for(j=0;j<total_interface_dofs;j++){
  //a interface dof should have at least one neighbour from other rank
    if(N_ranks_per_interface_dofs[j] == 1){
	printf("......................This shudnt happen.....j=%d............\n",j);
	MPI_Finalize();
	exit(0);
    }
  }
 
  memset(N_allocated_masters, 0, size*sizeof(int));
  for(i=0; i<N_Dof; i++){
    if(Verify[i]=='z'){
      N_allocated_masters[Master[i]]++;
    }
  }
  
  
//   for(aa=0;aa<size;aa++){
//      if(rank==aa){
//        printf("\nrank=%d\n",rank);
//        for(i=0;i<size;i++)
// 	 printf("Before Redistribution::             N_allocated_masters[%d] = %d\n",i,N_allocated_masters[i]);
//      }
//      MPI_Barrier(Comm);
//      //break;
//   }
  
  memset(N_allocated_masters, 0, size*sizeof(int));
  int min = 0;
  //a interface dof should have at least one neighbour from other rank
  for(i=2; i<=max_n_ranks_interface_dofs;i++){
    for(j=0;j<total_interface_dofs;j++){
      if(N_ranks_per_interface_dofs[j] == i){
	
	for(aa=0;aa<size;aa++){
	  tempc[aa] = Master_Table[aa][j];
	}
	
	min = find_min(N_allocated_masters,size,tempc);
	N_allocated_masters[min]++;
	
	for(aa=0;aa<N_Dof;aa++){
	  if(GlobalDofNo[aa] == j){
	    Master[aa] = min;
	  }
	}
	N_ranks_per_interface_dofs[j] = -1;
      }
    }
  }
  
//   for(aa=0;aa<size;aa++){
//      if(rank==aa){
//        printf("\nrank=%d\n",rank);
//        for(i=0;i<size;i++)
// 	 printf("After Redistribution::             N_allocated_masters[%d] = %d\n",i,N_allocated_masters[i]);
//      }
//      MPI_Barrier(Comm);
//   }
  
  memset(N_allocated_masters, 0, size*sizeof(int));
  for(i=0; i<N_Dof; i++){
    if(Verify[i]=='z'){
      N_allocated_masters[Master[i]]++;
      
      if(GlobalDofNo[i] == -1){
	printf("......................This shudnt happen.................\n");
	MPI_Finalize();
	exit(0);
      }
    }
  }
  
//   for(aa=0;aa<size;aa++){
//      if(rank==aa){
//        printf("\nrank=%d\n",rank);
//        for(i=0;i<size;i++)
// 	 printf("N_allocated_masters[%d] = %d\n",i,N_allocated_masters[i]);
//      }
//      MPI_Barrier(Comm);
//      //break;
//   }
  
  delete [] all_GlobalDofNo;                  all_GlobalDofNo            = NULL;
  for(i=0;i<size;i++) 
    delete [] Master_Table[i];                Master_Table[i]            = NULL;
  delete [] N_ranks_per_interface_dofs;       N_ranks_per_interface_dofs = NULL;     
  delete [] N_allocated_masters;              N_allocated_masters        = NULL;
  delete [] temp_arr;                         temp_arr                   = NULL;
  
  end_time = MPI_Wtime();
  if(rank == 0)
    printf("Total Time Taken for Redistribution of master dofs = %lf\n",end_time-start_time);
  //*/
//   MPI_Finalize();
//   exit(0);
//################################################################# Redistribution of interface dofs #######################################################################//
    }

//#################################################################  Marking The Dofs ######################################################################################//     
 //------------------------------------------------------------------//
 /** Gather information about other DOFs ::                          */
 /**    ## dofs marked as 'i'             -->independent             */
 /**    ## dofs marked as 's'             -->slave                   */
 /**    ## dofs marked as 'D'             -->dependent_type1         */
 /**    ## dofs marked as 'd'             -->dependent_type2         */
 /**    ## dofs marked as 'H'             -->halo_type1              */
 /**    ## dofs marked as 'h'             -->halo_type2              */
 //------------------------------------------------------------------//
 start_time = MPI_Wtime();

 memset(Verify, 'i', N_Dof*sizeof(char));		//all dofs marked as independent
 //all dofs in halo cell marked as halo_type2(unused)
 for(i=N_OwnCells;i<N_Cells;i++){
   cell     = Coll->GetCell(i);
   DOF      = GlobalNumbers + BeginIndex[i];
   N_LocDof = BeginIndex[i+1] - BeginIndex[i];
   
   for(j=0;j<N_LocDof;j++){
     N = DOF[j];
     Verify[N]='h';
   }
 }
 
 for(i=0;i<N_OwnCells;i++){
   cell = Coll->GetCell(i);
   //now mark dofs in dependent cells
   if(cell->IsDependentCell()){
     DOF      = GlobalNumbers + BeginIndex[i];
     N_LocDof = BeginIndex[i+1] - BeginIndex[i];
     
     for(j=0;j<N_LocDof;j++){
       N = DOF[j];
       if(Verify[N]!='h')	//if dof is not marked halo then it is dependent dof
	 Verify[N]='d';
     }
   }
 }
 
 for(i=0;i<N_OwnCells;i++){
   cell = Coll->GetCell(i);
   //now mark dofs in dependent cells
   if(cell->IsDependentCell()){
     DOF      = GlobalNumbers + BeginIndex[i];
     N_LocDof = BeginIndex[i+1] - BeginIndex[i];
     
     for(j=0;j<N_LocDof;j++){
       N = DOF[j];
       if(Verify[N] == 'h'){
	 if(Master[N]!=rank)
	   Verify[N]='s';	//if dof is marked halo and master of it is some other proc then it is a slave
	 else
	   Verify[N]='m';	//else it is a master
       }
     }
   }
 }
 
 //mark dependent type1 & type2
 bool flag = false;
 for(i=0;i<N_OwnCells;i++){
   cell = Coll->GetCell(i);
   //now mark dofs in dependent cells
   if(cell->IsDependentCell()){
     DOF      = GlobalNumbers + BeginIndex[i];
     N_LocDof = BeginIndex[i+1] - BeginIndex[i];
     
     for(j=0;j<N_LocDof;j++){
       N = DOF[j];
       if(Master[N]!=rank){
	 flag = true;
	 break;
       }
     }
     //dependent dofs connected to slave dofs are marked as type1(D) else type2(d)
     //type1 dofs must be smoothed first as they are the halo dofs (type1 i.e. useful) to other procs
     if(flag==true){
       for(j=0;j<N_LocDof;j++){
	 N = DOF[j];
	 if(Verify[N]=='d')
	   Verify[N]='D';
       }
       flag=false;
     }
   }
 }
 
 //mark halo type1(H)-->useful & type2(h)
 flag = false;
 for(i=N_OwnCells;i<N_Cells;i++){
   DOF      = GlobalNumbers + BeginIndex[i];
   N_LocDof = BeginIndex[i+1] - BeginIndex[i];
   
   for(j=0;j<N_LocDof;j++){
     N = DOF[j];
     if(Master[N]==rank){
       flag = true;
       break;
     }
   }
   //halo dofs connected to master dofs are marked as type1(H) else type2(h)
   //type2 halo dofs are not required in smoothing operations
   if(flag==true){
     for(j=0;j<N_LocDof;j++){
       N = DOF[j];
       if(Verify[N]=='h')
	 Verify[N]='H';  
     }
     flag=false;
   }
 }
 
 //identify dependent3 dofs
 //mark dependent3 type dofs as 'x'
 flag = false;
 for(i=0;i<N_OwnCells;i++){
   DOF      = GlobalNumbers + BeginIndex[i];
   N_LocDof = BeginIndex[i+1] - BeginIndex[i];
   
   cell = Coll->GetCell(i);
   if(cell->IsDependentCell()) continue;
   
   for(j=0;j<N_LocDof;j++){
     N = DOF[j];
     if(Verify[N] == 'd' || Verify[N] == 'D'){
       flag = true;
       break;
     }
   }
   //interior dofs connected to dependent dofs are marked as dependent3 type i.e. 'x'
   if(flag==true){
     for(j=0;j<N_LocDof;j++){
       N = DOF[j];
       if(Verify[N]=='i')
	 Verify[N]='x';  
     }
     flag=false;
   }
 }
 
 DofMarker = Verify;
 
 end_time = MPI_Wtime();
 if(rank == 0)
    printf("Total Time Taken for marking the dofs = %lf",end_time-start_time);
 
//#################################################################  Marking The Dofs ######################################################################################// 

 
//#############################################################  Mapper for Master Dof are set ##############################################################################//  
 start_time = MPI_Wtime();
 
 int **SlaveBuf,**MasterBuf;
 int *temp_arrH1,*temp_arrH2;
 
 N_InterfaceM = 0;      N_InterfaceS  = 0;
 N_Slave      = 0;     
 N_OwnDof     = 0; 
 N_Master     = 0;
 N_Int        = 0; 
 N_Dept       = 0;      N_Dept1 = 0;    N_Dept2 = 0;    N_Dept3 = 0;
 N_Halo       = 0;      N_Halo1 = 0;    N_Halo2 = 0;    
 
 N_DofSend    = new int[size];
 N_DofSendMS  = new int[size]; 
 N_DofSendH2  = new int[size];
 N_DofSendH1  = new int[size];
 
 N_DofRecv    = new int[size];
 N_DofRecvMS  = new int[size]; 
 N_DofRecvH2  = new int[size];
 N_DofRecvH1  = new int[size];
 
 memset(N_DofRecvMS  ,0,size*SizeOfInt);
 memset(N_DofRecvH1,0,size*SizeOfInt);
 memset(N_DofRecvH2,0,size*SizeOfInt);
 
 temp_arr   = new int[size];
 temp_arrH1 = new int[size];
 temp_arrH2 = new int[size];
 
 SlaveBuf  = new int*[2];
 MasterBuf = new int*[2];
 
 for(N=0;N<N_Dof;N++){
   
   if(Master[N] != rank){
     N_Slave++;
     N_DofRecv[Master[N]]++;
     
     if(Verify[N] == 's'){
       N_InterfaceS++;
       N_DofRecvMS[Master[N]]++;
     }
     else if(Verify[N] == 'H'){
       N_Halo1++;       
       N_DofRecvH1[Master[N]]++;
     }
     else{
       N_Halo2++;
       N_DofRecvH2[Master[N]]++;
     }
   }
   else{
     N_OwnDof++;
     N_Master++;
     if(Verify[N] == 'm')
       N_InterfaceM++;
     else if(Verify[N] == 'D')
       N_Dept1++;
     else if(Verify[N] == 'd')
       N_Dept2++;
     else if(Verify[N] == 'x')
       N_Dept3++;
     else
       N_Int++;
   }
 }
 
 N_Halo = N_Halo1 + N_Halo2;
 N_Dept = N_Dept1 + N_Dept2 + N_Dept3;
 
  //DEBUG
  if(size<5);
  for(aa=0;aa<size;aa++){
     if(rank==aa){
          printf("\nRank::%d\n",rank);
	  printf("N_Dof               = %d\t N_OwnDof           = %d\n", N_Dof, N_OwnDof);
	  printf("N_Master(Total)     = %d\t N_Slave(Total)     = %d\n", N_Master, N_Slave);
	  printf("N_Master(Interface) = %d\t N_Slave(Interface) = %d\n", N_InterfaceM, N_InterfaceS);
	  printf("N_Halo              = %d\t N_Halo1            = %d\t   N_Halo2  = %d\n",N_Halo,N_Halo1,N_Halo2);
	  printf("N_Depndt            = %d\t N_Dept1            = %d\t   N_Dept2  = %d\t  N_Dept3  = %d\n",N_Dept,N_Dept1,N_Dept2,N_Dept3);
	  printf("N_Indpdt            = %d\n",N_Int); 
     }
     MPI_Barrier(MPI_COMM_WORLD);
  }//verified

 rdispl[0] = 0; sdispl[0] = 0; 
 sdisplMS  = new int[size];   sdisplMS[0] = 0;  
 rdisplMS  = new int[size];   rdisplMS[0] = 0;
 sdisplH1  = new int[size];   sdisplH1[0] = 0;  
 rdisplH1  = new int[size];   rdisplH1[0] = 0;
 sdisplH2  = new int[size];   sdisplH2[0] = 0;
 rdisplH2  = new int[size];   rdisplH2[0] = 0;
 
 OwnDofs = new int[N_OwnDof];
 
 N_SendDofMS = 0; 
 N_SendDofH1 = 0; 
 N_SendDofH2 = 0;

 MPI_Alltoall(N_DofRecv  , 1, MPI_INT,N_DofSend   , 1, MPI_INT, Comm);
 MPI_Alltoall(N_DofRecvMS, 1, MPI_INT,N_DofSendMS , 1, MPI_INT, Comm);
 MPI_Alltoall(N_DofRecvH1, 1, MPI_INT,N_DofSendH1 , 1, MPI_INT, Comm);
 MPI_Alltoall(N_DofRecvH2, 1, MPI_INT,N_DofSendH2 , 1, MPI_INT, Comm);
 
 for(i=1;i<size;i++){
   rdispl[i]   = rdispl[i-1]   + N_DofRecv[i-1];
   rdisplMS[i] = rdisplMS[i-1] + N_DofRecvMS[i-1];
   rdisplH1[i] = rdisplH1[i-1] + N_DofRecvH1[i-1];
   rdisplH2[i] = rdisplH2[i-1] + N_DofRecvH2[i-1];
   
   sdispl[i]   = sdispl[i-1]   + N_DofSend[i-1];
   sdisplMS[i] = sdisplMS[i-1] + N_DofSendMS[i-1];
   sdisplH1[i] = sdisplH1[i-1] + N_DofSendH1[i-1];
   sdisplH2[i] = sdisplH2[i-1] + N_DofSendH2[i-1];
 }
 
 for(i=0;i<size;i++){
   N_SendDofMS += N_DofSendMS[i];
   N_SendDofH1 += N_DofSendH1[i];
   N_SendDofH2 += N_DofSendH2[i];
 }
 N_SendDof = N_SendDofMS + N_SendDofH1 + N_SendDofH2;
 
 memcpy (temp_arr,   rdisplMS, size*SizeOfInt );
 memcpy (temp_arrH1, rdisplH1, size*SizeOfInt );
 memcpy (temp_arrH2, rdisplH2, size*SizeOfInt );
 
 DofSend    = new int[N_SendDofMS+N_SendDofH1+N_SendDofH2];
 DofSendMS  = DofSend;
 DofSendH1  = DofSend + N_SendDofMS;
 DofSendH2  = DofSend + N_SendDofMS + N_SendDofH1;
 
 DofRecv    = new int[N_InterfaceS+N_Halo1+N_Halo2];
 DofRecvMS  = DofRecv;
 DofRecvH1  = DofRecv + N_InterfaceS;
 DofRecvH2  = DofRecv + N_InterfaceS + N_Halo1;
 
 for(i=0;i<2;i++){
   if(N_InterfaceS>0)
     SlaveBuf[i]  = new int[N_InterfaceS];
   if(N_SendDofMS>0)
     MasterBuf[i] = new int[N_SendDofMS];
   //memset (SlaveBuf[i] , 0, size*SizeOfInt);
   //memset (MasterBuf[i], 0, size*SizeOfInt);
 }
 
 Reorder = new int[N_Dof];
 NewGN = new int[N_Dof];
 int m = 0;

 int Mstr  = 0;                  Reorder_M  = Reorder;
 int Indpt = N_InterfaceM;       Reorder_I  = Reorder + Indpt;
 int Dept1 = Indpt + N_Int;      Reorder_D1 = Reorder + Dept1;
 int Dept2 = Dept1 + N_Dept1;    Reorder_D2 = Reorder + Dept2;
 int Dept3 = Dept2 + N_Dept2;    Reorder_D3 = Reorder + Dept3;
 int Slv   = Dept3 + N_Dept3;
 int Hl1   = Slv   + N_InterfaceS;
 int Hl2   = Hl1   + N_Halo1;
 int ts = 0, th1 = 0, th2 = 0, ti = 0, tm = 0, td1 = 0, td2 = 0, td3 = 0;
 
 for(i=0;i<N_Dof;i++){
   //slave dofs
   if(Master[i] != rank){
     if(Verify[i] == 's'){
       SlaveBuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
       SlaveBuf[DOF_NO]   [temp_arr[Master[i]]] = MappingData[DOF_NO]   [i];
   
       if(SlaveBuf[0][temp_arr[Master[i]]]<0 || Master[i]>=size || Master[i]<0 || temp_arr[Master[i]]>N_Dof){
	printf("\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"); 
	printf("\n.................This shudnt happen!!!(rank=%d,cell=%d,dof=%d)\n",rank,MappingData[GLOBAL_NO][i],i);
	exit(0);
       }
   
       DofRecvMS[temp_arr[Master[i]]] = i;
       temp_arr[Master[i]]++;
       
       {
         Reorder[Slv] = i;
	 NewGN[i]     = ts++;//Slv;
         Slv++;
       }
     }
     else if(Verify[i] == 'H'){
       {
	 Reorder[Hl1] = i;
	 NewGN[i]     = th1++;//Hl1; 
	 Hl1++;
       }
     }
     else{
       {
	 Reorder[Hl2] = i;
	 NewGN[i]     = th2++;//Hl2;
	 Hl2++;
       }
     }
   }
   else{
     OwnDofs[m++]=i;
       {
       if(Verify[i]=='m'){
	 Reorder[Mstr] = i;
	 NewGN[i]      = tm++;//Mstr;
	 Mstr++;
       }
       else if(Verify[i]=='D'){
	 Reorder[Dept1] = i;
	 NewGN[i]       = td1++;//Dept1;
	 Dept1++;
       }
       else if(Verify[i]=='d'){
	 Reorder[Dept2] = i;
	 NewGN[i]       = td2++;//Dept2;
	 Dept2++;
       }
       else if(Verify[i]=='x'){
	 Reorder[Dept3] = i;
	 NewGN[i]       = td3++;//Dept3;
	 Dept3++;
       }
       else{
	 Reorder[Indpt] = i;
	 NewGN[i]       = ti++;//Indpt;
	 Indpt++;
       }
     }
   }
 }
  
 MPI_Alltoallv(SlaveBuf[GLOBAL_NO], N_DofRecvMS, rdisplMS, MPI_INT, MasterBuf[GLOBAL_NO], N_DofSendMS, sdisplMS, MPI_INT, Comm);	
 MPI_Alltoallv(SlaveBuf[DOF_NO],    N_DofRecvMS, rdisplMS, MPI_INT, MasterBuf[DOF_NO],    N_DofSendMS, sdisplMS, MPI_INT, Comm);
 
 for(i=0;i<N_SendDofMS;i++){
   temp_globalno = MasterBuf[GLOBAL_NO][i];
   temp_dofno    = MasterBuf[DOF_NO][i];
   temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
   temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
   DofSendMS[i] = temp; 
 }

 for(i=0;i<2;i++){
    if(N_InterfaceS>0){
      delete [] SlaveBuf[i];
      SlaveBuf[i] = NULL;
    }
    if(N_SendDofMS>0){
      delete [] MasterBuf[i];
      MasterBuf[i] = NULL;
    }
 }
 delete [] SlaveBuf;    SlaveBuf  = NULL;
 delete [] MasterBuf;   MasterBuf = NULL;
 delete [] temp_arr;    temp_arr  = NULL;  

 //*********************************************************************************************************************************//
 //                                                  Mapper for Halo Dof are set                                              //
 //*********************************************************************************************************************************//
 int **SlaveBufH1,**MasterBufH1;
 int **SlaveBufH2,**MasterBufH2;
 
 SlaveBufH1  = new int*[2];
 MasterBufH1 = new int*[2];
 
 SlaveBufH2  = new int*[2];
 MasterBufH2 = new int*[2];
 
 for(i=0;i<2;i++){
   if(N_Halo1>0)
     SlaveBufH1[i]  = new int[N_Halo1];
   if(N_SendDofH1>0)
     MasterBufH1[i] = new int[N_SendDofH1];
   //memset (SlaveBufH1[i] , 0, size*SizeOfInt);
   //memset (MasterBufH1[i], 0, size*SizeOfInt);
   if(N_Halo2>0)
     SlaveBufH2[i]  = new int[N_Halo2];
   if(N_SendDofH2>0)
     MasterBufH2[i] = new int[N_SendDofH2];
   //memset (SlaveBufH2[i] , 0, size*SizeOfInt);
   //memset (MasterBufH2[i], 0, size*SizeOfInt);
 }
 
 for(i=0;i<N_Dof;i++){
   if(Verify[i] == 'H'){
     if(Master[i] == rank){
       printf("\nThis should be a HALO_1 dof\n");
       exit(0);
     }
     
     SlaveBufH1[GLOBAL_NO][temp_arrH1[Master[i]]] = MappingData[GLOBAL_NO][i];
     SlaveBufH1[DOF_NO]   [temp_arrH1[Master[i]]] = MappingData[DOF_NO]   [i];
	  
     DofRecvH1[temp_arrH1[Master[i]]] = i;
     temp_arrH1[Master[i]]++;
    }
    else if(Verify[i] == 'h'){
     if(Master[i] == rank){
       printf("\nThis should be a HALO_2 dof\n");
       exit(0);
     }
     
     SlaveBufH2[GLOBAL_NO][temp_arrH2[Master[i]]] = MappingData[GLOBAL_NO][i];
     SlaveBufH2[DOF_NO]   [temp_arrH2[Master[i]]] = MappingData[DOF_NO]   [i];
	  
     DofRecvH2[temp_arrH2[Master[i]]] = i;
     temp_arrH2[Master[i]]++;
    }
 }
 
 MPI_Alltoallv(SlaveBufH1[GLOBAL_NO], N_DofRecvH1, rdisplH1, MPI_INT, MasterBufH1[GLOBAL_NO], N_DofSendH1, sdisplH1, MPI_INT, Comm);	
 MPI_Alltoallv(SlaveBufH1[DOF_NO],    N_DofRecvH1, rdisplH1, MPI_INT, MasterBufH1[DOF_NO],    N_DofSendH1, sdisplH1, MPI_INT, Comm);
 
 for(i=0;i<N_SendDofH1;i++){
   temp_globalno = MasterBufH1[GLOBAL_NO][i];
   temp_dofno    = MasterBufH1[DOF_NO][i];
   temp          = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
   temp          = GlobalNumbers[temp*N_LocDof + temp_dofno];
   DofSendH1[i]  = temp; 
 }
 
 for(i=0;i<2;i++){
   if(N_Halo1>0){
    delete [] SlaveBufH1[i];
    SlaveBufH1[i] = NULL;
   }
   if(N_SendDofH1>0){
     delete [] MasterBufH1[i];
     MasterBufH1[i] = NULL;
   }
 }
 delete [] SlaveBufH1;    SlaveBufH1  = NULL;
 delete [] MasterBufH1;   MasterBufH1 = NULL;
 delete [] temp_arrH1;    temp_arrH1  = NULL;


 MPI_Alltoallv(SlaveBufH2[GLOBAL_NO], N_DofRecvH2, rdisplH2, MPI_INT, MasterBufH2[GLOBAL_NO], N_DofSendH2, sdisplH2, MPI_INT, Comm);	
 MPI_Alltoallv(SlaveBufH2[DOF_NO],    N_DofRecvH2, rdisplH2, MPI_INT, MasterBufH2[DOF_NO],    N_DofSendH2, sdisplH2, MPI_INT, Comm);

 for(i=0;i<N_SendDofH2;i++){
   temp_globalno = MasterBufH2[GLOBAL_NO][i];
   temp_dofno    = MasterBufH2[DOF_NO][i];
   temp          = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
   temp          = GlobalNumbers[temp*N_LocDof + temp_dofno];
   DofSendH2[i]  = temp; 
 }
 
 for(i=0;i<2;i++){
   if(N_Halo2>0){
    delete [] SlaveBufH2[i];
    SlaveBufH2[i] = NULL;
   }
   if(N_SendDofH2>0){
     delete [] MasterBufH2[i];
     MasterBufH2[i] = NULL;
   }
 }
 delete [] SlaveBufH2;    SlaveBufH2  = NULL;
 delete [] MasterBufH2;   MasterBufH2 = NULL;
 delete [] temp_arrH2;    temp_arrH2  = NULL;
 
 if(N_SendDof>0)        Send_Info   = new double[N_SendDof];
 if(N_SendDofMS>0)      Send_InfoMS = Send_Info;
 if(N_SendDofH1>0)      Send_InfoH1 = Send_Info + N_SendDofMS;
 if(N_SendDofH2>0)      Send_InfoH2 = Send_Info + N_SendDofMS + N_SendDofH1;
 
 if(N_Slave>0)          Recv_Info   = new double[N_Slave];
 if(N_Slave>0)          Recv_InfoMS = Recv_Info;
 if(N_Halo1>0)          Recv_InfoH1 = Recv_Info + N_InterfaceS;
 if(N_Halo2>0)          Recv_InfoH2 = Recv_Info + N_InterfaceS + N_Halo1;
 
 end_time = MPI_Wtime();
 if(rank == 0)
    printf("Total Time Taken for creating the mapper = %lf",end_time-start_time);
 
 //if(TDatabase::ParamDB->SC_VERBOSE>2)
    if(rank==TDatabase::ParamDB->Par_P0)
      printf("\n################       Mapping for slave-master dofs and halo_1, halo_2 dofs done !!!    ################\n");
  
  if(rank == 0)
      printf("total time taken by the ConstructDofMap_light() = %lf\n",MPI_Wtime()-global_start_time);  
      
 //only for checking purpose   
 if(0){
	  int Total_Own = 0;   
	  for(i=0;i<N_Dof;i++)
	  {
	    if(Master[i]!=rank){
	      if(!(Verify[i]=='s' || Verify[i]=='h' || Verify[i]=='H')){
		printf("\n................1.This shudnt happen...................\n");
		MPI_Finalize();
		exit(0);
	      }
	    }
	    else if(Master[i]==rank){
	      Total_Own++;
	      if(!(Verify[i]=='m' || Verify[i]=='d' || Verify[i]=='D' || Verify[i]=='i' || Verify[i]=='x')){
		printf("\n................2.This shudnt happen...........Verify[%d]=%c........\n",i,Verify[i]);
		MPI_Finalize();
		exit(0);
	      }
	    }
	    else{
		printf("\n................3.This shudnt happen...................\n");
		MPI_Finalize();
		exit(0);
	    } 
	  }
	  
// 	  //DEBUG
// 	  //if(size<5)
// 	  for(aa=0;aa<size;aa++){
// 	      if(rank==aa){
// 		    printf("\nRank::%d\n",rank);
// 		    printf("N_OwnDof     = %d\t Not_Own = %d\n",Total_Own,N_Dof-Total_Own);
// 	      }
// 	      MPI_Barrier(MPI_COMM_WORLD);
// 	    }//verified
 }
    MPI_Barrier(MPI_COMM_WORLD);
    
//  MPI_Finalize();   
//  exit(0);
}

void TParFECommunicator3D::ConstructDofMap()
{
   
 int rank, size, i,jj,ii, j, k, l, m, m1, P, N_Cells, N_U,N_Active, N_LocDof, ID;
 int M, N, N_Vert, N_Joints, N_JointDOF, Neib_ID;
 int *DOF, *GlobalNumbers, *BeginIndex, *JointDof, Disp, N_EdgeDOF, *EdgeDof,*LocalIndex,*Verify;
 int N_CrossEdgeNeibs, *CrossEdgeNeibsRank, N_Edges, N_VertInCell,N_OwnCells;
 int N_VertCrossNeibs, *VertCrossNeibs, VertDof;
 int test, N_Dof;

 double x,y,z;

 bool UPDATE;

 TCollection *Coll;
 TBaseCell *cell, *Neib_cell;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_U = FESpace->GetN_DegreesOfFreedom();
  N_Active = FESpace->GetN_ActiveDegrees();
  N_OwnCells = Coll->GetN_OwnCells();
  
 // printf("N_U is %d----Rank %d\n",N_U,rank);
  
  
  for(i=0;i<N_OwnCells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->IsHaloCell()) printf("This shudnt happen---------------------------------------------\n");
  }
  for(i=N_OwnCells;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(!cell->IsHaloCell()) printf("This shudnt happen---------------------------------------------\n");
  }

  /** *************************************************************************************/
  /** ARRAY CONTAINING GLOBAL DOF numbers (GLOBAL means LOCALLY GLOBAL for the subdomain) */ 
  /** BeginIndex: START ID IN GlobalNumbers ARRAY Cellwise*********************************/
  BeginIndex = FESpace->GetBeginIndex();		
  GlobalNumbers = FESpace->GetGlobalNumbers();          
   
  /** Array containing number of ranks(process ID) associated (surrounding) with each DOF (including halo DOF) */ 
  N_DofRankIndex = new int[N_U];
 
  memset(N_DofRankIndex, 1, N_U*SizeOfInt);
  /** Array containing global number of all Local cells */
  LocalIndex = new int[N_Cells]; 
  
  int **MappingData;
  MappingData = new int*[2];
  for(i=0;i<2;i++)
  {
    MappingData[i] = new int[N_U];
    memset(MappingData[i], 0, N_U*SizeOfInt);
  }
  
  Master = new int[N_U];
  for(i=0;i<N_U;i++) Master[i]=rank;
    
  Verify = new int[N_U];
  memset(Verify, 0, N_U*SizeOfInt);
  
 /** START ---> [ DofRankIndex ------- N_DofRankIndex -------- N_DependentCells ] **/ 
 for(i=0; i<N_Cells; i++)
  {
     cell = Coll->GetCell(i);
     LocalIndex[i] = cell->GetGlobalCellNo();
     
     if(cell->IsDependentCell() && !cell->IsHaloCell())
     {
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	Verify[N] = -1;
      }  // for(j=0; j<N_DOF; j++)
      
     }
  } //  for(i=0; i<N_Cells
    

 for(i=N_OwnCells;i<N_Cells;i++)
 {
      cell = Coll->GetCell(i);
      ID = cell->GetSubDomainNo();
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	//if(N>=N_Active) continue;
	if(Verify[N] == 0) 
	{
	  Verify[N] = 1;
	  Master[N] = ID;
	}
	if(ID < Master[N]) Master[N] = ID;
	MappingData[GLOBAL_NO][N]    = cell->GetGlobalCellNo();
	MappingData[DOF_NO][N]       = j; 
      }
      
  }
 
   
    int temp,temp_globalno,temp_dofno;
  
  
   /** END ---> [ Master ------------ Verify] **/ 
  
    int *N_DOFtobeverified;	  /** Array containing how many DOF's need to be verified from other processors(ranks) */
    int *N_DOFverifiedbythisrank; /** Array containing how many DOF's need to be verified by this rank for other processors(ranks) */
    int VerifiedDOF = 0;	  /** Total no of DOF's that need to be verified from other ranks */
    int Verifiedbythisrank =0;    /** Total no of DOF verified by this rank */
    int **sendbuf,**recvbuf,*verrecvbuf,*versendbuf;
     
    /** PROTECTED VARIABLES  * int *sdispl,*rdispl;*/
    int *temp_arr;
	
    sendbuf 	= new int*[2];
    recvbuf 	= new int*[2];
    sdispl  	= new int[size];
    rdispl  	= new int[size];
    temp_arr	= new int[size];
	
    N_DOFtobeverified = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
    N_DOFverifiedbythisrank = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
	
    for(i=0;i<N_U;i++)
    {
      if(Verify[i] == 1)
      {
	N_DOFtobeverified[Master[i]]++;
	VerifiedDOF++;
      }
    }
		
    MPI_Alltoall(N_DOFtobeverified,1, MPI_INT,N_DOFverifiedbythisrank ,1, MPI_INT, Comm);
					
    for(i=0;i<size;i++)
      Verifiedbythisrank +=  N_DOFverifiedbythisrank[i];
		
    
    for(i=0;i<2;i++)
    {
      sendbuf[i] = new int[VerifiedDOF];
      recvbuf[i] = new int[Verifiedbythisrank];
    }
		
    //printf("Rank F %d----------%d\n",rank,Verifiedbythisrank);
    sdispl[0] = 0;
    for(i=1;i<size;i++)
    {
      sdispl[i] = N_DOFtobeverified[i-1] + sdispl[i-1];
      if(rank == 0);
      //printf(" %d \n",sdispl[i]);
    }
					
    rdispl[0] = 0;
    for(i=1;i<size;i++)
      rdispl[i] = N_DOFverifiedbythisrank[i-1] + rdispl[i-1];
					
					
    memcpy (temp_arr,sdispl, size*SizeOfInt );
					
    for(i=0;i<N_U;i++)
     {
      if(Verify[i] == 1)
      {
	//printf("%d ",temp_arr[Master[i]]);
	sendbuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	sendbuf[DOF_NO][temp_arr[Master[i]]] = MappingData[DOF_NO][i];
	temp_arr[Master[i]]++;
      }
    }

    MPI_Alltoallv(sendbuf[GLOBAL_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[GLOBAL_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
    MPI_Alltoallv(sendbuf[DOF_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[DOF_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
	      
    verrecvbuf  	= new int[Verifiedbythisrank];
    versendbuf  	= new int[VerifiedDOF];

    for(i=0;i<Verifiedbythisrank;i++)
    {
      temp = GetLocalIndex(N_Cells,LocalIndex,recvbuf[GLOBAL_NO][i]);
      verrecvbuf[i] = Master[GlobalNumbers[temp * N_LocDof + recvbuf[DOF_NO][i]]];
    }
	
    MPI_Alltoallv(verrecvbuf, N_DOFverifiedbythisrank, rdispl, MPI_INT,versendbuf,N_DOFtobeverified, sdispl, MPI_INT, Comm);
					
    for(i=0;i<VerifiedDOF;i++)
    {
      temp_globalno = sendbuf[GLOBAL_NO][i];
      temp_dofno    = sendbuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      if(Verify[temp] != 1)
	printf("Error : This degree of Freedom (%d) didn't require verification\n",temp);
      Master[temp] = versendbuf[i];
    }
					
/**===================================== END : MASTER FOR ALL DEGREES OF FREEDOM HAS BEEN DECIDED ========================*/

    delete [] verrecvbuf;
    delete [] versendbuf;
    delete [] N_DOFtobeverified;
    delete [] N_DOFverifiedbythisrank;
    delete [] temp_arr;

		
    for(i=0;i<2;i++)
    {
      delete [] sendbuf[i];
      delete [] recvbuf[i];
    }

    delete [] sendbuf;
    delete [] recvbuf;
		
				
/**==================================START : NOW FIRST GATHERING INFORMATION REGARDING WHICH 
 *                                    DOF's have to be received from OTHER PROCESSORS==================*/
					
    int **SlaveBuf,**MasterBuf,SizeDofSend = 0;
    /** PROTECTED VARIABLES * int *N_DofSend,*N_DofRecv,*DofSend,*DofRecv;*/
    N_Slave = 0;
    N_OwnDof = 0;
					
    N_DofSend = new int[size];
    N_DofRecv = new int[size];
    temp_arr  = new int[size];
    SlaveBuf  = new int*[2];
    MasterBuf = new int*[2];
	
	//for(i=N_Active;i<N_U;i++) Master[i]=rank;
	
    memset(N_DofRecv,0,size*SizeOfInt);
			
    for(i=0;i<N_U;i++)
    {
     if(Master[i] != rank)
     {  
      N_Slave++;
      N_DofRecv[Master[i]]++;
     } 
     else 
       N_OwnDof++;
    }
  //  printf("Rank %d      OwnDofs is %d---------\n",rank,N_OwnDof);

    OwnDofs = new int[N_OwnDof];
    rdispl[0] = 0; 

    MPI_Alltoall(N_DofRecv,1, MPI_INT,N_DofSend ,1, MPI_INT, Comm);

    for(i=1;i<size;i++)
      rdispl[i] = rdispl[i-1] + N_DofRecv[i-1];
	
    sdispl[0] = 0;
    for(i=0;i<size;i++)
    {
      SizeDofSend += N_DofSend[i];
      sdispl[i] = sdispl[i-1] + N_DofSend[i-1];
    }


    DofSend  = new int[SizeDofSend];
    DofRecv  = new int[N_Slave];
    N_SendDof = SizeDofSend;
	
    for(i=0;i<2;i++)
    {
      SlaveBuf[i]  = new int[N_Slave];
      MasterBuf[i] = new int[SizeDofSend];
    }
	
    memcpy (temp_arr,rdispl, size*SizeOfInt );
	
    m=0;
    for(i=0;i<N_U;i++)
    {
      if(Master[i] != rank)
      {
	SlaveBuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	SlaveBuf[DOF_NO]   [temp_arr[Master[i]]] = MappingData[DOF_NO]   [i];
	DofRecv[temp_arr[Master[i]]] = i;
	temp_arr[Master[i]]++;
      }
      else 
      {
	OwnDofs[m]=i;
	m++;
      }
    }
			
    MPI_Alltoallv(SlaveBuf[GLOBAL_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[GLOBAL_NO], N_DofSend, sdispl, MPI_INT, Comm);
    MPI_Alltoallv(SlaveBuf[DOF_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[DOF_NO], N_DofSend, sdispl, MPI_INT, Comm);
		
    for(i=0;i<SizeDofSend;i++)
    {
      temp_globalno = MasterBuf[GLOBAL_NO][i];
      temp_dofno    = MasterBuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      DofSend[i] = temp;  
    }
  
  //  if(TDatabase::ParamDB->SC_VERBOSE>1)   
	//printf(" Rank %d ------ NUMBER OF DOF's to be sent = %d -------- NUMBER OF DOF's to be recv = %d\n",rank,N_SendDof,N_Slave); 
    for(i=0;i<2;i++)
    {
      delete [] SlaveBuf[i];
      delete [] MasterBuf[i];
    }

    delete [] SlaveBuf;
    delete [] MasterBuf;
	
    DofSendPos = new int[SizeDofSend];
    DofRecvPos = new int[N_Slave];
    for(i=0;i<SizeDofSend;i++) DofSendPos[i]=i;
    for(i=0;i<N_Slave;i++) DofRecvPos[i]=i;
     //SortPos(DofSend,DofSendPos,SizeDofSend);
    //quickSort(DofSend,DofSendPos,0,SizeDofSend);
  if(N_SendDof>0)
   Send_Info = new double[N_SendDof];
  if(N_Slave>0)
   Recv_Info = new double[N_Slave];
    
    if(TDatabase::ParamDB->SC_VERBOSE>2)
      if(rank==TDatabase::ParamDB->Par_P0)
	printf("ConstructDofMap done !!!\n");
      
       //DEBUG
 int aa;
 for(aa=0;aa<size;aa++){
     if(rank==aa){
          printf("\nRank::%d\n",rank);
	  printf("N_OwnDof     = %d\t Not_Own = %d\n",N_OwnDof,N_Slave);
    }
     MPI_Barrier(MPI_COMM_WORLD);
  }//verified
}

void TParFECommunicator3D::SetFENeibCommunicationSteps()
{
  int rank,size,i, j, k, l, m, N, P, Q;
  int *dofNeibIDs,N_DofNeibs_MaxAll;
  int *pos, *N_DofNeibs_Array, *DofNeibIDs_Array, N_Entries;
  bool UPDATE;
  
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  N_DofNeibs=0;
  dofNeibIDs= new int[size];
  m=0;
  for(i=0;i<size;i++)
  {
    if(N_DofRecv[i]>0)
    {
      N_DofNeibs++;
      dofNeibIDs[m]=i;
      m++;
    }
    else if(N_DofSend[i]>0)
    {
      N_DofNeibs++;
      dofNeibIDs[m]=i;
      m++;
    }
  }
  DofNeibIDs = new int[N_DofNeibs];
  memcpy(DofNeibIDs, dofNeibIDs, N_DofNeibs*SizeOfInt);
  delete [] dofNeibIDs;
  
  IndexOfNeibRank = new int[size];
  for(i=0; i<size; i++)
   IndexOfNeibRank[i] = -1;
  for(i=0; i<N_DofNeibs; i++)
   IndexOfNeibRank[DofNeibIDs[i]] = i;
  
  MPI_Allreduce(&N_DofNeibs, &N_DofNeibs_MaxAll, 1, MPI_INT, MPI_MAX, Comm); 
  if(rank==0)
  {
   pos = new int[size];
   for(i=0;i<size;i++)
     pos[i] = i*N_DofNeibs_MaxAll;
   N_DofNeibs_Array = new int[size];
   DofNeibIDs_Array= new int[size*N_DofNeibs_MaxAll];
  }

  MPI_Gather(&N_DofNeibs, 1, MPI_INT, N_DofNeibs_Array, 1, MPI_INT, 0, Comm);
  MPI_Gatherv(DofNeibIDs, N_DofNeibs, MPI_INT, DofNeibIDs_Array, N_DofNeibs_Array, pos,MPI_INT, 0, Comm);
 
  if(rank==0)
  {
   int step, *Entries, *ColInd, *RowPtr, *Collection, MaxPossibleSteps;
   N_Entries = 0;

   for(i=0;i<size;i++)
      N_Entries +=N_DofNeibs_Array[i];
   printf("No.of entries are %d\n",N_Entries);
   Entries = new int[N_Entries];         // adjacency matrix memory allocation
   ColInd     = new int[N_Entries];
   RowPtr = new int[size+1];

   RowPtr[0] = 0; m = 0;
   for(i=0;i<size;i++)
   {
     N =  N_DofNeibs_Array[i];
     RowPtr[i+1] = RowPtr[i] + N;
     for(j=0;j<N;j++)
     {
       ColInd[m] = DofNeibIDs_Array[i*N_DofNeibs_MaxAll + j];  // assmeble the matrix
       Entries[m] =-1;
       m++;
     }
   }
   
    // print adjacency matrix
   /*  for(i=0;i<size;i++)
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
        printf("a( %d , %d )=  %d \n",i, ColInd[j], Entries[j]);
    MPI_Finalize();exit(0);*/
    
   // find number of steps to communicate between neibs (edge coloring alogrithm)
   Collection = new int[size];
   Max_CommunicationSteps = -1;
   MaxPossibleSteps = (N_DofNeibs_MaxAll+1)*10;
   for(step=0;step<MaxPossibleSteps;step++)
   {
    m=0;
    for(k=0;k<size;k++)
    {
     UPDATE = TRUE;
     for(l=0;l<m;l++)
     if(Collection[l] == k)
     {
      UPDATE = FALSE;
       break;
     }
     if(UPDATE)
     {
       Collection[m] = k;
        m++;
     }
     else
     {
      continue;
     }
 
     for(l=RowPtr[k];l<RowPtr[k+1];l++)
     { 
        // only upper tirangular martix is enough (due to symmetry)
         if( ColInd[l]<k) continue;

         UPDATE = TRUE;
         for(P=0;P<m;P++)
         if(Collection[P] == ColInd[l])
         {
           // ColInd[l] is already in the collection
          UPDATE = FALSE;
          break;
         } 
       
         if(UPDATE)
         {
          if(Entries[l]<0)
          {
            Entries[l] = step;
            if(Max_CommunicationSteps<step+1) Max_CommunicationSteps=step+1;
            Collection[m] = ColInd[l];
            m++;
            break;
           }
         }
        } //  for(l=RowPtr[k];
      } // for(k=0;
    } //  for(step=0;step<

    printf("Total number of steps needed to communicate is %d\n",Max_CommunicationSteps);
   //  print adjacency matrix
//     for(i=0;i<size;i++)
//      for(j=RowPtr[i];j<RowPtr[i+1];j++)
//        printf("a( %d , %d )=  %d \n",i, ColInd[j], Entries[j]);
    

    N_CommunicationProcesses = new int[Max_CommunicationSteps];
    SendID_Index= new int[Max_CommunicationSteps*size];
    ReceiveID_Index = new int[Max_CommunicationSteps*size];

    for(j=0;j<Max_CommunicationSteps;j++)
     N_CommunicationProcesses[j] =0;

    for(i=0;i<Max_CommunicationSteps*size;i++)
     {
      SendID_Index[i] = -1;
      ReceiveID_Index[i] = -1;
     }

    for(i=0;i<size;i++)
     for(j=RowPtr[i];j<RowPtr[i+1];j++)
       if(i<ColInd[j])
        {
         if( Entries[j]<0)
          { 
           printf("Error in finding DOF communication steps %d\n",Entries[j]);
           MPI_Abort(MPI_COMM_WORLD,  0);
          }
         k =  Entries[j];
	// if(k==1) printf("i is %d j is %d\n",i,ColInd[j]);
         SendID_Index[k*size + N_CommunicationProcesses[k]] = i;
         ReceiveID_Index[k*size + N_CommunicationProcesses[k]] =  ColInd[j];
         N_CommunicationProcesses[k]++;
        }
   /*
     for(i=0;i<Max_CommunicationSteps;i++)
     {
       N= N_CommunicationProcesses[i];
       printf("Step %d: No.of processes is %d\n",i,N);
       for(j=0;j<N;j++)
        printf("Send from process %d to %d -------- and onto\n",
                   SendID_Index[i*size +j], ReceiveID_Index[i*size +j]);
       printf("\n");
     }   */

    delete [] Collection;
  }
   
   MPI_Bcast(&Max_CommunicationSteps, 1, MPI_INT, 0, Comm);
   MPI_Barrier(Comm);

   if(rank!=0)
    {
     N_CommunicationProcesses = new int[Max_CommunicationSteps];
     SendID_Index= new int[Max_CommunicationSteps*size];
     ReceiveID_Index = new int[Max_CommunicationSteps*size];
    }

   MPI_Bcast(N_CommunicationProcesses, Max_CommunicationSteps, MPI_INT, 0, Comm);
   MPI_Bcast(SendID_Index, Max_CommunicationSteps*size, MPI_INT, 0, Comm);
   MPI_Bcast(ReceiveID_Index, Max_CommunicationSteps*size, MPI_INT,  0, Comm);

 //  printf("Parallel SetIntraCommunicator ranks  %d, N_DofNeibs %d\n", rank, N_DofNeibs);
       //     MPI_Finalize();
       //     exit(0);
}

// #ifdef _HYBRID

void TParFECommunicator3D::Color(int &numColors, int *&ptrColors, char type)
{
  int i,j,k;
  int *RowPtr, *KCol;
  
  int *myReorder;
  int ndof;   

  //find the type of dof that needs to be colored
  if(type == 'm'){
    myReorder = Reorder_M;	ndof = N_InterfaceM; }
  else if(type == 'i'){
    myReorder = Reorder_I;	ndof = N_Int;        }
  else if(type == 'D'){
    myReorder = Reorder_D1;	ndof = N_Dept1;      }
  else if(type == 'd'){
    myReorder = Reorder_D2;	ndof = N_Dept2;      }
  else if(type == 'x'){
    myReorder = Reorder_D3;	ndof = N_Dept3;      }
  else{
    printf("wrong type!!!\n"); exit(0);              }
  
  //this stores the total number of colrs used for coloring
  numColors = 0;
  if(!ndof)	return;
  
  int max, temp, t;
  //this stores the color number for dofs
  int *allocatedColor = new int[ndof];
  //initialize all colors with default
  for(i=0;i<ndof;i++)
    allocatedColor[i] = -1;
  
  //get the system matrix  
  RowPtr = sqstruct->GetRowPtr();
  KCol   = sqstruct->GetKCol();
  
  //now we start coloring the dofs
  for(i=0;i<ndof;i++)
  {
    temp = -1;
    k = myReorder[i];
    for(j=RowPtr[k];j<RowPtr[k+1];j++)
    {
      if(KCol[j] >= k || DofMarker[KCol[j]] != type)	continue;
      
      t = NewGN[KCol[j]];		    //locate the pos of the dof in myreorder to identify its color in allocatedColor
      if(temp < allocatedColor[t])
	temp = allocatedColor[t];
    }
    
    allocatedColor[i] = temp+1;
    
    if(numColors < allocatedColor[i])
      numColors = allocatedColor[i];
  }
  //colors were numbered from 0, hece total will be 1 more
  numColors++;
  
  ptrColors = new int[numColors+1];
  temp = 0; k = 0;
  //arrange the dofs, such that same color dofs are together
  for(i=0;i<numColors;i++)
  {
    ptrColors[i] = k;
    for(j=0;j<ndof;j++)
    {
      if(allocatedColor[j] == i)
      {
	temp = myReorder[j];
	myReorder[j] = myReorder[k];
	myReorder[k] = temp;
	k++;
	
	allocatedColor[j] = -1;
      }
    }
  }
  
  ptrColors[numColors] = k;
  printf("numcolors (type = %c):: %d\t total dof = %d\n",type,numColors,ndof);
//   exit(0);
}

// #endif


void TParFECommunicator3D::CommUpdate(double *sol)
{
  int i, j, k, l, N,iter, rank, size;
  int sendID, recvID, *val_send, *val_recv;
  double t1,t2=0.0,t3,tSum=0.;
 
  MPI_Status status;
  MPI_Comm_rank(Comm,&rank);	
  MPI_Comm_size(Comm, &size);
  
  t1=MPI_Wtime();
  
  //ConstructDofMap_light
  if(TDatabase::ParamDB->SC_SMOOTHER_SCALAR==6){
#ifdef _HYBRID
    int t;
    printf("not yet implemented\n");
    MPI_Finalize();
    exit(0);
#else
  CommUpdateMS(sol);
  CommUpdateH1(sol);
  CommUpdateH2(sol);
#endif
  }
  else{
    for(i=0;i<N_SendDof;i++)
      Send_Info[i]=sol[DofSend[i]];
    
    for(i=0; i<Max_CommunicationSteps; i++)
    {
      N = N_CommunicationProcesses[i];
      for(j=0;j<N;j++)
      {
	sendID = SendID_Index[i*size +j];
	recvID = ReceiveID_Index[i*size +j];

	if((rank!=sendID) && (rank!=recvID))  continue;
	
	if(rank==sendID && N_DofSend[recvID]>0)
	  MPI_Send(&Send_Info[sdispl[recvID]], N_DofSend[recvID], MPI_DOUBLE,  recvID, 100, Comm);
	
	else if(rank==recvID && N_DofRecv[sendID]>0)
	  MPI_Recv(&Recv_Info[rdispl[sendID]], N_DofRecv[sendID], MPI_DOUBLE,  sendID, 100, Comm, &status);

	if(rank==recvID && N_DofSend[sendID])
	  MPI_Send(&Send_Info[sdispl[sendID]], N_DofSend[sendID], MPI_DOUBLE,  sendID, 200, Comm);
	
	else if(rank==sendID && N_DofRecv[recvID])
	  MPI_Recv(&Recv_Info[rdispl[recvID]], N_DofRecv[recvID], MPI_DOUBLE,  recvID, 200, Comm, &status);
      } //for(j=0;j<N
    } //for(i=0; i<Max_CommunicationSteps

    for(i=0;i<N_Slave;i++)
      sol[DofRecv[i]] = Recv_Info[i];
  
    MPI_Barrier(MPI_COMM_WORLD);
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

void TParFECommunicator3D::CommUpdate_M_H1(double *sol)
{
  double t1,t2;
  t1=MPI_Wtime();
  
  CommUpdateMS(sol);
  CommUpdateH1(sol);
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

void TParFECommunicator3D::CommUpdateH1(double *sol)
{
  int i;
  double t1,t2;
  t1=MPI_Wtime();
  
  for(i=0;i<N_SendDofH1;i++)
      Send_InfoH1[i]=sol[DofSendH1[i]];

  MPI_Alltoallv(Send_InfoH1,N_DofSendH1,sdisplH1,MPI_DOUBLE,Recv_InfoH1,N_DofRecvH1,rdisplH1,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_Halo1;i++)  
    sol[DofRecvH1[i]] = Recv_InfoH1[i];
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

void TParFECommunicator3D::CommUpdateH2(double *sol)
{
  int i;
  double t1,t2;
  t1=MPI_Wtime();
  
  for(i=0;i<N_SendDofH2;i++)
      Send_InfoH2[i]=sol[DofSendH2[i]];
  
  MPI_Alltoallv(Send_InfoH2,N_DofSendH2,sdisplH2,MPI_DOUBLE,Recv_InfoH2,N_DofRecvH2,rdisplH2,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_Halo2;i++)  
    sol[DofRecvH2[i]] = Recv_InfoH2[i];
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

void TParFECommunicator3D::CommUpdateMS(double *sol)
{
  int i;
  double t1,t2;
  t1=MPI_Wtime();
  
  for(i=0;i<N_SendDofMS;i++)
      Send_InfoMS[i]=sol[DofSendMS[i]];
  
  MPI_Alltoallv(Send_InfoMS,N_DofSendMS,sdisplMS,MPI_DOUBLE,Recv_InfoMS,N_DofRecvMS,rdisplMS,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_InterfaceS;i++)  
    sol[DofRecvMS[i]] = Recv_InfoMS[i];
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

void TParFECommunicator3D::CommUpdateAlltoAllv(double *sol)
{
  int i,rank;
  double t1,t2;
  t1=MPI_Wtime();
   
  MPI_Comm_rank(Comm,&rank);		

  Send_Info = new double[N_SendDof];
  Recv_Info = new double[N_Slave];

  for(i=0;i<N_SendDof;i++)
  {
      Send_Info[i]=sol[DofSend[i]];
  }
  MPI_Alltoallv(Send_Info,N_DofSend,sdispl,MPI_DOUBLE,Recv_Info,N_DofRecv,rdispl,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_Slave;i++)
  {  
    sol[DofRecv[i]] = Recv_Info[i];
  }
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

void TParFECommunicator3D::CommUpdateReduce(double *rhs)
{
  int i,j,rank,size,sendID, recvID,N;
  MPI_Status status;
  MPI_Comm_rank(Comm,&rank);		
  MPI_Comm_size(Comm, &size);
  double t1,t2;
  t1=MPI_Wtime();
  
  for(i=0;i<N_Slave;i++)
  {
      Recv_Info[i]=rhs[DofRecv[i]];
  }
 
//  MPI_Alltoallv(Recv_Info,N_DofRecv,rdispl,MPI_DOUBLE,Send_Info,N_DofSend,sdispl,MPI_DOUBLE,Comm);
  
    
  for(i=0; i<Max_CommunicationSteps; i++)
  {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
    {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];
     // if(rank==0) printf("rank %d step i is %d send ID %d recvID %d\n",rank,i,sendID,recvID);

      if((rank!=sendID) && (rank!=recvID))  continue;
     
      if(rank==sendID && N_DofRecv[recvID]>0)
      {
	//printf("rank is %d\n",rank);
        MPI_Send(&Recv_Info[rdispl[recvID]], N_DofRecv[recvID], MPI_DOUBLE,  recvID, 100, Comm);
      }
      
      else if(rank==recvID && N_DofSend[sendID]>0)
        MPI_Recv(&Send_Info[sdispl[sendID]], N_DofSend[sendID], MPI_DOUBLE,  sendID, 100, Comm, &status);


      if(rank==recvID && N_DofRecv[sendID])
         MPI_Send(&Recv_Info[rdispl[sendID]], N_DofRecv[sendID], MPI_DOUBLE,  sendID, 200, Comm);
       
      else if(rank==sendID && N_DofSend[recvID])
         MPI_Recv(&Send_Info[sdispl[recvID]], N_DofSend[recvID], MPI_DOUBLE,  recvID, 200, Comm, &status);
      
     } //for(j=0;j<N
   } //for(i=0; i<Max_CommunicationSteps*/
 
  for(i=0;i<N_SendDof;i++)
  {  
      rhs[DofSend[i]] += Send_Info[i];
  }

/** MASTER HAS FINAL VALUE **/  
  for(i=0;i<N_SendDof;i++)
  {
      Send_Info[i]=rhs[DofSend[i]];
  }
  
 // MPI_Alltoallv(Send_Info,N_DofSend,sdispl,MPI_DOUBLE,Recv_Info,N_DofRecv,rdispl,MPI_DOUBLE,Comm);
  
    
  for(i=0; i<Max_CommunicationSteps; i++)
  {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
    {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];
     // if(rank==0) printf("rank %d step i is %d send ID %d recvID %d\n",rank,i,sendID,recvID);

      if((rank!=sendID) && (rank!=recvID))  continue;
     
      if(rank==sendID && N_DofSend[recvID]>0)
      {
	//printf("rank is %d\n",rank);
        MPI_Send(&Send_Info[sdispl[recvID]], N_DofSend[recvID], MPI_DOUBLE,  recvID, 100, Comm);
      }
      
      else if(rank==recvID && N_DofRecv[sendID]>0)
        MPI_Recv(&Recv_Info[rdispl[sendID]], N_DofRecv[sendID], MPI_DOUBLE,  sendID, 100, Comm, &status);


      if(rank==recvID && N_DofSend[sendID])
         MPI_Send(&Send_Info[sdispl[sendID]], N_DofSend[sendID], MPI_DOUBLE,  sendID, 200, Comm);
       
      else if(rank==sendID && N_DofRecv[recvID])
         MPI_Recv(&Recv_Info[rdispl[recvID]], N_DofRecv[recvID], MPI_DOUBLE,  recvID, 200, Comm, &status);
      
     } //for(j=0;j<N
   } //for(i=0; i<Max_CommunicationSteps*/
  
  for(i=0;i<N_Slave;i++)
  {  
      rhs[DofRecv[i]] = Recv_Info[i];
  }
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

void TParFECommunicator3D::CommUpdateReduceMS(double *rhs)
{
  double t1,t2;
  t1=MPI_Wtime();
  int i,j,rank,size;
  MPI_Status status;
  MPI_Comm_rank(Comm,&rank);		
  MPI_Comm_size(Comm, &size);
  
  for(i=0;i<N_InterfaceS;i++)
      Recv_InfoMS[i]=rhs[DofRecvMS[i]];
 
  MPI_Alltoallv(Recv_InfoMS,N_DofRecvMS,rdisplMS,MPI_DOUBLE,Send_InfoMS,N_DofSendMS,sdisplMS,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_SendDofMS;i++) 
      rhs[DofSendMS[i]] += Send_InfoMS[i];

  /** MASTER HAS FINAL VALUE **/  
  for(i=0;i<N_SendDofMS;i++)
      Send_InfoMS[i]=rhs[DofSendMS[i]];
  
  MPI_Alltoallv(Send_InfoMS,N_DofSendMS,sdisplMS,MPI_DOUBLE,Recv_InfoMS,N_DofRecvMS,rdisplMS,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_InterfaceS;i++)
      rhs[DofRecvMS[i]] = Recv_InfoMS[i];
  
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
}

/** Communication between processor based on defined schedule */
/** arrays (send & recev) with same (senddisp) size */
int TParFECommunicator3D::FECommunicateNeib(int **sendbuf, int senddisp, int *sendlen,
                                              int **recvbuf, int *recevlen, int N_array)
{
 int i, j, k, l, m, M, N, rank, size, sendID, recvID;
 int totaldisp, *buffer, *tmp, pos, *Sendbuf, *Recvbuf, max;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  max = N_array*senddisp;
  MPI_Allreduce(&max, &totaldisp, 1, MPI_INT, MPI_MAX, Comm);

//  if(TDatabase::ParamDB->Par_P6==5)
//    printf(" rank %d FECommunicateNeib    totaldisp %d\n", rank, totaldisp);

//   if(max==0)  
//    totaldisp=0;

  if(totaldisp!=0)
  {
   buffer = new int[totaldisp];
   tmp = new int[totaldisp];

   for(i=0;i<Max_CommunicationSteps;i++)
   {
    N = N_ActiveCommuPross[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank
        pos = 0;

        if(sendlen[M]>0)
         {
         // pack all arrays and send 

         for(l=0;l<N_array; l++)
          {
           Sendbuf = sendbuf[l];
           MPI_Pack(Sendbuf+M*senddisp, sendlen[M], MPI_INT, buffer, totaldisp*SizeOfInt,  &pos, Comm);
          }

          MPI_Send(buffer, pos, MPI_PACKED,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; 

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfInt, MPI_PACKED,  sendID, 100, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfInt, &pos, tmp, recevlen[M], MPI_INT, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*senddisp + m] = tmp[m];
           } //  for(l=0;l<N_ar
         } // if(recevlen[M]>0)
       } //  else if(rank==recvID

      if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank
        pos = 0;

        if(sendlen[M]>0)
         {
          // pack all arrays and send 
          for(l=0;l<N_array; l++)
           {
            Sendbuf = sendbuf[l];
            MPI_Pack(Sendbuf+M*senddisp, sendlen[M], MPI_INT, buffer, totaldisp*SizeOfInt,  &pos, Comm);
           }

          MPI_Send(buffer, pos, MPI_PACKED,  sendID, 200, Comm);
         } // if(sendlen[M]>0)
       }
      else if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfInt, MPI_PACKED,  recvID, 200, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfInt, &pos, tmp, recevlen[M], MPI_INT, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*senddisp + m] = tmp[m];
           } //  for(l=0;l<N_ar
          } // if(recevlen[M]>0
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

   delete [] buffer;
   delete [] tmp;
  }

//   printf(" rank %d FECommunicateNeib senddisp %d totaldisp %d\n", rank, senddisp, totaldisp);
//   MPI_Finalize();
//   exit(0);  

  return 0;
} // FECommunicateNeib

/** Communication between processor based on defined schedule */
/** send and recev arrays with same size */
int TParFECommunicator3D::FECommunicateNeib(double **sendbuf, int senddisp, int *sendlen,
                                            double **recvbuf, int *recevlen, int N_array)
{
 int i, j, k, l, m, M, N, rank, size, sendID, recvID;
 int totaldisp,  pos,  max;

 double *buffer, *tmp, *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  max = N_array*senddisp;
  MPI_Allreduce(&max, &totaldisp, 1, MPI_INT, MPI_MAX, Comm);

  if(totaldisp!=0)
  {
   buffer = new double[totaldisp];
   tmp = new double[totaldisp];

   for(i=0;i<Max_CommunicationSteps;i++)
   {
    N = N_ActiveCommuPross[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank
        pos = 0;

        if(sendlen[M]>0)
         {
         // pack all arrays and send 

         for(l=0;l<N_array; l++)
          {
           Sendbuf = sendbuf[l];
           MPI_Pack(Sendbuf+M*senddisp, sendlen[M], MPI_DOUBLE, buffer, totaldisp*SizeOfDouble,  &pos, Comm);
          }

          MPI_Send(buffer, pos, MPI_PACKED,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; 

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfDouble, MPI_PACKED,  sendID, 100, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfDouble, &pos, tmp, recevlen[M], MPI_DOUBLE, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*senddisp + m] = tmp[m];
           } //  for(l=0;l<N_ar
         } // if(recevlen[M]>0)
       } //  else if(rank==recvID

      if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank
        pos = 0;

        if(sendlen[M]>0)
         {
          // pack all arrays and send 
          for(l=0;l<N_array; l++)
           {
            Sendbuf = sendbuf[l];
            MPI_Pack(Sendbuf+M*senddisp, sendlen[M], MPI_DOUBLE, buffer, totaldisp*SizeOfDouble,  &pos, Comm);
           }

          MPI_Send(buffer, pos, MPI_PACKED,  sendID, 200, Comm);
         } // if(sendlen[M]>0)
       }
      else if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfDouble, MPI_PACKED,  recvID, 200, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfDouble, &pos, tmp, recevlen[M], MPI_DOUBLE, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*senddisp + m] = tmp[m];
           } //  for(l=0;l<N_ar
          } // if(recevlen[M]>0
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

   delete [] buffer;
   delete [] tmp;
  }

//   printf(" rank %d FECommunicateNeib senddisp %d totaldisp %d\n", rank, senddisp, totaldisp);
//   MPI_Finalize();
//   exit(0);  

  return 0;
} // FECommunicateNeib

int TParFECommunicator3D::FECommunicateNeib(double *sendbuf, int senddisp, int *sendlen,
                                            double *recvbuf, int N_dim)
{
 int i, j, M, N, rank, size, sendID, recvID;
 double *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  for(i=0;i<Max_CommunicationSteps;i++)
   {
    N = N_ActiveCommuPross[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Sendbuf = sendbuf+N_dim*M*senddisp;
//           MPI_Send(Sendbuf, N_dim*sendlen[M], MPI_DOUBLE,  recvID, 100, Comm);
          MPI_Send(Sendbuf, N_dim*senddisp, MPI_DOUBLE,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Recvbuf = recvbuf+N_dim*M*senddisp;
//           MPI_Recv(Recvbuf, N_dim*sendlen[M], MPI_DOUBLE,  sendID, 100, Comm, &status);
          MPI_Recv(Recvbuf, N_dim*senddisp, MPI_DOUBLE,  sendID, 100, Comm, &status);
         } // if(sendlen[M]>0)
       } //  else if(rank==recvID

      if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Sendbuf = sendbuf+N_dim*M*senddisp;
//           MPI_Send(Sendbuf, N_dim*sendlen[M], MPI_DOUBLE,  sendID, 100, Comm);
          MPI_Send(Sendbuf, N_dim*senddisp, MPI_DOUBLE,  sendID, 100, Comm);
         } // if(sendlen[M]>0)
       }
      else if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Recvbuf = recvbuf+N_dim*M*senddisp;
//           MPI_Recv(Recvbuf, N_dim*sendlen[M], MPI_DOUBLE,  recvID, 100, Comm, &status);
          MPI_Recv(Recvbuf, N_dim*senddisp, MPI_DOUBLE,  recvID, 100, Comm, &status);
          } // if(sendlen[M]>0
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

  return 0;
} // FECommunicateNeib

int TParFECommunicator3D::FECommunicateNeib(double *sendbuf, int senddisp, int *sendlen,
                                            double *recvbuf, int *recevlen, int N_dim)
{
 int i, j, M, N, rank, size, sendID, recvID;
 double *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  for(i=0;i<Max_CommunicationSteps;i++)
   {
    N = N_ActiveCommuPross[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Sendbuf = sendbuf+N_dim*M*senddisp;
          MPI_Send(Sendbuf, N_dim*sendlen[M], MPI_DOUBLE,  recvID, 100, Comm);
//           MPI_Send(Sendbuf, N_dim*senddisp, MPI_DOUBLE,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          Recvbuf = recvbuf+N_dim*M*senddisp;
          MPI_Recv(Recvbuf, N_dim*recevlen[M], MPI_DOUBLE,  sendID, 100, Comm, &status);
//           MPI_Recv(Recvbuf, N_dim*senddisp, MPI_DOUBLE,  sendID, 100, Comm, &status);
         } // if(recevlen[M]>0)
       } //  else if(rank==recvID

      if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Sendbuf = sendbuf+N_dim*M*senddisp;
          MPI_Send(Sendbuf, N_dim*sendlen[M], MPI_DOUBLE,  sendID, 100, Comm);
//           MPI_Send(Sendbuf, N_dim*senddisp, MPI_DOUBLE,  sendID, 100, Comm);
         } // if(sendlen[M]>0)
       }
      else if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          Recvbuf = recvbuf+N_dim*M*senddisp;
          MPI_Recv(Recvbuf, N_dim*recevlen[M], MPI_DOUBLE,  recvID, 100, Comm, &status);
//           MPI_Recv(Recvbuf, N_dim*senddisp, MPI_DOUBLE,  recvID, 100, Comm, &status);
          } // if(recevlen[M]>0
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

  return 0;
} // FECommunicateNeib

int TParFECommunicator3D::FECommunicateOneWay(int **sendbuf, int disp, int *sendlen,
                                            int **recvbuf, int *recevlen, int N_array, int type)
{
 int i, j, k, l, m, M, N, rank, size, sendID, recvID;
 int totaldisp, *buffer, *tmp, pos, *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  totaldisp = N_array*disp;

// printf(" rank %d FECommunicateNeib    totaldisp %d\n", rank, totaldisp);

 if(totaldisp!=0)
  {
   buffer = new int[totaldisp];
   tmp = new int[totaldisp];

   for(i=0;i<Max_CommunicationSteps;i++)
    {
     N = N_ActiveCommuPross[i];
     for(j=0;j<N;j++)
      {
       if(type==0)
        {
         sendID = SendID_Index[i*size +j];
         recvID = ReceiveID_Index[i*size +j];
        }
       else
        {
         recvID = SendID_Index[i*size +j];
         sendID = ReceiveID_Index[i*size +j];
        }

       if(rank==sendID)
        {
         M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank
         pos = 0;

        if(sendlen[M]>0)
         {
         // pack all arrays and send 
         for(l=0;l<N_array; l++)
          {
           Sendbuf = sendbuf[l];
           MPI_Pack(Sendbuf+M*disp, sendlen[M], MPI_INT, buffer, totaldisp*SizeOfInt,  &pos, Comm);
          }

          MPI_Send(buffer, pos, MPI_PACKED,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfInt, MPI_PACKED,  sendID, 100, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfInt, &pos, tmp, recevlen[M], MPI_INT, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*disp + m] = tmp[m];
           } //  for(l=0;l<N_ar
         } // if(sendlen[M]>0)
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

   delete [] buffer;
   delete [] tmp;
  } // if(tot

// //   printf(" rank %d FECommunicateNeib senddisp %d totaldisp %d\n", rank, senddisp, totaldisp);
// //   MPI_Finalize();
// //   exit(0);  
// 
//   return 0;
} // FECommunicateNeib

/** send array of disp size between neibs */
int TParFECommunicator3D::FECommunicateOneWay(double *sendbuf, int disp, double *recvbuf,  int type)
{
 int i, j, k, l, m, M, N, rank, size, sendID, recvID;
 double *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

 if(disp!=0)
  {
   for(i=0;i<Max_CommunicationSteps;i++)
    {
     N = N_ActiveCommuPross[i];
     for(j=0;j<N;j++)
      {
       if(type==0)
        {
         sendID = SendID_Index[i*size +j];
         recvID = ReceiveID_Index[i*size +j];
        }
       else
        {
         recvID = SendID_Index[i*size +j];
         sendID = ReceiveID_Index[i*size +j];
        }
       if(rank==sendID)
        {
         M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank         
         Sendbuf = sendbuf+M*disp;
//printf("rank %d i %d j %d  recvID %d\n", rank, i, j,  recvID);
         MPI_Send(Sendbuf, disp, MPI_DOUBLE,  recvID , 100, Comm);
       }
      else if(rank==recvID)
       { 
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank
        Recvbuf = recvbuf+M*disp;
//printf("rank %d  i %d j %d   sendID %d \n", rank, i, j, sendID);
        MPI_Recv(Recvbuf, disp, MPI_DOUBLE,  sendID, 100, Comm, &status);
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max
  } // if(tot

// //   printf(" rank %d FECommunicateNeib senddisp %d totaldisp %d\n", rank, senddisp, totaldisp);
// //   MPI_Finalize();
// //   exit(0);  
// 
//   return 0;
} // FECommunicateNeib

#endif


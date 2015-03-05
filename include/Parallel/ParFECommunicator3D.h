// =======================================================================
// @(#)ParFECommunicator3D.h
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

#ifndef __PARFECOMMUNICATOR3D__
#define __PARFECOMMUNICATOR3D__

#include "mpi.h"

#include <FESpace3D.h>
#include <SquareStructure.h>
#include <SquareStructure3D.h>

class TParFECommunicator3D
{
  protected:
    
     int *N_DofSend,*N_DofRecv,*DofSend,*DofRecv,*sdispl,*rdispl,*DofSendPos,*DofRecvPos;  
     
     int *N_DofSendH,*N_DofRecvH,*DofSendH,*DofRecvH,*sdisplH,*rdisplH;
     
     int N_SendDof, N_SendDofH, *Master, *Reorder, *NewGN; 
     
     int N_Slave, N_Halo, N_Master, N_Dept, N_Int, N_ActiveSlave, N_ActiveHalo;
     
     int N_CMaster, N_CDept, N_CInt, *ptrCMaster, *ptrCDept, *ptrCInt;

     double *Send_Info,*Recv_Info; 
     
     double *Send_InfoM,*Recv_InfoM, *Send_InfoH,*Recv_InfoH;

    /** MPI_Comm for which the fespace communications are needed */
     MPI_Comm Comm;

    /** fespace for which the communications are needed */
    TFESpace3D *FESpace;
    
    TSquareStructure3D *sqstruct;

    /**number of dependent Cells*/
    int N_DependentCells;

    /**dependent Cells indices*/
    int *DependentCellIndex;

    /**dependent Cells indices*/
    int *DeptCellGlobalNo;

    /**maximum number(possible) of subdomains among all dofs */
    int MaxSubDomainPerDof;

    /**array containing list of all rank(process ID) associated with each dof (including Halo dof) */
    int *N_DofRankIndex;

    /**array rank(process) ID of each dof (including Halo dof) */
    int *DofRankIndex;

    /**number of own (self + master) dofs*/
    int N_OwnDof;

    /** dof index of own dofs*/
    int *OwnDofs;

    /**number of dependent dofs*/
    int N_DeptDofs;

    /** dof index of dependent dofs*/
    int *DeptDofs;

    /** list of owner of the dependent dofs*/
    int *DeptDofMaster;

    /**array containing list of all rank(process ID) associated with each dependent dof */
    int *N_DeptDofNeibs;

    /**Max number of Dept dofs among all process*/
    int MaxN_DeptDofs_All;

    /**array rank(process) ID of each dependent dof */
    int *DeptDofNeibRanks;

    /**array of local dof index of neib process of each dependent dof */
    int *DeptDofNeibRanksLocalDOF;

    /**array of local dof index of neib process of each dependent dof */
    int *DeptDofLocalIndexOfNeibRanks;

    /** number of global degrees of freedom among all processors*/
    int N_GlobalDegreesOfFreedom;

    /**array containing the corresponding global dof of local dof */
    int *GlobalDofOfLocalDof;

    /**total number of dependent dofs neibs */
    int N_Neibs;

    /**DOF neib ranks' IDs of this rank */
    int *NeibsRank;

    /**index in the neib list (NeibsRank) for the given neib rank */
    int *IndexOfNeibRank;

    /** Max number of Communication Steps for all ranks*/
    int Max_CommunicationSteps;

    /** Number of active Processes during Communication  in each step of  communication*/
    int *N_ActiveCommuPross;

    /** Receive rank ID in each Communication Process*/
    int *ReceiveID_Index;

    /** send rank ID in each Communication Process*/
    int *SendID_Index;

    /**request tags for non-blocking communication */
    MPI_Request request001, request002, request003, request004, request005, request006;

    /**Max number of master dofs among all process*/
    int MaxN_MasterDofs_All;

    /**Max number of slave dofs among all process */
    int MaxN_SlaveDofs_All;

    /**array of size [size] contains No. of distributed dof in each processor*/
    int *N_DistDofAll;

    /**number of local dof among all ranks*/
    int *N_LocalDofAllRank;

    /**maximum number of local dof among all ranks, except root */
    int MaxN_LocalDofAllRank;

    /**array containing the correponding global dof of local dof for all ranks*/
    int *GlobalDofOFLocalDofAllRank;

    /**array containing the correponding Dept. Dof-Index Of LocDof*/
    int *DepDofIndexOfLocDof;    
    
    /**total number of DOF neib ranks for this rank */
    int N_DofNeibs;

    /**DOF neib ranks' IDs of this rank */
    int *DofNeibIDs;
    
       /** Number of Communication Processes in each step of  communication*/
    int *N_CommunicationProcesses;
    
  public:
    TParFECommunicator3D(MPI_Comm comm, TFESpace3D *fespace, TSquareStructure3D *Sqstruct);

    ~TParFECommunicator3D();
    
    MPI_Comm GetComm()
     { return Comm;}
    /** Construct DofRankIndex based on used FESpace */
    void ConstructDofMap();
    
    /** Construct DofRankIndex based and find re-ordering */
    void ConstructDofMapRe();
    /**Set new Global numbers */
    void SetNewGN();
    
    /** Reorder matrix after coloring */
    void ColorAndReorder(int start, int end, int &numColors, int *&ptrColors);
    
    void ConstructDofRankIndex();
    
    void ConstructDofRankIndex_old();

    /** Construct DofRankIndex based on used FESpace */
    void ConstructDofRankIndex_New();

    /** Scheduling for communication between processor based on used FESpace */
    void ScheduleParFEComm3D();

    /** Mapping of DOFs between subdomains */
    void MapDofFromNeib3D();

    /** Mapping of DOFs between subdomains */
    void ConstructGlobalDofFromNeib3D();
    
    void SetFENeibCommunicationSteps();
    

    /** Communication between processor based on defined schedule */
    /** arrays with same size */
    int FECommunicateNeib(int **sendbuf, int senddisp, int *sendlen,
                          int **recvbuf, int *recevlen, int N_array);

     int FECommunicateNeib(double **sendbuf, int senddisp, int *sendlen,
                           double **recvbuf, int *recevlen, int N_array);

    int FECommunicateNeib(double *sendbuf, int senddisp, int *sendlen, double *recvbuf, int N_dim);

    int FECommunicateNeib(double *sendbuf, int senddisp, int *sendlen,
                          double *recvbuf, int *recevlen, int N_dim);

    int FECommunicateOneWay(int **sendbuf, int disp, int *sendlen,
                            int **recvbuf, int *recevlen, int N_array, int type);
    
    int FECommunicateOneWay(double *sendbuf, int disp, double *recvbuf,  int type);
    
    /** Mapping of global dof from root to all subdomains */
    void MapDofFromRoot3D();

    TFESpace3D *Getfespace()
      {return FESpace; }

    int WaitForGlobalDofFromNeib3D();

    int WaitForMapDofFromRoot3D();

    void GetDofNeibInfo(int &maxrankpDof, int &n_Neibs, int *&neibsRank, int &n_DeptDofs, 
                        int *&deptDofs, int *&n_DeptDofNeibs, 
                        int &maxN_DeptDofs_All, int *&deptDofNeibRanks, int *&indexOfNeibRank,
                        int *&deptDofLocalIndexOfNeibRanks, int *&deptDofNeibRanksLocalDOF)
     {
       maxrankpDof = MaxSubDomainPerDof;
       n_Neibs = N_Neibs;
       neibsRank = NeibsRank;
       n_DeptDofs = N_DeptDofs;
       deptDofs = DeptDofs;
       n_DeptDofNeibs = N_DeptDofNeibs;
       maxN_DeptDofs_All = MaxN_DeptDofs_All;
       deptDofNeibRanks = DeptDofNeibRanks;
       indexOfNeibRank = IndexOfNeibRank;
       deptDofLocalIndexOfNeibRanks = DeptDofLocalIndexOfNeibRanks;
       deptDofNeibRanksLocalDOF = DeptDofNeibRanksLocalDOF;
      }

    void GetNeibInfo(int &n_Neibs, int *&neibsRank)
     {
       n_Neibs = N_Neibs;
       neibsRank = NeibsRank;
     }

     void GetOwnDofs(int &n_OwnDof, int *&ownDofs)
      {
       n_OwnDof = N_OwnDof;
       ownDofs = OwnDofs;
      }

    int *GetGlobalDofOfLocalDof()
        { return GlobalDofOfLocalDof; }

    void GetLocalDofAllRankInfo(int *&n_LocalDofAllRank, int &maxN_LocalDofAllRank,
                                int *&globalDofOFLocalDofAllRank)
       {
        n_LocalDofAllRank = N_LocalDofAllRank;
        maxN_LocalDofAllRank =MaxN_LocalDofAllRank;
        globalDofOFLocalDofAllRank = GlobalDofOFLocalDofAllRank;
       }

    /** return  N_GlobalDegreesOfFreedom */
    int GetN_GlobalDegreesOfFreedom()
    { return N_GlobalDegreesOfFreedom; }

    /** return  N_GlobalDegreesOfFreedom */
    void GetGlobalDistDofInfo(int &a, int *&b, int *&c)
    {
     a = N_GlobalDegreesOfFreedom;
     b = GlobalDofOfLocalDof;
     c = N_DistDofAll;
    }

    /** return  DeptDofMaster */
    void GetDeptDofMasterInfo(int &a, int &b, int *&c)
     {
      a = MaxN_MasterDofs_All;
      b = MaxN_SlaveDofs_All;
      c = DeptDofMaster; 
     }

    /** return  DeptDofMaster */
    int *GetDeptDofMaster()
    { return DeptDofMaster; }

    /** return  N_DofRankIndex */
    int *GetN_DofRankIndex()
    { return N_DofRankIndex; }

      
     /** return  N_DofRankIndex */
    int *GetDepDofIndexOfLocDof()
    { return DepDofIndexOfLocDof; }
    
    /*=========================== MODOFIED PAR FE COMMUNICATOR using Alltoall =================================*/
     void SetSlaveDofRows(int *Row,int *KCol,double *Values,double *rhs);

    void CommUpdate(double *sol,double *rhs);
    
    void CommUpdateH(double *sol,double *rhs);
    
    void CommUpdateM(double *sol,double *rhs);
   // void CommUpdate(double *sol);   
    void CommUpdateAlltoAllv(double *sol,double *rhs);  
    
    void CommUpdateReduce(double *sol,double *rhs);
    
    void CommUpdate(double *sol);
    
    int *GetMaster()
    {return Master;}
    
    int* GetReorder()
    {return Reorder;}
    
    int GetN_Master()
    {return N_Master;}
    int GetN_Int()
    {return N_Master+N_Int;}
    int GetN_Dept()
    {return N_Master+N_Int+N_Dept;}
    
    int GetN_CMaster()
    {return N_CMaster;}
    int GetN_CDept()
    {return N_CDept;}
    int GetN_CInt()
    {return N_CInt;}
    
    int* GetptrCMaster()
    {return ptrCMaster;}
    int* GetptrCDept()
    {return ptrCDept;}
    int* GetptrCInt()
    {return ptrCInt;}
    
    
};

#endif
#endif

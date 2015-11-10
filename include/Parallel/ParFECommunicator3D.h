/** ***************************************************************************
 *
 * @name       TParFEMapper3D
 * @brief      Class containing communication routines between processes which all
 * 			   administer their own subdomain and corresponding degrees of freedom of
 * 			   an FE system.
 *
 *             Objects of this type are set up after the domain is split, so
 *             there exists a particular one in each process. The data members
 *             refer to the data of the processes' Domain and FESpace.
 *
 *			   All the information the class needs is known to an object of class
 *			   TParFEMapper3D.
 *
 * @author     Sashikumaar Ganesan
 * @date       2015/04/24
 *
 ******************************************************************************/

#ifdef _MPI

#ifndef __PARFECOMMUNICATOR3D__
#define __PARFECOMMUNICATOR3D__

#include "mpi.h"

#include <ParFEMapper3D.h>
#include <FESpace3D.h>
#include <SquareStructure.h>
#include <SquareStructure3D.h>

class TParFECommunicator3D
{
  protected:
    TParFEMapper3D *Mapper;
    
    int N_Dim;
    
    int N_Dof;
    
    int N_SendDof, N_SendDofMS, N_SendDofH1, N_SendDofH2;
    
    double *Send_Info, *Send_InfoMS, *Send_InfoH1, *Send_InfoH2;
    
    double *Recv_Info, *Recv_InfoMS, *Recv_InfoH1, *Recv_InfoH2;
    
    int *N_DofSend, *N_DofSendMS, *N_DofSendH1, *N_DofSendH2; 
    
    int *N_DofRecv, *N_DofRecvMS, *N_DofRecvH1, *N_DofRecvH2;
    
    int *sdispl, *sdisplMS, *sdisplH1, *sdisplH2;
    
    int *rdispl, *rdisplMS, *rdisplH1, *rdisplH2; 
    
    int *DofSend, *DofSendMS, *DofSendH1, *DofSendH2;
    
    int *DofRecv, *DofRecvMS, *DofRecvH1, *DofRecvH2;
    
    int N_Slave, N_InterfaceS, N_Halo1, N_Halo2;
    
    MPI_Comm Comm;
    
    int *Master;
    
  public:
    //TODO Comment the usage of this!
    TParFECommunicator3D(TParFEMapper3D *mapper);
    
    //TODO Comment the usage of this!
    TParFECommunicator3D();
    
    //TODO Comment the usage of this!
    void CommUpdateMS(double *sol);
    
    //TODO Comment the usage of this!
    void CommUpdateH1(double *sol);
    
    //TODO Comment the usage of this!
    void CommUpdateH2(double *sol);
    
    //TODO Comment the usage of this!
    void CommUpdate_M_H1(double *sol);
    
    //TODO Comment the usage of this!
    void CommUpdate(double *sol);
    
    //TODO Comment the usage of this!
    void CommUpdateReduce(double *rhs);
    
    //TODO Comment the usage of this!
    void CommUpdate(double *sol, double *rhs);
    
    //TODO Comment the usage of this!
    void GatherToRoot(double *&GlobalArray, int &GolabalSize, double *LocalArray, int LocalSize, int root);
    
    //TODO Comment the usage of this!
    void ScatterFromRoot(double *GlobalArray, int &GlobalSize, double *LocalArray, int LocalSize, int root);
    
    //TODO Comment the usage of this!
    int *GetMaster()
    {return Master;}
    
    //TODO Comment the usage of this!
    int GetNDof()
    {return N_Dof;}
    
    //TODO Comment the usage of this!
    int* Get_Local2Global()
    { return Mapper->Get_Local2Global();}
    
    //TODO Comment the usage of this!
    int GetN_Master()
    { return Mapper->GetN_Master();}
    
    //Special member functions. Rule of zero applied (class does not manage ressources).

    //Declaration of special member functions - rule of zero

    //! Default copy constructor.
    TParFECommunicator3D(const TParFECommunicator3D&) = default;

    //! Default move constructor.
    TParFECommunicator3D(TParFECommunicator3D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    TParFECommunicator3D& operator=(const TParFECommunicator3D&) = default;

    //! Default move assignment operator
    TParFECommunicator3D& operator=(TParFECommunicator3D&&) = default;

    //! Default destructor.
    ~TParFECommunicator3D() = default;
};

#endif
#endif























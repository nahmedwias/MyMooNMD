// =======================================================================
// @(#)ParVectorNSE3D.h
//
// Class:      TParVectorNSE3D
// Purpose:    Class containing all info needed for communication between subdomains
//             for solution/rhs vectors of NSE 2D/3D (including halo dof if any)
//
// Author:     Sashikumaar Ganesan (12.10.09)
//
// History:    Start of implementation 12.10.09 (Sashikumaar Ganesan)
//
// =======================================================================
#ifdef _MPI

#include <ParVector.h>

#ifndef __ParVectorNSE3D__
#define __ParVectorNSE3D__

class TParVectorNSE3D : public TParVector
{
  protected:

    /**  length of a part of the array (velocity) (including Halo dof) */
    int NU;

    /**  length of a part of the array (pressure) (including Halo dof) */
    int NP;

    /**  length of a part of the array (velocity) (excluding Halo dof) */
    int NU_own;

    /**  length of a part of the array (pressure) (excluding Halo dof) */
    int NP_own;

    /** fespace for which the vector is defined */
    TFESpace3D *VelocityFESpace;

    /** fespace for which the vector is defined */
    TFESpace3D *PressureFESpace;

    /** ParCommunicator of the Velocity vector's Fespace*/
    TParFECommunicator3D *VelocityCommunicator;

    /** ParCommunicator of the Pressure vector's Fespace*/
    TParFECommunicator3D *PressureCommunicator;

  public:
    /** constructor with given values of velocity and pressure*/
    TParVectorNSE3D(MPI_Comm comm, double *U, int n, int np, int N_Dim,  TParFECommunicator3D *velocityCommunicator,
                  TParFECommunicator3D *pressureCommunicator);


   /** Set the values of the vector */
   void SetOwnArray(int nu_own, int np_own, double *values_own)
    { NU_own=nu_own;  NP_own=np_own; Values_own = values_own;  }

   /** return the MPIComm of the vector */
   MPI_Comm GetMPIComm()
    { return Comm; }

   /** return the VeloCommunicator of the vector */
   TParFECommunicator3D *GetVeloCommunicator()
    { return VelocityCommunicator; }

   /** return the Pressure Communicator of the vector */
   TParFECommunicator3D *GetPressureCommunicator()
    { return PressureCommunicator; }

   /** return the values of the vector */
   double *GetValues()
    { return Values; }


   /** calculate norm based of the given method */
   void ParDdot(int assembletype, double &residual, double &impuls_residual);

   /** assemble the vector by adding the subdomain dof values*/
   void AssembleByADD();

   /** assemble the vector (ONLY VELOCITY) by adding the subdomain dof values*/
   void AssembleVeloAtRootByADD(double *GlobalVal);

   /** scatter the global vector from root to all sub domains (including halo dof)*/  
   void ScatterFromRoot(double *GlobalVal);
   
   /** Copy the master dof value to all slave dofs */   
   void AssembleWithMaster();
   
   /** project pressure into L20 */
   void IntoL20Pressure();
   
    /** destructor */
    ~TParVectorNSE3D();

};
#endif
#endif

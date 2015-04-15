// =======================================================================
// @(#)ParVector3D.h
//
// Class:      TParVector3D
// Purpose:    Class containing all info needed for communication between subdomains
//             for solution/rhs vectors
//
// Author:     Sashikumaar Ganesan (05.10.09)
//
// History:    Start of implementation 05.10.09 (Sashikumaar Ganesan)
//
// =======================================================================


#ifdef _MPI
#  include "mpi.h"
#include <ParVector.h>

#ifndef __PARVECTOR3D__
#define __PARVECTOR3D__

#include <ParFECommunicator3D.h>

class TParVector3D : public TParVector
{
  protected:

    /** length of the vector array (without Halo dof) */
    int N_own;
    
    /** fespace for which the vector is defined */
    TFESpace3D *FESpace;

    /** ParCommunicator of the vector's Fespace*/
    TParFECommunicator3D *Communicator;

  public:
    TParVector3D(MPI_Comm comm, double *U, int ndof, int N_Dim,  TParFECommunicator3D *communicator);

   /** assemble the Global array (slave dep. dofs and Hallo dofs only) with the Master dof value for iterative solver */    
    void UpdateDepAndHalloDof(double *GlobalVal);
    
    /** L2 norm of the vector */
    void ParDdot(int assembletype, double *residual);

   /** assemble the vector by adding  dofs values on subdomain interface */
    void AssembleByADD(double *tempval);

   /** Copy the master dof value to all slave dofs */
   void AssembleWithMaster();
   
   /** assemble the vector by averaging dof values on subdomain interface */   
   void AssembleByAverage();

   /** assemble the global vector at root using values from all sub domains */
   void AssembleAtRootByADD(double *GlobalVal);
    
   /** scatter the global vector from root to all sub domains */
   void ScatterFromRoot(double *GlobalVal);

   void SetMastDofToHalloNeib(int N_Neibs, int *N_MastDofToHalloNeib, 
                              int N_MastDofToHalloNeib_All, int *MastDofToHalloNeib);
      
   /** destructor */
    ~TParVector3D();

};

#endif

#endif

// =======================================================================
// @(#)ParVector.h
//
// Class:      TParVector
// Purpose:    Class containing all info needed for communication between subdomains
//             for solution/rhs vectors (including halo dof if any)
//
// Author:     Sashikumaar Ganesan (05.10.09)
//
// History:    Start of implementation 05.10.09 (Sashikumaar Ganesan)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"

#ifndef __PARVECTOR__
#define __PARVECTOR__

class TParVector
{
  protected:

    /** MPI_Comm for which the fespace communications are needed */
     MPI_Comm Comm;

    /** vector array (including Halo dof) */
    double *Values;

    /** vector array (excluding Halo dof) */
    double *Values_own;
        
    /** length of the vector array */
    int Len;
    
    /** dimension of the vector  */
    int N_DIM;

    /** leading dimension  of the vector  */
    int LDIM;    
    
   /** No. of Master Dofs which is needed in other processors in iterative solver in each neib */
   /** not only neib dof but also hallo dofs for the entire matrix row */
    int *TCLS_N_MastDofToHalloNeib;
    
   /** No. of Master Dofs which is needed in other processors in iterative solver in each neib*/
    int TCLS_N_MastDofToHalloNeib_All;
    
   /** list of Master Dofs which is needed in other processors in iterative solver in each neib*/
    int *TCLS_MastDofToHalloNeib;    

    MPI_Request requestAssembleAtRoot;
    
  public:

    /** constructor with given values*/
    TParVector(MPI_Comm comm, double *U, int n, int N_Dim);

   /** return the values of the vector */
   int GetDimension()
    { return N_DIM; }

   /** return the values of the vector */
   double *GetValues()
    { return Values; }   
    
    /** destructor */
    ~TParVector();
};
#endif
#endif

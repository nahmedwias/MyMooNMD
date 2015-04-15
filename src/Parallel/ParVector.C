// =======================================================================
// @(#)ParVector.C
//
// Class:      TParVector
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ParVector.h>


/** constructor with given values*/
TParVector::TParVector(MPI_Comm comm, double *U, int n, int N_Dim)
{
 Comm = comm;
 Values = U; 
 Len = n;
 N_DIM= N_Dim;
 LDIM = n/N_Dim;
 TCLS_N_MastDofToHalloNeib = NULL;
 TCLS_N_MastDofToHalloNeib_All = -1;
 TCLS_MastDofToHalloNeib = NULL; 
}



TParVector::~TParVector()
{
 if(Values)  delete [] Values;
}
#endif

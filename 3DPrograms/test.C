/* -------------------------------------------------------------------- */
/*      Example program to show the use of the "PARDISO" routine        */
/*      on for unsymmetric linear systems                               */
/* -------------------------------------------------------------------- */
/*      This program can be downloaded from the following site:         */
/*      http://www.pardiso-project.org                                  */
/*                                                                      */
/*  (C) Olaf Schenk, Institute of Computational Science                 */
/*      Universita della Svizzera italiana, Lugano, Switzerland.        */
/*      Email: olaf.schenk@usi.ch                                       */
/* -------------------------------------------------------------------- */

// C++ compatible

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <mpi.h>

#ifdef _HYBRID
#define _AnyParallel_
#endif
#ifdef _MPI
#define _AnyParallel_
#endif
#ifdef _OPENMP
#define _AnyParallel_
#endif

int main()
{
  int i;
  
#ifdef _HYBRID
  printf("IN HYBRID\n");
#endif

#ifdef _MPI
  printf("IN MPI\n");
#endif
  
#ifdef _OPENMP
  printf("IN OPENMP\n");
#endif

#ifdef _AnyParallel_  
  printf("#############\n");
#endif
//   int a =3;
//   printf("iparam(3)::%d\nNumber of OpenMP tasks per host",a);
//     printf("iparam(11)::%d\n", iparam(11), "Scaling");
  
  return 0;
}
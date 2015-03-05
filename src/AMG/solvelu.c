#define Int long

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "clapack.h"

#include "memc.h"
#include "dmem.h"
#include "utils.h"
#define ZSIZE 1
#include "zahl.h"
#include "ilu.h"
#include "zutils.h"

    #include <ctype.h>

    __const unsigned short int *__ctype_b;
    __const __int32_t *__ctype_tolower;
    __const __int32_t *__ctype_toupper;

    void ctSetup()
    {
    __ctype_b = *(__ctype_b_loc());
    __ctype_toupper = *(__ctype_toupper_loc());
    __ctype_tolower = *(__ctype_tolower_loc());
    }


// user parameters
 double iluEps1=0.01;  // should be ca. 100 times smaller then iluEps1
 double iluEps2=0.5;    // the main eps for preconditioner
 double iluEps3=0.0;    // for experts
 double iluShift=0.0;   // unitary shift
 double IterEps=1.0e-10; // the main eps for iterative method

 int iluMem=700*1024*1024; // memory usage - should be smaller that the total amount of main memory
 int iluCond=0;             // if your problem is too hard, then set it as 1
 char MatrixName[]="/tmp/matrix.mat";
 char PrecName[]="/tmp/prec.mat";
 char oocName[]="/tmp/data.1";
 int gmr=1;                // gmr=1 - use gmres, gmr=0 - use bicgstab
 int nIter=200;            // the total amount of iterative steps allowed

 LUPrec  LUP;
 LUMat   LUM;

 void MultA(void *X, void *Y)
 { SET(LUM.N, 0., Y, 1);
   BegTim(33);
   iluMV(&LUM, X, Y);
   EndTim(33);
   return;
 }


 void MultP(void *X, void *Y)
 { COPY(LUM.N, X, 1, Y, 1);
   BegTim(34);
   iluPV(&LUP, Y);
   EndTim(34);
   return;
 }


 Int PrintMsg(Int ID, Int Iter, double Val, double MFlop)
 { fprintf(stdout,"[%ld], RES=%g\n", Iter, Val);
   fflush(stdout);
   return 0;
 }


 void SolveLU(int aN, int NNZ, int *IA, int *JA, double *A, double *RHS, double *SOL)
 { Int i, j, k, N=aN, *tmp;
   int fd;

   ctSetup();
   printf("RHS NORM %lg\n", NRM2(N, RHS, 1));
   OOCName=oocName;
   mSetBuf(0,  iluMem, 1);
   unlink(MatrixName);
   fd=open(MatrixName, O_CREAT | O_WRONLY, 0644);
   tmp=mNew(sizeof(Int), N);
   write(fd, &N, sizeof(Int));
   i=0;
   write(fd, &i, sizeof(Int));
   for(i=0; i<N; i++)
   { k=IA[i+1]-IA[i];
     write(fd, &k, sizeof(Int));
     if(k)
     { for(j=0; j<k; j++)
         tmp[j]=JA[j+IA[i]];
       write(fd, tmp, sizeof(Int)*k);
       write(fd, A+IA[i], sizeof(double)*k);
     }
   }
   mFree(tmp);
   close(fd);

   ClrTim(22);
   BegTim(22);
   iluComputeSP(MatrixName, PrecName, iluShift, iluEps1, iluEps2, iluEps3, iluCond, OOCName);
   EndTim(22);
   PrnTim(22, "ILU ");
   iluLoadM(&LUM, MatrixName);
   iluLoadP(&LUP, PrecName);
   ClrTim(22);
   BegTim(22);
   SET(LUM.N, 0., SOL, 1);
   if(gmr)
     mgmres(LUM.N, MultA, MultP, PrintMsg, (void*)SOL, (void*)RHS, 1, nIter, 0, 1, IterEps, 0);
   else
     bicgstab(LUM.N, MultA, MultP, PrintMsg, (void*)SOL, (void*)RHS, nIter, 0, IterEps, 0);
   EndTim(22);
   PrnTim(22, "ITER ");
   mSetBuf(0, 0, -1);
   return;
 }


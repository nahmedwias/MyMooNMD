// =======================================================================
// @(#)ParDirectSolver.h
//
// Class:      TParDirectSolver
// Purpose:    Solve equation system by OpenMP or MPI based direct solvers
//
// Author:     Sashikumaar Ganesan (30.06.09)
//
// History:    Start of implementation 30.06.09 (Sashikumaar Ganesan)
//
// =======================================================================
# ifdef _OMP
#include <omp.h>
# endif

#include <stdlib.h>
#include <MainUtilities.h>
#include <ParDirectSolver.h>
#include <Database.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __OMPONLY__
extern "C"
{
  int pardisoinit_(void **, int *, int *);
}

extern "C"
{
   int pardiso_(void **, int *, int *, int *, int *, int *,
                      double *, int *, int *, int *, int *, int *,
                      int *, double *, double *, int *);

}
#endif
int TParDirectSolver::Init_SMP_Default(int InN_Eqn, int *InRowPtr, int *InKCol, int InN_Entries)
{
 int phase, error;

 char *var;

  mtype=11; 
  maxfct = 1;         /* Maximum number of numerical factorizations.  */
  mnum   = 1;         /* Which factorization to use. */

  if( TDatabase::ParamDB->SC_VERBOSE )
   msglvl = 1;         /* Print statistical information  */
  else
   msglvl = 0; 

  error  = 0;         /* Initialize error flag */
  nrhs   = 1;        /* number of rhs */

  N_Eqn = InN_Eqn; 
  N_Entries = InN_Entries;

  RowPtr = new int[N_Eqn+1];
  KCol = new int[N_Entries];

  memcpy(RowPtr, InRowPtr, (N_Eqn+1)*SizeOfInt);
  memcpy(KCol, InKCol, N_Entries*SizeOfInt);

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
     sscanf( var, "%d", &num_threads );
    else 
     {
        cout<<"Set environment OMP_NUM_THREADS to 1"<<endl;
        exit(1);
     }

    pt = new void*[64];
    iparm = new int [64];

    memset(iparm, 0, 64*SizeOfInt);
    iparm[2]  = num_threads; /* wrong in pardiso manual  */

/*  initialize the solvers      */
#ifdef __OMPONLY__
    pardisoinit_(pt,  &mtype, iparm);
# else
    OutPut("Error: Define __OMPONLY__  during compiling the main prg !!!!" <<endl)
    exit(0);
#endif
//     iparm[1] = 2; //Metis package for Fill-In reduction
//     iparm[3] = 0;
//     iparm[4] = 0;
//     iparm[7] = 0; /* Max numbers of iterative refinement steps. */
//     iparm[8] = 0; // unused
// 
// 
// 
//     iparm[10] = 1; //
//     iparm[11] = 0; // unused
//     iparm[12] = 1; //
//     iparm[17] = 0; //
//     iparm[18] = 0; //
//     iparm[20] = 0; //
    cout<<" Parallel Direct solver initialize completed "<<endl;

   return 0;
}

void TParDirectSolver::ParDirectSolver_SMP_Terminate(double *Values)
{
 int phase, error;

   phase = -1;                 /* Release internal memory. */
   error  = 0;
#ifdef __OMPONLY__
   pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eqn, Values, RowPtr, KCol,
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
# else
    OutPut("Error: Define __OMPONLY__  during compiling the main prg !!!!" <<endl)
    exit(0);
#endif
}

TParDirectSolver::~TParDirectSolver()
{

   delete [] RowPtr;
   delete [] KCol;
}


TParDirectSolver::TParDirectSolver(int InN_Eqn, int *InRowPtr, int *InKCol,
                                   int InN_Entries)
{
  Init_SMP_Default(InN_Eqn, InRowPtr, InKCol, InN_Entries);
}

void TParDirectSolver::SymbolicFactorize(double *Entries)
{
 int phase, error;

    error  = 0;
    phase = 11; /* .Analysis and Numerical factorization.   */
#ifdef __OMPONLY__
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eqn, Entries, RowPtr, KCol,
             &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
# else
    OutPut("Error: Define __OMPONLY__  during compiling the main prg !!!!" <<endl)
    exit(0);
#endif
    if(error != 0)
     {
      OutPut("ERROR during numerical factorization: "<< error<<endl);
      exit(0);
     }

    cout<<" Analysis and Numerical Factorization completed "<<endl;

}


void TParDirectSolver::FactorizeAndSolve(TSquareMatrix *matrix, double *rhs, double *sol)
{
  int phase, error=0;

  double *MatValues;

  MatValues = matrix->GetEntries();


/* -------------------------------------------------------------------- */    
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */    
    phase = 22;
#ifdef __OMPONLY__
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eqn, MatValues, RowPtr, KCol,
             &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
# else
    OutPut("Error: Define __OMPONLY__  during compiling the main prg !!!!" <<endl)
    exit(0);
#endif
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    printf("\nFactorization completed ...\n ");

/* .Back substitution and iterative refinement. */
    phase = 33;
    error  = 0;
#ifdef __OMPONLY__
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eqn, MatValues, RowPtr, KCol,
             &idum, &nrhs, iparm, &msglvl, rhs, sol, &error);
# else
    OutPut("Error: Define __OMPONLY__  during compiling the main prg !!!!" <<endl)
    exit(0);
#endif
    if(error != 0)
     {
      OutPut("ERROR during solution:  "<< error<<endl);
      exit(0);
     }


    if(iparm[19]<0)
     {
      cout<<"Solve completed ... iparm[19] "<< iparm[19] <<endl;
//       exit(0);
     }

}


void TParDirectSolver::Solve(TSquareMatrix *matrix, double *rhs, double *sol)
{
  int phase, error;

  double *MatValues;

  MatValues = matrix->GetEntries();

/* .Back substitution and iterative refinement. */
    iparm[3] = 61;
    phase = 23;
    error  = 0;
#ifdef __OMPONLY__
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eqn, MatValues, RowPtr, KCol,
             &idum, &nrhs, iparm, &msglvl, rhs, sol, &error);
# else
    OutPut("Error: Define __OMPONLY__  during compiling the main prg !!!!" <<endl)
    exit(0);
#endif
    if(error != 0)
     {
      OutPut("ERROR during solution:  "<< error<<endl);
      exit(0);
     }


    if(iparm[19]<0)
     {
      cout<<"Solve completed ... iparm[19] "<< iparm[19] <<endl;
//       exit(0);
     }

}


#ifdef __3D__
void TParDirectSolver::FactorizeAndSolve(TSquareMatrix3D *sqmatrixA11,TSquareMatrix3D *sqmatrixA12,
                         TSquareMatrix3D *sqmatrixA13, TSquareMatrix3D *sqmatrixA21,
                         TSquareMatrix3D *sqmatrixA22,TSquareMatrix3D *sqmatrixA23,
                         TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
                         TSquareMatrix3D *sqmatrixA33,
                         double *sol, double *rhs)
{
  int i,j,k, N_U, *KColA, *RowPtrA, th_id;
  int begin, end, pos, l, phase, len;
  int N_MatEntries, pos1, pos2, pos3, error;

  double *EntriesA11, *EntriesA12, *EntriesA13;
  double *EntriesA21, *EntriesA22, *EntriesA23;
  double *EntriesA31, *EntriesA32, *EntriesA33;
  double *Values;

  N_U = sqmatrixA11->GetN_Rows();
  KColA = sqmatrixA11->GetKCol();
  RowPtrA = sqmatrixA11->GetRowPtr();
  N_MatEntries = RowPtrA[N_U];

  Values = new double[N_Entries];

  EntriesA11 = sqmatrixA11->GetEntries();
  EntriesA12 = sqmatrixA12->GetEntries();
  EntriesA13 = sqmatrixA13->GetEntries();
  EntriesA21 = sqmatrixA21->GetEntries();
  EntriesA22 = sqmatrixA22->GetEntries();
  EntriesA23 = sqmatrixA23->GetEntries();
  EntriesA31 = sqmatrixA31->GetEntries();
  EntriesA32 = sqmatrixA32->GetEntries();
  EntriesA33 = sqmatrixA33->GetEntries();

//  #pragma omp parallel for shared(N_U, RowPtrA, Values, EntriesA11, EntriesA12, EntriesA13, \
//                           EntriesA21, EntriesA22, EntriesA23, EntriesA31, EntriesA32, EntriesA33) \
//                           private(i, j, begin, end, len, pos1, pos2, pos3, th_id)
  for(i=0;i<N_U;i++)
   {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    len = end - begin;

    pos1 = 3*begin;
    for(j=begin;j<end;j++)
     {
      Values[pos1] = EntriesA11[j];
      Values[len + pos1] = EntriesA12[j];
      Values[2*len + pos1] = EntriesA13[j];
      pos1++;
     }

    pos2 = 3*N_MatEntries + 3*begin;
    for(j=begin;j<end;j++)
     {
      Values[pos2] = EntriesA21[j];
      Values[len + pos2] = EntriesA22[j];
      Values[2*len + pos2] = EntriesA23[j];
      pos2++;
     }

    pos3 = 6*N_MatEntries + 3*begin;
    for(j=begin;j<end;j++)
     {
      Values[pos3] = EntriesA31[j];
      Values[len + pos3] = EntriesA32[j];
      Values[2*len + pos3] = EntriesA33[j];
      pos3++;
     }
   } // for(i=0;i<N_U;i++)


/** test **/
//     cout<<" iparm[3] Solve completed ... "<< iparm[3] <<endl;
//     cout<< " N_Eqn   "<< N_Eqn<<endl;
//     cout<< " N_Entries   "<< N_Entries<<endl;
//     cout<< " N_Entries   "<< 9*RowPtrA[N_U]<<endl;
//     cout<< " mtype   "<<mtype<<endl;
// 

/* -------------------------------------------------------------------- */    
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */    
    phase = 22;
    error  = 0;
#ifdef __OMPONLY__
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eqn, Values, RowPtr, KCol,
             &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
# else
    OutPut("Error: Define __OMPONLY__  during compiling the main prg !!!!" <<endl)
    exit(0);
#endif
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    printf("\nFactorization completed ...\n ");


/* .Back substitution and iterative refinement. */
    phase = 33;
    error  = 0;
#ifdef __OMPONLY__
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eqn, Values, RowPtr, KCol,
             &idum, &nrhs, iparm, &msglvl, rhs, sol, &error);
# else
    OutPut("Error: Define __OMPONLY__  during compiling the main prg !!!!" <<endl)
    exit(0);
#endif
    if(error != 0)
     {
      OutPut("ERROR during solution:  "<< error<<endl);
      exit(0);
     }

    if(iparm[19]<0)
     {
      cout<<"Solve completed ... iparm[19] "<< iparm[19] <<endl;
//       exit(0);
     }
 delete [] Values;
}




void TParDirectSolver::Solve(TSquareMatrix3D *sqmatrixA11,TSquareMatrix3D *sqmatrixA12,
                         TSquareMatrix3D *sqmatrixA13, TSquareMatrix3D *sqmatrixA21,
                         TSquareMatrix3D *sqmatrixA22,TSquareMatrix3D *sqmatrixA23,
                         TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
                         TSquareMatrix3D *sqmatrixA33,
                         double *sol, double *rhs)
{
  int i,j,k, N_U, *KColA, *RowPtrA, th_id;
  int begin, end, pos, l, phase, len;
  int N_MatEntries, pos1, pos2, pos3, error;

  double *EntriesA11, *EntriesA12, *EntriesA13;
  double *EntriesA21, *EntriesA22, *EntriesA23;
  double *EntriesA31, *EntriesA32, *EntriesA33;
  double *Values;

  N_U = sqmatrixA11->GetN_Rows();
  KColA = sqmatrixA11->GetKCol();
  RowPtrA = sqmatrixA11->GetRowPtr();
  N_MatEntries = RowPtrA[N_U];

  Values = new double[N_Entries];

  EntriesA11 = sqmatrixA11->GetEntries();
  EntriesA12 = sqmatrixA12->GetEntries();
  EntriesA13 = sqmatrixA13->GetEntries();
  EntriesA21 = sqmatrixA21->GetEntries();
  EntriesA22 = sqmatrixA22->GetEntries();
  EntriesA23 = sqmatrixA23->GetEntries();
  EntriesA31 = sqmatrixA31->GetEntries();
  EntriesA32 = sqmatrixA32->GetEntries();
  EntriesA33 = sqmatrixA33->GetEntries();

//  #pragma omp parallel for shared(N_U, RowPtrA, Values, EntriesA11, EntriesA12, EntriesA13, \
//                           EntriesA21, EntriesA22, EntriesA23, EntriesA31, EntriesA32, EntriesA33) \
//                           private(i, j, begin, end, len, pos1, pos2, pos3, th_id)
  for(i=0;i<N_U;i++)
   {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    len = end - begin;

    pos1 = 3*begin;
    for(j=begin;j<end;j++)
     {
      Values[pos1] = EntriesA11[j];
      Values[len + pos1] = EntriesA12[j];
      Values[2*len + pos1] = EntriesA13[j];
      pos1++;
     }

    pos2 = 3*N_MatEntries + 3*begin;
    for(j=begin;j<end;j++)
     {
      Values[pos2] = EntriesA21[j];
      Values[len + pos2] = EntriesA22[j];
      Values[2*len + pos2] = EntriesA23[j];
      pos2++;
     }

    pos3 = 6*N_MatEntries + 3*begin;
    for(j=begin;j<end;j++)
     {
      Values[pos3] = EntriesA31[j];
      Values[len + pos3] = EntriesA32[j];
      Values[2*len + pos3] = EntriesA33[j];
      pos3++;
     }
   } // for(i=0;i<N_U;i++)


/** test **/


/* .Back substitution and iterative refinement. */
    iparm[3] = 61;
    phase = 23;
    error  = 0;
#ifdef __OMPONLY__
    pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N_Eqn, Values, RowPtr, KCol,
             &idum, &nrhs, iparm, &msglvl, rhs, sol, &error);
# else
    OutPut("Error: Define __OMPONLY__  during compiling the main prg !!!!" <<endl)
    exit(0);
#endif
    if(error != 0)
     {
      OutPut("ERROR during solution:  "<< error<<endl);
      exit(0);
     }

    if(iparm[19]<0)
     {
      cout<<"Solve completed ... iparm[19] "<< iparm[19] <<endl;
//       exit(0);
     }
 delete [] Values;
}


#endif


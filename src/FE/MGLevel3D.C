// =======================================================================
// %W% %G%
//
// Class:       TMGLevel3D
// Purpose:     store all data for one level in a multi grid method in 3d
//
// Author:      Gunar Matthies 26.06.2000
//
// History:     26.06.2000 start of implementation
//
// =======================================================================

#include <MGLevel3D.h>
#include <FESpace3D.h>
#include <Database.h>
#include <LinAlg.h>
#include <MooNMD_Io.h>

#include <stdlib.h>
#include <string.h>
#include <ParFECommunicator3D.h>
#include <FEFunction3D.h>

/** constructor */
TMGLevel3D::TMGLevel3D(int level, TSquareMatrix3D *a,
                       double *rhs, double *sol, int n_aux,
                       int *permutation)
{
  int i;
  double *aux;

  Level = level;

  FESpace = a->GetFESpace();

  N_Active = FESpace->GetN_ActiveDegrees();
  HangingNodeBound = FESpace->GetHangingBound();
  N_Dirichlet = FESpace->GetN_Dirichlet();
  N_DOF = FESpace->GetN_DegreesOfFreedom();

/*
  cout << "N_Active: " << N_Active << endl;
  cout << "HangingNodeBound: " << HangingNodeBound << endl;
  cout << "N_Dirichlet: " << N_Dirichlet << endl;
  cout << "N_DOF: " << N_DOF << endl;
*/

  A = a;
  MatrixStructure = a->GetMatrixStructure();
  RowPtr = a->GetRowPtr();
  KCol = a->GetKCol();
  Entries = a->GetEntries();

  Rhs = rhs;
  X = sol;

  N_Aux = n_aux;
  Aux = new double* [N_Aux]; 
  aux = new double[N_Aux*N_DOF];
  for(i=0;i<N_Aux;i++)
    Aux[i] = aux+i*N_DOF;

  Additional = NULL;

  Permutation = permutation;
}


#ifdef _MPI
/** constructor for parallel */
TMGLevel3D::TMGLevel3D(int level, TSquareMatrix3D *a, double *rhs, double *sol, 
                       TFEFunction3D *c, TParFECommunicator3D *parComm,TFESpace3D *ownScalarSpace, int n_aux,
                       int *permutation)
{
  int i;
  double *aux;
  char UString[] = "u";

  Level = level;

  FESpace = a->GetFESpace();

  N_Active = FESpace->GetN_ActiveDegrees();
  HangingNodeBound = FESpace->GetHangingBound();
  
 // if(N_Active == HangingNodeBound) printf("actve hangbnd equal\n");
  N_Dirichlet = FESpace->GetN_Dirichlet();
  N_DOF = FESpace->GetN_DegreesOfFreedom();
  
  Reorder = parComm->GetReorder();
  N_Master = parComm->GetN_Master();
  N_Int = parComm->GetN_Int();
  N_Dept = parComm->GetN_Dept();
  
  N_CMaster = parComm->GetN_CMaster();
  ptrCMaster = parComm->GetptrCMaster();
  N_CInt = parComm->GetN_CInt();
  ptrCInt = parComm->GetptrCInt();
  N_CDept = parComm->GetN_CDept();
  ptrCDept = parComm->GetptrCDept();
/*
  cout << "N_Active: " << N_Active << endl;
  cout << "HangingNodeBound: " << HangingNodeBound << endl;
  cout << "N_Dirichlet: " << N_Dirichlet << endl;
  cout << "N_DOF: " << N_DOF << endl;
*/

  A = a;
  C = c;
  MatrixStructure = a->GetMatrixStructure();
  RowPtr = a->GetRowPtr();
  KCol = a->GetKCol();
  Entries = a->GetEntries();

  Rhs = rhs;
  X = sol;

  N_Aux = n_aux;
  Aux = new double* [N_Aux]; 
  aux = new double[N_Aux*N_DOF];
  for(i=0;i<N_Aux;i++)
    Aux[i] = aux+i*N_DOF;

  Additional = NULL;

  Permutation = permutation;
  
  
  ParComm = parComm;
  OwnScalarSpace = ownScalarSpace;  
  OwnN_DOF = OwnScalarSpace->GetN_DegreesOfFreedom();
  OwnSolArray = new double[OwnN_DOF];
 
  OwnC = new TFEFunction3D(OwnScalarSpace, UString, UString, OwnSolArray, OwnN_DOF);
}
#endif


/** destructor */
TMGLevel3D::~TMGLevel3D()
{
  delete Aux[0];
  delete Aux;

  if(Additional)
    delete Additional;
} // ~TMGLevel3D

/** return i-th auxiliary vector */
double *TMGLevel3D::GetAuxVector(int i)
{
  double *ret;

  if(i<N_Aux)
    ret = Aux[i];
  else
  {
    cerr << "Not enough aux vectors in MGLevel3D.C!" << endl;
    exit(-1);
  }
 
  return ret;
} // GetAuxVector


double tD=0.0,tS=0.0;
void TMGLevel3D::Defect(double *sol, double *f, double *d, double &res)
{
  double t1,t2;
#ifdef _MPI
  t1 = MPI_Wtime();
#endif
  ScalarDefect(A, sol, f, d, res);
  #ifdef _MPI
  t2 = MPI_Wtime();
  tS += (t2-t1);
  t1 = MPI_Wtime();
#endif
#ifdef _MPI  
  int i, rank, *MasterOfDof, dof,numThreads= TDatabase::ParamDB->OMPNUMTHREADS;
  double res_global,res1;
  
  MPI_Comm_rank(ParComm->GetComm(), &rank); 
  MasterOfDof = ParComm->GetMaster();
  
#ifdef _HYBRID
  omp_set_num_threads(numThreads);
#pragma omp parallel default(shared) private(i)
{
#pragma omp for schedule(guided) nowait reduction(+:res1)
#endif
  for(i=0; i<N_DOF; i++)
    if(MasterOfDof[i] == rank)
      res1 += d[i]*d[i];
#ifdef _HYBRID
}
#endif

  MPI_Allreduce(&res1, &res_global, 1, MPI_DOUBLE, MPI_SUM, ParComm->GetComm());
  res = sqrt(res_global); 
#endif 
#ifdef _MPI
t2 = MPI_Wtime();
tD +=(t2-t1);
#endif
} // end Defect

// SOR smoother
void TMGLevel3D::SOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int ii, i,j,k,l,index;
  double s, t, diag;
  double omega;

  omega = Parameters[0];
 
  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,ii,s,k,j,index,diag,t)
{
#pragma omp for schedule(dynamic) nowait 
#endif
 
  for(ii=0;ii<N_Active;ii++)
  {
    i = ii;
    // i = Permutation[ii];
    // cout << "row: " << i << "   " << endl;
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == i)
      {
        diag = Entries[j];
      }
      else
      {
        s -= Entries[j] * sol[index];
      }
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i


  // set hanging nodes 
  j = RowPtr[N_Active];
#ifdef _HYBRID
#pragma omp for schedule(dynamic) nowait 
#endif
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i
#ifdef _HYBRID
}
#endif

} // SOR

#ifdef _MPI
void TMGLevel3D::SOR_Re(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int ii, i,j,k,l,index,rank;
  double s, t, diag;
  double omega;

  omega = Parameters[0];
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);
  
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,ii,s,k,j,index,diag,t)
{
#pragma omp for schedule(dynamic) nowait 
#endif

  for(ii=0;ii<N_Master;ii++)
  {
    i = ii;
    // i = Permutation[ii];
    // cout << "row: " << i << "   " << endl;
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == i)
      {
        diag = Entries[j];
      }
      else
      {
        s -= Entries[j] * sol[index];
      }
    } // endfor j
 
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i
  
//  MPI_Barrier(MPI_COMM_WORLD);
  ParComm->CommUpdateM(sol,f);
//  MPI_Barrier(MPI_COMM_WORLD);
     for(ii=N_Master;ii<N_Int;ii++)
  {
    i = ii;
    // i = Permutation[ii];
    // cout << "row: " << i << "   " << endl;
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == i)
      {
        diag = Entries[j];
      }
      else
      {
        s -= Entries[j] * sol[index];
      }
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i
//   MPI_Barrier(MPI_COMM_WORLD);
 
  for(ii=N_Int;ii<N_Dept;ii++)
  {
    i = ii;
    // i = Permutation[ii];
    // cout << "row: " << i << "   " << endl;
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == i)
      {
        diag = Entries[j];
      }
      else
      {
        s -= Entries[j] * sol[index];
      }
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i
//   MPI_Barrier(MPI_COMM_WORLD);
  ParComm->CommUpdateH(sol,f);
//   MPI_Barrier(MPI_COMM_WORLD);
  

 /* 
  for(ii=N_Dept;ii<N_Active;ii++)
  {
    i = Reorder[ii];
    // i = Permutation[ii];
    // cout << "row: " << i << "   " << endl;
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == i)
      {
        diag = Entries[j];
      }
      else
      {
        s -= Entries[j] * sol[index];
      }
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i*/
// MPI_Barrier(MPI_COMM_WORLD);
  // set hanging nodes 
  j = RowPtr[N_Active];
#ifdef _HYBRID
#pragma omp for schedule(dynamic) nowait 
#endif
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i
#ifdef _HYBRID
}
#endif

} // SOR_Re

double tSor=0.0;
//double tM=0.0,tH=0.0;
void TMGLevel3D::SOR_Color(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int ii, i,j,jj,k,l,index,rank,tid,nrows=0,numThreads;
  double s, t, diag,t1,t2;
  double omega;
  numThreads = TDatabase::ParamDB->OMPNUMTHREADS;
  #ifdef _MPI
  t1 = MPI_Wtime();
#endif
  omega = Parameters[0];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

#ifdef _HYBRID
  omp_set_num_threads(numThreads);
#pragma omp parallel default(shared) private(i,ii,s,k,j,jj,tid,index,diag,t) 
{
  tid = omp_get_thread_num();
#endif
  for(ii=0;ii<N_CMaster;ii++)
  {
#ifdef _HYBRID
    #pragma omp for schedule(guided) 
#endif
    for(jj=ptrCMaster[ii];jj<ptrCMaster[ii+1];jj++)
    {
//       nrows+=1;
     i = jj;
     // i = Permutation[ii];
     // cout << "row: " << i << "   " << endl;
      s = f[i];
      k = RowPtr[i+1];
      for(j=RowPtr[i];j<k;j++)
      {
       index = KCol[j];
       if(index == i)
       {
          diag = Entries[j];
       }
       else
       {
        s -= Entries[j] * sol[index];
       }
      }  // endfor j
 
      t = sol[i];
     sol[i] = omega*(s/diag-t) + t;
    }
//      #pragma omp barrier
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i
  
//   MPI_Barrier(MPI_COMM_WORLD);
//  if(tid==0)
#pragma omp master
{
    ParComm->CommUpdateM(sol,f);
}
 //cout<<"01"<<"\n";
 //   MPI_Barrier(MPI_COMM_WORLD);
 //cout<<"02"<<"\n";
 #ifdef _HYBRID
  #pragma omp for schedule(guided) 
#endif
    for(jj=ptrCInt[0];jj<ptrCInt[1];jj++)
    {
//      nrows+=1;
     i = jj;
     
     // i = Permutation[ii];
     // cout << "row: " << i << "   " << endl;
     s = f[i];
     k = RowPtr[i+1];
     for(j=RowPtr[i];j<k;j++)
     {
      index = KCol[j];
      if(index == i)
      {
        diag = Entries[j];
      }
      else
      {
        s -= Entries[j] * sol[index];
      }
     } // endfor j
     t = sol[i];
     sol[i] = omega*(s/diag-t) + t;
    }
 #pragma omp barrier
 
  for(ii=0;ii<N_CDept;ii++)
  {
  
#ifdef _HYBRID
  #pragma omp for schedule(guided) 
#endif
    for(jj=ptrCDept[ii];jj<ptrCDept[ii+1];jj++)
    {
//         nrows+=1;
      i = jj;
      // i = Permutation[ii];
      // cout << "row: " << i << "   " << endl;
       s = f[i];
       k = RowPtr[i+1];
      for(j=RowPtr[i];j<k;j++)
      {
       index = KCol[j];
       if(index == i)
       {
        diag = Entries[j];
       }
       else
       {
         s -= Entries[j] * sol[index];
       }
      } // endfor j
     t = sol[i];
     sol[i] = omega*(s/diag-t) + t;
    }
    
//     #pragma omp barrier
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i
//   MPI_Barrier(MPI_COMM_WORLD);
//  if(tid==0)
 
#pragma omp master
{
   ParComm->CommUpdateH(sol,f);
}

   
  for(ii=1;ii<N_CInt;ii++)
  {
   
#ifdef _HYBRID
  #pragma omp for schedule(guided) 
#endif
    for(jj=ptrCInt[ii];jj<ptrCInt[ii+1];jj++)
    {
//      nrows+=1;
     i = jj;
     
     // i = Permutation[ii];
     // cout << "row: " << i << "   " << endl;
     s = f[i];
     k = RowPtr[i+1];
     for(j=RowPtr[i];j<k;j++)
     {
      index = KCol[j];
      if(index == i)
      {
        diag = Entries[j];
      }
      else
      {
        s -= Entries[j] * sol[index];
      }
     } // endfor j
     t = sol[i];
     sol[i] = omega*(s/diag-t) + t;
    }
//     #pragma omp barrier
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i
   
 //  cout<<"03"<<"\n";
//   MPI_Barrier(MPI_COMM_WORLD);
 /* 
  for(ii=N_Dept;ii<N_Active;ii++)
  {
    i = Reorder[ii];
    // i = Permutation[ii];
    // cout << "row: " << i << "   " << endl;
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == i)
      {
        diag = Entries[j];
      }
      else
      {
        s -= Entries[j] * sol[index];
      }
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i*/
// MPI_Barrier(MPI_COMM_WORLD);
  // set hanging nodes 
  j = RowPtr[N_Active];
#ifdef _HYBRID
#pragma omp for schedule(static) 
#endif
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i

#ifdef _HYBRID
}
#endif
#ifdef _MPI
 t2 = MPI_Wtime();
 tSor += (t2-t1);
#endif
} // SOR_Coloring
#endif
// SSOR smoother
void TMGLevel3D::SSOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k,l,index;
  double s, t, diag;
  double omega;

  omega = Parameters[0];
 
  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, N_Dirichlet*SizeOfDouble);

  // set active nodes
  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index == i)
        diag = Entries[j];
      else
        s -= Entries[j] * sol[index];
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
  } // endfor i

  // set hanging nodes 
  j = RowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i

  // set active nodes
  j = RowPtr[N_Active]-1;
  for(i=N_Active-1;i>=0;i--)
  {
    s = f[i];
    k = RowPtr[i];
    for(;j>=k;j--)
    {
      index = KCol[j];
      if(index == i)
        diag = Entries[j];
      else
        s -= Entries[j] * sol[index];
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
  } // endfor i

  // set hanging nodes 
  j = RowPtr[HangingNodeBound]-1;
  for(i=HangingNodeBound-1;i>=N_Active;i--)
  {
    s = f[i];
    k = RowPtr[i];
    for(;j>=k;j--)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i
} // SSOR

// Jacobi smoother
void TMGLevel3D::Jacobi(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k,l,index;
  double t, s, diag, omega;

  omega = Parameters[0];

  memcpy(aux, sol, N_DOF*SizeOfDouble);

  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

  // set active nodes
  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index == i)
        diag = Entries[j];
      else
        s -= Entries[j] * aux[index];
    } // endfor j
    t = aux[i];
    sol[i] = (1-omega)*t + omega*s/diag;
  } // endfor i

  // set hanging nodes 
  j = RowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i

} // Jacobi

void TMGLevel3D::Update(double *sol, double *upd)
{
  int i;
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i)
 printf("num thrds is %d---\n",omp_get_num_threads());
#pragma omp for schedule(static) nowait 
#endif
  for(i=0;i<N_DOF;i++)
  {
    sol[i] += TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR*upd[i];
  }
}

void TMGLevel3D::Reset(double *vect)
{
  memset(vect, 0, N_DOF*SizeOfDouble);
}

void TMGLevel3D::CorrectNodes(double *vect)
{
  int i, j, k, index;
  double s;

  memset(vect+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);

  // set hanging nodes 
  j = RowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = 0;
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
      {
        s -= Entries[j] * vect[index];
      }
    } // endfor j
    vect[i] = s;
  } // endfor i
} // CorrectNodes

/** correct defect */
void TMGLevel3D::CorrectDefect(double *vect)
{
  memset(vect+N_Active, 0, SizeOfDouble*(N_DOF-N_Active));
}

// block 2x2 smoother
void TMGLevel3D::Block2x2(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k,l,index,maxindex;
  double s1, s2, t, mt, m11, m12, m21, m22, max;
  double omega;

  omega = Parameters[0];

  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

  // set active nodes
  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = RowPtr[i+1];
    max = 0;
    maxindex = -1;
    m11 = m12 = m21 = m22 = 0;
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == i)
      {
        m11 = Entries[j];
      }
      else
      {
        t = Entries[j];
        if( t > max ) // consider only negative entries
        {
          max = t;
          maxindex = index;
        }
      }
    } // endfor j

    s1 = f[i];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == maxindex)
        m12 = Entries[j];
      else 
        if(index != i)
          s1 -= Entries[j] * sol[index];
    } // endfor j

    if(maxindex != -1)
    {
      s2 = f[maxindex];
      k=RowPtr[maxindex+1];
      for(j=RowPtr[maxindex];j<k;j++)
      {
        index = KCol[j];
        if(index == maxindex)
          m22 = Entries[j];
        else
          if(index == i)
            m21 = Entries[j];
          else
            s2 -= Entries[j] *sol[index];
      } // endfor j

      t = (m11*m22 - m21*m12);
      sol[i] = (s1*m22 - s2*m12) / t;
      sol[maxindex] = (s2*m11 - s1*m21) / t;
    }
    else
    {
      // no nondiagonal matrix entry found
      sol[i] = s1/m11;
    }
  } // endfor i

  // set hanging nodes 
  j = RowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s1 = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
      {
        s1 -= Entries[j] * sol[index];
      }
    } // endfor j
    sol[i] = s1;
  } // endfor i

} // Block2x2

// generate ILU decomposition
void TMGLevel3D::ILUDecomposition()
{
  int i,j,jj,k,l,N_;
  double diag, pivot, update;
  int begin, end, beginJ, endJ, found;
  static double beta_ilu = TDatabase::ParamDB->SC_ILU_BETA;

  N_=RowPtr[N_DOF];
  Additional = new double[N_];
  memcpy(Additional, Entries, N_*SizeOfDouble);

  for(i=0;i<N_DOF;i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    diag = Additional[begin];
    if(fabs(diag)<1e-8)
    {
      cerr << "ILU decomposition failed" << endl;
      return;
    }

    for(j=begin+1;j<end;j++)
    {
      if( (jj=KCol[j]) > i)
      {
        beginJ = RowPtr[jj];
        endJ = RowPtr[jj+1];
        found = 0;
        for(k=beginJ+1;k<endJ;k++)
        {
          if(KCol[k] != i) continue;
          pivot = Additional[k]/diag;
          Additional[k] = pivot;
          found = 1;
          break;
        } // endfor k
        if(!found) continue;
        Additional[beginJ] -= pivot*Additional[j];
        for(k=begin+1;k<end;k++)
        {
          if( KCol[k] <= KCol[j] ) continue;
          update = Additional[k] * pivot;
          for(l=beginJ+1;l<endJ;l++)
          {
            if(KCol[k] == KCol[l])
            {
              Additional[l] -= update;
              update = 0;
              break;
            } // endif
          } // endfor l
          Additional[beginJ] += beta_ilu * fabs(update);
        } // endfor k
      } //endif
    } // endfor j
  } // endfor i

/*
  cout << "Matrix A" << endl;
  for(i=0;i<N_DOF;i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin;j<end;j++)
      cout << i << "   " << KCol[j] << "  " << Entries[j] << endl;
  }
  cout << "-----------------" << endl;

  cout << "ILU decomposition" << endl;
  for(i=0;i<N_DOF;i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin;j<end;j++)
      cout << i << "   " << KCol[j] << "  " << Additional[j] << endl;
  }
  cout << "-----------------" << endl;
*/

} // ILUDecomposition

// ILU smoother
void TMGLevel3D::ILU(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k, begin, end;

  if (Additional==NULL)
  {
    if (TDatabase::ParamDB->SC_VERBOSE>1)    
      OutPut("do ILU decomposition" << endl);
    ILUDecomposition();
  }

  // do ILU smoothing

  // invert lower
  for(i=1; i<N_DOF; i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin;j<end;j++)
      if( (k=KCol[j]) < i)
        aux[i] -= Additional[j]*aux[k];
  }

  // invert upper
  for (i=N_DOF-1; i>=0; i--)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin;j<end;j++)
      if( (k=KCol[j]) > i)
        aux[i] -= Additional[j]*aux[k];
    aux[i] /= Additional[begin];
  }

  for(i=0;i<N_DOF;i++)
  {
    sol[i] += 1.0*aux[i];
  }
} // ILU
/** solve exact on this level */
void TMGLevel3D::SolveExact(double *u1, double *rhs1)
{  
  double *a, *b;
  int i,j,k,l,index, N_DOF2 = N_DOF*N_DOF, end, begin, m;
  double value;

  // arrays for matrix and one vector
  a = new double[N_DOF2];
  b = new double[N_DOF];
  // initialize
  memset(a, 0, N_DOF2*SizeOfDouble);
  // fill matrix columnwise
  j = RowPtr[0];
  for(i=0;i<N_DOF;i++)
  {
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      value = Entries[j];
      a[index*N_DOF+i] = value;
    }
  } // endfor i

  // copy into local data
  memcpy(b, rhs1, N_DOF*SizeOfDouble);
  
  // for (i=0;i<N_DOF;i++)
  // cout<< "b("<<i+1<<") = " << b[i] << ";"<<endl;

  SolveLinearSystemLapack(a, b, N_DOF, N_DOF);

  // for (i=0;i<N_DOF;i++)
  //  cout<< "c("<<i+1<<") = " << b[i] << ";"<<endl;
  //exit(1);
  // copy from local data
  memcpy(u1, b, N_DOF*SizeOfDouble);

  delete a;
  delete b;
}

/** step length control multigrid cycle */
double TMGLevel3D::StepLengthControl(double *u, 
                                     double *uold, 
                                     double *def,
                                     int N_Parameters, 
                                     double *Parameters)
{
  double *x,*y,omega,numerator,nominator;
  int i;

  // allocate auxiliary array
  x = new double[2*N_DOF];
  y = x+N_DOF;

  for (i=0;i<N_DOF;i++)
    x[i] = u[i]-uold[i];
  memset(y,0,N_DOF*SizeOfDouble);
  // if (Ddot(N_DOF-N_Active,x+N_Active,x+N_Active)>0)
  //   cout << "Update in non-active nodes found !!" << 
  //     Ddot(N_DOF-N_Active,x+N_Active,x+N_Active) << 
  //     __FILE__ << endl;
  if (Ddot(N_DOF-HangingNodeBound,x+HangingNodeBound,x+HangingNodeBound)>0)
    cout << "Update in non-active nodes found !!" << 
      Ddot(N_DOF-HangingNodeBound,x+HangingNodeBound,x+HangingNodeBound) <<
      __FILE__ << endl;

  
  // compute matrix times update
  MatVect(A,x,y);

  numerator = Ddot(N_DOF,def,y);
  nominator = Ddot(N_DOF,y,y);

  if (nominator > 0)
    omega = numerator/nominator;
  else
    {
      if(N_Parameters>0)
        omega = Parameters[0];
      else
        omega = 0.5;
      cout << "MESSAGE : Step length control failed. Set omega = " << omega << endl;
    }
  if (fabs(omega)<0.1)
    {
       if(N_Parameters>0)
        omega = Parameters[0];
      else
        omega = 0.9;
    }     
  delete x;
  //cout << sqrt(Ddot(N_DOF,def,def)) << " " << numerator << " " << nominator<< " omega " << omega << endl;
  return(omega);
}

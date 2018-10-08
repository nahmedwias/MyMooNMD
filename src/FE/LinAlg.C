// =======================================================================
// @(#)LinAlg.C        1.18 07/03/00
//
// Purpose:     basic routines for linear algebra
//
// Author:      Gunar Matthies          27.01.1999
//              Volker John             27.10.1999  
//              Sashikumaar Ganesan     08.10.2009 (eigen values)
// =======================================================================

#include <Database.h>
#include <MooNMD_Io.h>
#ifdef __2D__
  #include <FEDatabase2D.h>
#else
  #include <FEDatabase3D.h>
#endif  
#include <LinAlg.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

extern "C" {

void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info);

void dgetrs_(char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv,
             double *B, int *ldb, int *info);

void dsptrd_(char *uplo, int *n,  double *AP, double *D, double *E, 
             double *TAU, int *info);

void dopgtr_(char *uplo, int *n,  double *AP, double *TAU, double *Q, int *LDQ, double *work,int *info);

void  dsteqr_(char *compz, int *N, double *D,  double *E, double *Z, int *LDZ, double *work, int *info);
 
}

// ==========================================================================
// general routines independent of dimension
// ==========================================================================
#define AT(i,j) (a[j*LDA+i])
#define A(i,j) (a[i*LDA+j])

void MatVectFull(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  int dof =  TDatabase::ParamDB->INTERNAL_LOCAL_DOF; 
  double *a = (double *)A[0];
  int i,j;

  memset(y,0,dof*SizeOfDouble);
  // matrix is stored as transposed;
  for (i=0;i<dof;i++)
    for (j=0;j<dof;j++)
      y[j] += a[i*dof+j] * x[i];


  return;

}
 
void DefectFull(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
{
  int dof =  TDatabase::ParamDB->INTERNAL_LOCAL_DOF; 
  double *a = (double *)A[0];
  int i,j;

  memcpy(r,b,dof*SizeOfDouble);
  // matrix is stored as transposed;
  for (i=0;i<dof;i++)
    for (j=0;j<dof;j++)
      r[j] -= a[i*dof+j] * x[i];

  return;
}
void SolveLinearSystemLapack(double *a, double *b, int N_Eqn, int LDA)
{
// Arguments:
//    a         double array which contains the matrix columnwise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
  int m, n, nrhs, lda, ldb;
  int *ipivot, info;
  char t='n';

  m = N_Eqn;
  n = N_Eqn;
  lda = N_Eqn;
  ldb = N_Eqn;
  nrhs = 1;

  ipivot = new int[n];
  
  dgetrf_(&m, &n, a, &lda, ipivot, &info);
  dgetrs_(&t, &n, &nrhs, a, &lda, ipivot, b, &ldb, &info);

  delete ipivot;
}
void SolveLinearSystemTranspose(double *a, double *b, int N_Eqn, int LDA)
{
// Arguments:
//    a         double array which contains the matrix columnwise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
  int m, n, nrhs, lda, ldb;
  int *ipivot, info;
  char t='t';

  m = N_Eqn;
  n = N_Eqn;
  lda = N_Eqn;
  ldb = N_Eqn;
  nrhs = 1;

  ipivot = new int[n];
  
  dgetrf_(&m, &n, a, &lda, ipivot, &info);
  dgetrs_(&t, &n, &nrhs, a, &lda, ipivot, b, &ldb, &info);

  delete [] ipivot;
}

/* subroutine for solving a system of linear equations */
void SolveLinearSystem(double *a, double *b, int N_Eqn, int LDA)
// Arguments:
//    a         double array which contains the matrix columnwise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
{
  int i,l,m, row;
  double pivot, tmp, eps=1e-12;


  OutPut("Use SolveLinearSystemNew !!!"<< endl);
  exit(4711);
//  int ii, jj;
/*
  for(ii=0;ii<N_Eqn;ii++)
  {
    for(jj=0;jj<N_Eqn;jj++)
    {
      cout << "a(" << ii+1 << "," << jj+1 << ") = " << setw(8);
      cout << AT(ii,jj) << ";" << endl;
    }
  }
  cout << endl;
*/
  for(i=0;i<N_Eqn-1;i++)
  {
    pivot = 0;
    row = i;
    for(l=i;l<N_Eqn;l++)
    {
      tmp = fabs(AT(l,i));
      if(tmp > pivot)
      {
        pivot = tmp;
        row = l;
      } // endif
    } // endfor l
    //cout << i << " " << pivot << endl;
    if(pivot < eps) 
    {  
      OutPut("pivot " << pivot << " < " << eps);
      OutPut(" Error in solving linear System " << __FILE__ << endl);
      OutPut("equation: "<< i << endl);
      exit(4711);
    }
    if(i<row)
    {
      for(l=0;l<N_Eqn;l++)
      {
        pivot = AT(i,l);
        AT(i,l) = AT(row,l);
        AT(row, l) = pivot;
      }
      tmp = b[i];
      b[i] = b[row];
      b[row] = tmp;
    } // endif

    tmp = AT(i,i);
    for(l=i+1;l<N_Eqn;l++)
      AT(l,i) /= tmp;

    for(l=i+1;l<N_Eqn;l++)
    {
      tmp = AT(i,l);
      for(m=i+1;m<N_Eqn;m++)
        AT(m,l) -= AT(m,i) * tmp;
    }

  } // endfor i

  for(i=0;i<N_Eqn;i++)
  {
    tmp = b[i];
    for(l=i+1;l<N_Eqn;l++)
      b[l] -= AT(l,i)*tmp;
  }

  for(i=N_Eqn-1;i>=0;i--)
  {
    b[i] /= AT(i,i);
    tmp = b[i];
    for(l=0;l<i;l++)
      b[l] -= AT(l,i)*tmp;
  }
}

void SolveMultipleSystemsLapack(double *A, double *B, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs)
// Arguments:
//    A         double array which contains the matrix row wise
//    B         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
//
// This works for a -- stored row wise
//                b -- stored row wise
// The LAPACK Routines for the transposed matrices are used               
//
// The result will be stored row wise
{
  int info, *ipiv;
  char t='t';

  ipiv = new int[N_Eqn];

  dgetrf_(&N_Eqn, &N_Eqn, A, &LDA, ipiv, &info);
  dgetrs_(&t, &N_Eqn, &N_Rhs, A, &LDA, ipiv, B, &LDB, &info);

  delete ipiv;
}

void SolveMultipleSystemsLapack_old(double *A, double *B, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs)
// Arguments:
//    A         double array which contains the matrix columnwise
//              a[i,j] = a[i*LDA+j]
//    B         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
// This works for a -- stored row wise
//                b -- stored row wise
//                     since the LAPACK routines work for
//                     column wise stored right hand sides,
//                     the matrix b will be transposed and
//                     the output matrix containing the solution
//                     will be transposed back
//                LDB = N_Rhs
// The result will be stored column wise
{
  int info, *ipiv;
//  int i, j;
  double *transp;
  char t='t';
  transp = new double[N_Eqn*N_Rhs];
  ipiv = new int[N_Eqn];

  dgetrf_(&N_Eqn, &N_Eqn, A, &LDA, ipiv, &info);

  /*for(i=0;i<N_Eqn;i++)
    for(j=0;j<N_Rhs;j++)
      transp[j+i*N_Rhs]=B[j*N_Eqn+i];
*/
 /*for(i=0;i<N_Eqn*N_Rhs;i++)
      OutPut("B1("<<i<<")= "<<B[i]<<"         T1("<<i<<")= "<<transp[i]<<endl);*/

  dgetrs_(&t, &N_Eqn, &N_Rhs, A, &LDA, ipiv, B, &LDB, &info);

 /*for(i=0;i<N_Eqn*N_Rhs;i++)
      OutPut("B2("<<i<<")= "<<B[i]<<"         T2("<<i<<")= "<<transp[i]<<endl);*/


/*  for(j=0;j<N_Rhs;j++)
    for(i=0;i<N_Eqn;i++)
      B[j*N_Eqn+i]=transp[j+i*N_Rhs];
*/
  /*for(i=0;i<N_Eqn*N_Rhs;i++)
      OutPut("B3("<<i<<")= "<<B[i]<<"         T3("<<i<<")= "<<transp[i]<<endl);*/


//exit(1);

  delete transp;
  delete ipiv;
}


/* subroutine for solving a multiple systems of linear equations */
void SolveMultipleSystemsNew(double *a, double *b, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs)
// Arguments:
//    a         double array which contains the matrix row wise
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
// This works for a -- stored row wise
//                b -- stored column wise, b = (b1_1,..., b_NRhs_1, b1_2, ...)
//                LDB = N_Rhs
// The result will be stored column wise
{
  int i,j,l,m, row;  // k, info, *ipiv;
  double pivot, tmp, *f, *frow,eps=0.0;

//  int ii, jj;
/*
  for(ii=0;ii<N_Eqn;ii++)
  {
    for(jj=0;jj<N_Eqn;jj++)
    {
      cout << "a(" << ii+1 << "," << jj+1 << ") = " << setw(8);
      cout << A(ii,jj) << ";" << endl;
    }
  }
  cout << endl;
  for(jj=0;jj<N_Rhs;jj++)
  {
    cout << jj+1 << ": ";
    for(ii=0;ii<N_Eqn;ii++)
      cout << setw(15) << b[ii*LDB+jj];
    cout << endl;
  }
*/

  // LU decomposition of matrix A with pivot search
  for(i=0;i<N_Eqn-1;i++)
  {
    pivot = 0;
    row = i;
    // find pivot
    for(l=i;l<N_Eqn;l++)
    {
      tmp = fabs(A(l,i));
      if(tmp > pivot)
      {
        pivot = tmp;
        row = l;
      } // endif
    } // endfor l
    if(pivot <= eps)
    { 
      OutPut("Error in solving multiple Systems " << __FILE__ << endl);
      OutPut("equation: " << i << endl);
      Error("Error in solving multiple Systems " << __FILE__ << endl);
      exit(4711);
    }
    // change rows if necessary
    if(i<row)
    {
      for(l=0;l<N_Eqn;l++)
      {
        pivot = A(i,l);
        A(i,l) = A(row,l);
        A(row, l) = pivot;
      }
      // row-th row of rhs
      frow = b+row*LDB;
      // i-th row of rhs
      f = b+i*LDB;          
      for(j=0;j<N_Rhs;j++)
      {
        tmp = f[j];
        f[j] = frow[j];
        frow[j] = tmp;
      }
    } // endif

    tmp = A(i,i);
    for(l=i+1;l<N_Eqn;l++)
      A(l,i) /= tmp;

    for(l=i+1;l<N_Eqn;l++)
    {
      tmp = A(i,l);
      for(m=i+1;m<N_Eqn;m++)
        A(m,l) -= A(m,i) * tmp;
    }

  } // endfor i

  frow = b;
   // solve left lower system
  for(i=0;i<N_Eqn;i++)
  { 
    // i-th row of rhs
    f = b+i*LDB;
    // for all rhs
    for(j=0;j<N_Rhs;j++)
    {    
      tmp = f[j];
      // for remaining rows
      for(l=i+1;l<N_Eqn;l++)
        {
          frow[l*LDB+j] -= A(l,i)*tmp; 
        }
    }
  }
  // solve right upper system
  for(i=N_Eqn-1;i>=0;i--)
  {
    // i-th row of rhs
    f = b+i*LDB;    
    // for all rhs
    for(j=0;j<N_Rhs;j++)
     {
       f[j] /= A(i,i);
       tmp = f[j];
       for(l=0;l<i;l++)
         {
           frow[l*LDB+j] -= A(l,i)*tmp;
         }          
    }
  }
}

void SolveMultipleSystems(double *a, double *b, int N_Eqn, 
                       int LDA, int LDB, int N_Rhs)
// Arguments:
//    a         double array which contains the matrix row wise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
//
// This works for a -- stored row wise
//                b -- stored row wise, b = (b1,b2,b3, ...)
//                LDA = LDB
// The result will be stored row wise
{
  int i,j,l,m, row;
  double pivot, tmp, *f;

//  int ii, jj;

/*
  for(ii=0;ii<N_Eqn;ii++)
  {
    cout << ii;
    for(jj=0;jj<N_Eqn;jj++)
      cout << setw(8) << A(ii,jj);
    cout << endl;
  }
  cout << endl;

  cout << "number of Rhs: " << N_Rhs << endl;
*/

  // for all columns
  for(i=0;i<N_Eqn-1;i++)
  {
    // compute pivot element
    pivot = 0;
    row = i;
    // check current column from diagonal element
    for(l=i;l<N_Eqn;l++)
    {
      // a_li 
      tmp = fabs(A(l,i));
      if(tmp > pivot)
      {
        pivot = tmp;
        row = l;
      } // endif
    } // endfor l
    if(pivot == 0.0)
    { 
      OutPut("Error in solving multiple Systems" << __FILE__ << endl);
      Error("Error in solving multiple Systems" << __FILE__ << endl);
      exit(4711);
    }
    // change rows i and 'row' if necessary
    if(i<row)
    {
      for(l=0;l<N_Eqn;l++)
      {
	// N_Eqn subsequent entries since a is stored row wise
        pivot = A(i,l);
        A(i,l) = A(row,l);
        A(row, l) = pivot;
      }
      // same for rhs
      for(j=0;j<N_Rhs;j++)
      { 
	// since rhs is stored row wise, find corresponding index in each rhs
        f = b+j*LDB;
        tmp = f[i];
        f[i] = f[row];
        f[row] = tmp;
      }
    } // endif

    // apply pivoting
    tmp = A(i,i);
    // current column
    for(l=i+1;l<N_Eqn;l++)
      A(l,i) /= tmp;
    // remainder of the matrix
    for(l=i+1;l<N_Eqn;l++)
    {
      tmp = A(i,l);
      for(m=i+1;m<N_Eqn;m++)
        A(m,l) -= A(m,i) * tmp;
    }
  } // endfor i

  for(i=0;i<N_Eqn;i++)
  {
    for(j=0;j<N_Rhs;j++)
    {
      f = b+j*LDB;
      tmp = f[i];
      for(l=i+1;l<N_Eqn;l++)
        f[l] -= A(l,i)*tmp;
    }
  }

  for(i=N_Eqn-1;i>=0;i--)
  {
    for(j=0;j<N_Rhs;j++)
    {
      f = b+j*LDB;
      f[i] /= A(i,i);
      tmp = f[i];
      for(l=0;l<i;l++)
        f[l] -= A(l,i)*tmp;
    }
  }
}

/** calculate the eigenvalue of the system using Lapack routines*/
void FindEigenValues(double *ap, int N_Eqn, char &COMPZ, double *d, double *z)
{
// Arguments:
//  ap         double precision array which contains the packed upper triangular matrix column wise
//              a[i,j] = a[i +(j-1)*j/2]
//  UPLO    (input) CHARACTER*1
//          = 'U':  Upper triangle of A is packed;
//          = 'L':  Lower triangle of A is packed./**/
//  N_Eqn    : order of the matrix ap 
//  COMPZ   (input) CHARACTER*1
//           = 'N':  Compute eigenvalues only.
//           = 'V':  Compute eigenvalues and eigenvectors of the original
//                   symmetric matrix.  On entry, Z must contain the
//                   orthogonal matrix used to reduce the original matrix
//                   to tridiagonal form.
//           = 'I':  Compute eigenvalues and eigenvectors of the
//                   tridiagonal matrix.  Z is initialized to the identity
//                   matrix.
//  d       (output) DOUBLE PRECISION array, dimension (LDZ, N_Eqn) eigen values
//  z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
//           COMPZ = 'V', Z contains the
//           orthonormal eigenvectors of the original symmetric matrix,
//           and if COMPZ = 'I', Z contains the orthonormal eigenvectors
//           of the symmetric tridiagonal matrix.
//           If COMPZ = 'N', then Z is not referenced.


 int info, i;
 char uplo='U';
// char compz='N';
//  double *d = new double[N_Eqn];
 double *e = new double[N_Eqn-1];
 double *tau = new double[N_Eqn-1];
 double *q = new double[N_Eqn*N_Eqn];
 double *work =  new double[2*N_Eqn-1];
 
 dsptrd_(&uplo, &N_Eqn,  ap, d, e, tau, &info);
  
 dopgtr_(&uplo, &N_Eqn,  ap, tau, q, &N_Eqn, work, &info);

 dsteqr_(&COMPZ, &N_Eqn, d, e, q, &N_Eqn, work, &info);

 
  for(i=0; i<N_Eqn; i++)
   cout<< " Eigen values " << d[i] << endl;
  
  for (i=0;i<N_Eqn*N_Eqn;i++)
    cout<< "q["<<i<<"]="<<q[i]<<endl;
//  
}

/* **************************************
 * The following functions are still needed somewhere in the code and where
 * moved here form deprecated "MComponents.C"
 * **************************************/

#ifdef __3D__
// determine L2 and H1 error
void L1Int3D(int N_Points, double *X, double *Y, double *Z,
             double *AbsDetjk,
             double *Weights, double hK,
             double **Der, double **Exact,
             double **coeffs, double *Loc)
{
  int i;
  double *deriv, w, t;

  Loc[0] = 0.0;
  Loc[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    w = Weights[i]*AbsDetjk[i];

    // int(f)
    t = deriv[0];
    Loc[0] += w*t;

    // int(1)
    Loc[1] += w;
  } // endfor i
}

// FIXME: This does not work in MPI case and should therefore be removed.
void IntoL20Vector3D(double *v, int Length, int order)
{
  double s;
  int i;

#ifdef _MPI
  int totallength;
  double temp;
#endif

  switch(order)
  {
    case -11:
      s=0;
      for(i=0;i<Length;i+=4)
        s += v[i];

#ifdef _MPI
      MPI_Allreduce(&Length, &totallength, 1, MPI_INT, MPI_SUM, TDatabase::ParamDB->Comm);
      MPI_Allreduce(&s, &temp, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);

      s = temp/(double)totallength;
#else
      s /= Length;
#endif

      s *= 4;

      for(i=0;i<Length;i+=4)
        v[i] -= s;
      break;

    case -12:
      s=0;
      for(i=0;i<Length;i+=10)
        s += v[i];

#ifdef _MPI
      MPI_Allreduce(&Length, &totallength, 1, MPI_INT, MPI_SUM, TDatabase::ParamDB->Comm);
      MPI_Allreduce(&s, &temp, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
      s = temp/(double)totallength;
#else
      s /= Length;
#endif

      s *= 10;

      for(i=0;i<Length;i+=10)
        v[i] -= s;

      break;

    default :
      s=0;
      for(i=0;i<Length;i++)
        s += v[i];

#ifdef _MPI
      MPI_Allreduce(&Length, &totallength, 1, MPI_INT, MPI_SUM, TDatabase::ParamDB->Comm);
      MPI_Allreduce(&s, &temp, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
      s = temp/(double)totallength;
#else
      s /= Length;
#endif

      for(i=0;i<Length;i++)
        v[i] -= s;
      break;

  }
} // IntoL20


/** project vector v into L20 */
void IntoL20FEFunction3D(double *v, int Length, TFESpace3D *FESpace)
{
  double s;
  int i,j,l,N_LocalUsedElements;
  int N_Cells, N_Points, N_;
  int *N_BaseFunct; //Used[N_FEs3D];
//  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
//  TFE3D *ele;
  double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
//  RefTrans3D RefTrans;
//  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D], der[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF; //ActiveBound, DirichletBound, end, last, number;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool SecondDer[1];
  double error0, error1;
  double *interpol;
  TNodalFunctional3D *nf;
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];

// #ifdef _MPI
//   int rank, ID;
//   MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
// #endif

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  for(i=0;i<MaxN_QuadPoints_3D;i++)
    Derivatives[i] = der+i;

  SecondDer[0] = false;

  GlobalNumbers = FESpace->GetGlobalNumbers();
  BeginIndex = FESpace->GetBeginIndex();

  interpol = new double[Length];

  error0 = 0.0;
  error1 = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = FESpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 1;
    CurrentElement = FESpace->GetFE3D(i, cell);
    LocalUsedElements[0] = CurrentElement;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements,
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    // calculate all needed derivatives of this FE function
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = v[DOF[l]];

    OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
      {
        value += FEFunctValues[l] * Orig[l];
      } // endfor l
      Derivatives[j][0] = value;
    } // endfor j

    L1Int3D(N_Points, X, Y, Z, AbsDetjk, weights, hK, Derivatives,
              ExactVal, AuxArray, LocError);

    error0 += LocError[0];
    error1 += LocError[1];

  } // endfor i

#ifdef _MPI
  double temp;
  temp = error0;
  MPI_Allreduce(&temp, &error0, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
  temp = error1;
  MPI_Allreduce(&temp, &error1, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
#endif

  s = error0/error1;

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    CurrentElement = FESpace->GetFE3D(i, cell);
    N_ = N_BaseFunct[CurrentElement];
    nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(CurrentElement);
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    for(j=0;j<N_Points;j++)
      PointValues[j] = s;
    nf->GetAllFunctionals(Coll, cell, PointValues, FunctionalValues);
    DOF = GlobalNumbers+BeginIndex[i];
    for(j=0;j<N_;j++)
      interpol[DOF[j]] = FunctionalValues[j];
  } // endfor i

  for(i=0;i<Length;i++)
    v[i] -= interpol[i];

  delete interpol;
} // IntoL20Function

/** project function v into L20 */
void IntoL20FEFunction3D(double *v, int Length, TFESpace3D *FESpace,
                       int velocity_space, int pressure_space)
{
  double s;
  int i,j,l,N_LocalUsedElements;
  int N_Cells, N_Points, N_;
  int *N_BaseFunct;
//  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
//  TFE3D *ele;
  double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
//  RefTrans3D RefTrans;
//  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D], der[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF;// ActiveBound, DirichletBound, end, last, number;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool SecondDer[1];
  double error0, error1;
  double *interpol;
  TNodalFunctional3D *nf;
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];

#ifdef _MPI
//   MPI_Comm Comm;
//   int rank, ID;
  double temp;
//
//   Comm = TDatabase::ParamDB->Comm;
//   MPI_Comm_rank(Comm, &rank);
#endif

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  for(i=0;i<MaxN_QuadPoints_3D;i++)
    Derivatives[i] = der+i;

  SecondDer[0] = false;

  GlobalNumbers = FESpace->GetGlobalNumbers();
  BeginIndex = FESpace->GetBeginIndex();

  interpol = new double[Length];

  error0 = 0.0;
  error1 = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = FESpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

// #ifdef _MPI
//     ID = cell->GetSubDomainNo();
//     if(ID!=rank)
//       continue; // this cell not belongs to this processor
// #endif

    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 1;
    CurrentElement = FESpace->GetFE3D(i, cell);
    LocalUsedElements[0] = CurrentElement;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements,
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    // calculate all needed derivatives of this FE function
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = v[DOF[l]];

    OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
      {
        value += FEFunctValues[l] * Orig[l];
      } // endfor l
      Derivatives[j][0] = value;
    } // endfor j

    L1Int3D(N_Points, X, Y, Z, AbsDetjk, weights, hK, Derivatives,
              ExactVal, AuxArray, LocError);

    error0 += LocError[0];
    error1 += LocError[1];

  } // endfor i

#ifdef _MPI
  temp = error0;
  MPI_Allreduce(&temp, &error0, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
  temp = error1;
  MPI_Allreduce(&temp, &error1, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
#endif


  s = error0/error1;

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    CurrentElement = FESpace->GetFE3D(i, cell);
    N_ = N_BaseFunct[CurrentElement];
    nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(CurrentElement);
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    for(j=0;j<N_Points;j++)
      PointValues[j] = s;
    nf->GetAllFunctionals(Coll, cell, PointValues, FunctionalValues);
    DOF = GlobalNumbers+BeginIndex[i];
    for(j=0;j<N_;j++)
      interpol[DOF[j]] = FunctionalValues[j];
  } // endfor i

  for(i=0;i<Length;i++)
    v[i] -= interpol[i];

  delete [] interpol;
} // IntoL20Function

#endif

#ifdef __2D__
// determine L2 and H1 error
void L1Int(int N_Points, double *X, double *Y, double *AbsDetjk,
                double *Weights, double hK,
                double **Der, double **Exact,
                double **coeffs, double *Loc)
{
  int i;
  double *deriv, w, t;

  Loc[0] = 0.0;
  Loc[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    w = Weights[i]*AbsDetjk[i];

    // int(f)
    t = deriv[0];
    Loc[0] += w*t;

    // int(1)
    Loc[1] += w;
  } // endfor i
}


/** project vector v into L20 */
void IntoL20Vector2D(double *v, int Length, int order)
{
  double s;
  int i;

  switch(order)
  {
    case -11:
      s=0;
      for(i=0;i<Length;i+=3)
        s += v[i];

      s /= Length;
      s *= 3;

      for(i=0;i<Length;i+=3)
        v[i] -= s;
      break;

    case -12:
      s=0;
      for(i=0;i<Length;i+=6)
        s += v[i];

      s /= Length;
      s *= 6;

      for(i=0;i<Length;i+=6)
        v[i] -= s;

      break;

    case -13:
      s=0;
      for(i=0;i<Length;i+=10)
        s += v[i];

      s /= Length;
      s *= 10;

      for(i=0;i<Length;i+=10)
        v[i] -= s;

      break;

    case -14:
      s=0;
      for(i=0;i<Length;i+=15)
        s += v[i];

      s /= Length;
      s *= 15;

      for(i=0;i<Length;i+=15)
        v[i] -= s;

      break;

    case 201:
      break;

    default :
      s=0;
      for(i=0;i<Length;i++)
        s += v[i];

      s /= Length;
//      OutPut("VECTOR " << Length << " " << s << endl);
      for(i=0;i<Length;i++)
        v[i] -= s;
      break;

  }
} // IntoL20

/** project function v into L20 */
void IntoL20FEFunction_OLD(double *v, int Length, TFESpace2D *FESpace,
                       int velocity_space, int pressure_space)
{
  double s;
  int i,j,l, N_LocalUsedElements;
  int N_Cells, N_Points, N_;
  int *N_BaseFunct;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Derivatives[MaxN_QuadPoints_2D], der[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, number;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool SecondDer[1];
  double error0, error1;

  if (pressure_space==-4711)
  {
    switch (velocity_space)
    {
        case 1:
          number = 0;
          break;
        case 2:
        case 3:
        case 4:
        case 5:
          number = velocity_space;
          break;
        case 12:
        case 13:
        case 14:
          number = -velocity_space+1;
          break;
        case 22:
        case 23:
          number = -velocity_space+11;
          break;

        case -1:
          number = -TDatabase::ParamDB->VELOCITY_SPACE -1;
        break;

        case -2:
        case -3:
        case -4:
        case -5:
          number = TDatabase::ParamDB->VELOCITY_SPACE - 9;
        break;
    }
  }
  else
    number = pressure_space;

//  OutPut("NUMBER " << number << " " << Length << " ");

  switch(number)
  {
    // pw constant
    case 0:
      s=0;
      for(i=0;i<Length;i++)
        s += v[i];
//      OutPut(s << endl);
      s /= Length;
      for(i=0;i<Length;i++)
        v[i] -= s;
      break;

    // pw linear discontinuous
    case -11:
      s=0;
      for(i=0;i<Length;i+=3)
        s += v[i];

      s /= Length;
      s *= 3;

      for(i=0;i<Length;i+=3)
        v[i] -= s;
      break;

    // pw quadratics discontinuous
    case -12:
      s=0;
      for(i=0;i<Length;i+=6)
        s += v[i];

      s /= Length;
      s *= 6;

      for(i=0;i<Length;i+=6)
        v[i] -= s;
      break;

    // pw cubics discontinuous
    case -13:
      s=0;
      for(i=0;i<Length;i+=10)
        s += v[i];

      s /= Length;
      s *= 10;

      for(i=0;i<Length;i+=10)
        v[i] -= s;

      break;

    // discontinuous, pw P4
    case -14:
      s=0;
      for(i=0;i<Length;i+=15)
        s += v[i];

      s /= Length;
      s *= 15;

      for(i=0;i<Length;i+=15)
        v[i] -= s;

      break;

    // pw polynomials, continuous
    // with and without bubbles
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 12:
    case 13:
    case 14:
    case 22:
    case 23:
      BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
      N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

      for(i=0;i<MaxN_QuadPoints_2D;i++)
        Derivatives[i] = der+i;

      SecondDer[0] = false;

      GlobalNumbers = FESpace->GetGlobalNumbers();
      BeginIndex = FESpace->GetBeginIndex();

      error0 = 0.0;
      error1 = 0.0;

    // ########################################################################
    // loop over all cells
    // ########################################################################
      Coll = FESpace->GetCollection(); // all spaces use same Coll
      N_Cells = Coll->GetN_Cells();
      for(i=0;i<N_Cells;i++)
      {
        cell = Coll->GetCell(i);

        hK = cell->GetDiameter();

        // ####################################################################
        // find local used elements on this cell
        // ####################################################################
        N_LocalUsedElements = 1;
        CurrentElement = FESpace->GetFE2D(i, cell);
        LocalUsedElements[0] = CurrentElement;

        // ####################################################################
        // calculate values on original element
        // ####################################################################
        TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
                             Coll, cell, SecondDer,
                             N_Points, xi, eta, weights, X, Y, AbsDetjk);

        // calculate all needed derivatives of this FE function
        BaseFunct = BaseFuncts[CurrentElement];
        N_ = N_BaseFunct[CurrentElement];

        DOF = GlobalNumbers + BeginIndex[i];
        for(l=0;l<N_;l++)
          FEFunctValues[l] = v[DOF[l]];

        OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
        for(j=0;j<N_Points;j++)
        {
          Orig = OrigFEValues[j];
          value = 0;
          for(l=0;l<N_;l++)
          {
            value += FEFunctValues[l] * Orig[l];
          } // endfor l
          Derivatives[j][0] = value;
        } // endfor j

        L1Int(N_Points, X, Y, AbsDetjk, weights, hK, Derivatives,
                  ExactVal, AuxArray, LocError);

        error0 += LocError[0];
        error1 += LocError[1];

      } // endfor i

      s = error0/error1;

      for(i=0;i<Length;i++)
        v[i] -= s;

    break;

    default:
      cout << "The L^2_0 projection does not ";
      cout << "work for this type of elements" << endl;
      exit(-1);
  }
} // IntoL20Function

/** project function v into L20 */
void IntoL20FEFunction(double *v, int Length, const TFESpace2D *FESpace,
                       int velocity_space, int pressure_space
#ifdef _MPI
                        , MPI_Comm comm
#endif
                      )
{
  double s;
  int i,j,l, N_LocalUsedElements;
  int N_Cells, N_Points, N_;
  int *N_BaseFunct;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Derivatives[MaxN_QuadPoints_2D], der[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool SecondDer[1];
  double error0, error1;
  double *interpol;
  TNodalFunctional2D *nf;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];

  #ifdef _MPI
  int ID, rank;
  MPI_Comm_rank(comm, &rank);
  #endif

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  for(i=0;i<MaxN_QuadPoints_2D;i++)
    Derivatives[i] = der+i;

  SecondDer[0] = false;

  GlobalNumbers = FESpace->GetGlobalNumbers();
  BeginIndex = FESpace->GetBeginIndex();

  interpol = new double[Length];

  error0 = 0.0;
  error1 = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = FESpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    #ifdef _MPI
    ID = cell->GetSubDomainNo();
    if(ID!=rank)
      continue; // this cell is a halo cell
    #endif
    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 1;
    CurrentElement = FESpace->GetFE2D(i, cell);
    LocalUsedElements[0] = CurrentElement;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // calculate all needed derivatives of this FE function
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = v[DOF[l]];

    OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
      {
        value += FEFunctValues[l] * Orig[l];
      } // endfor l
      Derivatives[j][0] = value;
    } // endfor j

    L1Int(N_Points, X, Y, AbsDetjk, weights, hK, Derivatives,
              ExactVal, AuxArray, LocError);

    error0 += LocError[0];
    error1 += LocError[1];

  } // endfor i

  #ifdef _MPI
  temp = error0;
  MPI_Allreduce(&temp, &error0, 1, MPI_DOUBLE, MPI_SUM, comm);
  temp = error1;
  MPI_Allreduce(&temp, &error1, 1, MPI_DOUBLE, MPI_SUM, comm);
  #endif

  s = error0/error1;

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    CurrentElement = FESpace->GetFE2D(i, cell);
    N_ = N_BaseFunct[CurrentElement];
    nf = TFEDatabase2D::GetNodalFunctional2DFromFE2D(CurrentElement);
    nf->GetPointsForAll(N_Points, xi, eta);
    for(j=0;j<N_Points;j++)
      PointValues[j] = s;
    nf->GetAllFunctionals(Coll, cell, PointValues, FunctionalValues);
    DOF = GlobalNumbers+BeginIndex[i];
    for(j=0;j<N_;j++)
      interpol[DOF[j]] = FunctionalValues[j];
  } // endfor i

  for(i=0;i<Length;i++)
    v[i] -= interpol[i];

  delete [] interpol;
} // IntoL20Function
#endif




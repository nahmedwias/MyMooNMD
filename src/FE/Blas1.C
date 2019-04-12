// =======================================================================
// Blas1.C
//
// Purpose:     basic routines for linear algebra
//
// Author:      Gunar Matthies          17.01.2003
//
// =======================================================================

#include <string.h>
#include <math.h>

/** return inner product (x,y) */
double Ddot(int n, const double *x, const double *y)
{
  double r;
  int i;
  const double *a, *b;

  r = 0.0;
  a = x;
  b = y;
  for(i=0;i<n;i++)
  {
    r += *a * *b;
    a++;
    b++;
  }

  return r;
}

/** y := alpha*x + y */
void Daxpy(int n, double alpha, const double *x, double *y)
{
  int i;
  const double *a;
  double *b;
  double scal;

  a = x;
  b = y;
  scal = alpha;
  for(i=0;i<n;i++)
  {
    *b += scal * *a;
    a++;
    b++;
  }

}

/** z := alpha*x + beta*y */
void Dsum(int n, double alpha, double beta, double *x, double *y, double *z)
{
  int i;
  double *a, *b, *c;
  double scal1, scal2;

  a = x;
  b = y;
  c = z;
  scal1 = alpha;
  scal2 = beta;
  for(i=0;i<n;i++)
  {
    *c = scal1 * *a + scal2 * *b;
    a++;
    b++;
    c++;
  }

}

/** b := a */
void Dcopy(int n, const double *a, double *b)
{
  memcpy(b, a, n*sizeof(double));
}

/** x := alpha*x */
void Dscal(int n, double alpha, double *x)
{
  int i;
  double scal;
  double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a *= scal;
    a++;
  }
}

/** return Euclidian norm of x */
double Dnorm(int n, const double *x)
{
  double r;
  int i;
  const double *a;
  double z;

  a = x;
  r = 0.0;
  for(i=0; i<n; i++)
  {
    z = *a; 
    r += z*z;
    a++;
  }

  return sqrt(r);
}



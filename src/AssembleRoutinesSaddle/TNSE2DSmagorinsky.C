#include <TNSE2DSmagorinsky.h>
#include <CommonRoutineTNSE2D.h>
#include <MooNMD_Io.h>

// ======================================================================
// Type 1, Smagorinsky
// the nonlinear viscosity is treated implicitly
// ======================================================================
void TimeNSType1Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  
  mu = turbulentviscosity(hK,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}



// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void TimeNSType2Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  mu = turbulentviscosity(hK,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}



// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// Type 2, Smagorinsky, only nonlinear part
// ======================================================================
void TimeNSType1_2NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2, mu;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  mu = turbulentviscosity(hK,&param[2],&param[0],&param[6]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}



// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; //  **MatrixA21, **MatrixA12;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  mu = turbulentviscosity(hK,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}



// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row; // *Matrix21Row, *Matrix12Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  mu = turbulentviscosity(hK,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v), only nonlinear part
// Type 4, Smagorinsky, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  mu = turbulentviscosity(hK,&param[2],&param[0],&param[6]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, viscosity;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  mu = turbulentviscosity(hK,&param[2],&param[0],&param[0]);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}



// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType4SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  
  
  mu = turbulentviscosity(hK,&param[2],&param[0],&param[0]);
  mu = mu/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(2*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0, viscosity;
  double u1, u2, mu;
  // cout << "Sma" << endl;
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  mu = turbulentviscosity(hK,&param[2],&param[0],&param[6]);
  mu = mu/2.0;
  viscosity = mu+c0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


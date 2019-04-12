#include <TNSE3DSmagorinsky.h>
#include <CommonRoutineTNSE3D.h>

// ======================================================================
// Type 1, Smagorinsky
// the nonlinear viscosity is treated implicitly
// ======================================================================
void TimeNSType1Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void TimeNSType2Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];
  MatrixB1T = LocMatrices[5];
  MatrixB2T = LocMatrices[6];
  MatrixB3T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);

      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, mu;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);

      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, viscosity;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);

      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);

      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);

      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T  = LocMatrices[15];
  MatrixB2T  = LocMatrices[16];
  MatrixB3T  = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);

      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType4SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, viscosity;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T = LocMatrices[15];
  MatrixB2T = LocMatrices[16];
  MatrixB3T = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);

      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// Type 2, Smagorinsky, only nonlinear part
// ======================================================================
void TimeNSType1_2NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0, mu;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);


  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);

      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v), only nonlinear part
// Type 4, Smagorinsky, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0, mu;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);

      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1,val2, val3, val4;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0, viscosity, dummy_param[6] = {0, 0, 0, 0, 0, 0};
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &dummy_param[0];

  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  viscosity = Mult*(mu/2.0+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = viscosity*Orig0[i];
    test010 = viscosity*Orig1[i];
    test001 = viscosity*Orig2[i];
    test000 = Mult*Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;
      val1 += val2+val3+val4;
      Matrix11Row[j] += val2+val1;
      Matrix12Row[j] += test010*ansatz100;
      Matrix13Row[j] += test001*ansatz100;
      Matrix21Row[j] += test100*ansatz010;
      Matrix22Row[j] += val3+val1;
      Matrix23Row[j] += test001*ansatz010;
      Matrix31Row[j] += test100*ansatz001;
      Matrix32Row[j] += test010*ansatz001;
      Matrix33Row[j] += val4+val1;
    } // endfor j
  } // endfor i
}


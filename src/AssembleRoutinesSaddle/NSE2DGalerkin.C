#include <../../include/AssembleRoutinesSaddle/NSE2DGalerkin.h>
#include <Database.h>

#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!
Hotfixglobal_AssembleNSE assemble_nse(Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION);

// =================================================================================

void NSType1Galerkin(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

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

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);

      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
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

// =================================================================================
void NSType2Galerkin(double Mult, double *coeff, double *param, double hK, 
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

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

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
     //HOTFIX: Check the documentation!
     if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
       val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
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

// =================================================================================
void NSType1_2NLGalerkin(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                         double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

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

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}
// =================================================================================
void NSType3Galerkin(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row; // *Matrix21Row, *Matrix12Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

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

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

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

// =================================================================================
void NSType4Galerkin(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;  // *Matrix21Row, *Matrix12Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1,u2;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  int offset = TDatabase::ParamDB->NSTYPE == 14 ? 1 : 0;
  MatrixB1 = LocMatrices[4+offset];
  MatrixB2 = LocMatrices[5+offset];
  MatrixB1T = LocMatrices[6+offset];
  MatrixB2T = LocMatrices[7+offset];

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

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

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
// =================================================================================
void NSType3_4NLGalerkin(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                         double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

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

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}
// =================================================================================

void NSType3GalerkinDD(double Mult, double *coeff, double *param, double hK,
                       double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                       double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

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

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

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

// =================================================================================
void NSType4GalerkinDD(double Mult, double *coeff, double *param, double hK,
                       double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                       double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  int offset = TDatabase::ParamDB->NSTYPE == 14 ? 1 : 0;
  MatrixB1 = LocMatrices[4+offset];
  MatrixB2 = LocMatrices[5+offset];
  MatrixB1T = LocMatrices[6+offset];
  MatrixB2T = LocMatrices[7+offset];

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

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

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

// =================================================================================
void NSType3_4NLGalerkinDD(double Mult, double *coeff, double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                           double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

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

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}

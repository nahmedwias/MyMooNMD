#include <NSE3DGalerkin.h>
#include <Hotfixglobal_AssembleNSE.h>
// this file gets only compiled #ifdef __3D__ - this is why it is no problem
// to define the global variable assemble_nse here AND in NSE2D_FixPo (which
// gets compiled only #ifdef __2D__)
Hotfixglobal_AssembleNSE assemble_nse(Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION);



void NSLaplaceGradGradSingle(double Mult, double *coeff, double *param,
                             double hK, double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs)
{
  double ansatz100, ansatz010, ansatz001, test100, test010, test001;
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = OrigValues[4];
  
  double nu = coeff[0]; // = 1/reynolds_number
  
  for(int i = 0; i < N_U; i++)
  {
    test100 = u_x[i];
    test010 = u_y[i];
    test001 = u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz100 = u_x[j];
      ansatz010 = u_y[j];
      ansatz001 = u_z[j];
      MatrixA[i][j] += Mult * nu*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
    }
  }
}

void NSLaplaceGradGrad(double Mult, double *coeff, double *param,
                       double hK, double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double ansatz100, ansatz010, ansatz001, test100, test010, test001;
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[4];
  double ** MatrixA33 = LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = OrigValues[4];
  double nu = coeff[0]; // = 1/reynolds_number
  
  for(int i = 0; i < N_U; i++)
  {
    test100 = u_x[i];
    test010 = u_y[i];
    test001 = u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz100 = u_x[j];
      ansatz010 = u_y[j];
      ansatz001 = u_z[j];
      double val = Mult * nu*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      MatrixA33[i][j] += val;
    }
  }
}

void NSLaplaceDeformation(double Mult, double *coeff, double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs)
{
  double ansatz100, ansatz010, ansatz001;
  double test100, test010, test001;
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA12 = LocMatrices[1];
  double ** MatrixA13 = LocMatrices[2];
  double ** MatrixA21 = LocMatrices[3];
  double ** MatrixA22 = LocMatrices[4];
  double ** MatrixA23 = LocMatrices[5];
  double ** MatrixA31 = LocMatrices[6];
  double ** MatrixA32 = LocMatrices[7];
  double ** MatrixA33 = LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = OrigValues[4];
  double nu = coeff[0]; 
  for(int i = 0; i < N_U; i++)
  {
    test100 = u_x[i];
    test010 = u_y[i];
    test001 = u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz100 = u_x[j];
      ansatz010 = u_y[j];
      ansatz001 = u_z[j];
      MatrixA11[i][j] += Mult * 2*nu*(test100*ansatz100+0.5*test010*ansatz010
                                      +0.5*test001*ansatz001);
      MatrixA12[i][j] += Mult * nu*(test010*ansatz100);
      MatrixA13[i][j] += Mult * nu*(test001*ansatz100);
      MatrixA21[i][j] += Mult * nu*(test100*ansatz010);
      MatrixA22[i][j] += Mult * 2*nu*(0.5*test100*ansatz100+test010*ansatz010
                                      +0.5*test001*ansatz001);
      MatrixA23[i][j] += Mult * nu*(test001*ansatz010);
      MatrixA31[i][j] += Mult * nu*(test100*ansatz001);
      MatrixA32[i][j] += Mult * nu*(test010*ansatz001);
      MatrixA33[i][j] += Mult * 2*nu*(0.5*test100*ansatz100+0.5*test010*ansatz010
                                      +test001*ansatz001);
    } 
  } 
}

void NSDivergenceBlocks(double Mult, double *coeff, double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs)
{
  double ansatz100, ansatz010, ansatz001, test000;
  double** MatrixB1 = LocMatrices[10];
  double** MatrixB2 = LocMatrices[11];
  double** MatrixB3 = LocMatrices[12];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * p = OrigValues[1];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = OrigValues[4];
  for(int i=0;i<N_P;i++)
  {
    test000 = p[i];
    for(int j=0;j<N_U;j++)
    {
      ansatz100 = u_x[j];
      ansatz010 = u_y[j];
      ansatz001 = u_z[j];

      MatrixB1[i][j] += -Mult*test000*ansatz100;
      MatrixB2[i][j] += -Mult*test000*ansatz010;
      MatrixB3[i][j] += -Mult*test000*ansatz001;
    }
  }
}


void NSGradientBlocks(double Mult, double *coeff, double *param, double hK,
                      double **OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs)
{
  double test100, test010, test001, ansatz000;
  double** MatrixB1T = LocMatrices[13];
  double** MatrixB2T = LocMatrices[14];
  double** MatrixB3T = LocMatrices[15];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * p = OrigValues[1];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = OrigValues[4];
  for(int i=0;i<N_U;i++)
  {
    test100 = u_x[i];
    test010 = u_y[i];
    test001 = u_z[i];
    for(int j=0;j<N_P;j++)
    {
      ansatz000 = p[j];
      MatrixB1T[i][j] += -Mult*ansatz000*test100;
      MatrixB2T[i][j] += -Mult*ansatz000*test010;
      MatrixB3T[i][j] += -Mult*ansatz000*test001;
    }
  }
}

void NSRightHandSide(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs)
{
  double test000;
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs3 = LocRhs[2];
  double * Rhs4 = LocRhs[3];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * u = OrigValues[0];
  double * p = OrigValues[1];
  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = coeff[3];
  double g = coeff[4]; // divergence
  for(int i = 0; i < N_U; i++)
  {
    test000 = u[i];
    Rhs1[i] += Mult*test000*f1;
    Rhs2[i] += Mult*test000*f2;
    Rhs3[i] += Mult*test000*f3;
  }
  for(int i = 0; i < N_P; i++)
  {
    test000 = p[i];
    Rhs4[i] += -Mult*test000*g;
  }
}


void NSNonlinearTermSingle(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                         double **LocRhs)
{
  double test000, ansatz100, ansatz010, ansatz001;
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  double * u = OrigValues[0]; 
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = OrigValues[4];
  double u1 = param[0];
  double u2 = param[1];
  double u3 = param[2];
  for(int i = 0; i < N_U; i++)
  {
    test000 = u[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz100 = u_x[j];
      ansatz010 = u_y[j];
      ansatz001 = u_z[j];
      MatrixA[i][j] += Mult * (u1*ansatz100 + u2*ansatz010 + u3*ansatz001) * test000;
    }
  }
}

void NSNonlinearTerm(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                         double **LocRhs)
{
  double test000, ansatz100, ansatz010, ansatz001;
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[4];
  double ** MatrixA33 = LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  double * u = OrigValues[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = OrigValues[4];
  double u1 = param[0];
  double u2 = param[1];
  double u3 = param[2];
  for(int i = 0; i < N_U; i++)
  {
    test000 = u[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz100 = u_x[j];
      ansatz010 = u_y[j];
      ansatz001 = u_z[j];
      MatrixA11[i][j] += Mult * (u1*ansatz100 + u2*ansatz010 + u3*ansatz001)*test000;
      MatrixA22[i][j] += Mult * (u1*ansatz100 + u2*ansatz010 + u3*ansatz001)*test000;
      MatrixA33[i][j] += Mult * (u1*ansatz100 + u2*ansatz010 + u3*ansatz001)*test000;
    }
  }
}
void NSType1Galerkin(double Mult, double *coeff, double *param, double hK, 
                     double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs)
{
  
  double** MatrixA = LocMatrices[0];
  double** MatrixB1 = LocMatrices[1];
  double** MatrixB2 = LocMatrices[2];
  double** MatrixB3 = LocMatrices[3];

  double* Rhs1 = LocRhs[0];
  double* Rhs2 = LocRhs[1];
  double* Rhs3 = LocRhs[2];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double* Orig0 = OrigValues[0]; // u_x
  double* Orig1 = OrigValues[1]; // u_y
  double* Orig2 = OrigValues[2]; // u_y
  double* Orig3 = OrigValues[3]; // u
  double* Orig4 = OrigValues[4]; // p

  double c0 = coeff[0]; // nu
  double c1 = coeff[1]; // f1
  double c2 = coeff[2]; // f2
  double c3 = coeff[3]; // f3

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old
  
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double val, val1;
  
    
  for(int i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(int i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;
    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

void NSType2Galerkin(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs)
{  
  double** MatrixA = LocMatrices[0];
  double** MatrixB1 = LocMatrices[1];
  double** MatrixB2 = LocMatrices[2];
  double** MatrixB3 = LocMatrices[3];
  double** MatrixB1T = LocMatrices[4];
  double** MatrixB2T = LocMatrices[5];
  double** MatrixB3T = LocMatrices[6];

  double* Rhs1 = LocRhs[0];
  double* Rhs2 = LocRhs[1];
  double* Rhs3 = LocRhs[2];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double* Orig0 = OrigValues[0]; // u_x
  double* Orig1 = OrigValues[1]; // u_y
  double* Orig2 = OrigValues[2]; // u_z
  double* Orig3 = OrigValues[3]; // u
  double* Orig4 = OrigValues[4]; // p

  double c0 = coeff[0]; // nu
  double c1 = coeff[1]; // f1
  double c2 = coeff[2]; // f2
  double c3 = coeff[3]; // f3

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old

  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double val, val1;
  for(int i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      
      MatrixRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(int j=0;j<N_P;j++)
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

  for(int i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

void NSType1_2NLGalerkin(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs)
{
  double **MatrixA, val, *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3;
  
  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
    
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
     
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

//===========================================================================================
void NSType3Galerkin(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs)
{  
  double** MatrixA11 = LocMatrices[0];
  double** MatrixA22 = LocMatrices[4];
  double** MatrixA33 = LocMatrices[8];
  double** MatrixB1  = LocMatrices[9];
  double** MatrixB2  = LocMatrices[10];
  double** MatrixB3  = LocMatrices[11];

  double* Rhs1 = LocRhs[0];
  double* Rhs2 = LocRhs[1];
  double* Rhs3 = LocRhs[2];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double* Orig0 = OrigValues[0]; // u_x
  double* Orig1 = OrigValues[1]; // u_y
  double* Orig2 = OrigValues[2]; // u_y
  double* Orig3 = OrigValues[3]; // u
  double* Orig4 = OrigValues[4]; // p

  double c0 = coeff[0]; // nu
  double c1 = coeff[1]; // f1
  double c2 = coeff[2]; // f2
  double c3 = coeff[3]; // f3

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old
  
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double val, val1;

  for(int i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(int i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

void NSType4Galerkin(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                       double **LocRhs)
{
  double** MatrixA11 = LocMatrices[0];
//  double** MatrixA12 = LocMatrices[1];
//  double** MatrixA13 = LocMatrices[2];
//  double** MatrixA21 = LocMatrices[3];
  double** MatrixA22 = LocMatrices[4];
//  double** MatrixA23 = LocMatrices[5];
//  double** MatrixA31 = LocMatrices[6];
//  double** MatrixA32 = LocMatrices[7];
  double** MatrixA33 = LocMatrices[8];
  double** MatrixB1  = LocMatrices[9];
  double** MatrixB2  = LocMatrices[10];
  double** MatrixB3  = LocMatrices[11];
  double** MatrixB1T = LocMatrices[12];
  double** MatrixB2T = LocMatrices[13];
  double** MatrixB3T = LocMatrices[14];

  double* Rhs1 = LocRhs[0];
  double* Rhs2 = LocRhs[1];
  double* Rhs3 = LocRhs[2];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double* Orig0 = OrigValues[0]; // u_x
  double* Orig1 = OrigValues[1]; // u_y
  double* Orig2 = OrigValues[2]; // u_y
  double* Orig3 = OrigValues[3]; // u
  double* Orig4 = OrigValues[4]; // p

  double c0 = coeff[0]; // nu
  double c1 = coeff[1]; // f1
  double c2 = coeff[2]; // f2
  double c3 = coeff[3]; // f3

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old
  
  double *Matrix11Row;  // *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row;  // *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double val, val1;

  for(int i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];  // not used in the rest of this function
//    Matrix13Row = MatrixA13[i];  // not used in the rest of this function
//    Matrix21Row = MatrixA21[i];  // not used in the rest of this function
    Matrix22Row = MatrixA22[i];
//    Matrix23Row = MatrixA23[i];  // not used in the rest of this function
//    Matrix31Row = MatrixA31[i];  // not used in the rest of this function
//    Matrix32Row = MatrixA32[i];  // not used in the rest of this function
    Matrix33Row = MatrixA33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(int j=0;j<N_P;j++)
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

  for(int i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];
    
    test000 = Orig4[i];
    val1 = Mult*test000;
    
    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;
      
      val = -val1*ansatz010;
      MatrixRow2[j] += val;
      
      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
    
  } // endfor i

}


void NSType3_4NLGalerkin(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
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
      
      val  = c0*(test100*ansatz100+test010*ansatz010
        + test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

    } // endfor j
  } // endfor i
}

void NSType3GalerkinDD(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                       double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixB1  = LocMatrices[9];
  MatrixB2  = LocMatrices[10];
  MatrixB3  = LocMatrices[11];

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

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
    
      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

void NSType4GalerkinDD(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixB1  = LocMatrices[9];
  MatrixB2  = LocMatrices[10];
  MatrixB3  = LocMatrices[11];
  MatrixB1T = LocMatrices[12];
  MatrixB2T = LocMatrices[13];
  MatrixB3T = LocMatrices[14];

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

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;
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
    val1 = Mult*test000;
    
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;
      
      val = -val1*ansatz010;
      MatrixRow2[j] += val;
      
      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
    
  } // endfor i
}

void NSType3_4NLGalerkinDD(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs)
{
  double **MatrixA11, **MatrixA22,  **MatrixA33;
  double val, val1, val2, val3, val4;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, u3;
 
  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // 2*nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

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
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;

      val  = c0*(2*val2+val3+val4)+val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(val2+2*val3+val4)+val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(val2+val3+2*val4)+val1;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ========================================================================
void NSParamsVelo3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old  
  out[2] = in[5]; // u3old  
}

void mat_coriolis(double Mult, double *coeff, double *param, double hK, 
                  double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                  double **LocRhs)
{
  double** MatrixA12 = LocMatrices[1];
  double** MatrixA13 = LocMatrices[2];
  double** MatrixA21 = LocMatrices[3];
  double** MatrixA23 = LocMatrices[5];
  double** MatrixA31 = LocMatrices[6];
  double** MatrixA32 = LocMatrices[7];

  int N_U = N_BaseFuncts[0];

  double* Orig3 = OrigValues[3]; // u

  double Omega1 = coeff[5];
  double Omega2 = coeff[6];
  double Omega3 = coeff[7];

  for(int i=0;i<N_U;i++)
  {
    const double test000 = Orig3[i];
    for(int j=0;j<N_U;j++)
    {
      const double ansatz000 = Orig3[j];
      const double val1 = Mult * ansatz000 * test000;
      MatrixA12[i][j] -= Omega3 * val1;
      MatrixA13[i][j] += Omega2 * val1;
      MatrixA21[i][j] += Omega3 * val1;
      MatrixA23[i][j] -= Omega1 * val1;
      MatrixA31[i][j] -= Omega2 * val1;
      MatrixA32[i][j] += Omega1 * val1;
    } // endfor j
  } // endfor i
}

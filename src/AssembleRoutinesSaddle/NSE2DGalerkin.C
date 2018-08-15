#include <../../include/AssembleRoutinesSaddle/NSE2DGalerkin.h>
#include <Database.h>

#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!
Hotfixglobal_AssembleNSE assemble_nse(Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION);

// =================================================================================
void NSLaplaceGradGradSingle(double Mult, double *coeff, double *param,
                             double hK, double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs)
{
  double ansatz10, ansatz01, test10, test01;
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double nu = coeff[0]; // = 1/reynolds_number
  
  for(int i = 0; i < N_U; i++)
  {
    test10 = u_x[i];
    test01 = u_y[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz10 = u_x[j];
      ansatz01 = u_y[j];
      MatrixA[i][j] += Mult * nu*(test10*ansatz10+test01*ansatz01);
    }
  }
}

void NSLaplaceGradGrad(double Mult, double *coeff, double *param,
                       double hK, double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double ansatz10, ansatz01, test10, test01;
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[3];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double nu = coeff[0]; // = 1/reynolds_number
  
  for(int i = 0; i < N_U; i++)
  {
    test10 = u_x[i];
    test01 = u_y[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz10 = u_x[j];
      ansatz01 = u_y[j];
      double val = Mult * nu*(test10*ansatz10+test01*ansatz01);
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
    }
  }
}

void NSLaplaceDeformation(double Mult, double *coeff, double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs)
{
  double ansatz10, ansatz01, test10, test01;
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA21 = LocMatrices[2];
  double **MatrixA22 = LocMatrices[3];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double nu = coeff[0];
  for(int i = 0; i < N_U; i++)
  {
    test10 = u_x[i];
    test01 = u_y[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz10 = u_x[j];
      ansatz01 = u_y[j];
      MatrixA11[i][j] += Mult * 2*nu*(test10*ansatz10+0.5*test01*ansatz01);
      MatrixA12[i][j] += Mult * nu*(test01*ansatz10);
      MatrixA21[i][j] += Mult * nu*(test10*ansatz01);
      MatrixA22[i][j] += Mult * 2*nu*(0.5*test10*ansatz10+test01*ansatz01);
    }
  }
}

void NSDivergenceBlocks(double Mult, double *coeff, double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs)
{
  double ansatz10, ansatz01, test00;
  double ** MatrixB1 = LocMatrices[5];
  double ** MatrixB2 = LocMatrices[6];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * Orig1 = OrigValues[1];         // p
  double * Orig2 = OrigValues[2];         // u_x
  double * Orig3 = OrigValues[3];         // u_y
  for(int i = 0; i < N_P; i++)
  {
    test00 = Orig1[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      MatrixB1[i][j] += -Mult*test00*ansatz10;
      MatrixB2[i][j] += -Mult*test00*ansatz01;
    }
  }
}

void NSGradientBlocks(double Mult, double *coeff, double *param, double hK,
                      double **OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs)
{
  double ansatz00, test10, test01;
  double ** MatrixB1T = LocMatrices[7];
  double ** MatrixB2T = LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * Orig1 = OrigValues[1];         // p
  double * Orig2 = OrigValues[2];         // u_x
  double * Orig3 = OrigValues[3];         // u_y
  for(int i = 0; i < N_U; i++)
  {
    test10 = Orig2[i];
    test01 = Orig3[i];
    for(int j = 0; j < N_P; j++)
    {
      ansatz00 = Orig1[j];
      MatrixB1T[i][j] += -Mult*ansatz00*test10;
      MatrixB2T[i][j] += -Mult*ansatz00*test01;
    }
  }
}

void NSRightHandSide(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs)
{
  double test00;
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs3 = LocRhs[2];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * u = OrigValues[0];
  double * p = OrigValues[1];
  double f1 = coeff[1];
  double f2 = coeff[2];
  double g = coeff[3]; // divergence
  for(int i = 0; i < N_U; i++)
  {
    test00 = u[i];
    Rhs1[i] += Mult*test00*f1;
    Rhs2[i] += Mult*test00*f2;
  }
  for(int i = 0; i < N_P; i++)
  {
    test00 = p[i];
    Rhs3[i] += -Mult*test00*g;
  }
}

void NSNonlinearTermSingle(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                         double **LocRhs)
{
  double test00, ansatz10, ansatz01;
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  double * Orig0 = OrigValues[0];         // u
  double * Orig2 = OrigValues[2];         // u_x
  double * Orig3 = OrigValues[3];         // u_y
  double u1 = param[0];                 // u1old
  double u2 = param[1];                 // u2old
  for(int i = 0; i < N_U; i++)
  {
    test00 = Orig0[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      MatrixA[i][j] += Mult * (u1*ansatz10 + u2*ansatz01) * test00;
    }                            // endfor j
  }                              // endfor i
}

void NSNonlinearTerm(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                         double **LocRhs)
{
  double test00, ansatz10, ansatz01;
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[3];
  int N_U = N_BaseFuncts[0];
  double * Orig0 = OrigValues[0];         // u
  double * Orig2 = OrigValues[2];         // u_x
  double * Orig3 = OrigValues[3];         // u_y
  double u1 = param[0];                 // u1old
  double u2 = param[1];                 // u2old
  for(int i = 0; i < N_U; i++)
  {
    test00 = Orig0[i];
    for(int j = 0; j < N_U; j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      MatrixA11[i][j] += Mult * (u1*ansatz10 + u2*ansatz01)*test00;
      MatrixA22[i][j] += Mult * (u1*ansatz10 + u2*ansatz01)*test00;
    }
  }
}


void NSParamsVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
}

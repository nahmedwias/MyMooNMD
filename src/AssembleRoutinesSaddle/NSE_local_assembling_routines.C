#include "NSE_local_assembling_routines.h"

template<int d>
void NSLaplaceGradGradSingle(double Mult, double *coeff, double *param,
                             double hK, double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs)
{
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = d == 2 ? nullptr : OrigValues[4];
  double nu = coeff[0]; // = 1/reynolds_number
  
  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixA[i][j] += Mult * nu*(test_x*ansatz_x + test_y*ansatz_y
                                 +test_z * ansatz_z);
    }
  }
}

template <int d>
void NSLaplaceGradGrad(double Mult, double *coeff, double *param,
                       double hK, double**OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = d == 2 ? nullptr : OrigValues[4];
  double nu = coeff[0]; // = 1/reynolds_number
  
  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = Mult * nu*( test_x * ansatz_x + test_y * ansatz_y 
                             + test_z * ansatz_z);
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if(d == 3)
        MatrixA33[i][j] += val;
    }
  }
}

template <int d>
void NSLaplaceDeformation(double Mult, double *coeff, double *param, double hK,
                          double**OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices,  double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double ** MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d+1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d+2];
  double ** MatrixA31 = d == 2 ? nullptr : LocMatrices[6];
  double ** MatrixA32 = d == 2 ? nullptr : LocMatrices[7];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = d == 2 ? nullptr : OrigValues[4];
  double nu = coeff[0];
  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixA11[i][j] += Mult * 2*nu*(test_x*ansatz_x + 0.5*test_y*ansatz_y 
                                     +0.5*test_z*ansatz_z);
      MatrixA12[i][j] += Mult * nu*(test_y*ansatz_x);
      if(d == 3)
        MatrixA13[i][j] += Mult * nu*(test_z*ansatz_x);
      MatrixA21[i][j] += Mult * nu*(test_x*ansatz_y);
      MatrixA22[i][j] += Mult * 2*nu*(0.5*test_x*ansatz_x+test_y*ansatz_y
                                     +0.5*test_z*ansatz_z);
      if(d == 3)
      {
        MatrixA23[i][j] += Mult * nu*(test_z*ansatz_y);
        MatrixA31[i][j] += Mult * nu*(test_x*ansatz_z);
        MatrixA32[i][j] += Mult * nu*(test_y*ansatz_z);
        MatrixA33[i][j] += Mult * 2*nu*(0.5*test_x*ansatz_x+0.5*test_y*ansatz_y
                                       +test_z*ansatz_z);
      }
    }
  }
}

template <int d>
void NSDivergenceBlocks(double Mult, double *coeff, double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs)
{
  double ** MatrixB1 = LocMatrices[d == 2 ? 5 : 10];
  double ** MatrixB2 = LocMatrices[d == 2 ? 6 : 11];
  double ** MatrixB3 = d == 2 ? nullptr : LocMatrices[12];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * p = OrigValues[1];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = d == 2 ? nullptr : OrigValues[4];
  for(int i = 0; i < N_P; i++)
  {
    double test = p[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixB1[i][j] -= Mult*test*ansatz_x;
      MatrixB2[i][j] -= Mult*test*ansatz_y;
      if(d == 3)
        MatrixB3[i][j] -= Mult*test*ansatz_z;
    }
  }
}

template <int d>
void NSGradientBlocks(double Mult, double *coeff, double *param, double hK,
                      double**OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs)
{
  double ** MatrixB1T = LocMatrices[d == 2 ? 7 : 13];
  double ** MatrixB2T = LocMatrices[d == 2 ? 8 : 14];
  double ** MatrixB3T = d == 2 ? nullptr : LocMatrices[15];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * p = OrigValues[1];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = d == 2 ? nullptr : OrigValues[4];
  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_P; j++)
    {
      double ansatz = p[j];
      MatrixB1T[i][j] -= Mult*ansatz*test_x;
      MatrixB2T[i][j] -= Mult*ansatz*test_y;
      if(d == 3)
        MatrixB3T[i][j] -= Mult*ansatz*test_z;
    }
  }
}

template <int d>
void NSRightHandSide(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs)
{
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs3 = d == 2 ? nullptr : LocRhs[2];
  double * Rhs_div = LocRhs[d];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * u = OrigValues[0];
  double * p = OrigValues[1];
  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = d == 2 ? 0. : coeff[3];
  double g = coeff[d+1]; // divergence
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    Rhs1[i] += Mult*test*f1;
    Rhs2[i] += Mult*test*f2;
    if(d == 3)
      Rhs3[i] += Mult*test*f3;
  }
  for(int i = 0; i < N_P; i++)
  {
    double test = p[i];
    Rhs_div[i] -= Mult*test*g;
  }
}

template <int d>
void NSNonlinearTermSingle(double Mult, double *coeff, double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  double * u   = OrigValues[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = d == 2 ? nullptr : OrigValues[4];
  double u1 = param[0];        
  double u2 = param[1];       
  double u3 = d == 2 ? 0. : param[2];
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixA[i][j] += Mult * (u1*ansatz_x + u2*ansatz_y + u3*ansatz_z) * test;
    }
  }
}

template <int d>
void NSNonlinearTerm(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs)
{
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  double * u = OrigValues[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = d == 2 ? nullptr : OrigValues[4];
  double u1 = param[0];
  double u2 = param[1];
  double u3 = d == 2 ? 0. : param[2];
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = Mult * (u1*ansatz_x + u2*ansatz_y + u3*ansatz_z) * test;
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if(d == 3)
        MatrixA33[i][j] += val;
    }
  }
}

template <int d>
void NSCoriolis(double Mult, double *coeff, double *param, double hK, 
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

double compute_PSPG_delta(double hK, double nu)
{
  return 0.01 * hK * hK / nu;
}

template <int d>
void NSPSPG(double Mult, double *coeff, double *param, double hK,
            double **OrigValues, int *N_BaseFuncts,
            double ***LocMatrices, double **LocRhs)
{
  double ** MatrixB1 = LocMatrices[d == 2 ? 5 : 10];
  double ** MatrixB2 = LocMatrices[d == 2 ? 6 : 11];
  double ** MatrixB3 = d == 2 ? nullptr : LocMatrices[12];
  double ** MatrixC = LocMatrices[d == 2 ? 4 : 9];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  double * p_x = OrigValues[2+d];
  double * p_y = OrigValues[3+d];
  double * p_z = d == 2 ? nullptr : OrigValues[4+d];
  double * u_xx = OrigValues[2+2*d];
  double * u_yy = OrigValues[2+2*d+d];
  double * u_zz = d == 2 ? nullptr : OrigValues[13];
  double nu = coeff[0]; // = 1/reynolds_number
  double delta = compute_PSPG_delta(hK, nu);
  for(int i = 0; i < N_P; i++)
  {
    double test_x = delta * p_x[i];
    double test_y = delta * p_y[i];
    double test_z = d == 2 ? 0. : delta * p_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double laplace = nu * (u_xx[j] + u_yy[j] + (d == 2 ? 0. : u_zz[j]));
      MatrixB1[i][j] += Mult*test_x*laplace;
      MatrixB2[i][j] += Mult*test_y*laplace;
      if(d == 3)
        MatrixB3[i][j] += Mult*test_z*laplace;
    }
    for(int j = 0; j < N_P; j++)
    {
      double ansatz_x = p_x[j];
      double ansatz_y = p_y[j];
      double ansatz_z = d == 2 ? 0. : p_z[j];
      MatrixC[i][j] -= Mult * (test_x*ansatz_x + test_y*ansatz_y 
                              +test_z*ansatz_z);
    }
  }
}

template <int d>
void NSPSPG_RightHandSide(double Mult, double *coeff, double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs)
{
  double * Rhs_div = LocRhs[d];
  int N_P = N_BaseFuncts[1];
  double * p_x = OrigValues[2+d];
  double * p_y = OrigValues[3+d];
  double * p_z = d == 2 ? nullptr : OrigValues[4+d];
  double nu = coeff[0]; // = 1/reynolds_number
  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = d == 2 ? 0. : coeff[3];
  double delta = compute_PSPG_delta(hK, nu);
  for(int i = 0; i < N_P; i++)
  {
    double test_x = p_x[i];
    double test_y = p_y[i];
    double test_z = d == 2 ? 0. : p_z[i];
    Rhs_div[i] -= Mult * delta * (test_x*f1 + test_y*f2 + test_z*f3);
  }
}

template <int d>
void NSParamsVelocity(double *in, double *out)
{
  out[0] = in[d]; // u1old
  out[1] = in[d+1]; // u2old
  if (d == 3)
    out[2] = in[d+2]; // u3old
}

// explicit instatiations below:
#ifdef __2D__
template void NSLaplaceGradGradSingle<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSLaplaceGradGrad<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSLaplaceDeformation<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSDivergenceBlocks<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSGradientBlocks<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSRightHandSide<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSNonlinearTermSingle<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSNonlinearTerm<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSPSPG<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSPSPG_RightHandSide<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
// not yet available
// template void NSCoriolis<2>(
//   double Mult, double *coeff, double *param, double hK, double **OrigValues,
//   int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSParamsVelocity<2>(double *in, double *out);
#endif // 2D
#ifdef __3D__
template void NSLaplaceGradGradSingle<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSLaplaceGradGrad<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSLaplaceDeformation<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSDivergenceBlocks<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSGradientBlocks<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSRightHandSide<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSNonlinearTermSingle<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSNonlinearTerm<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSCoriolis<2>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSParamsVelocity<3>(double *in, double *out);
template void NSPSPG<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void NSPSPG_RightHandSide<3>(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
#endif // 3D

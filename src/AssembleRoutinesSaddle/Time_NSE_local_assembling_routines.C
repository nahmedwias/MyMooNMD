#include "Time_NSE_local_assembling_routines.h"
#include "CommonRoutineTNSE3D.h"
#include <MooNMD_Io.h>

template <int d>
void NSMassMatrixSingle(double Mult, double *, double *, double,
  double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **)
{
  double ** MatrixM = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  double * u = OrigValues[0];  
  
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz = u[j];
      MatrixM[i][j] += Mult * (test * ansatz);
    }
  }
}

template<int d>
void NSMassMatrix(double Mult, double *, double *, double, 
  double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **)
{
  double ** MatrixM11 = LocMatrices[0];
  double ** MatrixM22 = LocMatrices[1];
  double ** MatrixM33 = d == 2 ? nullptr : LocMatrices[2];
  int N_U = N_BaseFuncts[0];
  double * u = OrigValues[0];
  
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz = u[j];
      double val = Mult * (test * ansatz);
      MatrixM11[i][j] += val;
      MatrixM22[i][j] += val;
      if(d == 3)
        MatrixM33[i][j] += val;
    }
  }
}

template<int d> 
void NSLaplaceGradGradSingleSmagorinsky(double Mult, double *, double *param, double hK,
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **)
{
  if(d==2)
    ErrThrow("smagorinsky is not supported in 2D yet");
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = d == 2 ? nullptr : OrigValues[4];
  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  
  double mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
  
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
      MatrixA[i][j] += Mult * mu * (test_x * ansatz_x + test_y * ansatz_y
                                 + test_z * ansatz_z);
    }
  }
}

template <int d>
void NSLaplaceGradGradSmagorinsky(double Mult, double *, double *param, double hK,
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **)
{
  if(d==2)
    ErrThrow("smagorinsky is not supported in 2D yet");
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  double * u_x = OrigValues[2];
  double * u_y = OrigValues[3];
  double * u_z = d == 2 ? nullptr : OrigValues[4];
  
  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  double mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
  
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
      double val = Mult * mu* ( test_x * ansatz_x + test_y * ansatz_y
                             + test_z * ansatz_z);
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if(d == 3)
        MatrixA33[i][j] += val;
    }
  }
}

template <int d>
void NSLaplaceDeformationSmagorinsky(double Mult, double *, double *param, double hK,
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **)
{
  if(d==2)
    ErrThrow("smagorinsky is not supported in 2D yet");
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
  
  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  double nu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711)/2.;
  
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

template<>
void NSParamVelGradSmagorinsky<2>(const double *in, double *out)
{
  ErrThrow("Smagorinsky Model in 2-dimensional case is not supported so far!!");
}

template<>
void NSParamVelGradSmagorinsky<3>(const double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old
  
  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3
  
  out[12] = in[0]; // x - coordinate for van Driest damping
  out[13] = in[1]; // y - coordinate for van Driest damping
  out[14] = in[2]; // z - coordinate for van Driest damping
}



#ifdef __2D__
template void NSMassMatrixSingle<2>(
  double Mult, double *coeff,double *param,double hK,
  double **OrigValues,int *N_BaseFuncts,double ***LocMatrices,
  double **LocRhs);
template void NSMassMatrix<2>(double Mult, double *, double *, double, 
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **);
#endif
#ifdef __3D__
template void NSMassMatrixSingle<3>(
  double Mult, double *coeff,double *param,double hK,
  double **OrigValues,int *N_BaseFuncts,double ***LocMatrices,
  double **LocRhs);
template void NSMassMatrix<3>(
  double Mult, double *, double *, double, double **OrigValues, int *N_BaseFuncts,
  double ***LocMatrices, double **);
template void NSLaplaceGradGradSingleSmagorinsky<3>(
  double Mult, double *coeff, double *, double, double **OrigValues, int *N_BaseFuncts,
  double ***LocMatrices, double **);
template void NSLaplaceGradGradSmagorinsky<3>(double Mult, double *coeff, double *param, double hK,
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);
template void NSLaplaceDeformationSmagorinsky<3>(double Mult, double *coeff, double *param, double hK,
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);
#endif

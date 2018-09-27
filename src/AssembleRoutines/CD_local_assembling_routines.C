#include "CD_local_assembling_routines.h"

template <int d>
void TCDStiff(double Mult, double *coeff, double *param,
                    double hK, double**OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix = LocMatrices[0];
  int N_ = N_BaseFuncts[0];
  double * u  =OrigValues[0];
  double * u_x=OrigValues[1];
  double * u_y=OrigValues[2];
  double * u_z = d == 2 ? nullptr : OrigValues[3];
  
  double eps = coeff[0];
  double b1 = coeff[1];
  double b2 = coeff[2];
  double b3 = d == 2 ? 0. : coeff[3];
  double c = coeff[d+1];
  
  for(int i=0; i<N_; i++)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j=0; j<N_; j++)
    {
      double ansatz=u[j];
      double ansatz_x=u_x[j];
      double ansatz_y=u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = (eps*(test_x*ansatz_x + test_y*ansatz_y + test_z*ansatz_z));
      val += (b1*ansatz_x + b2*ansatz_y+b3*ansatz_z)*test;
      val += c*test*ansatz;
      Matrix[i][j] += Mult * val;
    }
  }
}

template<int d>
void TCDMass(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double **Matrix = LocMatrices[1];
  int N_ = N_BaseFuncts[0];
  double * u  =OrigValues[0];
  
  for(int i=0; i<N_; i++)
    for(int j=0; j<N_; j++)
      Matrix[i][j] += Mult * u[i] * u[j];
}

template<int d>
void TCDRhs(double Mult, double* coeff, double* param, double hK, double ** OrigValues, 
            int* N_BaseFuncts, double *** LocMatrices, double ** LocRhs)
{
  int N_ = N_BaseFuncts[0];
  double *u = OrigValues[0];
  double f = coeff[d+2];
  double *Rhs = LocRhs[0];
  for(int i=0; i<N_; i++)
    Rhs[i] += Mult * u[i] * f;
}

#ifdef __2D__
template void TCDStiff<2>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDMass<2>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDRhs<2>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
#endif // 2D
#ifdef __3D__
template void TCDStiff<3>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDMass<3>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDRhs<3>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
#endif // 3D

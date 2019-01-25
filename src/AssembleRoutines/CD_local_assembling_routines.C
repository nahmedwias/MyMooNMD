#include "CD_local_assembling_routines.h"
#include "Database.h"
#include "ConvDiff.h"
#include <MooNMD_Io.h>

template <int d>
void TCDStiff(double Mult, double *coeff, double *, double , double**OrigValues,
              int *N_BaseFuncts, double ***LocMatrices, double **)
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
void TCDMass(double Mult, double *, double *, double, double**OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **)
{
  double **Matrix = LocMatrices[1];
  int N_ = N_BaseFuncts[0];
  double * u  =OrigValues[0];
  
  for(int i=0; i<N_; i++)
    for(int j=0; j<N_; j++)
      Matrix[i][j] += Mult * u[i] * u[j];
}

// this function is the same as "TCDMass", but assemble the mass matrix in LocMatrices[0]
// needed for POD
template<int d>
void TCDMassPOD(double Mult, double *, double *, double, double**OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **)
{
  double **Matrix = LocMatrices[0];
  int N_ = N_BaseFuncts[0];
  double * u  =OrigValues[0];
  
  for(int i=0; i<N_; i++)
    for(int j=0; j<N_; j++)
      Matrix[i][j] += Mult * u[i] * u[j];
}

template<int d>
void TCDRhs(double Mult, double* coeff, double*, double, double ** OrigValues, 
            int* N_BaseFuncts, double ***, double ** LocRhs)
{
  int N_ = N_BaseFuncts[0];
  double *u = OrigValues[0];
  double f = coeff[d+2];
  double *Rhs = LocRhs[0];
  for(int i=0; i<N_; i++)
    Rhs[i] += Mult * u[i] * f;
}

template <int d>
void TCDStiffSUPG(double Mult, double *coeff, double *, double hK, double**OrigValues, 
              int *N_BaseFuncts, double ***LocMatrices, double **)
{
  double **Matrix = LocMatrices[0];
  int N_ = N_BaseFuncts[0];
  double * u  =OrigValues[0];
  double * u_x=OrigValues[1];
  double * u_y=OrigValues[2];
  double * u_z = d == 2 ? nullptr : OrigValues[3];
  double * u_xx = d== 2 ? OrigValues[3] : OrigValues[4];
  double * u_yy = d== 2 ? OrigValues[5] : OrigValues[7];
  double * u_zz = d== 2 ? nullptr : OrigValues[9];
  
  double eps = coeff[0];
  double b1 = coeff[1];
  double b2 = coeff[2];
  double b3 = d == 2 ? 0. : coeff[3];
  double c = coeff[d+1];
  double delta = 0.;
  if(d==2)
  {
    // stabilization parameter
    double t1 = TDatabase::TimeDB->THETA1;
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    eps = coeff[0] * t1 * tau;
    
    b1  = coeff[1] * t1 * tau;
    b2  = coeff[2] * t1 * tau;
    b3  = d==2 ? 0.0 : coeff[3];
    int sdfem_type = TDatabase::ParamDB->SDFEM_TYPE;
    c = 1. + coeff[d+1] * t1 *tau;
    if(sdfem_type == 8)
      c = coeff[3] * t1 * tau;
    double bb = 0.;
    if(fabs(b1) > fabs(b2) ) 
       bb = fabs(b1);
    else
      bb = fabs(b2);
    delta = Compute_SDFEM_delta<2>(hK, eps, {{b1, b2}}, c, bb);
  }
  else 
  {
    double b_norm = coeff[6];
    delta = Compute_SDFEM_delta<3>(hK, eps, {{b1, b2, b3}}, c, b_norm);
  }  
  
  for(int i=0; i<N_; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    
    double bgradv = delta * (b1 * test_x + b2 * test_y + b3 * test_z);
    
    for(int j=0; j<N_; j++)
    {
      double ansatz=u[j];
      double ansatz_x=u_x[j];
      double ansatz_y=u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double ansatz_xx = u_xx[j];
      double ansatz_yy = u_yy[j];
      double ansatz_zz = d==2 ? 0. : u_zz[j];
      double val = (b1*ansatz_x + b2*ansatz_y+b3*ansatz_z) * bgradv;
      val += c * bgradv * ansatz;
      val -= eps * (ansatz_xx + ansatz_yy + ansatz_zz) * bgradv;
      Matrix[i][j] += Mult * val;
    }
  }
}
template <int d>
void TCDMassSUPG(double Mult, double *coeff, double *, double hK, double**OrigValues, 
             int *N_BaseFuncts, double ***LocMatrices, double **)
{
  double **Matrix = LocMatrices[1];
  int N_ = N_BaseFuncts[0];
  double * u  = OrigValues[0];
  double * ux = OrigValues[1];
  double * uy = OrigValues[2];
  double * uz = d==2 ? nullptr : OrigValues[3];
  
  double eps = coeff[0];
  double b1 = coeff[1];
  double b2 = coeff[2];
  double b3 = d == 2 ? 0. : coeff[3];
  double c = coeff[d+1];
  double delta = 0.;
  // stabilization parameter
  int sdfem_type = TDatabase::ParamDB->SDFEM_TYPE;
  double t1 = TDatabase::TimeDB->THETA1;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  if(d==2)
  {
    // stabilization parameter
    eps = coeff[0] * t1 * tau;
    
    b1  = coeff[1] * t1 * tau;
    b2  = coeff[2] * t1 * tau;
    b3  = d==2 ? 0.0 : coeff[3];
    c = 1. + coeff[d+1] * t1 *tau;
    if(sdfem_type == 8)
      c = coeff[3] * t1 * tau;
    double bb = 0.;
    if(fabs(b1) > fabs(b2) ) 
       bb = fabs(b1);
    else
      bb = fabs(b2);
    delta = Compute_SDFEM_delta<2>(hK, eps, {{b1, b2}}, c, bb);
  }
  else
  {
    double b_norm = coeff[6];
    delta = Compute_SDFEM_delta<3>(hK, eps, {{b1, b2, b3}}, c, b_norm);
  }  
  
  if(sdfem_type !=9 && sdfem_type != 10 && sdfem_type != 11)
    delta *= t1 * tau;
  
  for(int i=0; i<N_; i++)
  {
    double test_x = ux[i];
    double test_y = uy[i];
    double test_z = d==2 ? 0.0 : uz[i];
    double bgradv = delta*(b1 * test_x + b2 * test_y + b3 * test_z);
    
    for(int j=0; j<N_; j++)
      Matrix[i][j] += Mult * u[j] * bgradv;
  }
}

template <int d>
void TCDRhsSUPG(double Mult, double *coeff, double *, double hK,
                double**OrigValues,  int *N_BaseFuncts, double ***,
                double **LocRhs)
{
  int N_ = N_BaseFuncts[0];
  
  double *u_x = OrigValues[1];
  double *u_y = OrigValues[2];
  double *u_z = d==2 ? nullptr : OrigValues[3];
  double f = coeff[d+2];
  double *Rhs = LocRhs[0];
  double eps = coeff[0];
  double b1 = coeff[1];
  double b2 = coeff[2];
  double b3 = d==2 ? 0.0 : coeff[3];
  double c = coeff[d+1];
  double delta;
  if(d==2)
  {
    // stabilization parameter
    double t1 = TDatabase::TimeDB->THETA1;
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    eps = coeff[0] * t1 * tau;
    
    b1  = coeff[1] * t1 * tau;
    b2  = coeff[2] * t1 * tau;
    b3  = d==2 ? 0.0 : coeff[3];
    int sdfem_type = TDatabase::ParamDB->SDFEM_TYPE;
    c = 1. + coeff[d+1] * t1 *tau;
    if(sdfem_type == 8)
      c = coeff[3] * t1 * tau;
    double bb = 0.;
    if(fabs(b1) > fabs(b2) ) 
       bb = fabs(b1);
    else
      bb = fabs(b2);
    delta = Compute_SDFEM_delta<2>(hK, eps, {{b1, b2}}, c, bb);
  }
  else 
  {
    double b_norm = coeff[6];
    delta = Compute_SDFEM_delta<3>(hK, eps, {{b1, b2, b3}}, c, b_norm);
  }  
  
  for(int i=0; i<N_; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d==2 ? 0.0 : u_z[i];
    double bgradv = delta * (b1 * test_x + b2 * test_y + b3 * test_z);
    Rhs[i] += Mult * f * bgradv;
  }
}
#ifdef __2D__
template void TCDStiff<2>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDMass<2>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDMassPOD<2>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDRhs<2>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDStiffSUPG<2>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDMassSUPG<2>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDRhsSUPG<2>(
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
template void TCDMassPOD<3>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDRhs<3>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDStiffSUPG<3>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDMassSUPG<3>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template void TCDRhsSUPG<3>(
  double Mult, double *coeff, double *param, double hK, double**OrigValues, 
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
#endif // 3D

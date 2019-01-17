#include <Coupled_NS_Stress_Local_routines.h>
#include <Database.h>

template<int d>
void Stress_Stress(double Mult, double* coeff, double* param, double hK, double** OrigValues, 
                   int* N_BaseFuncts, double*** LocMatrices, double** LocRhs, double eta)
{
  double **Matrix1 = LocMatrices[0];
  double **Matrix2 = LocMatrices[1];
  double **Matrix3 = LocMatrices[2];
  int N_ = N_BaseFuncts[0];
  double * s  =OrigValues[0];  
  
  for(int i = 0; i < N_; i++)
  {
    double test = s[i];
    for(int j = 0; j < N_; j++)
    {
      double ansatz = s[j];
      double val = Mult * eta * (test * ansatz);
      Matrix1[i][j] += val;
      Matrix2[i][j] += val;
      Matrix3[i][j] += val;
    }
  } 
}


template<int d>
void Stress_velocity(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
              int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double **MatrixST00 = LocMatrices[3];
  //double **MatrixST01 = LocMatrices[4];
  double **MatrixST10 = LocMatrices[5];
  double **MatrixST11 = LocMatrices[6];
  //double **MatrixST20 = LocMatrices[7];
  double **MatrixST21 = LocMatrices[8];
  
  int N_S = N_BaseFuncts[0];
  int N_U = N_BaseFuncts[1];
  //int N_P = N_BaseFuncts[2];
  
  
  double * s  = OrigValues[0];
  //double * u  = OrigValues[1];
  //double * p = OrigValues[2];
  //double * sx = OrigValues[3];
  //double * sy = OrigValues[4];
  double * ux = OrigValues[5];
  double * uy = OrigValues[6];
  
  for(int i=0; i<N_U; i++)
  {
    double testx = ux[i];
    double testy = uy[i];
    for(int j=0; j<N_S; j++)
    {
      double ansatz = s[j];
      
      MatrixST00[i][j] += -Mult * testx * ansatz;
      MatrixST10[i][j] += -Mult * testy * ansatz;
      MatrixST11[i][j] += -Mult * testx * ansatz;
      MatrixST21[i][j] += -Mult * testy * ansatz;
    }
  }
}


template<int d>
void Velocity_stress(double Mult, double *coeff, double *param, double hK, double**OrigValues, 
                     int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double **MatrixS00 = LocMatrices[9];
  double **MatrixS01 = LocMatrices[10];
  //double **MatrixS02 = LocMatrices[11];
  //double **MatrixS10 = LocMatrices[12];
  double **MatrixS11 = LocMatrices[13];
  double **MatrixS12 = LocMatrices[14];
  
  int N_S = N_BaseFuncts[0];
  int N_U = N_BaseFuncts[1];
  //int N_P = N_BaseFuncts[2];
  
  double * s  = OrigValues[0];
  // double * u  = OrigValues[1];
  //double * p = OrigValues[2];
  //double * sx = OrigValues[3];
  //double * sy = OrigValues[4];
  double * ux = OrigValues[5];
  double * uy = OrigValues[6];
  
  for(int i=0; i<N_S; i++)
  {
    double test = s[i];
    for(int j=0; j<N_U; j++)
    {
      double ansatzx = ux[j];
      double ansatzy = uy[j];
      
      MatrixS00[i][j] += Mult * ansatzx * test;
      MatrixS01[i][j] += Mult * ansatzy * test;
      MatrixS11[i][j] += Mult * ansatzx * test;
      MatrixS12[i][j] += Mult * ansatzy * test;
    }
  }
}

template<int d>
void Stress_rhs(double Mult, double* coeff, double* param, double hK, double** OrigValues, 
                int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
{
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs3 = LocRhs[2];
  
  int N_S = N_BaseFuncts[0];
  double f1 = coeff[4];
  double f2 = coeff[5];
  double f3 = coeff[6];
  
  double * s = OrigValues[0];
  for(int i = 0; i < N_S; i++)
  {
    double test = s[i];
    Rhs1[i] += Mult * test * f1;
    Rhs2[i] += Mult * test * f2;
    Rhs3[i] += Mult * test * f3;
  }  
}
#ifdef __2D__
template void Stress_Stress<2>(
  double Mult, double* coeff, double* param, double hK, double** OrigValues, 
  int* N_BaseFuncts, double*** LocMatrices, double** LocRhs, double eta);
template void Stress_velocity<2>(
  double Mult, double* coeff, double* param, double hK, double** OrigValues, 
  int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);
template void Velocity_stress<2>(
  double Mult, double* coeff, double* param, double hK, double** OrigValues, 
  int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);
template void Stress_rhs<2>(
  double Mult, double* coeff, double* param, double hK, double** OrigValues, 
  int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);
#else
template void Stress_Stress<3>(
  double Mult, double* coeff, double* param, double hK, double** OrigValues, 
  int* N_BaseFuncts, double*** LocMatrices, double** LocRhs, double eta);
template void Stress_velocity<3>(
  double Mult, double* coeff, double* param, double hK, double** OrigValues, 
  int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);
template void Velocity_stress<3>(
  double Mult, double* coeff, double* param, double hK, double** OrigValues, 
  int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);
template void Stress_rhs<3>(
  double Mult, double* coeff, double* param, double hK, double** OrigValues, 
  int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);
#endif

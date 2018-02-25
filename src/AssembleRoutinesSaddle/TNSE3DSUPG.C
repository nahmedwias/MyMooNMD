#include <TNSE3DSUPG.h>
#include <CommonRoutineTNSE3D.h>
#include <Database.h>


//================================================================================
void TimeNSType4SUPGDD3D(double Mult, double *coeff, double *param, double hK, 
  double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = LocMatrices[2];
  double **MatrixA21 = LocMatrices[3];
  double **MatrixA22 = LocMatrices[4];
  double **MatrixA23 = LocMatrices[5];
  double **MatrixA31 = LocMatrices[6];
  double **MatrixA32 = LocMatrices[7];
  double **MatrixA33 = LocMatrices[8];
  double **MatrixM11 = LocMatrices[9];
  double **MatrixM22 = LocMatrices[10];
  double **MatrixM33 = LocMatrices[11];
  double **MatrixB1  = LocMatrices[12];
  double **MatrixB2  = LocMatrices[13];
  double **MatrixB3  = LocMatrices[14];
  double **MatrixB1T = LocMatrices[15];
  double **MatrixB2T = LocMatrices[16];
  double **MatrixB3T = LocMatrices[17];
  
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[1];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_z
  double *Orig3 = OrigValues[3]; // u
  double *Orig4 = OrigValues[4]; // p
  double *Orig5 = OrigValues[5]; // px
  double *Orig6 = OrigValues[6]; // py
  double *Orig7 = OrigValues[7]; // pz
    
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  
  double c0=coeff[0];
  double c1=coeff[1];
  double c2=coeff[2];
  double c3=coeff[3];

  double u1=param[0];
  double u2=param[1];
  double u3=param[2];

  double test100, test010, test001, test000;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  // stabilization parameters
  double tau_m, tau_c;
  
  tau_c = TDatabase::ParamDB->DELTA1;
  tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row; 
  double *Matrix21Row, *Matrix22Row, *Matrix23Row;
  double *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixM11Row,*MatrixM22Row,*MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  
  for(int i=0; i<N_U; ++i)
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
    
    double ugradv = tau_m * (u1*test100+u2*test010+u3*test001);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*(test000+ugradv)*c1;
    Rhs2[i] += Mult*(test000+ugradv)*c2;
    Rhs3[i] += Mult*(test000+ugradv)*c3;
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      // Galerkin part
      double val  = c0*(2*test100*ansatz100+test010*ansatz010
                        +test001*ansatz001); // diffusion term

      double conv = (u1*ansatz100 + u2*ansatz010 + u3*ansatz001)*test000; // convective term
      val += conv;
      // supg contribution
      val +=  ((u1*ansatz100 + u2*ansatz010 + u3*ansatz001) ) *ugradv;
      // grad div contribution
      val += tau_c*test100*ansatz100;
      Matrix11Row[j] += Mult * val;// A11 block
            
      val  = c0*(test010*ansatz100);
      // grad div contribution
      val += tau_c * test100 * ansatz010;
      Matrix12Row[j] += Mult * val;
      
      val  = c0*(test001*ansatz100);
      // grad div contribution
      val += tau_c * test100 * ansatz001;
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      // grad div contribution
      val += tau_c * test010 * ansatz100;
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      // convective
      val += conv;
      // supg contribution
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugradv; 
      // grad div contribution
      val += tau_c * test010 * ansatz010;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010) + tau_c * test010 * ansatz001;
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      // grad div contribution
      val += tau_c * test001 * ansatz100;
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      // grad div contribution
      val += tau_c * test001 * ansatz010;
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      // convective term
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      // supg contribution
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugradv;
      // grad div contribution
      val += tau_c * test001 * ansatz001;
      Matrix33Row[j] += Mult * val;
      
      // weighted mass matrix
      val = ansatz000 * (test000 + ugradv);
      MatrixM11Row[j] += Mult * val;
      MatrixM22Row[j] += Mult * val;
      MatrixM33Row[j] += Mult * val;
    }// endfor j<N_U
    
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz functions
      ansatz000 = Orig4[j];
      ansatz100 = Orig5[j];
      ansatz010 = Orig6[j];
      ansatz001 = Orig7[j];
      
      // pressure term
      MatrixRow1[j] += Mult * (-ansatz000 * test100 + ansatz100 * ugradv);      
      MatrixRow2[j] += Mult * (-ansatz000 * test010 + ansatz010 * ugradv);      
      MatrixRow3[j] += Mult * (-ansatz000 * test001 + ansatz001 * ugradv);
    }
  }
  // pressure test functions
  for(int i=0; i<N_P; ++i)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];
    
    // pressure test
    test000 = Orig4[i];
    
    // velocity-pressure block
    for(int j=0;j<N_U;j++)
    {
      // velocity ansatz functions
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      // divergence constraint
      MatrixRow1[j] += Mult*(-test000*ansatz100);
      MatrixRow2[j] += Mult*(-test000*ansatz010);
      MatrixRow3[j] += Mult*(-test000*ansatz001);
    }    
  }
}
//================================================================================
void TimeNSType4NLSUPGDD3D(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = LocMatrices[2];
  double **MatrixA21 = LocMatrices[3];
  double **MatrixA22 = LocMatrices[4];
  double **MatrixA23 = LocMatrices[5];
  double **MatrixA31 = LocMatrices[6];
  double **MatrixA32 = LocMatrices[7];
  double **MatrixA33 = LocMatrices[8];
  double **MatrixM11 = LocMatrices[9];
  double **MatrixM22 = LocMatrices[10];
  double **MatrixM33 = LocMatrices[11];
  double **MatrixB1T  = LocMatrices[12];
  double **MatrixB2T  = LocMatrices[13];
  double **MatrixB3T  = LocMatrices[14];
  
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[1];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_z
  double *Orig3 = OrigValues[3]; // u
  double *Orig4 = OrigValues[4]; // p
  double *Orig5 = OrigValues[5]; // px
  double *Orig6 = OrigValues[6]; // py
  double *Orig7 = OrigValues[7]; // pz
    
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  
  double c0=coeff[0];
  double c1=coeff[1];
  double c2=coeff[2];
  double c3=coeff[3];

  double u1=param[0];
  double u2=param[1];
  double u3=param[2];

  double test100, test010, test001, test000;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  // stabilization parameters
  double tau_m, tau_c;
  
  tau_c = TDatabase::ParamDB->DELTA1;
  tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row; 
  double *Matrix21Row, *Matrix22Row, *Matrix23Row;
  double *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixM11Row,*MatrixM22Row,*MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  
  for(int i=0; i<N_U; ++i)
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
    
    double ugradv = tau_m * (u1*test100+u2*test010+u3*test001);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*(test000+ugradv)*c1;
    Rhs2[i] += Mult*(test000+ugradv)*c2;
    Rhs3[i] += Mult*(test000+ugradv)*c3;
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      // Galerkin part
      double val  = c0*(2*test100*ansatz100+test010*ansatz010
                        +test001*ansatz001); // diffusion term

      double conv = (u1*ansatz100 + u2*ansatz010 + u3*ansatz001)*test000; // convective term
      val += conv;
      // supg contribution
      val +=  ((u1*ansatz100 + u2*ansatz010 + u3*ansatz001) ) *ugradv;
      // grad div contribution
      val += tau_c*test100*ansatz100;
      Matrix11Row[j] += Mult * val;// A11 block
            
      val  = c0*(test010*ansatz100);
      // grad div contribution
      val += tau_c * test100 * ansatz010;
      Matrix12Row[j] += Mult * val;
      
      val  = c0*(test001*ansatz100);
      // grad div contribution
      val += tau_c * test100 * ansatz001;
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      // grad div contribution
      val += tau_c * test010 * ansatz100;
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      // convective
      val += conv;
      // supg contribution
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugradv; 
      // grad div contribution
      val += tau_c * test010 * ansatz010;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010) + tau_c * test010 * ansatz001;
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      // grad div contribution
      val += tau_c * test001 * ansatz100;
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      // grad div contribution
      val += tau_c * test001 * ansatz010;
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      // convective term
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      // supg contribution
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugradv;
      // grad div contribution
      val += tau_c * test001 * ansatz001;
      Matrix33Row[j] += Mult * val;
      
      // weighted mass matrix
      val = ansatz000 * (test000 + ugradv);
      MatrixM11Row[j] += Mult * val;
      MatrixM22Row[j] += Mult * val;
      MatrixM33Row[j] += Mult * val;
    }// endfor j<N_U
    
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz functions
      ansatz000 = Orig4[j];
      ansatz100 = Orig5[j];
      ansatz010 = Orig6[j];
      ansatz001 = Orig7[j];
      
      // pressure term
      MatrixRow1[j] += Mult * (-ansatz000 * test100 + ansatz100 * ugradv);      
      MatrixRow2[j] += Mult * (-ansatz000 * test010 + ansatz010 * ugradv);      
      MatrixRow3[j] += Mult * (-ansatz000 * test001 + ansatz001 * ugradv);
    }
  }
}
//================================================================================
void TimeNSType4RHSSUPG3D(double Mult, double* coeff, double* param, double hK, 
               double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, 
               double** LocRhs)
{  
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_z
  double *Orig3 = OrigValues[3]; // u
  
  double c1=coeff[1];
  double c2=coeff[2];
  double c3=coeff[3];

  double u1=param[0];
  double u2=param[1];
  double u3=param[2];

  double tau_m =  TDatabase::ParamDB->DELTA0*hK*hK;
  
  double test100, test010, test001, test000;
  int N_U = N_BaseFuncts[0];
  for(int i=0; i<N_U; ++i)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    
    double ugradv = tau_m * (u1*test100+u2*test010+u3*test001);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*c1*(test000 + ugradv);
    Rhs2[i] += Mult*c2*(test000 + ugradv);
    Rhs3[i] += Mult*c3*(test000 + ugradv);
  }
}

//================================================================================
void TimeNSType14SUPGDD3D(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  
}
//================================================================================
void TimeNSType14NLSUPGDD3D(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  
}

//================================================================================
void TimeNSType14RHSSUPGDD3D(double Mult, double* coeff, double* param, double hK, 
               double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, 
               double** LocRhs)
{
  
}

void TimeNSType4Params_SUPG(double *in, double *out)
{
  out[0] = in[3];
  out[1] = in[4];
  out[2] = in[5];
}

void TimeNSType14Params_SUPG(double *in, double *out)
{
  out[0] = in[3];
  out[1] = in[4];
  out[2] = in[5];
  // u1old, u2old, u3old previous time 
  out[3] = in[6]; 
  out[4] = in[7]; 
  out[5] = in[8];
}
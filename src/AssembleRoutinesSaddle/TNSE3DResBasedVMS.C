#include <TNSE3DResBasedVMS.h>
#include <Database.h>
#include <CommonRoutineTNSE3D.h>

//================================================================================
// Local assemble routines for the residual based VMS method:: It's implemented
// in the following way, the test nonlinearity is dealt using extrapolation 
// of velocity and pressure, whereas the ansatz are deal with the fix-point 
// iterations. 
//================================================================================
void TimeNSType4Residual_VMSDD3D(double Mult, double* coeff, double* param, double hK, 
  double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
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
  double **MatrixM12 = LocMatrices[10];
  double **MatrixM13 = LocMatrices[11];
  double **MatrixM21 = LocMatrices[12];
  double **MatrixM22 = LocMatrices[13];
  double **MatrixM23 = LocMatrices[14];
  double **MatrixM31 = LocMatrices[15];
  double **MatrixM32 = LocMatrices[16];
  double **MatrixM33 = LocMatrices[17];
  
  double **MatrixB1  = LocMatrices[18];
  double **MatrixB2  = LocMatrices[19];
  double **MatrixB3  = LocMatrices[20];
  double **MatrixB1T = LocMatrices[21];
  double **MatrixB2T = LocMatrices[22];
  double **MatrixB3T = LocMatrices[23];
  
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_z
  double *Orig3 = OrigValues[3]; // u
  double *Orig4 = OrigValues[4]; // p
  double *Orig5 = OrigValues[5]; // px
  double *Orig6 = OrigValues[6]; // py
  double *Orig7 = OrigValues[7]; // pz
  
  double *Orig8 = OrigValues[8];
  double *Orig9 = OrigValues[9];
  double *Orig10 = OrigValues[10];
    
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  
  double c0=coeff[0];
  double c1=coeff[1];
  double c2=coeff[2];
  double c3=coeff[3];

  // extrapolated velocity:
  // NOTE: In the residual based VMS, we only use the extraplation of the velocity and pressure 
  double u1=param[0];
  double u2=param[1];
  double u3=param[2];
  
  double u1x=param[3];
  double u2x=param[4];
  double u3x=param[5];
  
  double u1y=param[6];
  double u2y=param[7];
  double u3y=param[8];
  
  double u1z=param[9];
  double u2z=param[10];
  double u3z=param[11];
  
  double u1xx=param[12];
  double u2xx=param[13];
  double u3xx=param[14];
  
  double u1yy=param[15];
  double u2yy=param[16];
  double u3yy=param[17];
  
  double u1zz=param[18];
  double u2zz=param[19];
  double u3zz=param[20];
  
  //double p= param[21];
  double px=param[22];
  double py=param[23];
  double pz=param[24];
  
  double u1td=param[25];
  double u2td=param[26];
  double u3td=param[27];
  
  
  double test100, test010, test001, test000;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double ansatz200, ansatz020, ansatz002;
  // stabilization parameters
  double tau_m, tau_c;
  
  tau_c = TDatabase::ParamDB->DELTA1;
  tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  
  
  // old residual to be used Eq(51)
  double res1 = tau_m *(c1 - u1td + c0*(u1xx + u1yy + u1zz) - (u1*u1x + u2*u1y + u3*u1z) - px );
  double res2 = tau_m *(c2 - u2td + c0*(u2xx + u2yy + u2zz) - (u1*u2x + u2*u2y + u3*u2z) - py );
  double res3 = tau_m *(c3 - u3td + c0*(u3xx + u3yy + u3zz) - (u1*u3x + u2*u3y + u3*u3z) - pz );
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row; 
  double *Matrix21Row, *Matrix22Row, *Matrix23Row;
  double *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixM11Row,*MatrixM22Row,*MatrixM33Row;
  double *MatrixM12Row,*MatrixM13Row,*MatrixM23Row;
  double *MatrixM21Row,*MatrixM31Row,*MatrixM32Row;
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
    MatrixM12Row  = MatrixM12[i];
    MatrixM13Row  = MatrixM13[i];
    MatrixM21Row  = MatrixM21[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM23Row  = MatrixM23[i];
    MatrixM31Row  = MatrixM31[i];
    MatrixM32Row  = MatrixM32[i];
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
    // contribution from the second and third nonlinear terms
    Rhs1[i] += Mult*tau_m*(u1+res1)*(c1*test100 + c2*test010 + c3*test001);
    Rhs2[i] += Mult*tau_m*(u2+res2)*(c1*test100 + c2*test010 + c3*test001);
    Rhs3[i] += Mult*tau_m*(u3+res3)*(c1*test100 + c2*test010 + c3*test001);
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
      // second cross term and subgrid term
      ansatz200 = Orig8[j];
      ansatz020 = Orig9[j];
      ansatz002 = Orig10[j];
      double laplacian = -c0*(ansatz200 + ansatz020 + ansatz002);
      val += tau_m * (laplacian + conv)*(u1+res1)*test100;
      // grad div contribution
      val += tau_c*test100*ansatz100;
      Matrix11Row[j] += Mult * val;// A11 block
            
      val  = c0*(test010*ansatz100);
      // second cross term and subgrid term
      val += tau_m * conv*(u1+res1)*test010;
      // grad div contribution
      val += tau_c * test100 * ansatz010;
      Matrix12Row[j] += Mult * val;
      
      val  = c0*(test001*ansatz100);
      // second cross term and subgrid term
      val += tau_m * conv*(u1+res1)*test001;
      // grad div contribution
      val += tau_c * test100 * ansatz001;
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      // second cross term and subgrid term
      val += tau_m * conv * (u2+res2)*test100;
      // grad div contribution
      val += tau_c * test010 * ansatz100;
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      // convective
      val += conv;
      // supg contribution
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugradv; 
      // second cross term and subgrid term
      val += tau_m * (laplacian + conv)*(u2+res2)*test010;
      // grad div contribution
      val += tau_c * test010 * ansatz010;
      Matrix22Row[j] += Mult * val;

      //galerkin and grad-div contribution
      val  = c0*(test001*ansatz010) + tau_c * test010 * ansatz001;
      // second cross term and subgrid term
      val += tau_m * conv * (u2+res2) * test001;
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      // grad div contribution
      val += tau_c * test001 * ansatz100;
      // second cross term and subgrid term
      val += tau_m * conv * (u3+res3)* test100;
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      // grad div contribution
      val += tau_c * test001 * ansatz010;
      val += tau_m * conv * (u3+res3) * test010;
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      // convective term
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      // supg contribution
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugradv;
      // grad div contribution
      val += tau_c * test001 * ansatz001;
      // second cross term and subgrid term
      val += tau_m * (laplacian + conv) * u3 * test001;
      Matrix33Row[j] += Mult * val;
      
      // mass matrix blocks
      // Galerkin + SUPG term 
      val = ansatz000*(test000 + ugradv);
      // second cross term and subgrid term 
      val += tau_m * (u1 + res1)   * ansatz000 * test100;
      MatrixM11Row[j] += Mult * val;

      val = tau_m * (u1 + res1)  * ansatz000 * test010;
      MatrixM12Row[j] += Mult*val;
      
      val = tau_m * (u1 +res1) * ansatz000 * test001;
      MatrixM13Row[j] += Mult*val;
      
      val = tau_m * (u2 + res2) * ansatz000 * test100;
      MatrixM21Row[j] += Mult * val;
      
      val = ansatz000*(test000 + ugradv);
      val += tau_m * (u2+res2) * ansatz000 * test010;
      MatrixM22Row[j] += Mult * val;
      
      val = tau_m * (u2 + res2)   * ansatz000 * test001;
      MatrixM23Row[j] += Mult * val ;
      
      val = tau_m * (u3 +res3) * ansatz000 * test100;
      MatrixM31Row[j] += Mult * val;
      
      val = tau_m * (u3 +res3) * ansatz000 * test010;
      MatrixM32Row[j] += Mult * val;
      
      val = ansatz000*(test000 + ugradv);
      val += tau_m * (u3 + res3) * ansatz000 * test001;
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
      
      // galerkin and supg terms
      double val = (-ansatz000 * test100 + ansatz100 * ugradv);
      // second cross term and subgrid term 
      val += tau_m * (u1+res1) * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixRow1[j] += Mult * val;
      
      // galerkin and supg terms
      val  = (-ansatz000 * test010 + ansatz010 * ugradv);
      // second cross term and subgrid term 
      val += tau_m * (u2+res2) * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixRow2[j] += Mult * val;
      
      // galerkin and supg terms
      val  = -ansatz000 * test001 + ansatz001 * ugradv;
      // second cross term and subgrid term 
      val += tau_m * (u3+res3) * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixRow3[j] += Mult * val;
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

void TimeNSType4NLResidual_VMSDD3D(double Mult, double* coeff, double* param, double hK, 
  double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
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
  double **MatrixM12 = LocMatrices[10];
  double **MatrixM13 = LocMatrices[11];
  double **MatrixM21 = LocMatrices[12];
  double **MatrixM22 = LocMatrices[13];
  double **MatrixM23 = LocMatrices[14];
  double **MatrixM31 = LocMatrices[15];
  double **MatrixM32 = LocMatrices[16];
  double **MatrixM33 = LocMatrices[17];
  
  double **MatrixB1T  = LocMatrices[18];
  double **MatrixB2T  = LocMatrices[19];
  double **MatrixB3T  = LocMatrices[20];
  
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_z
  double *Orig3 = OrigValues[3]; // u
  double *Orig4 = OrigValues[4]; // p
  double *Orig5 = OrigValues[5]; // px
  double *Orig6 = OrigValues[6]; // py
  double *Orig7 = OrigValues[7]; // pz
  
  double *Orig8 = OrigValues[8];
  double *Orig9 = OrigValues[9];
  double *Orig10 = OrigValues[10];
    
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  
  double c0=coeff[0];
  double c1=coeff[1];
  double c2=coeff[2];
  double c3=coeff[3];

  // extrapolated velocity:
  // NOTE: In the residual based VMS, we only use the extraplation of the velocity and pressure 
  double u1=param[0];
  double u2=param[1];
  double u3=param[2];
  
  double u1x=param[3];
  double u2x=param[4];
  double u3x=param[5];
  
  double u1y=param[6];
  double u2y=param[7];
  double u3y=param[8];
  
  double u1z=param[9];
  double u2z=param[10];
  double u3z=param[11];
  
  double u1xx=param[12];
  double u2xx=param[13];
  double u3xx=param[14];
  
  double u1yy=param[15];
  double u2yy=param[16];
  double u3yy=param[17];
  
  double u1zz=param[18];
  double u2zz=param[19];
  double u3zz=param[20];
  
  //double p= param[21];
  double px=param[22];
  double py=param[23];
  double pz=param[24];
  
  double u1td=param[25];
  double u2td=param[26];
  double u3td=param[27];
  
  
  double test100, test010, test001, test000;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double ansatz200, ansatz020, ansatz002;
  // stabilization parameters
  double tau_m, tau_c;
  
  tau_c = TDatabase::ParamDB->DELTA1;
  tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  
  
  // old residual to be used Eq(51)
  double res1 = tau_m *(c1 - u1td + c0*(u1xx + u1yy + u1zz) - (u1*u1x + u2*u1y + u3*u1z) - px );
  double res2 = tau_m *(c2 - u2td + c0*(u2xx + u2yy + u2zz) - (u1*u2x + u2*u2y + u3*u2z) - py );
  double res3 = tau_m *(c3 - u3td + c0*(u3xx + u3yy + u3zz) - (u1*u3x + u2*u3y + u3*u3z) - pz );
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row; 
  double *Matrix21Row, *Matrix22Row, *Matrix23Row;
  double *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixM11Row,*MatrixM22Row,*MatrixM33Row;
  double *MatrixM12Row,*MatrixM13Row,*MatrixM23Row;
  double *MatrixM21Row,*MatrixM31Row,*MatrixM32Row;
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
    MatrixM12Row  = MatrixM12[i];
    MatrixM13Row  = MatrixM13[i];
    MatrixM21Row  = MatrixM21[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM23Row  = MatrixM23[i];
    MatrixM31Row  = MatrixM31[i];
    MatrixM32Row  = MatrixM32[i];
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
    // contribution from the second and third nonlinear terms
    Rhs1[i] += Mult*tau_m*(u1+res1)*(c1*test100 + c2*test010 + c3*test001);
    Rhs2[i] += Mult*tau_m*(u2+res2)*(c1*test100 + c2*test010 + c3*test001);
    Rhs3[i] += Mult*tau_m*(u3+res3)*(c1*test100 + c2*test010 + c3*test001);
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
      // second cross term and subgrid term
      ansatz200 = Orig8[j];
      ansatz020 = Orig9[j];
      ansatz002 = Orig10[j];
      double laplacian = -c0*(ansatz200 + ansatz020 + ansatz002);
      val += tau_m * (laplacian + conv)*(u1+res1)*test100;
      // grad div contribution
      val += tau_c*test100*ansatz100;
      Matrix11Row[j] += Mult * val;// A11 block
            
      val  = c0*(test010*ansatz100);
      // second cross term and subgrid term
      val += tau_m * conv*(u1+res1)*test010;
      // grad div contribution
      val += tau_c * test100 * ansatz010;
      Matrix12Row[j] += Mult * val;
      
      val  = c0*(test001*ansatz100);
      // second cross term and subgrid term
      val += tau_m * conv*(u1+res1)*test001;
      // grad div contribution
      val += tau_c * test100 * ansatz001;
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      // second cross term and subgrid term
      val += tau_m * conv * (u2+res2)*test100;
      // grad div contribution
      val += tau_c * test010 * ansatz100;
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      // convective
      val += conv;
      // supg contribution
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugradv; 
      // second cross term and subgrid term
      val += tau_m * (laplacian + conv)*(u2+res2)*test010;
      // grad div contribution
      val += tau_c * test010 * ansatz010;
      Matrix22Row[j] += Mult * val;

      //galerkin and grad-div contribution
      val  = c0*(test001*ansatz010) + tau_c * test010 * ansatz001;
      // second cross term and subgrid term
      val += tau_m * conv * (u2+res2) * test001;
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      // grad div contribution
      val += tau_c * test001 * ansatz100;
      // second cross term and subgrid term
      val += tau_m * conv * (u3+res3)* test100;
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      // grad div contribution
      val += tau_c * test001 * ansatz010;
      val += tau_m * conv * (u3+res3) * test010;
      Matrix32Row[j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      // convective term
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      // supg contribution
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*ugradv;
      // grad div contribution
      val += tau_c * test001 * ansatz001;
      // second cross term and subgrid term
      val += tau_m * (laplacian + conv) * u3 * test001;
      Matrix33Row[j] += Mult * val;
      
      // mass matrix blocks
      // Galerkin + SUPG term 
      val = ansatz000*(test000 + ugradv);
      // second cross term and subgrid term 
      val += tau_m * (u1 + res1)   * ansatz000 * test100;
      MatrixM11Row[j] += Mult * val;

      val = tau_m * (u1 + res1)  * ansatz000 * test010;
      MatrixM12Row[j] += Mult*val;
      
      val = tau_m * (u1 +res1) * ansatz000 * test001;
      MatrixM13Row[j] += Mult*val;
      
      val = tau_m * (u2 + res2) * ansatz000 * test100;
      MatrixM21Row[j] += Mult * val;
      
      val = ansatz000*(test000 + ugradv);
      val += tau_m * (u2+res2) * ansatz000 * test010;
      MatrixM22Row[j] += Mult * val;
      
      val = tau_m * (u2 + res2)   * ansatz000 * test001;
      MatrixM23Row[j] += Mult * val ;
      
      val = tau_m * (u3 +res3) * ansatz000 * test100;
      MatrixM31Row[j] += Mult * val;
      
      val = tau_m * (u3 +res3) * ansatz000 * test010;
      MatrixM32Row[j] += Mult * val;
      
      val = ansatz000*(test000 + ugradv);
      val += tau_m * (u3 + res3) * ansatz000 * test001;
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
      
      // galerkin and supg terms
      double val = (-ansatz000 * test100 + ansatz100 * ugradv);
      // second cross term and subgrid term 
      val += tau_m * (u1+res1) * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixRow1[j] += Mult * val;
      
      // galerkin and supg terms
      val  = (-ansatz000 * test010 + ansatz010 * ugradv);
      // second cross term and subgrid term 
      val += tau_m * (u2+res2) * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixRow2[j] += Mult * val;
      
      // galerkin and supg terms
      val  = -ansatz000 * test001 + ansatz001 * ugradv;
      // second cross term and subgrid term 
      val += tau_m * (u3+res3) * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixRow3[j] += Mult * val;
    }
  }
}

void TimeNSType4RHS_Residual_VMS(double Mult, double* coeff, double* param, double hK, 
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
  
  double c0=coeff[0];
  double c1=coeff[1];
  double c2=coeff[2];
  double c3=coeff[3];

  double u1=param[0];
  double u2=param[1];
  double u3=param[2];
  
  double u1x=param[3];
  double u2x=param[4];
  double u3x=param[5];
  
  double u1y=param[6];
  double u2y=param[7];
  double u3y=param[8];
  
  double u1z=param[9];
  double u2z=param[10];
  double u3z=param[11];
  
  double u1xx=param[12];
  double u2xx=param[13];
  double u3xx=param[14];
  
  double u1yy=param[15];
  double u2yy=param[16];
  double u3yy=param[17];
  
  double u1zz=param[18];
  double u2zz=param[19];
  double u3zz=param[20];
  
  //double p= param[21];
  double px=param[22];
  double py=param[23];
  double pz=param[24];
  
  double u1td=param[25];
  double u2td=param[26];
  double u3td=param[27];

  double tau_m =  TDatabase::ParamDB->DELTA0*hK*hK;
  
  double test100, test010, test001, test000;
  int N_U = N_BaseFuncts[0];
  
  
  double res1 = tau_m *(c1 - u1td + c0*(u1xx + u1yy + u1zz) - (u1*u1x + u2*u1y + u3*u1z) - px );
  double res2 = tau_m *(c2 - u2td + c0*(u2xx + u2yy + u2zz) - (u1*u2x + u2*u2y + u3*u2z) - py );
  double res3 = tau_m *(c3 - u3td + c0*(u3xx + u3yy + u3zz) - (u1*u3x + u2*u3y + u3*u3z) - pz );
  
  for(int i=0; i<N_U; ++i)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    
    // standard terms 
    double ugradv = tau_m * (u1*test100+u2*test010+u3*test001);
    // standard terms 
    Rhs1[i] += Mult*(test000+ugradv)*c1;
    Rhs2[i] += Mult*(test000+ugradv)*c2;
    Rhs3[i] += Mult*(test000+ugradv)*c3;
    // contribution from the second and third nonlinear terms
    Rhs1[i] += Mult*tau_m*(u1+res1)*(c1*test100 + c2*test010 + c3*test001);
    Rhs2[i] += Mult*tau_m*(u2+res2)*(c1*test100 + c2*test010 + c3*test001);
    Rhs3[i] += Mult*tau_m*(u3+res3)*(c1*test100 + c2*test010 + c3*test001);
  }
}

//===============================================================================
//================================================================================
void TimeNSType14L_NL_Residual_VMS3D(double Mult, double *coeff, double *param, double hK, 
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
  
  // pressure-pressure block
  double **MatrixC   = LocMatrices[12];
  
  double **MatrixB1  = LocMatrices[13];
  double **MatrixB2  = LocMatrices[14];
  double **MatrixB3  = LocMatrices[15];
  double **MatrixB1T = LocMatrices[16];
  double **MatrixB2T = LocMatrices[17];
  double **MatrixB3T = LocMatrices[18];
  
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];
  double *Rhs4 = LocRhs[3];
  
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
  
  double u1_old_time = param[3];
  double u2_old_time = param[4];
  double u3_old_time = param[5];

  double test100, test010, test001, test000;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  // stabilization parameters
  double stab_params[2];
  stabilization_parameters_equal_order3D(Mult, param, coeff, stab_params);
  double tau_m = stab_params[0];
  double tau_c = stab_params[1];
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row; 
  double *Matrix21Row, *Matrix22Row, *Matrix23Row;
  double *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixM11Row,*MatrixM22Row,*MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixRowC;
  double val;
  
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
      val  = c0*(2*test100*ansatz100+test010*ansatz010
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
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 
  double t1 = TDatabase::TimeDB->THETA1;
  double factor = 1./(dt*t1);
  // pressure test functions
  for(int i=0; i<N_P; ++i)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];
    MatrixRowC = MatrixC[i];
    
    // pressure test
    test000 = Orig4[i];
    test100 = Orig5[i];
    test010 = Orig6[i];
    test001 = Orig7[i];
    
    // 
    val  = (u1_old_time/dt + c1)*test100;
    val += (u2_old_time/dt+c2)*test010;
    val += (u3_old_time/dt+c3)*test001;
    Rhs4[i] += Mult * tau_m * val;
    
    // velocity-pressure block
    for(int j=0;j<N_U;j++)
    {
      // velocity ansatz functions
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      // divergence constraint
      val = -test000*ansatz100;
      // supg contribution
      val -= tau_m*(ansatz000 * factor + u1*ansatz100 + u2 * ansatz010 + u2*ansatz001)*test100;
      MatrixRow1[j] -= Mult*val;
      
      val = -test000*ansatz010;
      // supg contribution
      val -= tau_m*(ansatz000 * factor + u1*ansatz100 + u2 * ansatz010 + u2*ansatz001)*test010;
      MatrixRow2[j] -= Mult*val;
      
      val = -test000*ansatz001;
      // supg contribution
      val -= tau_m*(ansatz000 * factor + u1*ansatz100 + u2 * ansatz010 + u2*ansatz001)*test001;
      MatrixRow3[j] -= Mult*val;
    }
    // pressure-pressure block
    for(int j=0; j<N_P; j++)
    {
      // pressure ansatz
      ansatz100 = Orig5[j];
      ansatz010 = Orig6[j];
      ansatz001 = Orig7[j];
      
      val = tau_m * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixRowC[j] += Mult*val;
    }
  }
}

//===============================================================================
void TimeNSType14RHS_Residual_VMS3D(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];
  // pressure part
  double *Rhs4 = LocRhs[3];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_z
  double *Orig3 = OrigValues[3]; // u
  double *Orig4 = OrigValues[4]; // px
  double *Orig5 = OrigValues[5]; // py
  double *Orig6 = OrigValues[6]; // pz
    
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  
  double c1=coeff[1];
  double c2=coeff[2];
  double c3=coeff[3];

  double u1=param[0];
  double u2=param[1];
  double u3=param[2];
  
  double u1_old_time = param[3];
  double u2_old_time = param[4];
  double u3_old_time = param[5];

  double test100, test010, test001, test000;
  // stabilization parameters
  double stab_params[2];
  stabilization_parameters_equal_order3D(Mult, param, coeff, stab_params);
  double tau_m = stab_params[0];
  
  double val;
  
  for(int i=0; i<N_U; ++i)
  {
    
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
  }
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 
  // pressure test functions
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test100 = Orig4[i];
    test010 = Orig5[i];
    test001 = Orig6[i];
    
    // 
    val  = (u1_old_time/dt + c1)*test100;
    val += (u2_old_time/dt +c2)*test010;
    val += (u3_old_time/dt +c3)*test001;
    Rhs4[i] += Mult * tau_m * val;
  }
}
//===============================================================================
void TimeNSParams_Type4Residual_VMS3D(double* in, double* out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3]= in[6]; //  u1x
  out[4]= in[7]; //  u2x
  out[5]= in[8]; //  u3x
      
  out[6]= in[9];  // u1y
  out[7]= in[10]; // u2y
  out[8]= in[11]; // u3y
      
  out[9]= in[12];  //  u1z
  out[10]= in[13]; //  u2z
  out[11]= in[14]; //  u3z
      
  out[12]= in[15]; //  u1xx
  out[13]= in[16]; //  u2xx
  out[14]= in[17]; //  u3xx
      
  out[15]= in[18]; // u1yy
  out[16]= in[19]; // u2yy
  out[17]= in[20]; // u3yy
  
  out[18]= in[21]; //  u1zz
  out[19] = in[22]; // u2zz
  out[20] = in[23]; // u3zz
  
  out[21] = in[24]; // p      
  out[22] = in[25]; // px
  out[23] = in[26]; // py
  out[24] = in[27]; // pz
  // combined old time for bdf 2: this is used in the right hand side computation
  // for the right hand side computation      
  out[25] = in[28]; // u1'
  out[26] = in[29]; // u2'
  out[27] = in[30]; // u3'
}

void TimeNSParams_Type14Residual_VMS3D(double* in, double* out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old
  // combined old time for bdf 2: this is used in the right hand side computation
  // for the right hand side computation
  out[3]= in[6]; // combined old sols u1
  out[4]= in[7]; // combined old sols u2
  out[5]= in[8]; // combined old sols u3
      
  out[6]= in[9]; // u1x
  out[7]= in[10]; // u2x
  out[8]= in[11]; // u3x
      
  out[9]= in[12]; // u1y
  out[10]= in[13]; // u2y
  out[11]= in[14]; // u3y
      
  out[12]= in[15]; // u1z
  out[13]= in[16]; // u2z
  out[14]= in[17]; // u3z
      
  out[15]= in[18]; // px  
  out[16]= in[19]; // py  
  out[17]= in[20]; // pz  
  out[18]= in[21]; // p
                   
  out[19] = in[22]; // u1yy
  out[20] = in[23]; // u2yy
  out[21] = in[24]; // u3yy
      
  out[22] = in[25]; // u1zz
  out[23] = in[26]; // u2zz
  out[24] = in[27]; // u3zz
      
  out[25] = in[28]; // u1xx
  out[26] = in[29]; // u2xx
  out[27] = in[30]; // u3xx

  // time derivative from previous time   
  out[28] = in[31]; // u1 td
  out[29] = in[32]; // u2 td
  out[30] = in[33]; // u3 td
}

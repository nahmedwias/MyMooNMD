#include <TNSE2DResBasedVMS.h>
#include <Database.h>
#include <CommonRoutineTNSE2D.h>

//================================================================================
// Local assemble routines for the residual based VMS method:: It's implemented
// in the following way, the test nonlinearity is dealt using extrapolation 
// of velocity and pressure, whereas the ansatz are deal with the fix-point 
// iterations. 
//================================================================================
void TimeNSType4Residual_VMS(double Mult, double* coeff, double* param, double hK,
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
{
  double **sqMatrixA11 = LocMatrices[0];
  double **sqMatrixA12 = LocMatrices[1];
  double **sqMatrixA21 = LocMatrices[2];
  double **sqMatrixA22 = LocMatrices[3];
  // weighted mass matrix
  double **sqMatrixM11 = LocMatrices[4];
  double **sqMatrixM12 = LocMatrices[5];
  double **sqMatrixM21 = LocMatrices[6];
  double **sqMatrixM22 = LocMatrices[7];
  
  double **reMatrixB1 = LocMatrices[8];
  double **reMatrixB2 = LocMatrices[9];
  double **reMatrixB1T= LocMatrices[10];
  double **reMatrixB2T= LocMatrices[11];
  
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u
  
  double *Orig3 = OrigValues[3]; // p
  double *Orig4 = OrigValues[4]; // p_x
  double *Orig5 = OrigValues[5]; // p_y
  
  double *Orig6 = OrigValues[6]; // u_xx
  double *Orig7 = OrigValues[7]; // u_yy
  
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  
  double c0=coeff[0];
  double f1=coeff[1];
  double f2=coeff[2];

  double u1=param[0];
  double u2=param[1];
  double u1_pts=param[2]; // 2 u1_m-1 - 1/2 u1_m-2
  double u2_pts=param[3]; // 2 u2_m-1 - 1/2 u2_m-2
  // solution from the previous time step and derivatives
  // this can also be constant or linearly extrapolated
  double u1x = param[4]; // u1x_old
  double u2x = param[5]; // u2x_old
  double u1y = param[6]; // u1y_old
  double u2y = param[7]; // u2y_old
  double u1xx= param[8]; // u1xx_old
  double u2xx= param[9]; // u2xx_old
  double u1yy= param[10]; // u1yy_old
  double u2yy= param[11]; // u2yy_old
  double px  = param[12]; // p_x
  double py  = param[13]; // p_y
  double p   = param[14]; // not needed
  // time derivatives
  double u1_t= param[15]; // u1t
  double u2_t= param[16]; // u2t 

  double test10, test01, test00;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;
  double tau_m, tau_c;
  tau_m =  hK*hK*TDatabase::ParamDB->DELTA0;
  tau_c =  TDatabase::ParamDB->DELTA1;
  
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + c0*(u1xx + u1yy) - (u1*u1x + u2*u1y) - px);
  res2 = tau_m * (f2-u2_t + c0*(u2xx + u2yy) - (u1*u2x + u2*u2y) - py);
  
  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugradv = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*(test00+ugradv)*f1;
    Rhs2[i] += Mult*(test00+ugradv)*f2;
    // contribution from the second and third nonlinear terms
    Rhs1[i] += Mult*tau_m*(u1+res1)*(f1*test10 + f2*test01);
    Rhs2[i] += Mult*tau_m*(u2+res2)*(f1*test10 + f2*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];
      
      double laplacian = -c0*(ansatz20 + ansatz02);
      double ugradu = u1*ansatz10 + u2*ansatz01;
      
      // Galerkin part
      double val  = c0*(test10*ansatz10 + test01*ansatz01) // diffusion term
                     + ugradu*test00; // convective term
      // supg contribution
      val +=  (laplacian + ugradu ) *ugradv;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u1+res1)*test10;
      // grad div contribution
      sqMatrixA11[i][j] += Mult * (val + tau_c*test10*ansatz10); // A11 block
      
      val = tau_c * test10*ansatz01;
      // second cross term and subgrid term
      val += tau_m * ugradu*(u1+res1)*test01;
      sqMatrixA12[i][j] += Mult * val; // A12 block 
      
      val = tau_c * test01*ansatz10;
      // second cross term and subgrid term
      val += tau_m * ugradu * (u2+res2)*test10;
      sqMatrixA21[i][j] += Mult * val; // A21 block 
      
      // viscous term 
      val = c0*(test10*ansatz10 + test01*ansatz01)
            + ugradu * test00; 
      // supg contribution
      val += (laplacian + ugradu) * ugradv;
      // grad-div contribution
      val += tau_c*test01*ansatz01;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u2+res2)*test01;
      sqMatrixA22[i][j] += Mult * val;// + tau_c*test01*ansatz01); // A22 block     
      
      // weighted mass matrices
      val = ansatz00 * (test00 + ugradv);
      // second cross term and subgrid term
      val += tau_m * (u1 + res1) * ansatz00 * test10;
      sqMatrixM11[i][j] += Mult * val;
      
      val = tau_m * (u1 + res1) * ansatz00 * test01;
      sqMatrixM12[i][j] += Mult * val;
      
      val = tau_m * (u2 + res2) * ansatz00 * test10;
      sqMatrixM21[i][j] += Mult * val;
      
      val = ansatz00 * (test00 + ugradv);
      val += tau_m * (u2 + res2) * ansatz00 * test01;
      sqMatrixM22[i][j] += Mult * val;
    }// endfor j<N_U
    
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz functions
      ansatz00 = Orig3[j];
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      
      // pressure term
      double val = -ansatz00*test10;
      val += ansatz10 * ugradv;
      val += tau_m * (u1 + res1) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB1T[i][j] += Mult * val;
      
      val = -ansatz00*test01;
      val += ansatz01*ugradv;
      val += tau_m * (u2 + res2) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB2T[i][j] += Mult * val;
    }
  }
  
  // pressure test functions
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test00 = Orig3[i];

    // velocity-pressure block
    for(int j=0;j<N_U;j++)
    {
      // velocity ansatz functions
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      // divergence constraint
      reMatrixB1[i][j] +=Mult*(-test00*ansatz10);
      reMatrixB2[i][j] += Mult*(-test00*ansatz01);
    }    
  }
}

//================================================================================
void TimeNSType4NLResidual_VMS(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
{
  double **sqMatrixA11 = LocMatrices[0];
  double **sqMatrixA12 = LocMatrices[1];
  double **sqMatrixA21 = LocMatrices[2];
  double **sqMatrixA22 = LocMatrices[3];
  // weighted mass matrix
  double **sqMatrixM11 = LocMatrices[4];
  double **sqMatrixM12 = LocMatrices[5];
  double **sqMatrixM21 = LocMatrices[6];
  double **sqMatrixM22 = LocMatrices[7];
    
  double **reMatrixB1T= LocMatrices[8];
  double **reMatrixB2T= LocMatrices[9];
  
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u
  
  double *Orig3 = OrigValues[3]; // p
  double *Orig4 = OrigValues[4]; // p_x
  double *Orig5 = OrigValues[5]; // p_y
  
  double *Orig6 = OrigValues[6]; // u_xx
  double *Orig7 = OrigValues[7]; // u_yy
  
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  
  double c0=coeff[0];
  double f1=coeff[1];
  double f2=coeff[2];

  double u1=param[0];
  double u2=param[1];
  double u1_pts=param[2]; // 2 u1_m-1 - 1/2 u1_m-2
  double u2_pts=param[3]; // 2 u2_m-1 - 1/2 u2_m-2
  // solution from the previous time step and derivatives
  // this can also be constant or linearly extrapolated
  double u1x = param[4]; // u1x_old
  double u2x = param[5]; // u2x_old
  double u1y = param[6]; // u1y_old
  double u2y = param[7]; // u2y_old
  double u1xx= param[8]; // u1xx_old
  double u2xx= param[9]; // u2xx_old
  double u1yy= param[10]; // u1yy_old
  double u2yy= param[11]; // u2yy_old
  double px  = param[12]; // p_x
  double py  = param[13]; // p_y
  double p   = param[14]; // not needed
  // time derivatives
  double u1_t= param[15]; // u1t
  double u2_t= param[16]; // u2t 

  double test10, test01, test00;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;
  double tau_m, tau_c;
  tau_m =  hK*hK*TDatabase::ParamDB->DELTA0;
  tau_c =  TDatabase::ParamDB->DELTA1;
  
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + c0*(u1xx + u1yy) - (u1*u1x + u2*u1y) - px);
  res2 = tau_m * (f2-u2_t + c0*(u2xx + u2yy) - (u1*u2x + u2*u2y) - py);
  
  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugradv = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*(test00+ugradv)*f1;
    Rhs2[i] += Mult*(test00+ugradv)*f2;
    // contribution from the second and third nonlinear terms
    Rhs1[i] += Mult*tau_m*(u1+res1)*(f1*test10 + f2*test01);
    Rhs2[i] += Mult*tau_m*(u2+res2)*(f1*test10 + f2*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];
      
      double laplacian = -c0*(ansatz20 + ansatz02);
      double ugradu = u1*ansatz10 + u2*ansatz01;
      
      // Galerkin part
      double val  = c0*(test10*ansatz10 + test01*ansatz01) // diffusion term
                     + ugradu*test00; // convective term
      // supg contribution
      val +=  (laplacian + ugradu ) *ugradv;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u1+res1)*test10;
      // grad div contribution
      sqMatrixA11[i][j] += Mult * (val + tau_c*test10*ansatz10); // A11 block
      
      val = tau_c * test10*ansatz01;
      // second cross term and subgrid term
      val += tau_m * ugradu*(u1+res1)*test01;
      sqMatrixA12[i][j] += Mult * val; // A12 block 
      
      val = tau_c * test01*ansatz10;
      // second cross term and subgrid term
      val += tau_m * ugradu * (u2+res2)*test10;
      sqMatrixA21[i][j] += Mult * val; // A21 block 
      
      // viscous term 
      val = c0*(test10*ansatz10 + test01*ansatz01)
            + ugradu * test00; 
      // supg contribution
      val += (laplacian + ugradu) * ugradv;
      // grad-div contribution
      val += tau_c*test01*ansatz01;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u2+res2)*test01;
      sqMatrixA22[i][j] += Mult * val;// + tau_c*test01*ansatz01); // A22 block     
      
      // weighted mass matrices
      val = ansatz00 * (test00 + ugradv);
      // second cross term and subgrid term
      val += tau_m * (u1 + res1) * ansatz00 * test10;
      sqMatrixM11[i][j] += Mult * val;
      
      val = tau_m * (u1 + res1) * ansatz00 * test01;
      sqMatrixM12[i][j] += Mult * val;
      
      val = tau_m * (u2 + res2) * ansatz00 * test10;
      sqMatrixM21[i][j] += Mult * val;
      
      val = ansatz00 * (test00 + ugradv);
      val += tau_m * (u2 + res2) * ansatz00 * test01;
      sqMatrixM22[i][j] += Mult * val;
    }// endfor j<N_U
    
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz functions
      ansatz00 = Orig3[j];
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      
      // pressure term
      double val = -ansatz00*test10;
      val += ansatz10 * ugradv;
      val += tau_m * (u1 + res1) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB1T[i][j] += Mult * val;
      
      val = -ansatz00*test01;
      val += ansatz01*ugradv;
      val += tau_m * (u2 + res2) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB2T[i][j] += Mult * val;
    }
  }
}

//================================================================================
void TimeNSType4RHS_Residual_VMS(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
{
  
  // velocity part 
  double *rhs1 = LocRhs[0];
  double *rhs2 = LocRhs[1];
  
  double *Orig0 = OrigValues[0];
  double *Orig1 = OrigValues[1];
  double *Orig2 = OrigValues[2]; 
  
  double *Orig3 = OrigValues[3]; 
  double *Orig4 = OrigValues[4];
  
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  // problem parameters
  double nu=coeff[0];
  double f1=coeff[1];
  double f2=coeff[2];
  
  // solution from previous iterations
  // or the extrapolated solution from previous time step
  // in the case of extrapolations also in the nonlinear 
  // term on the left hand side of NS-equation
  double u1=param[0];
  double u2=param[1];
    // solution from previous time steps:
  // BDF1: u^n-1
  // BDF2: 2*u^n-1 - 0.5 u^n-2
  // double u1_pts = param[17];
  // double u2_pts = param[18];
  double u1_pts = param[2];
  double u2_pts = param[3];
  // solution from the previous time step and derivatives
  // this can also be constant or linearly extrapolated
  double u1x = param[4]; // u1x_old
  double u2x = param[5]; // u2x_old
  double u1y = param[6]; // u1y_old
  double u2y = param[7]; // u2y_old
  double u1xx= param[8]; // u1xx_old
  double u2xx= param[9]; // u2xx_old
  double u1yy= param[10]; // u1yy_old
  double u2yy= param[11]; // u2yy_old
  double px  = param[12]; // p_x
  double py  = param[13]; // p_y
  double p   = param[14]; // not needed
  // time derivatives
  double u1_t= param[15]; // u1t
  double u2_t= param[16]; // u2t 

  
  double test10, test01, test00, val;
  double tau_m =  hK*hK*TDatabase::ParamDB->DELTA0;
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1*u1x + u2*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1*u2x + u2*u2y) - py);

  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugradv = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    rhs1[i] += Mult*(test00+ugradv)*f1;
    rhs2[i] += Mult*(test00+ugradv)*f2;
    // contribution from the second and third nonlinear terms
    rhs1[i] += Mult*tau_m*(u1+res1)*(f1*test10 + f2*test01);
    rhs2[i] += Mult*tau_m*(u2+res2)*(f1*test10 + f2*test01);
  }  
}

//================================================================================
void TimeNSType14Residual_VMS(double Mult, double* coeff, double* param, double hK, 
  double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
{
  double **sqMatrixA11 = LocMatrices[0];
  double **sqMatrixA12 = LocMatrices[1];
  double **sqMatrixA21 = LocMatrices[2];
  double **sqMatrixA22 = LocMatrices[3];
  // weighted mass matrix
  double **sqMatrixM11 = LocMatrices[4];
  double **sqMatrixM12 = LocMatrices[5];
  double **sqMatrixM21 = LocMatrices[6];
  double **sqMatrixM22 = LocMatrices[7];
  
  // pressure pressure block
  double **sqMatrixC = LocMatrices[8];
  
  double **reMatrixB1 = LocMatrices[9];
  double **reMatrixB2 = LocMatrices[10];
  double **reMatrixB1T= LocMatrices[11];
  double **reMatrixB2T= LocMatrices[12];
  
  // velocity part 
  double *rhs1 = LocRhs[0];
  double *rhs2 = LocRhs[1];
  double *rhs3 = LocRhs[2];
  
  double *Orig0 = OrigValues[0];
  double *Orig1 = OrigValues[1];
  double *Orig2 = OrigValues[2]; 
  
  double *Orig3 = OrigValues[3]; 
  double *Orig4 = OrigValues[4];
  double *Orig5 = OrigValues[5];
  
  double *Orig6 = OrigValues[6];
  double *Orig7 = OrigValues[7];
  
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  // problem parameters
  double nu=coeff[0];
  double f1=coeff[1];
  double f2=coeff[2];
  
  double f1old = coeff[3];
  double f2old = coeff[4];
  
  // solution from previous iterations
  // or the extrapolated solution from previous time step
  // in the case of extrapolations also in the nonlinear 
  // term on the left hand side of NS-equation
  double u1=param[0]; // u1 
  double u2=param[1]; // u2
  double u1_pts=param[2]; // 2 u1_m-1 - 1/2 u1_m-2
  double u2_pts=param[3]; // 2 u2_m-1 - 1/2 u2_m-2
  // solution from the previous time step and derivatives
  // this can also be constant or linearly extrapolated
  double u1x = param[4]; // u1x_old
  double u2x = param[5]; // u2x_old
  double u1y = param[6]; // u1y_old
  double u2y = param[7]; // u2y_old
  double u1xx= param[8]; // u1xx_old
  double u2xx= param[9]; // u2xx_old
  double u1yy= param[10]; // u1yy_old
  double u2yy= param[11]; // u2yy_old
  double px  = param[12]; // p_x
  double py  = param[13]; // p_y
  double p   = param[14]; // not needed
  // time derivatives
  double u1_t= param[15]; // u1t
  double u2_t= param[16]; // u2t 
  
  
  double test10, test01, test00, val;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;
  
  // stabilization parameters
  double tau_m, tau_c;
  if(TDatabase::ParamDB->P1==100)
  {
    double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(dt*dt);
    double tmp = pow(0.5*hK,4);
    tau_m += 32.*nu*nu/tmp;
    tmp = 4.*(param[0]*param[0] + param[1]*param[1])/coeff[19];
    tau_m += tmp/(hK*hK/4.);
    
    tau_m = 1./sqrt(tau_m);
    tau_c = (hK*hK/4.0) / (8.*tau_m);
  }
  else if(TDatabase::ParamDB->P1==200)
  {
    double t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(t*t);
    double u = u1*u1+u2*u2;
    tau_m += 4.*u/(hK*hK);
    double hp4 = pow(hK,4);
    tau_m += 32.*nu*nu/hp4;
    tau_m = 1./sqrt(tau_m);
    tau_c = hK*hK/(8.*tau_m);
  }
  else
  {
    double stab_param[2];  
    stabilization_parameters_equal_order(Mult, param, coeff, stab_param);
    tau_m = stab_param[0];
    tau_c = stab_param[1];
  }
  TDatabase::ParamDB->P14 = tau_m;
  TDatabase::ParamDB->P15 = tau_c;
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1*u1x + u2*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1*u2x + u2*u2y) - py);

  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugradv = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    rhs1[i] += Mult*(test00+ugradv)*f1;
    rhs2[i] += Mult*(test00+ugradv)*f2;
    // contribution from the second and third nonlinear terms
    rhs1[i] += Mult*tau_m*(u1+res1)*(f1*test10 + f2*test01);
    rhs2[i] += Mult*tau_m*(u2+res2)*(f1*test10 + f2*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];
      
      double ugradu = u1*ansatz10 + u2*ansatz01;
      double laplacian = -nu*(ansatz20 + ansatz02);
      // viscous and convection term
      val = nu*(test10*ansatz10 + test01*ansatz01) + ugradu * test00; 
      // supg contribution
      val += (laplacian + ugradu) * ugradv;
      // grad-div contribution
      val += tau_c*test10*ansatz10;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u1+res1)*test10;
      sqMatrixA11[i][j] += Mult * val;
      
      val = tau_c * test10*ansatz01;
      // second cross term and subgrid term
      val += tau_m * ugradu*(u1+res1)*test01;
      sqMatrixA12[i][j] += Mult * val; 
      
      val = tau_c * test01*ansatz10;
      // second cross term and subgrid term
      val += tau_m * ugradu * (u2+res2)*test10;
      sqMatrixA21[i][j] += Mult * val; 
      
      // viscous term 
      val = nu*(test10*ansatz10 + test01*ansatz01)
            + ugradu * test00; 
      // supg contribution
      val += (-nu*(ansatz20 + ansatz02) + ugradu) * ugradv;
      // grad-div contribution
      val += tau_c*test01*ansatz01;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u2+res2)*test01;
      sqMatrixA22[i][j] += Mult * val;
      
      // mass matrix 
      val = ansatz00 * (test00 + ugradv);
      // second cross term and subgrid term
      val += tau_m * (u1 + res1) * ansatz00 * test10;
      sqMatrixM11[i][j] += Mult * val;
      
      val = tau_m * (u1 + res1) * ansatz00 * test01;
      sqMatrixM12[i][j] += Mult * val; 
      
      val = tau_m * (u2 + res2) * ansatz00 * test10;
      sqMatrixM21[i][j] += Mult * val; 
      
      
      val = ansatz00 * (test00 + ugradv);
      val += tau_m * (u2 + res2) * ansatz00 * test01;
      sqMatrixM22[i][j] += Mult * val; 
    }
    
    // pressure (ansatz) velocity(test) blocks
    for(int j=0; j<N_P; ++j)
    {
      ansatz00 = Orig3[j];
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      
      val = -ansatz00*test10;
      val += ansatz10 * ugradv;
      val += tau_m * (u1 + res1) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB1T[i][j] += Mult * val;
      
      val = -ansatz00*test01;
      val += ansatz01*ugradv;
      val += tau_m * (u2 + res2) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB2T[i][j] += Mult * val;
    }
  }
  
  double dt=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double t1 = TDatabase::TimeDB->THETA1;
  // velocity(ansatz) pressure (test) blocks
  for(int i=0; i<N_P; ++i)
  {
    test00 = Orig3[i];
    test10 = Orig4[i];
    test01 = Orig5[i];
    
    // u*_combined_old_time = (2 * u_{n-1} - 0.5 * u_{n-2} ) / dt
    rhs3[i] += Mult*tau_m*((u1_pts/dt + f1 )*test10 + 
                           (u2_pts/dt + f2 )*test01);
    for(int j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];

      val = -test00*ansatz10;
      // supg terms 
      val -=  tau_m * (ansatz00/(t1*dt) - nu*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test10;
      reMatrixB1[i][j] -= Mult * val;
      
      val = -test00*ansatz01;
      val -=  tau_m * (ansatz00/(t1*dt) - nu*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test01;
      reMatrixB2[i][j] -= Mult * val;
    }
    
    for(int j=0; j<N_P; ++j)
    {
      ansatz10=Orig4[j];
      ansatz01=Orig5[j];
      
      double val = tau_m * (ansatz10*test10+ansatz01*test01);
      sqMatrixC[i][j] += Mult * val;
    }
  }
}

//================================================================================
void TimeNSType14NLResidual_VMS(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
{
  double **sqMatrixA11 = LocMatrices[0];
  double **sqMatrixA12 = LocMatrices[1];
  double **sqMatrixA21 = LocMatrices[2];
  double **sqMatrixA22 = LocMatrices[3];
  
  double **reMatrixB1 = LocMatrices[4];
  double **reMatrixB2 = LocMatrices[5];
  
  double *Orig0 = OrigValues[0];
  double *Orig1 = OrigValues[1];
  double *Orig2 = OrigValues[2]; 
  
  double *Orig3 = OrigValues[3]; 
  double *Orig4 = OrigValues[4];
  double *Orig5 = OrigValues[5];
  
  double *Orig6 = OrigValues[6];
  double *Orig7 = OrigValues[7];
  
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  // problem parameters
  double nu=coeff[0];
  double f1=coeff[1];
  double f2=coeff[2];
  
  double f1old = coeff[3];
  double f2old = coeff[4];
  
  // solution from previous iterations
  // or the extrapolated solution from previous time step
  // in the case of extrapolations also in the nonlinear 
  // term on the left hand side of NS-equation
  double u1=param[0];
  double u2=param[1];
  // solution from the previous time step and derivatives
  // this can also be constant or linearly extrapolated
  double u1sigma=param[2];
  double u2sigma=param[3];
  double u1x = param[4]; // u1x_old
  double u2x = param[5]; // u2x_old
  double u1y = param[6]; // u1y_old
  double u2y = param[7]; // u2y_old
  double u1xx= param[8]; // u1xx_old
  double u2xx= param[9]; // u2xx_old
  double u1yy= param[10]; // u1yy_old
  double u2yy= param[11]; // u2yy_old
  double px  = param[12]; // p_x
  double py  = param[13]; // p_y
  double p   = param[14]; // not needed
  // time derivatives
  double u1_t= param[15]; // u1t
  double u2_t= param[16]; // u2t 
  // solution from previous time steps:
  // BDF1: u^n-1
  // BDF2: 2*u^n-1 - 0.5 u^n-2
  // double u1_pts = param[17]; // not needed
  // double u2_pts = param[18];
  
  double test10, test01, test00, val;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;

  // stabilization parameters
  double tau_m, tau_c;
  if(TDatabase::ParamDB->P1==100)
  {
    double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(dt*dt);
    double tmp = pow(0.5*hK,4);
    tau_m += 32.*nu*nu/tmp;
    tmp = 4.*(param[0]*param[0] + param[1]*param[1])/coeff[19];
    tau_m += tmp/(hK*hK/4.);
    
    tau_m = 1./sqrt(tau_m);
    tau_c = (hK*hK/4.0) / (8.*tau_m);
  }
  else if(TDatabase::ParamDB->P1==200)
  {
    double t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(t*t);
    double u = u1*u1+u2*u2;
    tau_m += 4.*u/(hK*hK);
    double hp4 = pow(hK,4);
    tau_m += 32.*nu*nu/hp4;
    tau_m = 1./sqrt(tau_m);
    tau_c = hK*hK/(8.*tau_m);
  }
  else
  {
    double stab_param[2];  
    stabilization_parameters_equal_order(Mult, param, coeff, stab_param);
    tau_m = stab_param[0];
    tau_c = stab_param[1];
  }
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1sigma*u1x + u2sigma*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1sigma*u2x + u2sigma*u2y) - py);
  
  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugradv = tau_m * (u1sigma*test10+u2sigma*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];
      
      double ugradu = u1*ansatz10 + u2*ansatz01;
      double laplacian = -nu*(ansatz20 + ansatz02);
      // viscous and convection term
      val = nu*(test10*ansatz10 + test01*ansatz01) + ugradu * test00; 
      // supg contribution
      val += (laplacian + ugradu) * ugradv;
      // grad-div contribution
      val += tau_c*test10*ansatz10;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u1sigma+res1)*test10;
      sqMatrixA11[i][j] += Mult * val;
      
      val = tau_c * test10*ansatz01;
      // second cross term and subgrid term
      val += tau_m * ugradu*(u1sigma+res1)*test01;
      sqMatrixA12[i][j] += Mult * val; 
      
      val = tau_c * test01*ansatz10;
      // second cross term and subgrid term
      val += tau_m * ugradu * (u2sigma+res2)*test10;
      sqMatrixA21[i][j] += Mult * val; 
      
      // viscous term 
      val = nu*(test10*ansatz10 + test01*ansatz01)
            + ugradu * test00; 
      // supg contribution
      val += (-nu*(ansatz20 + ansatz02) + ugradu) * ugradv;
      // grad-div contribution
      val += tau_c*test01*ansatz01;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u2sigma+res2)*test01;
      sqMatrixA22[i][j] += Mult * val;
    }
  }
  
  double dt=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double t1 = TDatabase::TimeDB->THETA1;
  
  // velocity(ansatz) pressure (test) blocks
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test00 = Orig3[i];
    test10 = Orig4[i];
    test01 = Orig5[i];
    
    for(int j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];

      // divergence constraint
      val = -test00*ansatz10;
      // supg terms 
      val -=  tau_m * (ansatz00/(t1*dt) - nu*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test10;
      reMatrixB1[i][j] -= Mult * val;
      
      val = -test00*ansatz01;
      val -=  tau_m * (ansatz00/(t1*dt) - nu*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test01;                       
      reMatrixB2[i][j] -= Mult * val;
    }
  }
}

//================================================================================
void TimeNSType14RHS_Residual_VMS(double Mult, double* coeff, double* param, 
 double hK, double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, 
 double** LocRhs)
{
  // velocity part 
  double *rhs1 = LocRhs[0];
  double *rhs2 = LocRhs[1];
  double *rhs3 = LocRhs[2];
  
  double *Orig0 = OrigValues[0];
  double *Orig1 = OrigValues[1];
  double *Orig2 = OrigValues[2]; 
  
  double *Orig3 = OrigValues[3]; 
  double *Orig4 = OrigValues[4];
  
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  // problem parameters
  double nu=coeff[0];
  double f1=coeff[1];
  double f2=coeff[2];
  
  // solution from previous iterations
  // or the extrapolated solution from previous time step
  // in the case of extrapolations also in the nonlinear 
  // term on the left hand side of NS-equation
  double u1=param[0];
  double u2=param[1];
    // solution from previous time steps:
  // BDF1: u^n-1
  // BDF2: 2*u^n-1 - 0.5 u^n-2
  // double u1_pts = param[17];
  // double u2_pts = param[18];
  double u1_pts = param[2];
  double u2_pts = param[3];
  // solution from the previous time step and derivatives
  // this can also be constant or linearly extrapolated
  double u1x = param[4]; // u1x_old
  double u2x = param[5]; // u2x_old
  double u1y = param[6]; // u1y_old
  double u2y = param[7]; // u2y_old
  double u1xx= param[8]; // u1xx_old
  double u2xx= param[9]; // u2xx_old
  double u1yy= param[10]; // u1yy_old
  double u2yy= param[11]; // u2yy_old
  double px  = param[12]; // p_x
  double py  = param[13]; // p_y
  double p   = param[14]; // not needed
  // time derivatives
  double u1_t= param[15]; // u1t
  double u2_t= param[16]; // u2t 

  
  double test10, test01, test00, val;
  double tau_m, tau_c;
  if(TDatabase::ParamDB->P1==100)
  {
    double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(dt*dt);
    double tmp = pow(0.5*hK,4);
    tau_m += 32.*nu*nu/tmp;
    tmp = 4.*(param[0]*param[0] + param[1]*param[1])/coeff[19];
    tau_m += tmp/(hK*hK/4.);
    
    tau_m = 1./sqrt(tau_m);
    tau_c = (hK*hK/4.0) / (8.*tau_m);
  }
  else if(TDatabase::ParamDB->P1==200)
  {
    double t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(t*t);
    double u = u1*u1+u2*u2;
    tau_m += 4.*u/(hK*hK);
    double hp4 = pow(hK,4);
    tau_m += 32.*nu*nu/hp4;
    tau_m = 1./sqrt(tau_m);
    tau_c = hK*hK/(8.*tau_m);
  }
  else
  {
    double stab_param[2];  
    stabilization_parameters_equal_order(Mult, param, coeff, stab_param);
    tau_m = stab_param[0];
    tau_c = stab_param[1];
  }
  TDatabase::ParamDB->P14 = tau_m;
  TDatabase::ParamDB->P15 = tau_c;
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1*u1x + u2*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1*u2x + u2*u2y) - py);

  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugradv = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    rhs1[i] += Mult*(test00+ugradv)*f1;
    rhs2[i] += Mult*(test00+ugradv)*f2;
    // contribution from the second and third nonlinear terms
    rhs1[i] += Mult*tau_m*(u1+res1)*(f1*test10 + f2*test01);
    rhs2[i] += Mult*tau_m*(u2+res2)*(f1*test10 + f2*test01);
  }
  double dt=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double t1 = TDatabase::TimeDB->THETA1;
  // velocity(ansatz) pressure (test) blocks
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test10 = Orig3[i];
    test01 = Orig4[i];
    
    // u*_combined_old_time = (2 * u_{n-1} - 0.5 * u_{n-2} ) / dt
    rhs3[i] += Mult*tau_m*((u1_pts/dt + f1 )*test10 + 
                           (u2_pts/dt + f2 )*test01);
  }
}

//================================================================================
void TimeNSParams_Residual_VMS(double* in, double* out)
{
  out[0] = in[2]; // u1old
  out[1] = in[3]; // u2old
  // combined old time for bdf 2: this is used in the right hand side computation
  // for the right hand side computation
  out[2] = in[4]; // combined old sols u1
  out[3] = in[5]; // combined old sols u2

  out[4] = in[6]; // u1x
  out[5] = in[7]; // u2x
  out[6] = in[8]; // u1y
  out[7] = in[9]; // u2y
  out[8] = in[10]; //u1xx
  out[9] = in[11]; // u2xx
  out[10] = in[12]; // u1yy
  out[11] = in[13]; // u2yy
  out[12] = in[14]; // px
  out[13] = in[15]; // py
  out[14] = in[16]; // p
// time derivative from previous time   
  out[15] = in[17]; // u1 td
  out[16] = in[18]; // u2 td
}

#include <assemble_routine_tnse2D_RVMS.h>
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
  double *rhs1 = LocRhs[0];
  double *rhs2 = LocRhs[1];
  
  double *uv_xorig = OrigValues[0];
  double *uv_yorig = OrigValues[1];
  double *uv_orig = OrigValues[2]; 
  
  double *pq_orig = OrigValues[3]; 
  double *pq_xorig = OrigValues[4];
  double *pq_yorig = OrigValues[5];
  
  double *uxx_orig = OrigValues[6];
  double *uyy_orig = OrigValues[7];
  
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
  double u1_pts=param[2];
  double u2_pts=param[3];
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
  double tau_m =  hK*hK*TDatabase::ParamDB->DELTA0;
  double tau_c =  TDatabase::ParamDB->DELTA1;
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1_pts*u1x + u2_pts*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1_pts*u2x + u2_pts*u2y) - py);

  for(int i=0; i<N_U; ++i)
  {
    test10 = uv_xorig[i];
    test01 = uv_yorig[i];
    test00 = uv_orig[i];
    
    double tau_m_uptsgradv = tau_m * (u1_pts*test10+u2_pts*test01);
    // right hand side 
    // standard terms 
    rhs1[i] += Mult*(test00+tau_m_uptsgradv)*f1;
    rhs2[i] += Mult*(test00+tau_m_uptsgradv)*f2;
    // contribution from the second and third nonlinear terms
    rhs1[i] += Mult*tau_m*(u1_pts+res1)*(f1*test10 + f2*test01);
    rhs2[i] += Mult*tau_m*(u2_pts+res2)*(f1*test10 + f2*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];
      ansatz00 = uv_orig[j];
      
      ansatz20 = uxx_orig[j];
      ansatz02 = uyy_orig[j];
      
      double ugradu = u1*ansatz10 + u2*ansatz01;
      double laplacian = -nu*(ansatz20 + ansatz02);
      // viscous and convection term
      val = nu*(test10*ansatz10 + test01*ansatz01) + ugradu * test00; 
      // supg contribution
      val += (laplacian + ugradu) * tau_m_uptsgradv;
      // grad-div contribution
      val += tau_c*test10*ansatz10;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u1_pts+res1)*test10;
      sqMatrixA11[i][j] += Mult * val;
      
      val = tau_c * test10*ansatz01;
      // second cross term and subgrid term
      val += tau_m * ugradu*(u1_pts+res1)*test01;
      sqMatrixA12[i][j] += Mult * val; 
      
      val = tau_c * test01*ansatz10;
      // second cross term and subgrid term
      val += tau_m * ugradu * (u2_pts+res2)*test10;
      sqMatrixA21[i][j] += Mult * val; 
      
      // viscous term 
      val = nu*(test10*ansatz10 + test01*ansatz01)
            + ugradu * test00; 
      // supg contribution
      val += (-nu*(ansatz20 + ansatz02) + ugradu) * tau_m_uptsgradv;
      // grad-div contribution
      val += tau_c*test01*ansatz01;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u2_pts+res2)*test01;
      sqMatrixA22[i][j] += Mult * val;
      
      // mass matrix 
      val = ansatz00 * (test00 + tau_m_uptsgradv);
      // second cross term and subgrid term
      val += tau_m * (u1_pts + res1) * ansatz00 * test10;
      sqMatrixM11[i][j] += Mult * val;
      
      val = tau_m * (u1_pts + res1) * ansatz00 * test01;
      sqMatrixM12[i][j] += Mult * val; 
      
      val = tau_m * (u2_pts + res2) * ansatz00 * test10;
      sqMatrixM21[i][j] += Mult * val; 
      
      
      val = ansatz00 * (test00 + tau_m_uptsgradv);
      val += tau_m * (u2_pts + res2) * ansatz00 * test01;
      sqMatrixM22[i][j] += Mult * val; 
    }
    
    // pressure (ansatz) velocity(test) blocks
    for(int j=0; j<N_P; ++j)
    {
      ansatz00 = pq_orig[j];
      ansatz10 = pq_xorig[j];
      ansatz01 = pq_yorig[j];
      
      val = -ansatz00*test10;
      val += ansatz10 * tau_m_uptsgradv;
      val += tau_m * (u1_pts + res1) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB1T[i][j] += Mult * val;
      
      val = -ansatz00*test01;
      val += ansatz01*tau_m_uptsgradv;
      val += tau_m * (u2_pts + res2) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB2T[i][j] += Mult * val;
    }
  }
  
  // velocity(ansatz) pressure (test) blocks
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test00 = pq_orig[i];
    for(int j=0;j<N_U;j++)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];

      val = -test00*ansatz10;
      reMatrixB1[i][j] += Mult * val;
      
      val = -test00*ansatz01;
      reMatrixB2[i][j] += Mult * val;
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
  
  double *uv_xorig = OrigValues[0];
  double *uv_yorig = OrigValues[1];
  double *uv_orig = OrigValues[2]; 
  
  double *pq_orig = OrigValues[3]; 
  double *pq_xorig = OrigValues[4];
  double *pq_yorig = OrigValues[5];
  
  double *uxx_orig = OrigValues[6];
  double *uyy_orig = OrigValues[7];
  
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
  double u1_pts=param[2];
  double u2_pts=param[3];
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
  double tau_m =  hK*hK*TDatabase::ParamDB->DELTA0;
  double tau_c =  TDatabase::ParamDB->DELTA1;
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1_pts*u1x + u2_pts*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1_pts*u2x + u2_pts*u2y) - py);
  
  for(int i=0; i<N_U; ++i)
  {
    test10 = uv_xorig[i];
    test01 = uv_yorig[i];
    test00 = uv_orig[i];
    
    double tau_m_uptsgradv = tau_m * (u1_pts*test10+u2_pts*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];
      ansatz00 = uv_orig[j];
      
      ansatz20 = uxx_orig[j];
      ansatz02 = uyy_orig[j];
      
      double ugradu = u1*ansatz10 + u2*ansatz01;
      double laplacian = -nu*(ansatz20 + ansatz02);
      // viscous and convection term
      val = nu*(test10*ansatz10 + test01*ansatz01) + ugradu * test00; 
      // supg contribution
      val += (laplacian + ugradu) * tau_m_uptsgradv;
      // grad-div contribution
      val += tau_c*test10*ansatz10;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u1_pts+res1)*test10;
      sqMatrixA11[i][j] += Mult * val;
      
      val = tau_c * test10*ansatz01;
      // second cross term and subgrid term
      val += tau_m * ugradu*(u1_pts+res1)*test01;
      sqMatrixA12[i][j] += Mult * val; 
      
      val = tau_c * test01*ansatz10;
      // second cross term and subgrid term
      val += tau_m * ugradu * (u2_pts+res2)*test10;
      sqMatrixA21[i][j] += Mult * val; 
      
      // viscous term 
      val = nu*(test10*ansatz10 + test01*ansatz01)
            + ugradu * test00; 
      // supg contribution
      val += (-nu*(ansatz20 + ansatz02) + ugradu) * tau_m_uptsgradv;
      // grad-div contribution
      val += tau_c*test01*ansatz01;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u2_pts+res2)*test01;
      sqMatrixA22[i][j] += Mult * val;
    }
  }
}

//================================================================================
void TimeNSType4RHS_Residual_VMS(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
{
  // weighted mass matrix
  double **sqMatrixM11 = LocMatrices[0];
  double **sqMatrixM12 = LocMatrices[1];
  double **sqMatrixM21 = LocMatrices[2];
  double **sqMatrixM22 = LocMatrices[3];
  double **reMatrixB1  = LocMatrices[4];
  double **reMatrixB2  = LocMatrices[5];
  double **reMatrixB1T = LocMatrices[6];
  double **reMatrixB2T = LocMatrices[7];
  
  // velocity part 
  double *rhs1 = LocRhs[0];
  double *rhs2 = LocRhs[1];
  
  double *uv_xorig = OrigValues[0];
  double *uv_yorig = OrigValues[1];
  double *uv_orig = OrigValues[2]; 
  
  double *pq_orig = OrigValues[3]; 
  double *pq_xorig = OrigValues[4];
  double *pq_yorig = OrigValues[5];
  
  double *uxx_orig = OrigValues[6];
  double *uyy_orig = OrigValues[7];
  
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
  double u1_pts=param[2];
  double u2_pts=param[3];
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
  double tau_m =  hK*hK*TDatabase::ParamDB->DELTA0;
  double tau_c =  TDatabase::ParamDB->DELTA1;
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1_pts*u1x + u2_pts*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1_pts*u2x + u2_pts*u2y) - py);

  for(int i=0; i<N_U; ++i)
  {
    test10 = uv_xorig[i];
    test01 = uv_yorig[i];
    test00 = uv_orig[i];
    
    double tau_m_uptsgradv = tau_m * (u1_pts*test10+u2_pts*test01);
    // right hand side 
    // standard terms 
    rhs1[i] += Mult*(test00+tau_m_uptsgradv)*f1;
    rhs2[i] += Mult*(test00+tau_m_uptsgradv)*f2;
    // contribution from the second and third nonlinear terms
    rhs1[i] += Mult*tau_m*(u1_pts+res1)*(f1*test10 + f2*test01);
    rhs2[i] += Mult*tau_m*(u2_pts+res2)*(f1*test10 + f2*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];
      ansatz00 = uv_orig[j];
      
      ansatz20 = uxx_orig[j];
      ansatz02 = uyy_orig[j];
      
      double ugradu = u1*ansatz10 + u2*ansatz01;
      double laplacian = -nu*(ansatz20 + ansatz02);      
      // mass matrix 
      val = ansatz00 * (test00 + tau_m_uptsgradv);
      // second cross term and subgrid term
      val += tau_m * (u1_pts + res1) * ansatz00 * test10;
      sqMatrixM11[i][j] += Mult * val;
      
      val = tau_m * (u1_pts + res1) * ansatz00 * test01;
      sqMatrixM12[i][j] += Mult * val; 
      
      val = tau_m * (u2_pts + res2) * ansatz00 * test10;
      sqMatrixM21[i][j] += Mult * val; 
      
      
      val = ansatz00 * (test00 + tau_m_uptsgradv);
      val += tau_m * (u2_pts + res2) * ansatz00 * test01;
      sqMatrixM22[i][j] += Mult * val; 
    }
    
    // pressure (ansatz) velocity(test) blocks
    for(int j=0; j<N_P; ++j)
    {
      ansatz00 = pq_orig[j];
      ansatz10 = pq_xorig[j];
      ansatz01 = pq_yorig[j];
      
      val = -ansatz00*test10;
      val += ansatz10 * tau_m_uptsgradv;
      val += tau_m * (u1_pts + res1) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB1T[i][j] += Mult * val;
      
      val = -ansatz00*test01;
      val += ansatz01*tau_m_uptsgradv;
      val += tau_m * (u2_pts + res2) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB2T[i][j] += Mult * val;
    }
  }
  
  // velocity(ansatz) pressure (test) blocks
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test00 = pq_orig[i];
    for(int j=0;j<N_U;j++)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];

      val = -test00*ansatz10;
      reMatrixB1[i][j] += Mult * val;
      
      val = -test00*ansatz01;
      reMatrixB2[i][j] += Mult * val;
    }    
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
  
  double *uv_xorig = OrigValues[0];
  double *uv_yorig = OrigValues[1];
  double *uv_orig = OrigValues[2]; 
  
  double *pq_orig = OrigValues[3]; 
  double *pq_xorig = OrigValues[4];
  double *pq_yorig = OrigValues[5];
  
  double *uxx_orig = OrigValues[6];
  double *uyy_orig = OrigValues[7];
  
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
  double u1_extr=param[2];
  double u2_extr=param[3];
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
  double u1_pts = param[17];
  double u2_pts = param[18];
  
  
  double test10, test01, test00, val;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;
  
  // stabilization parameters
  double u[2];
  u[0]=u1_extr; u[1]=u2_extr;
  double stab_param[2];
  stabilization_parameters_equal_order(Mult, u, coeff, stab_param);
  double tau_m =  stab_param[0];
  double tau_c =  stab_param[1];
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1_extr*u1x + u2_extr*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1_extr*u2x + u2_extr*u2y) - py);

  for(int i=0; i<N_U; ++i)
  {
    test10 = uv_xorig[i];
    test01 = uv_yorig[i];
    test00 = uv_orig[i];
    
    double tau_m_uptsgradv = tau_m * (u1_extr*test10+u2_extr*test01);
    // right hand side 
    // standard terms 
    rhs1[i] += Mult*(test00+tau_m_uptsgradv)*f1;
    rhs2[i] += Mult*(test00+tau_m_uptsgradv)*f2;
    // contribution from the second and third nonlinear terms
    rhs1[i] += Mult*tau_m*(u1_extr+res1)*(f1*test10 + f2*test01);
    rhs2[i] += Mult*tau_m*(u2_extr+res2)*(f1*test10 + f2*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];
      ansatz00 = uv_orig[j];
      
      ansatz20 = uxx_orig[j];
      ansatz02 = uyy_orig[j];
      
      double ugradu = u1*ansatz10 + u2*ansatz01;
      double laplacian = -nu*(ansatz20 + ansatz02);
      // viscous and convection term
      val = nu*(test10*ansatz10 + test01*ansatz01) + ugradu * test00; 
      // supg contribution
      val += (laplacian + ugradu) * tau_m_uptsgradv;
      // grad-div contribution
      val += tau_c*test10*ansatz10;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u1_extr+res1)*test10;
      sqMatrixA11[i][j] += Mult * val;
      
      val = tau_c * test10*ansatz01;
      // second cross term and subgrid term
      val += tau_m * ugradu*(u1_extr+res1)*test01;
      sqMatrixA12[i][j] += Mult * val; 
      
      val = tau_c * test01*ansatz10;
      // second cross term and subgrid term
      val += tau_m * ugradu * (u2_extr+res2)*test10;
      sqMatrixA21[i][j] += Mult * val; 
      
      // viscous term 
      val = nu*(test10*ansatz10 + test01*ansatz01)
            + ugradu * test00; 
      // supg contribution
      val += (-nu*(ansatz20 + ansatz02) + ugradu) * tau_m_uptsgradv;
      // grad-div contribution
      val += tau_c*test01*ansatz01;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u2_extr+res2)*test01;
      sqMatrixA22[i][j] += Mult * val;
      
      // mass matrix 
      val = ansatz00 * (test00 + tau_m_uptsgradv);
      // second cross term and subgrid term
      val += tau_m * (u1_extr + res1) * ansatz00 * test10;
      sqMatrixM11[i][j] += Mult * val;
      
      val = tau_m * (u1_extr + res1) * ansatz00 * test01;
      sqMatrixM12[i][j] += Mult * val; 
      
      val = tau_m * (u2_extr + res2) * ansatz00 * test10;
      sqMatrixM21[i][j] += Mult * val; 
      
      
      val = ansatz00 * (test00 + tau_m_uptsgradv);
      val += tau_m * (u2_extr + res2) * ansatz00 * test01;
      sqMatrixM22[i][j] += Mult * val; 
    }
    
    // pressure (ansatz) velocity(test) blocks
    for(int j=0; j<N_P; ++j)
    {
      ansatz00 = pq_orig[j];
      ansatz10 = pq_xorig[j];
      ansatz01 = pq_yorig[j];
      
      val = -ansatz00*test10;
      val += ansatz10 * tau_m_uptsgradv;
      val += tau_m * (u1_extr + res1) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB1T[i][j] += Mult * val;
      
      val = -ansatz00*test01;
      val += ansatz01*tau_m_uptsgradv;
      val += tau_m * (u2_extr + res2) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB2T[i][j] += Mult * val;
    }
  }
  
  double dt=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double t1 = TDatabase::TimeDB->THETA1;
  // velocity(ansatz) pressure (test) blocks
  for(int i=0; i<N_P; ++i)
  {
    test00 = pq_orig[i];
    test10 = pq_xorig[i];
    test01 = pq_yorig[i];
    
    // u*_combined_old_time = (2 * u_{n-1} - 0.5 * u_{n-2} ) / dt
    rhs3[i] += Mult*tau_m*((u1_pts/dt + f1 )*test10 + 
                           (u2_pts/dt + f2 )*test01);
    for(int j=0;j<N_U;j++)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];
      ansatz00 = uv_orig[j];
      ansatz20 = uxx_orig[j];
      ansatz02 = uyy_orig[j];

      val = -test00*ansatz10;
      // supg terms 
      val -=  tau_m * (ansatz00/(t1*dt) - nu*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test10;
      reMatrixB1[i][j] += Mult * val;
      
      val = -test00*ansatz01;
      val -=  tau_m * (ansatz00/(t1*dt) - nu*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test01;
      reMatrixB2[i][j] += Mult * val;
    }
    
    for(int j=0; j<N_P; ++j)
    {
      ansatz10=pq_xorig[j];
      ansatz01=pq_yorig[j];
      
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
  
  double *uv_xorig = OrigValues[0];
  double *uv_yorig = OrigValues[1];
  double *uv_orig = OrigValues[2]; 
  
  double *pq_orig = OrigValues[3]; 
  double *pq_xorig = OrigValues[4];
  double *pq_yorig = OrigValues[5];
  
  double *uxx_orig = OrigValues[6];
  double *uyy_orig = OrigValues[7];
  
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
  double u1_extr=param[2];
  double u2_extr=param[3];
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
  double u[2];
  u[0]=u1_extr; u[1]=u2_extr;
  double stab_param[2];
  stabilization_parameters_equal_order(Mult, u, coeff, stab_param);
  double tau_m =  stab_param[0];
  double tau_c =  stab_param[1];  
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1_extr*u1x + u2_extr*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1_extr*u2x + u2_extr*u2y) - py);
  
  for(int i=0; i<N_U; ++i)
  {
    test10 = uv_xorig[i];
    test01 = uv_yorig[i];
    test00 = uv_orig[i];
    
    double tau_m_uptsgradv = tau_m * (u1_extr*test10+u2_extr*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];
      ansatz00 = uv_orig[j];
      
      ansatz20 = uxx_orig[j];
      ansatz02 = uyy_orig[j];
      
      double ugradu = u1*ansatz10 + u2*ansatz01;
      double laplacian = -nu*(ansatz20 + ansatz02);
      // viscous and convection term
      val = nu*(test10*ansatz10 + test01*ansatz01) + ugradu * test00; 
      // supg contribution
      val += (laplacian + ugradu) * tau_m_uptsgradv;
      // grad-div contribution
      val += tau_c*test10*ansatz10;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u1_extr+res1)*test10;
      sqMatrixA11[i][j] += Mult * val;
      
      val = tau_c * test10*ansatz01;
      // second cross term and subgrid term
      val += tau_m * ugradu*(u1_extr+res1)*test01;
      sqMatrixA12[i][j] += Mult * val; 
      
      val = tau_c * test01*ansatz10;
      // second cross term and subgrid term
      val += tau_m * ugradu * (u2_extr+res2)*test10;
      sqMatrixA21[i][j] += Mult * val; 
      
      // viscous term 
      val = nu*(test10*ansatz10 + test01*ansatz01)
            + ugradu * test00; 
      // supg contribution
      val += (-nu*(ansatz20 + ansatz02) + ugradu) * tau_m_uptsgradv;
      // grad-div contribution
      val += tau_c*test01*ansatz01;
      // second cross term and subgrid term
      val += tau_m * (laplacian + ugradu)*(u2_extr+res2)*test01;
      sqMatrixA22[i][j] += Mult * val;
    }
  }
  
  double dt=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double t1 = TDatabase::TimeDB->THETA1;
  
  // velocity(ansatz) pressure (test) blocks
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test00 = pq_orig[i];
    test10 = pq_xorig[i];
    test01 = pq_yorig[i];
    
    for(int j=0;j<N_U;j++)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];
      ansatz00 = uv_orig[j];
      ansatz20 = uxx_orig[j];
      ansatz02 = uyy_orig[j];

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
  // weighted mass matrix
  double **sqMatrixM11 = LocMatrices[0];
  double **sqMatrixM12 = LocMatrices[1];
  double **sqMatrixM21 = LocMatrices[2];
  double **sqMatrixM22 = LocMatrices[3];
  double **sqMatrixC   = LocMatrices[4];
  double **reMatrixB1  = LocMatrices[5];
  double **reMatrixB2  = LocMatrices[6];
  double **reMatrixB1T = LocMatrices[7];
  double **reMatrixB2T = LocMatrices[8];
  
  // velocity part 
  double *rhs1 = LocRhs[0];
  double *rhs2 = LocRhs[1];
  double *rhs3 = LocRhs[2];
  
  double *uv_xorig = OrigValues[0];
  double *uv_yorig = OrigValues[1];
  double *uv_orig = OrigValues[2]; 
  
  double *pq_orig = OrigValues[3]; 
  double *pq_xorig = OrigValues[4];
  double *pq_yorig = OrigValues[5];
  
  double *uxx_orig = OrigValues[6];
  double *uyy_orig = OrigValues[7];
  
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
  double u1_extr=param[2];
  double u2_extr=param[3];
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
  double u1_pts = param[17];
  double u2_pts = param[18];
  
  double test10, test01, test00, val;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;
  
  // stabilization parameters
  double u[2];
  u[0]=u1_extr; u[1]=u2_extr;
  double stab_param[2];
  stabilization_parameters_equal_order(Mult, u, coeff, stab_param);
  double tau_m =  stab_param[0];
  double tau_c =  stab_param[1];
  
  double res1, res2;
  res1 = tau_m * (f1-u1_t + nu*(u1xx + u1yy) - (u1_extr*u1x + u2_extr*u1y) - px);
  res2 = tau_m * (f2-u2_t + nu*(u2xx + u2yy) - (u1_extr*u2x + u2_extr*u2y) - py);

  for(int i=0; i<N_U; ++i)
  {
    test10 = uv_xorig[i];
    test01 = uv_yorig[i];
    test00 = uv_orig[i];
    
    double tau_m_uptsgradv = tau_m * (u1_extr*test10+u2_extr*test01);
    // right hand side 
    // standard terms 
    rhs1[i] += Mult*(test00+tau_m_uptsgradv)*f1;
    rhs2[i] += Mult*(test00+tau_m_uptsgradv)*f2;
    // contribution from the second and third nonlinear terms
    rhs1[i] += Mult*tau_m*(u1_extr+res1)*(f1*test10 + f2*test01);
    rhs2[i] += Mult*tau_m*(u2_extr+res2)*(f1*test10 + f2*test01);
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];
      ansatz00 = uv_orig[j];
      
      ansatz20 = uxx_orig[j];
      ansatz02 = uyy_orig[j];
      
      double ugradu = u1*ansatz10 + u2*ansatz01;
      double laplacian = -nu*(ansatz20 + ansatz02);      
      // mass matrix 
      val = ansatz00 * (test00 + tau_m_uptsgradv);
      // second cross term and subgrid term
      val += tau_m * (u1_extr + res1) * ansatz00 * test10;
      sqMatrixM11[i][j] += Mult * val;
      
      val = tau_m * (u1_extr + res1) * ansatz00 * test01;
      sqMatrixM12[i][j] += Mult * val; 
      
      val = tau_m * (u2_extr + res2) * ansatz00 * test10;
      sqMatrixM21[i][j] += Mult * val; 
      
      
      val = ansatz00 * (test00 + tau_m_uptsgradv);
      val += tau_m * (u2_extr + res2) * ansatz00 * test01;
      sqMatrixM22[i][j] += Mult * val; 
    }
    
    // pressure (ansatz) velocity(test) blocks
    for(int j=0; j<N_P; ++j)
    {
      ansatz00 = pq_orig[j];
      ansatz10 = pq_xorig[j];
      ansatz01 = pq_yorig[j];
      
      val = -ansatz00*test10;
      val += ansatz10 * tau_m_uptsgradv;
      val += tau_m * (u1_extr + res1) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB1T[i][j] += Mult * val;
      
      val = -ansatz00*test01;
      val += ansatz01*tau_m_uptsgradv;
      val += tau_m * (u2_extr + res2) * (ansatz10*test10 + ansatz01*test01);
      reMatrixB2T[i][j] += Mult * val;
    }
  }
  double dt=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double t1 = TDatabase::TimeDB->THETA1;
  // velocity(ansatz) pressure (test) blocks
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test00 = pq_orig[i];
    test10 = pq_xorig[i];
    test01 = pq_yorig[i];
    
    // u*_combined_old_time = (2 * u_{n-1} - 0.5 * u_{n-2} ) / dt
    rhs3[i] += Mult*tau_m*((u1_pts/dt + f1 )*test10 + 
                           (u2_pts/dt + f2 )*test01);
    for(int j=0;j<N_U;j++)
    {
      ansatz10 = uv_xorig[j];
      ansatz01 = uv_yorig[j];
      ansatz00 = uv_orig[j];
      ansatz20 = uxx_orig[j];
      ansatz02 = uyy_orig[j];

      val = -test00*ansatz10;
      // supg terms 
      val -=  tau_m * (ansatz00/(t1*dt) - nu*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test10;
      reMatrixB1[i][j] += Mult * val;
      
      val = -test00*ansatz01;
      val -=  tau_m * (ansatz00/(t1*dt) - nu*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test01;                       
      reMatrixB2[i][j] += Mult * val;
    }
    
    for(int j=0; j<N_P; ++j)
    {
      ansatz10=pq_xorig[j];
      ansatz01=pq_yorig[j];
      
      double val = tau_m * (ansatz10*test10+ansatz01*test01);
      sqMatrixC[i][j] += Mult * val;
    }
  }
}

//================================================================================
void TimeNSParams_Residual_VMS(double* in, double* out)
{
  out[0] = in[2]; // u1old
  out[1] = in[3]; // u2old
  // extrapolated solution
  out[2] = in[4]; 
  out[3] = in[5]; 
  // derivatives of the extrapolated solution
  out[4] = in[6];
  out[5] = in[7];
  out[6] = in[8];
  out[7] = in[9];
  out[8] = in[10];
  out[9] = in[11];
  out[10] = in[12];
  out[11] = in[13];
  out[12] = in[14];
  out[13] = in[15];
  out[14] = in[16];
// time derivative from previous time   
  out[15] = in[17];
  out[16] = in[18];
  // previous time sols
  out[17] = in[19];
  out[18] = in[20];
}
#include <TNSE2DSUPG.h>
#include <CommonRoutineTNSE2D.h>
#include <Database.h>


//================================================================================
void TimeNSType4SUPG(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA21 = LocMatrices[2];
  double **MatrixA22 = LocMatrices[3];
  // weighted mass matrix
  double **MatrixM11 = LocMatrices[4];
  double **MatrixM22 = LocMatrices[5];
  
  double **MatrixB1 = LocMatrices[6];
  double **MatrixB2 = LocMatrices[7];
  double **MatrixB1T= LocMatrices[8];
  double **MatrixB2T= LocMatrices[9];
  
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
  double c1=coeff[1];
  double c2=coeff[2];

  double u1=param[0];
  double u2=param[1];

  double test10, test01, test00;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;
  // stabilization parameters
  double tau_m, tau_c;
  
  tau_c = TDatabase::ParamDB->DELTA1;
  tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  
  double skew_11 = 0.;
  double skew_22 = 0.;
  double skew_12 = 0.;
  double skew_21 = 0.;
  for(int i=0; i<N_U; ++i)
  {
    double *MatrixRowA11 = MatrixA11[i];
    double *MatrixRowA12 = MatrixA12[i];
    double *MatrixRowA21 = MatrixA21[i];
    double *MatrixRowA22 = MatrixA22[i];
    
    double *MatrixRowM11 = MatrixM11[i];
    double *MatrixRowM22 = MatrixM22[i];
    
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugrad = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];
      
      // Galerkin part
      double val  = c0*(test10*ansatz10 + test01*ansatz01); // diffusion term
      // NSE_NONLINEAR_FORM: global database parameter is deleted
      // uncomment after getting the idea
      // if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
      val += (u1*ansatz10 + u2*ansatz01)*test00; // convective term
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
      // {
        // skew symmetric form 
        // val += 0.5*(u1*ansatz10 + u2*ansatz01)*test00;
        // val -= 0.5*(u1*test10 + u2*test01)*ansatz00;
      // }
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 2)
      // {
      //   // another skew symmetric form n(u,v,w) = (u.grad v, w) + 0.5 * (div u, v.w)
      //   val += (u1*ansatz10 + u2*ansatz01)*test00;
      //   skew_11 = 0.5*(ansatz10*u1*test00);
      //   skew_22 = 0.5*(ansatz01*u2*test00);
      //   skew_12 = 0.5*ansatz01*u1*test00;
      //   skew_21 = 0.5*ansatz10*u2*test00;
      // }
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 3)
      // {
      //   ErrThrow("EMEC should be implemented");
      // }
      // supg contribution
      val +=  (-c0*(ansatz20 + ansatz02) + (u1*ansatz10 + u2*ansatz01) ) *ugrad;
      // grad div contribution
      MatrixRowA11[j] += Mult * (val + skew_11+tau_c*test10*ansatz10); // A11 block
      MatrixRowA22[j] += Mult * (val + skew_22+tau_c*test01*ansatz01); // A22 block
      
      val = tau_c * test10*ansatz01 + skew_12;
      MatrixRowA12[j] += Mult * val; // A12 block
      
      val = tau_c * test01*ansatz10 + skew_21;
      MatrixRowA21[j] += Mult * val; // A21 block 
      
      // weighted mass matrix
      val = ansatz00 * (test00 + ugrad);
      MatrixRowM11[j] += Mult * val;
      MatrixRowM22[j] += Mult * val;
    }// endfor j<N_U
    
    double *MatrixRow1 = MatrixB1T[i];
    double *MatrixRow2 = MatrixB2T[i];
    
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz functions
      ansatz00 = Orig3[j];
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      
      // pressure term
      MatrixRow1[j] += Mult * (-ansatz00 * test10 + ansatz10 * ugrad);
      
      MatrixRow2[j] += Mult * (-ansatz00 * test01 + ansatz01 * ugrad);
    }
  }
  
  // pressure test functions
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test00 = Orig3[i];
    double *MatrixRow1 = MatrixB1[i];
    double *MatrixRow2 = MatrixB2[i];
    // velocity-pressure block
    for(int j=0;j<N_U;j++)
    {
      // velocity ansatz functions
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      // divergence constraint
      MatrixRow1[j] +=Mult*(-test00*ansatz10);
      MatrixRow2[j] += Mult*(-test00*ansatz01);
    }    
  }
}
//================================================================================
void TimeNSType4NLSUPG(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA21 = LocMatrices[2];
  double **MatrixA22 = LocMatrices[3];
  // weighted mass matrix
  double **MatrixM11 = LocMatrices[4];
  double **MatrixM22 = LocMatrices[5];
  
  double **MatrixB1T = LocMatrices[6];
  double **MatrixB2T = LocMatrices[7];
  
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
  double c1=coeff[1];
  double c2=coeff[2];

  double u1=param[0];
  double u2=param[1];

  double test10, test01, test00;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;
  
  double tau_c = TDatabase::ParamDB->DELTA1;
  double tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  double skew_11 = 0.;
  double skew_22 = 0.;
  double skew_12 = 0.;
  double skew_21 = 0.;
  for(int i=0; i<N_U; ++i)
  {
    double *MatrixRowA11 = MatrixA11[i];
    double *MatrixRowA12 = MatrixA12[i];
    double *MatrixRowA21 = MatrixA21[i];
    double *MatrixRowA22 = MatrixA22[i];
    
    double *MatrixRowM11 = MatrixM11[i];
    double *MatrixRowM22 = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugrad = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];
      
      // Galerkin part
      double val  = c0*(test10*ansatz10 + test01*ansatz01); // diffusion term
      // if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
      val += (u1*ansatz10 + u2*ansatz01)*test00; // convective term
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
      // {
      //   // skew symmetric form 
      //   val += 0.5*(u1*ansatz10 + u2*ansatz01)*test00;
      //   val -= 0.5*(u1*test10 + u2*test01)*ansatz00;
      // }
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 2)
      // {
      //   // another skew symmetric form n(u,v,w) = (u.grad v, w) + 0.5 * (div u, v.w)
      //   val += (u1*ansatz10 + u2*ansatz01)*test00;
      //   skew_11 = 0.5*(ansatz10*u1*test00);
      //   skew_22 = 0.5*(ansatz01*u2*test00);
      //   skew_12 = 0.5*ansatz01*u1*test00;
      //   skew_21 = 0.5*ansatz10*u2*test00;
      // }
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 3)
      // {
      //   ErrThrow("EMEC should be implemented");
      // }
      // supg contribution
      val +=  (-c0*(ansatz20 + ansatz02) + (u1*ansatz10 + u2*ansatz01) ) *ugrad;
      // grad div contribution
      MatrixRowA11[j] += Mult * (val + skew_11+tau_c*test10*ansatz10); // A11 block
      MatrixRowA22[j] += Mult * (val + skew_22+tau_c*test01*ansatz01); // A22 block
      
      val = tau_c * test10*ansatz01 + skew_12;
      MatrixRowA12[j] += Mult * val; // A12 block
      
      val = tau_c * test01*ansatz10 + skew_21;
      MatrixRowA21[j] += Mult * val; // A21 block 
      
      // weighted mass matrix
      val = ansatz00 * (test00 + ugrad);
      MatrixRowM11[j] += Mult * val;
      MatrixRowM22[j] += Mult * val;
    }// endfor j<N_U
    double *MatrixRow1 = MatrixB1T[i];
    double *MatrixRow2 = MatrixB2T[i];
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz functions
      ansatz00 = Orig3[j];
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      
      // pressure term
      MatrixRow1[j] += Mult * (-ansatz00 * test10 + ansatz10 * ugrad);
      
      MatrixRow2[j] += Mult * (-ansatz00 * test01 + ansatz01 * ugrad);
    }
  }
}
//================================================================================
void TimeNSType4RHSSUPG(double Mult, double* coeff, double* param, double hK, 
               double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, 
               double** LocRhs)
{  
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u
  
  double c1=coeff[1];
  double c2=coeff[2];

  double u1=param[0];
  double u2=param[1];

  double tau_m =  TDatabase::ParamDB->DELTA0*hK*hK;
  
  double test10, test01, test00;  
  int N_U = N_BaseFuncts[0];
  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugradv = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*c1*(test00 + ugradv);
    Rhs2[i] += Mult*c2*(test00 + ugradv);    
  }
}

//================================================================================
void TimeNSType14SUPG(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA21 = LocMatrices[2];
  double **MatrixA22 = LocMatrices[3];
  // weighted mass matrix
  double **MatrixM11 = LocMatrices[4];
  double **MatrixM22 = LocMatrices[5];
  // pressure pressure block
  double **MatrixC = LocMatrices[6];
  
  double **MatrixB1 = LocMatrices[7];
  double **MatrixB2 = LocMatrices[8];
  double **MatrixB1T= LocMatrices[9];
  double **MatrixB2T= LocMatrices[10];
  
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2]; // pressure part  
  
  
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
  double c1=coeff[1];
  double c2=coeff[2];

  double u1=param[0];
  double u2=param[1];
  
  double u1_combined_old_time = param[2];
  double u2_combined_old_time = param[3];

  double test10, test01, test00;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;
  
  // stabilization parameters
  double tau_m, tau_c;
  if(TDatabase::ParamDB->P1==100)
  {
    double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(dt*dt);
    double tmp = pow(0.5*hK,4);
    tau_m += 32.*c0*c0/tmp;
    tmp = 4.*(param[0]*param[0] + param[1]*param[1])/coeff[19];
    tau_m += tmp/(hK*hK/4.);
    
    tau_m = 1./sqrt(tau_m);
    tau_c = (hK*hK/4.0) / (8.*tau_m);
    // ErrThrow("file is not copied which computes the stabilization parameters");
  }
  else if(TDatabase::ParamDB->P1==200)
  {
    double t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(t*t);
    double u = u1*u1+u2*u2;
    tau_m += 4.*u/(hK*hK);
    double hp4 = pow(hK,4);
    tau_m += 32.*c0*c0/hp4;
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
  
  double skew_11 = 0.;
  double skew_22 = 0.;
  double skew_12 = 0.;
  double skew_21 = 0.;
  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugradv = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*(test00+ugradv)*c1;
    Rhs2[i] += Mult*(test00+ugradv)*c2;
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];
      
      // Galerkin part
      double val  = c0*(test10*ansatz10 + test01*ansatz01); // diffusion term
      // if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 0)
      val += (u1*ansatz10 + u2*ansatz01)*test00; // convective term
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 1)
      // {
      //   // skew symmetric form 
      //   val += 0.5*(u1*ansatz10 + u2*ansatz01)*test00;
      //   val -= 0.5*(u1*test10 + u2*test01)*ansatz00;
      // }
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 2)
      // {
      //   // another skew symmetric form n(u,v,w) = (u.grad v, w) + 0.5 * (div u, v.w)
      //   val += (u1*ansatz10 + u2*ansatz01)*test00;
      //   skew_11 = 0.5*(ansatz10*u1*test00);
      //   skew_22 = 0.5*(ansatz01*u2*test00);
      //   skew_12 = 0.5*ansatz01*u1*test00;
      //   skew_21 = 0.5*ansatz10*u2*test00;
      // }
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 3)
      // {
      //   ErrThrow("EMEC or other forms should be implemented");
      // }

      //               + (u1*ansatz10 + u2*ansatz01)*test00; // convective term
      // supg contribution
      val +=  (-c0*(ansatz20 + ansatz02) + (u1*ansatz10 + u2*ansatz01) ) *ugradv;
      // grad div contribution
      MatrixA11[i][j] += Mult * (val + skew_11 + tau_c*test10*ansatz10); // A11 block
      MatrixA22[i][j] += Mult * (val + skew_22 + tau_c*test01*ansatz01); // A22 block
      
      val = tau_c * test10*ansatz01 +skew_12;
      MatrixA12[i][j] += Mult * val; // A12 block
      
      val = tau_c * test01*ansatz10 +skew_21;
      MatrixA21[i][j] += Mult * val; // A21 block 
      
      // weighted mass matrix
      val = ansatz00 * (test00 + ugradv);
      MatrixM11[i][j] += Mult * val;     
      MatrixM22[i][j] += Mult * val;
    }// endfor j<N_U
    
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz functions
      ansatz00 = Orig3[j];
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      
      // pressure term
      MatrixB1T[i][j] += Mult * (-ansatz00 * test10 + ansatz10 * ugradv);
      
      MatrixB2T[i][j] += Mult * (-ansatz00 * test01 + ansatz01 * ugradv);
    }
  }
  
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 
  double t1 = TDatabase::TimeDB->THETA1;
  
  // pressure test functions
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test00 = Orig3[i];
    test10 = Orig4[i];
    test01 = Orig5[i];
    
    // u*_combined_old_time = (2 * u_{n-1} - 0.5 * u_{n-2} ) / dt
    Rhs3[i] += Mult*tau_m*((u1_combined_old_time/dt + c1 )*test10 + 
                           (u2_combined_old_time/dt + c2 )*test01);
    
    for(int j=0;j<N_U;j++)
    {
      // velocity ansatz functions
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];

      // divergence constraint
      double val = -test00*ansatz10; // galerkin part
      // supg terms 
      val -=  tau_m * (ansatz00/(t1*dt) - c0*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test10;
      MatrixB1[i][j] -= Mult*val;

      val = -test00*ansatz01;
      val -=  tau_m * (ansatz00/(t1*dt) - c0*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test01;

      MatrixB2[i][j] -= Mult*val;      
    }
    // pressure-pressure block
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz
      ansatz10=Orig4[j];
      ansatz01=Orig5[j];
      
      double val = tau_m * (ansatz10*test10+ansatz01*test01);
      MatrixC[i][j] += Mult * val;
    }
  }
}
//================================================================================
void TimeNSType14NLSUPG(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA21 = LocMatrices[2];
  double **MatrixA22 = LocMatrices[3];
  // weighted mass matrix
  double **MatrixM11 = LocMatrices[4];
  double **MatrixM22 = LocMatrices[5];
  
  // pressure pressure block
  double **MatrixC = LocMatrices[6];
  
  double **MatrixB1 = LocMatrices[7];
  double **MatrixB2 = LocMatrices[8];
  double **MatrixB1T= LocMatrices[9];
  double **MatrixB2T= LocMatrices[10];
  
  // velocity part 
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2]; // pressure part  
  
  
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
  double c1=coeff[1];
  double c2=coeff[2];

  double u1=param[0];
  double u2=param[1];
  
  double u1_combined_old_time = param[2];
  double u2_combined_old_time = param[3];

  double test10, test01, test00;  
  double ansatz10, ansatz01, ansatz00, ansatz20, ansatz02;
  
  // stabilization parameters
  double tau_m, tau_c;
  if(TDatabase::ParamDB->P1==100)
  {
    double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(dt*dt);
    double tmp = pow(0.5*hK,4);
    tau_m += 32.*c0*c0/tmp;
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
    tau_m += 32.*c0*c0/hp4;
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
  double skew_11 = 0.;
  double skew_22 = 0.;
  double skew_12 = 0.;
  double skew_21 = 0.;
  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugradv = tau_m * (u1*test10+u2*test01);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*(test00+ugradv)*c1;
    Rhs2[i] += Mult*(test00+ugradv)*c2;
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];
      
      // Galerkin part
      double val  = c0*(test10*ansatz10 + test01*ansatz01); // diffusion term
      // if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
        val += (u1*ansatz10 + u2*ansatz01)*test00; // convective term
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==1)
      // {
      //   // skew symmetric form 
      //   val += 0.5*(u1*ansatz10 + u2*ansatz01)*test00;
      //   val -= 0.5*(u1*test10 + u2*test01)*ansatz00;
      // }
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 2)
      // {
      //   // another skew symmetric form n(u,v,w) = (u.grad v, w) + 0.5 * (div u, v.w)
      //   val += (u1*ansatz10 + u2*ansatz01)*test00;
      //   skew_11 = 0.5*(ansatz10*u1*test00);
      //   skew_22 = 0.5*(ansatz01*u2*test00);
      //   skew_12 = 0.5*ansatz01*u1*test00;
      //   skew_21 = 0.5*ansatz10*u2*test00;
      // }
      // else if(TDatabase::ParamDB->NSE_NONLINEAR_FORM == 3)
      // {
      //   ErrThrow("EMEC or other forms should be implemented");
      // }

      //               + (u1*ansatz10 + u2*ansatz01)*test00; // convective term
      // supg contribution
      val +=  (-c0*(ansatz20 + ansatz02) + (u1*ansatz10 + u2*ansatz01) ) *ugradv;
      // grad div contribution
      MatrixA11[i][j] += Mult * (val + skew_11 + tau_c*test10*ansatz10); // A11 block
      MatrixA22[i][j] += Mult * (val + skew_22 + tau_c*test01*ansatz01); // A22 block
      
      val = skew_12 + tau_c * test10*ansatz01;
      MatrixA12[i][j] += Mult * val; // A12 block
      
      val = skew_21 + tau_c * test01*ansatz10;
      MatrixA21[i][j] += Mult * val; // A21 block 
      
      // weighted mass matrix
      val = ansatz00 * (test00 + ugradv);
      MatrixM11[i][j] += Mult * val;     
      MatrixM22[i][j] += Mult * val;
    }// endfor j<N_U
    
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz functions
      ansatz00 = Orig3[j];
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      
      // pressure term
      MatrixB1T[i][j] += Mult * (-ansatz00 * test10 + ansatz10 * ugradv);
      
      MatrixB2T[i][j] += Mult * (-ansatz00 * test01 + ansatz01 * ugradv);
    }
  }
  
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 
  double t1 = TDatabase::TimeDB->THETA1;
  // pressure test functions
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    test00 = Orig3[i];
    test10 = Orig4[i];
    test01 = Orig5[i];
    
    // u*_combined_old_time = (2 * u_{n-1} - 0.5 * u_{n-2} )
    Rhs3[i] += Mult*tau_m*((u1_combined_old_time/dt + c1 )*test10 + 
                           (u2_combined_old_time/dt + c2 )*test01  );
    
    for(int j=0;j<N_U;j++)
    {
      // velocity ansatz functions
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig7[j];

      // divergence constraint
      // scaling with factor*step_length is done within the main class
      double val = -test00*ansatz10; // galerkin part
      // supg terms 
      val -=  tau_m * (ansatz00/(t1*dt) - c0*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test10;
      MatrixB1[i][j] -= Mult*val;

      val = -test00*ansatz01;
      val -=  tau_m * (ansatz00/(t1*dt) - c0*(ansatz20+ansatz02) 
                       + (u1*ansatz10+u2*ansatz01) ) * test01;

      MatrixB2[i][j] -= Mult*val;      
    }
    // pressure-pressure block
    for(int j=0; j<N_P; ++j)
    {
      // pressure ansatz
      ansatz10=Orig4[j];
      ansatz01=Orig5[j];
      
      double val = tau_m * (ansatz10*test10+ansatz01*test01);
      MatrixC[i][j] += Mult * val;
    }
  }
}


//================================================================================
void TimeNSType14RHSSUPG(double Mult, double* coeff, double* param, double hK, 
               double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, 
               double** LocRhs)
{
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2]; // pressure part 
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u
  
  double *Orig3 = OrigValues[3]; // p_x
  double *Orig4 = OrigValues[4]; // p_y
  
  double c0=coeff[0];
  double c1=coeff[1];
  double c2=coeff[2];

  double u1sigma=param[0];
  double u2sigma=param[1];
  double u1_combined_old_time = param[2];
  double u2_combined_old_time = param[3];

  // stabilization parameters
  double tau_m, tau_c;
  if(TDatabase::ParamDB->P1==100)
  {
    double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(dt*dt);
    double tmp = pow(0.5*hK,4);
    tau_m += 32.*c0*c0/tmp;
    tmp = 4.*(param[0]*param[0] + param[1]*param[1])/coeff[19];
    tau_m += tmp/(hK*hK/4.);
    
    tau_m = 1./sqrt(tau_m);
    tau_c = (hK*hK/4.0) / (8.*tau_m);
  }
  else if(TDatabase::ParamDB->P1==200)
  {
    double t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    tau_m = 4./(t*t);
    double u = u1sigma*u1sigma+u2sigma*u2sigma;
    tau_m += 4.*u/(hK*hK);
    double hp4 = pow(hK,4);
    tau_m += 32.*c0*c0/hp4;
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

  double test10, test01, test00;  
  int N_U = N_BaseFuncts[0];
 
  for(int i=0; i<N_U; ++i)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    
    double ugrad = tau_m * (u1sigma*test10+u2sigma*test01);
    // right hand side 
    // standard terms 
    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;
  }
  
  int N_P = N_BaseFuncts[1];
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; // time step length
  for(int i=0; i<N_P; ++i)
  {
    test10 = Orig3[i];
    test01 = Orig4[i];
    Rhs3[i] += Mult*tau_m*((u1_combined_old_time/dt + c1 )*test10 + 
                           (u2_combined_old_time/dt + c2 )*test01  );
  }
}

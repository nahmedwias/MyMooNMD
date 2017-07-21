// ======================================================================
// @(#)TNSE2D_FixPo.C        1.3 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#include <Database.h>
#include <Convolution.h>
#include <MooNMD_Io.h>
#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!

#include <stdlib.h>

// ======================================================================
// compute turbulent viscosity for LES
// ======================================================================
double TurbulentViscosity(double delta, double* gradU, double* u, double* uConv)
{
  int nu_type = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
  double nu_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  int nu_tensor =  TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR;
  double nu_power, nu_sigma;
  double frobenius_norm_tensor,nu,a,b,c;

  // compute turbulent viscosity
  switch(nu_type)
  {
      // no turbulent viscosity
      case 0:
	  nu = 0;
	  return(nu);
	  break;
      case 5:
	  nu = nu_constant;
	  return(nu);
	  break;
	  // defect correction  
      case 6:
	  if (TDatabase::ParamDB->DEFECT_CORRECTION_TYPE)	  
	      nu = delta;
	  else
	      nu = 0;  // bwe
	  return(nu);
	  break;
  }

  // compute square of the Frobenius norm of the tensor
  // use deformation tensor
  if (nu_tensor==0)
  {
    // compute (grad(u)+grad(u)^T)/2
    a = gradU[0]+gradU[0];
    b = gradU[1]+gradU[2];
    c = gradU[3]+gradU[3];
    frobenius_norm_tensor = (a*a+ 2*b*b + c*c)/4.0;
  }
  // use grad u
  else
    frobenius_norm_tensor =  gradU[0]*gradU[0] + gradU[1]*gradU[1] +
      gradU[2]*gradU[2] + gradU[3]*gradU[3];

  // compute turbulent viscosity
  switch(nu_type)
  {
    case 1:                      // Smagorinsky
    case 17:                     // Smagorinsky
      nu = nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
      break;
    case 2:                      // p laplacian
      nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
      nu = nu_constant * delta * delta * pow(frobenius_norm_tensor, nu_power/2.0);
      nu = nu_constant * delta * delta * pow(frobenius_norm_tensor, nu_power-2.0);
      break;
    case 3:                      // Layton SIAM J. Sci. Comput. 1996
      // nu = nu_0 |ln(h)|^(-0.5*(p-1)) h^sigma \|\nabla u\|^(p-2)
      nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
      nu_sigma = TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
      nu = nu_constant * pow(delta,nu_sigma) *
        pow(frobenius_norm_tensor,(nu_power-2)/2.0)
        / pow(fabs(log(delta)),0.5*(nu_power-1));
      break;
    case 4:                      // \|u-g_\delta \ast u\|_2
      nu = nu_constant * delta * sqrt((u[0]-uConv[0])*(u[0]-uConv[0])+
        (u[1]-uConv[1])*(u[1]-uConv[1]));
      //   OutPut(u[0] << " " << uConv[0] << "        "<< endl);
      //    OutPut(u[1] << " " << uConv[1] << endl);
      break;
    default:
      OutPut("This type of turbulent viscosity is not implemented !!!" << endl);
      exit(4711);
      break;
  }
  return(nu);
}

/******************************************************************************/
//
// computation of SUPG parameter following 
// Bazilevs, Calo, Cottrell, Hughes, Reali, Scovazzi
//
/******************************************************************************/

void SUPG_Param2D(double Mult, double* u, double* coeff, double* params)
{
    double x0, x1, x2, y0, y1, y2, g11, g12, g22;
    double d11, d12, d21, d22, nu, tau_c, tau_m;
    double u1, u2, rec_detjk;
    double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;    
    double C_I = TDatabase::ParamDB->DELTA0;

    nu = coeff[0];        
    rec_detjk = coeff[19];
    rec_detjk = 1/rec_detjk;
    u1 = u[0];
    u2 = u[1];
    
    x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
    y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
    x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
    y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
    x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
    y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
    
    // triangle
    if (TDatabase::ParamDB->INTERNAL_VERTEX_X[3]== -4711)
    {
	d11 = (y2-y1) * rec_detjk;  //dxi/dx
	d12 = (x1-x2) * rec_detjk;  //dxi/dy
	d21 = (y0-y1) * rec_detjk;  //deta/dx
	d22 = (x1-x0) * rec_detjk;  //deta/dy
    }
    else
    {
	// quadrilateral
	d11 = (y2-y1) * 0.5 * rec_detjk;  //dxi/dx
	d12 = (x1-x2) * 0.5 * rec_detjk;  //dxi/dy
	d21 = (y0-y1) * 0.5 * rec_detjk;  //deta/dx
	d22 = (x1-x0) * 0.5 * rec_detjk;  //deta/dy
    }
	
    g11 = d11*d11 + d21*d21;
    g12 = d11*d12 + d21*d22;
    g22 = d12*d12 + d22*d22;

    tau_m = g11*g11 + 2*g12*g12 + g22*g22; // G:G
    tau_m *= C_I*nu*nu;
    tau_m +=  4/(time_step*time_step); 
    tau_m += u1 * (g11*u1+g12*u2) + u2*(g12*u1+g22*u2);
    tau_m = 1/sqrt(tau_m); // this is the parameter for the momentum equation
 
    tau_c = (d11+d21)*(d11+d21)+(d12+d22)*(d12+d22);
    tau_c *= tau_m;
    tau_c = 1/tau_c;

    params[0] = tau_m;
    params[1] = tau_c;

/*    delta = (d11*d11+d21*d21)*(d11*d11+d21*d21)+2*(d11*d12+d21*d22)*(d11*d12+d21*d22)+  // G:G
	(d12*d12+d22*d22)*(d12*d12+d22*d22);
    delta *= C_I*nu*nu;         
    delta += 4/(time_step*time_step);  
    delta += u1*u1*(d11*d11+d21*d21)+2*u1*u2*(d11*d12+d21*d22)+u2*u2*(d12*d12+d22*d22);  // uGu
    delta = 1/sqrt(delta); // this is the parameter for the momentum equation
    
    delta *= ((d11+d21)*(d11+d21)+(d12+d22)*(d12+d22));   // gg
    delta = 1/delta; // this is the parameter for the mass equation
*/    
}


// ======================================================================
// Type 1, Standard Galerkin
// Type 1, Coletti
// Type 1, GL00Convolution
// ======================================================================



// ======================================================================
// for Type 1 is no SDFEM available
// ======================================================================

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 1, Smagorinsky
// the nonlinear viscosity is treated implicitly
// ======================================================================
void TimeNSType1Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 1, Galdi-Layton 98 Model auxiliary problem
// ======================================================================
void TimeNSType1GL00AuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM, **AuxMatrix;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2, mu2, delta, mu;
  double u1, u2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[3];
  MatrixB2 = LocMatrices[4];
  AuxMatrix = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 1, for group fem (the matrices for the convective term)
// ======================================================================
void TimeNSType1GroupFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixCx, **MatrixCy;
  double val;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;

  MatrixCx = LocMatrices[0];
  MatrixCy = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  for(i=0;i<N_U;i++)
  {
    MatrixRow1 = MatrixCx[i];
    MatrixRow2 = MatrixCy[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = ansatz10*test00;
      MatrixRow1[j] += Mult * val;

      val = ansatz01*test00;
      MatrixRow2[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 2, Standard Galerkin
// Type 2, Coletti
// Type 2, GL00Convolution
// ======================================================================



// ======================================================================
// for Type 2 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void TimeNSType2Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  
  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void TimeNSType2Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 2, GL00AuxProblem
// ======================================================================
void TimeNSType2GL00AuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM, **AuxMatrix;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow, *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[3];
  MatrixB2 = LocMatrices[4];
  MatrixB1T = LocMatrices[5];
  MatrixB2T = LocMatrices[6];
  AuxMatrix = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    AuxMatrixRow = AuxMatrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// Type 3, Coletti, (grad u, grad v)
// Type 3, GL00Convolution, (grad u, grad v)`



// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// Type 3, Coletti, D(u):D(v)
// Type 3, GL00Convolution, D(u):D(v)


// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = Mult*c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(2* test10*ansatz10+ test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; //  **MatrixA21, **MatrixA12;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, viscosity, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, GL00AuxProblem (grad u, grad v)
// ======================================================================
void TimeNSType3GL00AuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **AuxMatrix;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *AuxMatrixRow;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double delta, mu, mu2;
  double u1, u2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[7];
  MatrixB2  = LocMatrices[8];
  AuxMatrix = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    AuxMatrixRow = AuxMatrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType3GL00AuxProblemDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22, **AuxMatrix;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, viscosity, val1, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[7];
  MatrixB2  = LocMatrices[8];
  AuxMatrix = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, VMSProjection, D(u):D(v)
// ======================================================================
void TimeNSType3VMSProjectionDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixL, **Matrix_tilde_G11 ;
  double **Matrix_tilde_G22, **Matrix_G11, **Matrix_G22;
  double *Rhs1, *Rhs2, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P, N_L;
  double c0, c1, c2;
  double u1, u2, mu, viscosity, delta;

  //cout << "vms";

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixL   = LocMatrices[6];
  MatrixB1  = LocMatrices[7];
  MatrixB2  = LocMatrices[8];
  Matrix_tilde_G11  = LocMatrices[9];
  Matrix_tilde_G22  = LocMatrices[10];
  Matrix_G11  = LocMatrices[11];
  Matrix_G22  = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // l

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row  = Matrix_tilde_G22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      Matrix11Row[j] -= Mult * 2*mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * 2*mu * ansatz00 * test01;
    }
  }

  for(i=0;i<N_L;i++)
  {
    Matrix11Row = Matrix_G11[i];
    Matrix22Row = Matrix_G22[i];
    test00 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      Matrix11Row[j] -= Mult * ansatz10 * test00;
      Matrix22Row[j] -= Mult * ansatz01 * test00;
    }
  }

  for(i=0;i<N_L;i++)
  {
    test00 = Orig4[i];
    MatrixRow1 = MatrixL[i];
    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      MatrixRow1[j] += Mult * ansatz00 * test00;
    }
  }
}


// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// Type 4, Coletti, (grad u, grad v)
// Type 4, GL00Convolution, (grad u, grad v)


// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// Type 4, Coletti, D(u):D(v)
// Type 4, GL00Convolution, D(u):D(v)

// ======================================================================
// for Type 4 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void TimeNSType4Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row; // *Matrix21Row, *Matrix12Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

//  u1 = param[0];                 // u1old
//  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void TimeNSType4UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

//  u1 = param[0];                 // u1old
//  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row; // *Matrix21Row, *Matrix12Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
//void TimeNSType4SmagorinskyDD(double Mult, double *coeff,
//double *param, double hK,
//double **OrigValues, int *N_BaseFuncts,
//double ***LocMatrices, double **LocRhs)
//{
//  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
//  double **MatrixM11, **MatrixM22;
//  double **MatrixB1, **MatrixB2;
//  double **MatrixB1T, **MatrixB2T;
//  double *Rhs1, *Rhs2, val;
//  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
//  double *MatrixM11Row, *MatrixM22Row;
//  double *MatrixRow1, *MatrixRow2;
//  double ansatz00, ansatz10, ansatz01;
//  double test00, test10, test01;
//  double *Orig0, *Orig1, *Orig2, *Orig3;
//  int i,j, N_U, N_P;
//  double c0, c1, c2;
//  double u1, u2, mu, delta;
//
//  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
//  MatrixA22 = LocMatrices[3];
//  MatrixM11 = LocMatrices[4];
//  MatrixM22 = LocMatrices[5];
//  MatrixB1 = LocMatrices[6];
//  MatrixB2 = LocMatrices[7];
//  MatrixB1T = LocMatrices[8];
//  MatrixB2T = LocMatrices[9];
//
//  Rhs1 = LocRhs[0];
//  Rhs2 = LocRhs[1];
//
//  N_U = N_BaseFuncts[0];
//  N_P = N_BaseFuncts[1];
//
//  Orig0 = OrigValues[0];         // u_x
//  Orig1 = OrigValues[1];         // u_y
//  Orig2 = OrigValues[2];         // u
//  Orig3 = OrigValues[3];         // p
//
//  c0 = coeff[0];                 // nu
//  c1 = coeff[1];                 // f1
//  c2 = coeff[2];                 // f2
//
//  u1 = param[0];                 // u1old
//  u2 = param[1];                 // u2old
//
//  delta =  CharacteristicFilterWidth(hK);
//  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
//  mu = mu/2.0;
//
//  for(i=0;i<N_U;i++)
//  {
//    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
//    Matrix22Row = MatrixA22[i];
//    MatrixM11Row  = MatrixM11[i];
//    MatrixM22Row  = MatrixM22[i];
//    test10 = Orig0[i];
//    test01 = Orig1[i];
//    test00 = Orig2[i];
//
//    Rhs1[i] += Mult*test00*c1;
//    Rhs2[i] += Mult*test00*c2;
//
//    for(j=0;j<N_U;j++)
//    {
//      ansatz10 = Orig0[j];
//      ansatz01 = Orig1[j];
//      ansatz00 = Orig2[j];
//
//      val  = (c0+mu)*(2*test10*ansatz10+test01*ansatz01);
//      val += (u1*ansatz10+u2*ansatz01)*test00;
//      Matrix11Row[j] += Mult * val;
//
//      val  = (c0+mu)*(test01*ansatz10);
//      Matrix12Row[j] += Mult * val;
//
//      val  = (c0+mu)*(test10*ansatz01);
//      Matrix21Row[j] += Mult * val;
//
//      val  = (c0+mu)*(test10*ansatz10+2*test01*ansatz01);
//      val += (u1*ansatz10+u2*ansatz01)*test00;
//      Matrix22Row[j] += Mult * val;
//
//      val = Mult*(ansatz00*test00);
//      MatrixM11Row[j] += val;
//      MatrixM22Row[j] += val;
//    }                            // endfor j
//
//    MatrixRow1 = MatrixB1T[i];
//    MatrixRow2 = MatrixB2T[i];
//    for(j=0;j<N_P;j++)
//    {
//      ansatz00 = Orig3[j];
//
//      val = -Mult*ansatz00*test10;
//      MatrixRow1[j] += val;
//
//      val = -Mult*ansatz00*test01;
//      MatrixRow2[j] += val;
//    }
//  }                              // endfor i
//
//  for(i=0;i<N_P;i++)
//  {
//    MatrixRow1 = MatrixB1[i];
//    MatrixRow2 = MatrixB2[i];
//
//    test00 = Orig3[i];
//
//    for(j=0;j<N_U;j++)
//    {
//      ansatz10 = Orig0[j];
//      ansatz01 = Orig1[j];
//
//      val = -Mult*test00*ansatz10;
//      MatrixRow1[j] += val;
//
//      val = -Mult*test00*ansatz01;
//      MatrixRow2[j] += val;
//    }                            // endfor j
//
//  }                              // endfor i
//}


// ======================================================================
// Type 4, GL00AuxProblem, (grad u, grad v)
// ======================================================================
void TimeNSType4GL00AuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **AuxMatrix;  //**MatrixA21, **MatrixA2;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row; // *Matrix21Row, *Matrix12Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[7];
  MatrixB2 = LocMatrices[8];
  MatrixB1T = LocMatrices[9];
  MatrixB2T = LocMatrices[10];
  AuxMatrix = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    AuxMatrixRow = AuxMatrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType4GL00AuxProblemDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22, **AuxMatrix;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *AuxMatrixRow;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[7];
  MatrixB2 = LocMatrices[8];
  MatrixB1T = LocMatrices[9];
  MatrixB2T = LocMatrices[10];
  AuxMatrix = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(2*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, VMSProjection, D(u):D(v)
// ======================================================================
void TimeNSType4VMSProjectionDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixL, **Matrix_tilde_G11 ;
  double **Matrix_tilde_G22, **Matrix_G11, **Matrix_G22;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P, N_L;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixL   = LocMatrices[6];
  MatrixB1 = LocMatrices[7];
  MatrixB2 = LocMatrices[8];
  MatrixB1T = LocMatrices[9];
  MatrixB2T = LocMatrices[10];
  Matrix_tilde_G11  = LocMatrices[11];
  Matrix_tilde_G22  = LocMatrices[12];
  Matrix_G11  = LocMatrices[13];
  Matrix_G22  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // l

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(2*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row  = Matrix_tilde_G22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      Matrix11Row[j] -= Mult * 2*mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * 2*mu * ansatz00 * test01;
    }
  }

  for(i=0;i<N_L;i++)
  {
    Matrix11Row = Matrix_G11[i];
    Matrix22Row = Matrix_G22[i];
    test00 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      Matrix11Row[j] -= Mult * ansatz10 * test00;
      Matrix22Row[j] -= Mult * ansatz01 * test00;
    }
  }

  for(i=0;i<N_L;i++)
  {
    test00 = Orig4[i];
    MatrixRow1 = MatrixL[i];
    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      MatrixRow1[j] += Mult * ansatz00 * test00;
    }
  }
}


// ======================================================================
// assemble matrix for auxiliary problem
// ======================================================================

void MatrixAuxiliaryProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *AuxMatrixRow, **AuxMatrix;;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double delta, mu2, val;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  // solution does not need to be convolved
  AuxMatrix = LocMatrices[0];
  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  // filter width
  delta =  CharacteristicFilterWidth(hK);
  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  for(i=0;i<N_U;i++)
  {
    AuxMatrixRow = AuxMatrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Assembling routine for all nonlinear matrices
// ======================================================================

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// Type 2, Standard Galerkin, only nonlinear part
// Type 1, Coletti, only nonlinear part
// Type 2, Coletti, only nonlinear part
// Type 1, GL00Convolution, only nonlinear part
// Type 2, GL00Convolution, only nonlinear part
// Type 1, GL00AuxProblem, only nonlinear part
// Type 2, GL00AuxProblem, only nonlinear part
// ======================================================================
void TimeNSType1_2NLGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// Type 2, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1_2NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j, N_U;
  double c0;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = c0*(test10*ansatz10+test01*ansatz01);

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// Type 2, Smagorinsky, only nonlinear part
// ======================================================================
void TimeNSType1_2NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 3, Coletti, (grad u, grad v), only nonlinear part
// Type 4, Coletti, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 3, GL00AuxProblem, (grad u, grad v), only nonlinear part
// Type 4, GL00AuxProblem, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// + convection with higher order velocity
// ======================================================================
void TimeNSType3_4NLGalerkin_VMS_1_DD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2, ho_u1, ho_u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  ho_u1 = param[0];              // higher order u1old
  ho_u2 = param[1];              // higher order u2old
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = ((u1+ho_u1)*ansatz10+(u2+ho_u2)*ansatz01)*test00;
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// Type 4, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3_4NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j, N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = Mult*c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// Type 4, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3_4NLUpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j,N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v), only nonlinear part
// Type 4, Smagorinsky, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
//void TimeNSType3_4NLSmagorinskyDD(double Mult, double *coeff,
//double *param, double hK,
//double **OrigValues, int *N_BaseFuncts,
//double ***LocMatrices, double **LocRhs)
//{
//  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
//  double val, val1;
//  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
//  double ansatz10, ansatz01;
//  double test00, test10, test01;
//  double *Orig0, *Orig1, *Orig2;
//  int i,j,N_U;
//  double c0, viscosity, delta;
//  double u1, u2, mu;
//  // cout << "Sma" << endl;
//  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
//  MatrixA22 = LocMatrices[3];
//
//  N_U = N_BaseFuncts[0];
//
//  Orig0 = OrigValues[0];         // u_x
//  Orig1 = OrigValues[1];         // u_y
//  Orig2 = OrigValues[2];         // u
//
//  c0 = coeff[0];                 // nu
//
//  u1 = param[0];                 // u1old
//  u2 = param[1];                 // u2old
//
//  delta =  CharacteristicFilterWidth(hK);
//  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);
//  mu = mu/2.0;
//  viscosity = mu+c0;
//
//  for(i=0;i<N_U;i++)
//  {
//    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
//    Matrix22Row = MatrixA22[i];
//    test10 = Orig0[i];
//    test01 = Orig1[i];
//    test00 = Orig2[i];
//
//    for(j=0;j<N_U;j++)
//    {
//      ansatz10 = Orig0[j];
//      ansatz01 = Orig1[j];
//
//      val1 = (u1*ansatz10+u2*ansatz01)*test00;
//      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
//      val += val1;
//      Matrix11Row[j] += Mult * val;
//
//      val  = viscosity*(test01*ansatz10);
//      Matrix12Row[j] += Mult * val;
//
//      val  = viscosity*(test10*ansatz01);
//      Matrix21Row[j] += Mult * val;
//
//      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
//      val += val1;
//      Matrix22Row[j] += Mult * val;
//
//    }                            // endfor j
//  }                              // endfor i
//}


// ======================================================================
// Type 3, VMSProjection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMSProjection, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLVMSProjectionDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double **Matrix_tilde_G11,  **Matrix_tilde_G22;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_L;
  double c0, viscosity, delta;
  double u1, u2, mu;
  // cout << "Sma" << endl;
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  Matrix_tilde_G11  = LocMatrices[4];
  Matrix_tilde_G22  = LocMatrices[5];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);
  mu = mu/2.0;
  viscosity = mu+c0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row = Matrix_tilde_G22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig3[j];
      Matrix11Row[j] -= Mult * 2*mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * 2*mu * ansatz00 * test01;
    }
  }
}


// ======================================================================
// ROSENBROCK
// ======================================================================
void TimeNSType1GalerkinRHS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
//  double ansatz00, ansatz10, ansatz01;
  double test00; // test10, test01;
  double *Orig2;
  int i,N_U;
  double c1, c2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

//  Orig0 = OrigValues[0];         // u_x
//  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

//  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
//    test10 = Orig0[i];
//    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
  }                              // endfor i
}


void TimeNSType1GalerkinJ(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0; // c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];
  //N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  //Orig3 = OrigValues[3]; // p

  c0 = coeff[0];                 // nu
//  c1 = coeff[1];                 // f1
//  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
//      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      //val += (ansatz00)
      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


void TimeNSType1GalerkinC(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
//  double ansatz00, ansatz10, ansatz01;
  double test00;
  double *Orig2;
  int i,N_U;
  double c3, c4;

  //cout << "TimeNSType1GalerkinC" << endl;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig2 = OrigValues[2];         // u

  c3 = coeff[3];                 // dot f1
  c4 = coeff[4];                 // dot f2

  for(i=0;i<N_U;i++)
  {
    // test10 = Orig0[i];
    // test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c3;
    Rhs2[i] += Mult*test00*c4;
  }
}


void TimeNSType3GalerkinJ(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;// N_P;
  double c0;
  double u1, u2, u1_x, u1_y, u2_x, u2_y;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];
//  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1_x = param[2];               // u1old
  u2_x = param[3];               // u2old
  u1_y = param[4];               // u1old
  u2_y = param[5];               // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += u1_x*ansatz00*test00;
      Matrix11Row[j] += Mult * val;

      val  = u1_y*ansatz00*test00;
      Matrix12Row[j] += Mult * val;

      val  = u2_x*ansatz00*test00;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += u2_y*ansatz00*test00;
      Matrix22Row[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Assembling routine for right-hand sides only
// ======================================================================

// ======================================================================
// right-hand side ONLY, for NSE
// ======================================================================
void TimeNSRHS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i, N_U;
  double c1, c2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " ";
  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY for auxiliary problem applied to velocity
// ======================================================================
void TimeNSRHSAuxProblemU(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i, N_U;
  double c1, c2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = param[0];                 // f1
  c2 = param[1];                 // f2

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " ";
  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY, Coletti model
// ======================================================================
void TimeNSRHSColetti(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test10, test01;
  double *Orig0, *Orig1; // *Orig2;
  int i, N_U;
//  double c1, c2;

  double delta, val1, mu1;
  double D1u1, D2u1, D1u2, D2u2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
//  Orig2 = OrigValues[2];         // u

//  c1 = coeff[1];                 // f1
//  c2 = coeff[2];                 // f2

  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // (delta^2)/(2 gamma)
  mu1 = 0.5*delta*delta/gamma;

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
//    test00 = Orig2[i];

    // LES term
    val1 =  Mult* mu1;
    Rhs1[i] += val1*( test10* (D1u1*D1u1 + D2u1*D2u1)
      +test01* (D1u1*D1u2 + D2u1*D2u2));
    Rhs2[i] += val1*( test10* (D1u2*D1u1 + D2u2*D2u1)
      +test01* (D1u2*D1u2 + D2u2*D2u2));
  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY, Galdi-Layton model with convolution
// right-hand side ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSRHSLESModel(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test10, test01;
  double *Orig0, *Orig1;
  int i, N_U;
  double delta, val1, mu1;
  double gdT11, gdT12, gdT22;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  gdT11 = param[0];              // gDeltaT11 or 11 - component of auxiliary problem
  gdT12 = param[1];              // gDeltaT12 or 12 - and 21 - component
  gdT22 = param[2];              // gDeltaT22 or 22 - component of auxiliary problem

  // filter width
  delta =  CharacteristicFilterWidth(hK);
  mu1 = 0.5*delta*delta/gamma;
  val1 =  Mult* mu1;

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];

    // LES term
    Rhs1[i] += val1*( test10* gdT11 + test01* gdT12);
    Rhs2[i] += val1*( test10* gdT12 + test01* gdT22);
  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSRHSGL00AuxProblemPaper2(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;
  double c1, c2;

  double delta, val1, mu1;
//  double D1u1, D2u1, D1u2, D2u2;
  double AuxProblem11, AuxProblem12, AuxProblem22;
  double AuxProblem11Exact_x, AuxProblem12Exact_x, AuxProblem12Exact_y;
  double AuxProblem22Exact_y;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

//  D1u1 = param[0];               // D1u1
//  D1u2 = param[1];               // D1u2;
//  D2u1 = param[2];               // D2u1
//  D2u2 = param[3];               // D2u2;
  AuxProblem11 = param[4];       // 11 - component of auxiliary problem
  AuxProblem12 = param[5];       // 12 - and 21 - component
  AuxProblem22 = param[6];       // 22 - component of auxiliary problem
  AuxProblem11Exact_x = param[7];// 11 - component of auxiliary problem
  AuxProblem12Exact_x = param[8];// 12 - and 21 - component
  AuxProblem12Exact_y = param[9];// 12 - and 21 - component
                                 // 22 - component of auxiliary problem
  AuxProblem22Exact_y = param[10];

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // turbulent viscosity
//  mu = TurbulentViscosity(delta,&param[0],&param[7],&param[9]);

  // (delta^2)/(2 gamma)
  mu1 = 0.5*delta*delta/gamma;

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    // divergence of term with exact solution
    val1 =  Mult*mu1;
    Rhs1[i] -= val1*test00*(AuxProblem11Exact_x+AuxProblem12Exact_y);
    Rhs2[i] -= val1*test00*(AuxProblem12Exact_x+AuxProblem22Exact_y);

    // explicit treatment of artificial viscosity
    //val =  Mult *mu;
    //Rhs1[i] -= val*(D1u1*test10 + D2u1*test01);
    //Rhs2[i] -= val*(D1u2*test10 + D2u2*test01);

    val1 =  Mult*mu1;
    Rhs1[i] += val1*( test10* AuxProblem11 + test01* AuxProblem12);
    Rhs2[i] += val1*( test10* AuxProblem12 + test01* AuxProblem22);

  }                              // endfor i
}


// ======================================================================
// right-hand side for auxiliary problem
// ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSGL00AuxProblemRHS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, val;
  double test00;
  double *Orig0;
  int i, N_U;
  double D1u1, D2u1, D1u2, D2u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  D1u1 = param[0];               // D1u1
  D1u2 = param[1];               // D1u2;
  D2u1 = param[2];               // D2u1
  D2u2 = param[3];               // D2u2;

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    val = Mult*test00;

    Rhs1[i] += val*(D1u1*D1u1+D2u1*D2u1);
    Rhs2[i] += val*(D1u1*D1u2+D2u1*D2u2);
    Rhs3[i] += val*(D1u2*D1u2+D2u2*D2u2);

  }                              // endfor i
}


// ======================================================================
// right-hand side for auxiliary problem
// ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSGL00AuxProblemRHSPaper2(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, val;
  double test00;
  double *Orig0;
  int i, N_U;
//  double D1u1, D2u1, D1u2, D2u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    val = Mult*test00;

    Rhs1[i] += val*coeff[0];
    Rhs2[i] += val*coeff[1];
    Rhs3[i] += val*coeff[2];

  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY, Smagorinsky Explicit
// ======================================================================
//void TimeNSRHSSmagorinskyExplicit(double Mult, double *coeff,
//double *param, double hK,
//double **OrigValues, int *N_BaseFuncts,
//double ***LocMatrices, double **LocRhs)
//{
//  double *Rhs1, *Rhs2, val;
//  double test00, test10, test01;
//  double *Orig0, *Orig1, *Orig2;
//  int i, N_U;
//  double c1, c2, delta;
//
//  double mu;
//  double D1u1, D2u1, D1u2, D2u2;
//
//  Rhs1 = LocRhs[0];
//  Rhs2 = LocRhs[1];
//
//  N_U = N_BaseFuncts[0];
//
//  Orig0 = OrigValues[0];         // u_x
//  Orig1 = OrigValues[1];         // u_y
//  Orig2 = OrigValues[2];         // u
//
//  c1 = coeff[1];                 // f1
//  c2 = coeff[2];                 // f2
//
//  D1u1 = param[2];               // D1u1
//  D1u2 = param[3];               // D1u2;
//  D2u1 = param[4];               // D2u1
//  D2u2 = param[5];               // D2u2;
//
//  // turbulent viscosity
//  delta =  CharacteristicFilterWidth(hK);
//  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
//
//  for(i=0;i<N_U;i++)
//  {
//    test10 = Orig0[i];
//    test01 = Orig1[i];
//    test00 = Orig2[i];
//
//    Rhs1[i] += Mult*test00*c1;
//    Rhs2[i] += Mult*test00*c2;
//
//    val =  Mult * mu;
//    Rhs1[i] -= val*(D1u1*test10 + D2u1*test01);
//    Rhs2[i] -= val*(D1u2*test10 + D2u2*test01);
//
//  }                              // endfor i
//}


// ======================================================================
// right-hand side for additional terms in rhs of small scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_SmallRhs2D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double u1, u2, ho_u1, ho_u2;
  double D1u1, D2u1, D1u2, D2u2;
  double ho_D1u1, ho_D2u1, ho_D1u2, ho_D2u2;
  double p, c0;
  int i, N_U;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  // large scales
  u1 = param[0];
  u2 = param[1];
  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;
  // small scales
  ho_u1 = param[6];
  ho_u2 = param[7];
  ho_D1u1 = param[8];            // D1u1
  ho_D1u2 = param[9];            // D1u2;
  ho_D2u1 = param[10];           // D2u1
  ho_D2u2 = param[11];           // D2u2;
  // pressure
  p = param[12];

  c0 = coeff[0];

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += 2*c0*(D1u1*test10+(D1u2+D2u1)*test01/4);
    Rhs1[i] += (u1*D1u1+u2*D2u1) * test00;
    Rhs1[i] += (u1*ho_D1u1+u2*ho_D2u1) * test00;
    Rhs1[i] += (ho_u1*D1u1+ho_u2*D2u1) * test00;
    Rhs1[i] -= p * test10;
    Rhs1[i] *= Mult;
    Rhs2[i] += 2*c0*(D2u2*test01+(D1u2+D2u1)*test10/4);
    Rhs2[i] += (u1*D1u2+u2*D2u2) * test00;
    Rhs2[i] += (u1*ho_D1u2+u2*ho_D2u2) * test00;
    Rhs2[i] += (ho_u1*D1u2+ho_u2*D2u2) * test00;
    Rhs2[i] -= p * test01;
    Rhs2[i] *= Mult;
  }                              // endfor i
}


// ======================================================================
// right-hand side for additional terms in rhs of large scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_Large_0_Rhs2D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  // OutPut("l");
  double *Rhs1, *Rhs2;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double u1, u2, ho_u1, ho_u2;
  double D1u1, D2u1, D1u2, D2u2;
  double ho_D1u1, ho_D2u1, ho_D1u2, ho_D2u2;
  double p, c0;
  int i, N_U;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  // large scales
  u1 = param[0];
  u2 = param[1];
  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;
  // small scales
  ho_u1 = param[6];
  ho_u2 = param[7];
  ho_D1u1 = param[8];            // D1u1
  ho_D1u2 = param[9];            // D1u2;
  ho_D2u1 = param[10];           // D2u1
  ho_D2u2 = param[11];           // D2u2;
  // pressure (small scales)
  p = param[12];
  //OutPut(p);
  c0 = coeff[0];

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += 2*c0*(ho_D1u1*test10+(ho_D1u2+ho_D2u1)*test01/4);
    Rhs1[i] += (ho_u1*ho_D1u1+ho_u2*ho_D2u1) * test00;
    Rhs1[i] += (u1*ho_D1u1+u2*ho_D2u1) * test00;
    Rhs1[i] += (ho_u1*D1u1+ho_u2*D2u1) * test00;
    Rhs1[i] -= p * test10;
    Rhs1[i] *= Mult;
    Rhs2[i] += 2*c0*(ho_D2u2*test01+(ho_D1u2+ho_D2u1)*test10/4);
    Rhs2[i] += (ho_u1*ho_D1u2+ho_u2*ho_D2u2) * test00;
    Rhs2[i] += (u1*ho_D1u2+u2*ho_D2u2) * test00;
    Rhs2[i] += (ho_u1*D1u2+ho_u2*D2u2) * test00;
    Rhs2[i] -= p * test01;
    Rhs2[i] *= Mult;
  }                              // endfor i
}


// ======================================================================
// right-hand side for additional terms in rhs of large scale systems
// for VMS, Variant 1
// ======================================================================
void TimeNS_VMS_Large_1_Rhs2D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double u1, u2, ho_u1, ho_u2;
//  double D1u1, D2u1, D1u2, D2u2;
  double ho_D1u1, ho_D2u1, ho_D1u2, ho_D2u2;
  double p, c0;
  int i, N_U;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  // large scales
  u1 = param[0];
  u2 = param[1];
//  D1u1 = param[2];               // D1u1
//  D1u2 = param[3];               // D1u2;
//  D2u1 = param[4];               // D2u1
//  D2u2 = param[5];               // D2u2;
  // small scales
  ho_u1 = param[6];
  ho_u2 = param[7];
  ho_D1u1 = param[8];            // D1u1
  ho_D1u2 = param[9];            // D1u2;
  ho_D2u1 = param[10];           // D2u1
  ho_D2u2 = param[11];           // D2u2;
  // pressure (small scales)
  p = param[12];
  //OutPut(p);
  c0 = coeff[0];

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += 2*c0*(ho_D1u1*test10+(ho_D1u2+ho_D2u1)*test01/4);
    Rhs1[i] += (ho_u1*ho_D1u1+ho_u2*ho_D2u1) * test00;
    Rhs1[i] += (u1*ho_D1u1+u2*ho_D2u1) * test00;
    //    Rhs1[i] += (ho_u1*D1u1+ho_u2*D2u1) * test00;
    Rhs1[i] -= p * test10;
    Rhs1[i] *= Mult;
    Rhs2[i] += 2*c0*(ho_D2u2*test01+(ho_D1u2+ho_D2u1)*test10/4);
    Rhs2[i] += (ho_u1*ho_D1u2+ho_u2*ho_D2u2) * test00;
    Rhs2[i] += (u1*ho_D1u2+u2*ho_D2u2) * test00;
    //    Rhs2[i] += (ho_u1*D1u2+ho_u2*D2u2) * test00;
    Rhs2[i] -= p * test01;
    Rhs2[i] *= Mult;
  }                              // endfor i
}

// ======================================================================
// right-hand side ONLY, defect correction type 1, u2
// ======================================================================
void TimeNSRHSDefectCorrectionU2(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;
  double c1, c2;

  double D1u1, D2u1, D1u2, D2u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;
  
  //cout << D1u1 << " " << D1u2  << " " << D2u1  << " " << D2u2 <<endl;
  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs1[i] += Mult*hK*(D1u1*test10 + D2u1*test01);
    Rhs2[i] += Mult*test00*c2;
    Rhs2[i] += Mult*hK*(D1u2*test10 + D2u2*test01);
  }                              // endfor i
}
// ======================================================================
// right-hand side ONLY, defect correction type 2, u2
// ======================================================================
void TimeNSRHSDefectCorrectionU2_1(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;
  double c1, c2;
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double u1, u2, um11, um12, um21, um22;
  double D1u1, D2u1, D1u2, D2u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  
  u1 = param[0];
  u2 = param[1];
  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;
  um11 = param[6];  
  um12 = param[7];
  um21 = param[8];
  um22 = param[9];

  //OutPut( u1<< " " << u2 << " " << um11 << " " << um12  << " " << um21  << " " << um22 <<endl);
  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs1[i] += Mult*hK*(D1u1*test10 + D2u1*test01);
    Rhs1[i] -= Mult *(u1-2*um11+um21)*test00/(2*dt);
    Rhs2[i] += Mult*test00*c2;
    Rhs2[i] += Mult*hK*(D1u2*test10 + D2u2*test01);
    Rhs2[i] -= Mult *(u2-2*um12+um22)*test00/(2*dt);
  }                              // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v) for Axialsymmetric
// ======================================================================
void TimeNSType4GalerkinDD_Axial3D(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, r, x, sign;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old - grid VX   !!check aux parameters
  u2 = param[1]; // u2old - grid VY   !!check aux parameters

  x  = param[2]; // x
  r  = fabs(x);

  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }

// exit(0);

  if (x>=0)
     sign = 1;
  else
     sign = -1;
//      sign = 1;    // unles all integral oints are inside it is not feasible
// if(fabs(u2)>1e-2)
// cout << " u1 " << u1<<  " u2 " << u2<< endl;
// cout << "Mult" << Mult<< endl;
// cout << "f1 :" << c1<<  " f2 :" << c2<< endl;
// cout << "c0 :" << c0<<  " Mult :" << Mult<< endl;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i]*sign;
    test01 = Orig1[i];
    test00 = Orig2[i];


     Rhs1[i] += Mult*test00*c1*r;
     Rhs2[i] += Mult*test00*c2*r;


    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j]*sign;
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*((2.*test10*ansatz10+test01*ansatz01)*r   + 2.*ansatz00*test00/r) ;
      val += (u1*ansatz10+u2*ansatz01)*test00*r;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10)*r;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01)*r;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2.*test01*ansatz01)*r;
      val += (u1*ansatz10+u2*ansatz01)*test00*r;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00)*r;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*(ansatz00*test10*r + ansatz00*test00);
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01*r;
      MatrixRow2[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j]*sign;
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = -Mult*(test00*ansatz10*r + ansatz00*test00);
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01*r;
      MatrixRow2[j] += val;
    } // endfor j

  } // endfor i
}

// ======================================================================
//  for axial symmetric case
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD_Axial3D(double Mult, double *coeff,
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2, r, x, sign;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  x  = param[2]; // x
  r  = fabs(x);

  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }

  if (x>=0)
     sign = 1;
  else
     sign = -1;
//      sign = 1;    // unles all integral oints are inside it is not feasible
 // cout << "u1" <<u1 << endl;
 // cout << "u2" <<u2 << endl;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i]*sign;
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j]*sign;
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00*r;
      val  = c0*((2.*test10*ansatz10+test01*ansatz01)*r + 2.*ansatz00*test00/r );
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01)*r;
      val += val1;
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 4, Two phase Standard Galerkin, D(u):D(v)
// Type 4, Two phase Coletti, D(u):D(v)
// Type 4, Two phase GL00Convolution, D(u):D(v)
// ======================================================================
// /*
void TimeNSType4GalerkinDD_2PhaseAxial3D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, x, r;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // 1/Re_k
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  c3 = coeff[3]; // density ratio (inner/outer)

  u1 = param[0]; // u1old - grid VX   !!check aux parameters
  u2 = param[1]; // u2old - grid VY   !!check aux parameters

  x  = param[2]; // x
  r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }
// if(fabs(u2)>1e-2)
// cout << " 2phase u1 " << u1<<  " u2 " << u2<< endl;
// cout << "Mult" << Mult<< endl;
// cout << "f1 :" << c1<<  " f2 :" << c2<< endl;
// cout << "c3 :" << c3 <<   endl;
// exit(0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1* c3*r;
    Rhs2[i] += Mult*test00*c2* c3*r;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*((2*test10*ansatz10+test01*ansatz01)*r   + 2.*ansatz00*test00/r);
      val += (u1*ansatz10+u2*ansatz01)*test00*c3*r;
      Matrix11Row[j] += Mult *val;

      val  = c0*(test01*ansatz10)*r;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01)*r;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01)*r;
      val += (u1*ansatz10+u2*ansatz01)*test00*c3*r;
      Matrix22Row[j] += Mult * val;

      val = Mult* c3*(ansatz00*test00)*r;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*(ansatz00*test10*r + ansatz00*test00);
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01*r;
      MatrixRow2[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = -Mult*(test00*ansatz10*r + ansatz00*test00);
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01*r;
      MatrixRow2[j] += val;
    } // endfor j

  } // endfor i
}


void TimeNSType3_4NLGalerkinDD_2PhaseAxial3D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0, c3;
  double u1, u2, x, r;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c3 = coeff[3]; // density ratio (inner/outer)

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  x  = param[2]; // x
  r  = fabs(x);

  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }

 // cout << "u1" <<u1 << endl;
 // cout << "u2" <<u2 << endl;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00*c3*r;
      val  = c0*((2.*test10*ansatz10+test01*ansatz01)*r + 2.*ansatz00*test00/r);
      val += val1;
      Matrix11Row[j] += Mult* val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01)*r;
      val += val1;
      Matrix22Row[j] += Mult* val;

    } // endfor j
  } // endfor i
}


// ======================================================================
// Standard Galerkin
// ======================================================================
void GridAssemble4(double Mult, double *coeff,
                  double *param, double hK,
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, lambda=0;
  double *MatrixRow11, *MatrixRow12, *MatrixRow21, *MatrixRow22;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j,N_U;
  double c0=1., detjk;

//   double lame = TDatabase::ParamDB->LameC;

  if(TDatabase::ParamDB->P0)
   { detjk = coeff[19]; } // see DiscreteForm2D.C    
  else
   { detjk = 1; }

//   varying stiffness and divergence term values
//   lame *= detjk;


// cout<< "lame " <<lame<<endl;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y

  for(i=0;i<N_U;i++)
  {
    MatrixRow11 = MatrixA11[i];
    MatrixRow12 = MatrixA12[i];
    MatrixRow21 = MatrixA21[i];
    MatrixRow22 = MatrixA22[i];

    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

     // New part D(v):D(w)
  /*
      val = (2*test10*ansatz10+test01*ansatz01) + lame*test10*ansatz10;
      MatrixRow11[j] += Mult * val/detjk;

      val  = test01*ansatz10/detjk;
      val += lame*test10*ansatz01;
      MatrixRow12[j] += Mult * val;

      val = test10*ansatz01/detjk;
      val += lame*test01*ansatz10;
      MatrixRow21[j] += Mult * val;

      val = (test10*ansatz10+2*test01*ansatz01) + lame*test01*ansatz01;
      MatrixRow22[j] += Mult * val/detjk;

   */
   // old part gradv:gradv

      val  = c0*(test10*ansatz10+test01*ansatz01 + lambda*test10*ansatz10);
      MatrixRow11[j] += Mult * val/detjk;
      val  = c0*(lambda*test10*ansatz01);
      MatrixRow12[j] += Mult * val/detjk;
      val  = c0*(lambda*test01*ansatz10);
      MatrixRow21[j] += Mult * val/detjk;
      val  = c0*(test10*ansatz10+test01*ansatz01 + lambda*test01*ansatz01);
      MatrixRow22[j] += Mult * val/detjk;

    } // endfor j
  } // endfor i
 }


// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// used for : GALERKIN
//            COLETTI (non linear steps)
// ========================================================================
void TimeNSParams2(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
}


// comment: routine can be replaced by  TimeNSParamsVelo_GradVelo !!!
void TimeNSParams2RB(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2
}


// ========================================================================
// parameters: u1old, u2old,
// used for : COLETTI, Smagorinsky
// ========================================================================
void TimeNSParamsVelo_GradVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2
}

// ========================================================================
// parameters: u1old, u2old,
// used for : COLETTI, Smagorinsky
// ========================================================================
void TimeNSParamsVelo_GradVelo_ALE(double *in, double *out)
{
  out[0] = in[2] - in[8];       // u1old - w1old
  out[1] = in[3] - in[9];       // u2old - w2old

  out[2] = in[4];               // D1u1
  out[3] = in[5];               // D1u2
  out[4] = in[6];               // D2u1
  out[5] = in[7];               // D2u2
  
//   cout <<in[8] << " E " << in[9] << endl;
}

// ========================================================================
// parameters: u1old, u2old,
// used for : COLETTI, Smagorinsky
// ========================================================================
void TimeNSParamsVelo_GradVelo_ConvVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  out[6] = in[8];                // u1old
  out[7] = in[9];                // u2old
}


// ========================================================================
// parameters: u1old, u2old, Frobenius Norm(grad(u))
// used for : COLETTI, first iteration step, viscosity type = 4
// ========================================================================
void TimeNSParamsVelo_GradVeloNuT4(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  // convolution of the solution
  out[6] = in[8];                // g_\delta \ast u1
  out[7] = in[9];                // g_\delta \ast u2
}


// ========================================================================
// used for : GL00Convolution, first assembling
// ========================================================================
void TimeNSParamsGL00Convolution(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  // grad u
  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  // components of convoluted tensor
  out[6] = in[ 8];               // g_d\ast(D1u1*D1u1+D1u2*D1u2)
  out[7] = in[ 9];               // g_d\ast(D1u1*D1u2+D2u1*D2u2)
  out[8] = in[10];               // g_d\ast(D2u1*D2u1+D2u2*D2u2)

}


// ========================================================================
// used for : GL00Convolution, rhs assembling, viscosity type = 4
// ========================================================================
void TimeNSParamsRHSGL00ConvolutionNuT4(double *in, double *out)
{
  // grad u
  out[0] = in[4];                // D1u1
  out[1] = in[5];                // D1u2
  out[2] = in[6];                // D2u1
  out[3] = in[7];                // D2u2

  // components of convoluted tensor
  out[4] = in[ 8];               // g_d\ast(D1u1*D1u1+D1u2*D1u2)
  out[5] = in[ 9];               // g_d\ast(D1u1*D1u2+D2u1*D2u2)
  out[6] = in[10];               // g_d\ast(D2u1*D2u1+D2u2*D2u2)

  // velocity
  out[7] = in[2];                // u1
  out[8] = in[3];                // u2

  // convoluted velocity
  out[9] = in[11];               //  g_d\ast u1
  out[10] = in[12];              //  g_d\ast u2

}


// ========================================================================
// used for assembling of term coming from the LES model
// GL00AuxProb and GL00Conv
// parameters used in TimeNSRHSLESModel
// ========================================================================
void TimeNSParamsRHSLES(double *in, double *out)
{
  //  solution of auxiliary problem
  out[0] = in[2];
  out[1] = in[3];
  out[2] = in[4];
}


// ========================================================================
// parameters:
// used for : GL00AuxProblem, viscosity type = 4
// ========================================================================
void TimeNSParamsGL00AuxProblemNuT4(double *in, double *out)
{

  // \nabla u
  out[0] = in[4];                // D1u1
  out[1] = in[5];                // D1u2
  out[2] = in[6];                // D2u1
  out[3] = in[7];                // D2u2

  //  solution of auxiliary problem
  out[4] = in[ 8];               // sol11
  out[5] = in[ 9];               // sol12 = sol21
  out[6] = in[10];               // sol22

  // solution
  out[7] = in[2];                // u1
  out[8] = in[3];                // u2

  // convolution of the solution
  out[9] = in[11];               // g_\delta \ast u1
  out[10] = in[12];              // g_\delta \ast u2
}


// ========================================================================
// parameters:
// used for : GL00AuxProblem
// ========================================================================
void TimeNSParamsGL00AuxProblemPaper2(double *in, double *out)
{
  // \nabla u
  out[0] = in[4];                // D1u1
  out[1] = in[5];                // D1u2
  out[2] = in[6];                // D2u1
  out[3] = in[7];                // D2u2

  //  solution of auxiliary problem
  out[4] = in[ 8];
  out[5] = in[ 9];
  out[6] = in[10];

  // les term with exact solution
  out[7] = in[11];               // (t11)_x
  out[8] = in[12];               // (t12)_x
  out[9] = in[13];               // (t12)_y
  out[10] = in[14];              // (t22)_y

}


// ========================================================================
// parameters: gradient(u1), gradient(u2)
// ========================================================================
void TimeNSParamsGrad(double *in, double *out)
{
  // cout << "GRAD" << endl;
  out[0] = in[4];                // D10(u1old)
  out[1] = in[5];                // D10(u2old)
  out[2] = in[6];                // D01(u1old)
  out[3] = in[7];                // D01(u2old)
  // cout << in[4] << " E " << in[5] << " " << in[6] << " " << in[7] << endl;
}


// ========================================================================
// used for VMS, assembling of rhs for small scale equation
// ========================================================================
void TimeNSParams_VMS_SmallRhs2D(double *in, double *out)
{
  // large scales
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  // small scales
  out[6] = in[8];                // u1old
  out[7] = in[9];                // u2old
  out[8] = in[10];               // D1u1
  out[9] = in[11];               // D1u2
  out[10] = in[12];              // D2u1
  out[11] = in[13];              // D2u2

  // large pressure
  out[12] = in[14];              // p
}


// ========================================================================
// used for VMS, assembling of rhs for large scale equation
// ========================================================================
void TimeNSParams_VMS_LargeRhs2D(double *in, double *out)
{
  // large scales
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  // small scales
  out[6] = in[8];                // u1old
  out[7] = in[9];                // u2old
  out[8] = in[10];               // D1u1
  out[9] = in[11];               // D1u2
  out[10] = in[12];              // D2u1
  out[11] = in[13];              // D2u2

  // small pressure
  out[12] = in[14];              // p
}


// ========================================================================
// used for VMS, assembling of matrix, variant 1
// ========================================================================
void TimeNSParams_NLGalerkin_VMS_1_2D(double *in, double *out)
{
  // large scales
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  // small scales
  out[2] = in[4];                // u1current
  out[3] = in[5];                // u2current
}
// ========================================================================
// parameters: u1old, u2old, gradients; values of previous time steps
// defect correction type 2
// ========================================================================
void TimeNSParamsVelo_GradVeloOld2(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  out[6] = in[8];     //    u_1m11 
  out[7] = in[9];  // u_1m12
  out[8] = in[10]; // u_1m21
  out[9] = in[11]; // u_1m22
}

// ========================================================================
// parameters: x, y, u1old, u2old
// used for : SSMUM
// ========================================================================
void TimeNSParamsVeloPos(double *in, double *out)
{
  out[0] = in[0];                // x
  out[1] = in[1];                // y
  out[2] = in[2];                // u1old
  out[3] = in[3];                // u2old
}

// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void MovingTNSParams(double *in, double *out)
{
  out[0] = in[2] - in[4];
  out[1] = in[3] - in[5];

// cout << " out[0] " << in[4] << " out[1] " << in[4] <<endl;
}
// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void MovingTNSParams_Axial3D(double *in, double *out)
{
  out[0] = in[2] - in[4];
  out[1] = in[3] - in[5];
  out[2] = in[0];     // x value for axial symmetric

// cout<< "out[2]  " << out[2]<<endl;
}

// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void MovingTNSParams_Axial3D_HeatLine(double *in, double *out)
{
  out[0] = in[2]; 
  out[1] = in[3];
  out[2] = in[4]; 
  out[3] = in[5] - in[7];
  out[4] = in[6] - in[8]; 
  out[5] = in[9] - in[11]; 
  out[6] = in[10] - in[12]; 
  out[7] = in[0];  // x value for axial symmetric

// cout<< "out[0]  " << out[0]<<endl;
}



void TimeNSType2SUPG(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatA, **MatM, **MatK;
  double **MatB1, **MatB2, **MatB1T, **MatB2T;
  double *Rhs1, *Rhs2;
  double *MatRowA, *MatRowM, *MatRowK, *MatRow1, *MatRow2;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6, *Orig7;
  double test00, test10, test01;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double c0, c1, c2, u1, u2, val;
  int N_U, N_P;

  MatA=LocMatrices[0];
  MatM=LocMatrices[1];
  MatK=LocMatrices[2];
  
  MatB1=LocMatrices[3];
  MatB2=LocMatrices[4];
  MatB1T=LocMatrices[5];
  MatB2T=LocMatrices[6];

  Rhs1=LocRhs[0];
  Rhs2=LocRhs[1];

  N_U=N_BaseFuncts[0];
  N_P=N_BaseFuncts[1];
  

  Orig0=OrigValues[0]; // u_x
  Orig1=OrigValues[1]; // u_y 
  Orig2=OrigValues[2]; // u 
  Orig3=OrigValues[3]; // u_xx
  Orig4=OrigValues[4]; // u_yy
  Orig5=OrigValues[5]; // p_x 
  Orig6=OrigValues[6]; // p_y
  Orig7=OrigValues[7]; // p

  c0=coeff[0];
  c1=coeff[1];
  c2=coeff[2];

  u1=param[0];
  u2=param[1];

  // stabilization parameters have to defined
  // for initial test its setted to delta0*hK*hK
  double delta=TDatabase::ParamDB->DELTA0*hK*hK;
  double ugrad;

  for(int i=0;i<N_U;i++)
  {
    MatRowA=MatA[i];
    MatRowM=MatM[i];
    MatRowK=MatK[i];
    
    test10=Orig0[i];
    test01=Orig1[i];
    test00=Orig2[i];

    ugrad= delta*(u1*test10 + u2*test01);

    Rhs1[i] += Mult*(test00 + ugrad)*c1;
    Rhs2[i] += Mult*(test00 + ugrad)*c2;

    for(int j=0;j<N_U;j++)
    {
      ansatz10=Orig0[j];
      ansatz01=Orig1[j];
      ansatz00=Orig2[j];
      ansatz20=Orig3[j];
      ansatz02=Orig4[j];

      // Galerkin terms
      val =c0*(test10*ansatz10+test01*ansatz01); // diffusion term
      val +=(u1*ansatz10 + u2*ansatz01)*test00; // nonlinear (convective term)
      val +=(-c0*(ansatz20+ansatz02) + (u1*ansatz10 + u2*ansatz01))*ugrad; // SUPG terms
      
      MatRowA[j] += Mult*val; // A block
      MatRowM[j] += Mult*ansatz00*test00;

      MatRowK[j] += Mult*ansatz00*ugrad;
    }

    MatRow1=MatB1T[i];
    MatRow2=MatB2T[i];
    
    for(int j=0;j<N_P;j++)
    {
      ansatz10=Orig5[j];
      ansatz01=Orig6[j];
      ansatz00=Orig7[j];

      val = -ansatz00*test10; // Galerkin term
      val += ansatz10*ugrad;
      MatRow1[j] += Mult*val;
      
      val = -ansatz00*test01; // term due to stabilization
      val += ansatz01*ugrad;
      MatRow2[j] += Mult*val;
    }
  }

  for(int i=0;i<N_P; i++)
  {
    MatRow1=MatB1[i];
    MatRow2=MatB2[i];

    test00=Orig7[i];
    for(int j=0;j<N_U;j++)
    {
      ansatz10=Orig0[j];
      ansatz01=Orig1[j];

      val = -Mult*test00*ansatz10;
      MatRow1[j] += val;
      
      val = -Mult*test00*ansatz01;
      MatRow2[j] += val;
    }
  }
}


void TimeNSType2NLSUPG(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatA, **MatK, **MatB1T, **MatB2T;
  double *MatRowA, *MatRowK, *MatRow1, *MatRow2;
  double *Orig0, *Orig1, *Orig2, *Orig5, *Orig6, *Orig7;
  double test00, test10, test01, ansatz00, ansatz10, ansatz01; // ansatz20, ansatz02;
  double *Rhs1, *Rhs2;
  double c0, c1, c2, u1, u2, val, delta, ugrad;
  int N_U, N_P;

  MatA=LocMatrices[0];
  MatK=LocMatrices[1];
  MatB1T=LocMatrices[2];
  MatB2T=LocMatrices[3];


  Rhs1=LocRhs[0];
  Rhs2=LocRhs[1];

  N_U=N_BaseFuncts[0];
  N_P=N_BaseFuncts[1];
  

  Orig0=OrigValues[0]; // u_x
  Orig1=OrigValues[1]; // u_y 
  Orig2=OrigValues[2]; // u 
//  Orig3=OrigValues[3]; // u_xx
//  Orig4=OrigValues[4]; // u_yy
  Orig5=OrigValues[5]; // p_x 
  Orig6=OrigValues[6]; // p_y
  Orig7=OrigValues[7]; // p

  c0=coeff[0];
  c1=coeff[1];
  c2=coeff[2];

  u1=param[0];
  u2=param[1];
  
  if(c0<hK)
    delta=TDatabase::ParamDB->DELTA0*hK*hK;
  else 
    delta=TDatabase::ParamDB->DELTA1*hK*hK;

  for(int i=0;i<N_U;i++)
  {
    MatRowA=MatA[i];
    MatRowK=MatK[i];
    
    test10=Orig0[i]; 
    test01=Orig1[i];
    test00=Orig2[i];
    
    ugrad=delta*(u1*test10 + u2*test01);
    
    Rhs1[i] += Mult*ugrad*c1;
    Rhs2[i] += Mult*ugrad*c2;
    
    for(int j=0;j<N_U;j++)
    {
      ansatz10=Orig0[j];
      ansatz01=Orig1[j];
      ansatz00=Orig2[j];
//      ansatz20=Orig3[j];
//      ansatz02=Orig4[j];
      
      val =c0*(test10*ansatz10 + test01*ansatz01); // diffusion term
      val +=(u1*ansatz10+u2*ansatz01)*test00; // nonlinear term
      val +=(/*-c0*(ansatz20 + ansatz02) +*/ (u1*ansatz10 + u2*ansatz01))*ugrad; // SUPG terms
      
      MatRowA[j] += Mult*val;
      
      MatRowK[j] += Mult*ansatz00*ugrad;
    }
    
    MatRow1=MatB1T[i];
    MatRow2=MatB2T[i];
    for(int j=0;j<N_P;j++)
    {
      ansatz10=Orig5[j];
      ansatz01=Orig6[j];
      ansatz00=Orig7[j];
      
      val = -ansatz00*test10;
      val += ansatz10*ugrad;
      MatRow1[j] += Mult*val;
      
      val = -ansatz00*test01;
      val += ansatz01*ugrad;
      MatRow2[j] += Mult*val;
    }
  }
}

void TimeNSRHSSUPG(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double *Rhs3, *Rhs4;
  double *Orig0, *Orig1;  // *Orig2;
  double test10, test01;
  double u1, u2, c0, c1, c2, ugrad, delta;
  int N_U;
  
//   Rhs1=LocRhs[0];
//   Rhs2=LocRhs[1];
  Rhs3=LocRhs[0];
  Rhs4=LocRhs[1];

  N_U = N_BaseFuncts[0];
  
  Orig0=OrigValues[0]; // u_x
  Orig1=OrigValues[1]; // u_y 
//  Orig2=OrigValues[2]; // u

  c0=coeff[0];
  c1=coeff[1];
  c2=coeff[2];

  u1=param[0];
  u2=param[1];
  
  if(c0<hK)
    delta=TDatabase::ParamDB->DELTA0*hK*hK;
  else 
    delta=TDatabase::ParamDB->DELTA1*hK*hK;

  for(int i=0;i<N_U;i++)
  {
    test10=Orig0[i]; 
    test01=Orig1[i];
//    test00=Orig2[i];
    
    ugrad=delta*(u1*test10 + u2*test01);
    /*
    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;*/
    
    Rhs3[i] += Mult*ugrad*c1;
    Rhs4[i] += Mult*ugrad*c2;
  }   
}

void TimeNSParams4(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // u1, previous time
  out[3] = in[5];                // u2, previous time
}







/* **** BELOW THIS LINE USER SPECIFIC CODE*************** */
void TimeNSType1Galerkin_dimensional(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c1, c2;
  double u1, u2, u3, u4;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

//  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
//  c3 = coeff[3];                 // rho taken as a coefficient from examples
//  c4 = coeff[4];                 // mu  taken as a coefficient from examples

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // rho_field taken as a param from fe_function in local_assembling
  u4 = param[3];                 // mu_field taken as a param from fe_function in local_assembling

//        cout << " mu = " << u4 << " ";
//        cout << " rho= " << u3 << " ";
//
//      cout << "valeur de c4 = mu = " << c4 << " ";
//      cout << "valeur de c3 = rho= " << c3 << endl;

  /** NOTES: there are 2 ways to consider the property fields in the equations : take it
   * as input from example objects (as a coefficient, written in the example methods or given
   * as user input), use c3 and c4 to use this case and replace them in val, below. The other way
   * is to read them as values from fe_functions taken as a Param in a local Assembling. This is
   * with u3 and u4.
   */

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += u3*Mult*test00*c1;  // rhs should be multiplied by rho in the EXAMPLE CLASS
    Rhs2[i] += u3*Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      // hand test, reynolds number =2 ,rho = 4, mu=2
      val  = u4*(test10*ansatz10+test01*ansatz01);
      val += u3*(u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = u3*ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -u3*Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -u3*Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
}

void TimeNSType1_2NLGalerkin_dimensional(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double u1, u2, u3, u4;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

//  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u3 = param[2];                 // rho_field taken as a param from fe_function in local_assembling
  u4 = param[3];                 // mu_field taken as a param from fe_function in local_assembling

  /** NOTES: there are 2 ways to consider the property fields in the equations : take it
   * as input from example objects (as a coefficient, written in the example methods or given
   * as user input), use c3 and c4 to use this case and replace them in val, below. The other way
   * is to read them as values from fe_functions taken as a Param in a local Assembling. This is
   * with u3 and u4.
   */

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = u4*(test10*ansatz10+test01*ansatz01);
      val += u3*(u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


void TimeNSParamsVelo_dimensional(double *in, double *out)
{
  out[0] = in[2];   // u1old
  out[1] = in[3];   // u2old
  out[2] = in[4];   // rho_field
  out[3] = in[5];   // mu_field

  // below are derivatives from phase field Phi
  out[4] = in[6];   // derivative w.r.t. x Phix
  out[5] = in[7];   // Phi_y
  out[6] = in[8];   // Phi_xy = Phi_yx
  out[7] = in[9];   // Phi_xx
  out[8] = in[10];   // Phi_yy

//  for (int i=0; i<7;i++)
//  {
//    cout << "in i = " << i << " " << in[i] << endl;
//  }
}

void TimeNSParamsVelo_GradVelo_dimensional(double *in, double *out)
{
  out[0] = in[2];   // u1old
  out[1] = in[3];   // u2old

  out[2] = in[4];   // D1u1
  out[3] = in[5];   // D1u2
  out[4] = in[6];   // D2u1
  out[5] = in[7];   // D2u2

  out[6] = in[8];   // rho_field
  out[7] = in[9];   // mu_field

}

// the following function is the ParamFct for RHS_dimensional
// taking into account the gradient and 2nd derivatives of
// phase fraction (or density) to be used in the RHS routine to
// calculate the CSF.
void TimeNSParamsRhs_dimensional(double *in, double *out)
{
  out[0] = in[2];   // rho_denoted as R

  // below are derivatives from phase field Phi
  out[1] = in[3];   // derivative w.r.t. x Phix
  out[2] = in[4];   // Phi_y
  out[3] = in[5];   // Phi_xy = Phi_yx
  out[4] = in[6];   // Phi_xx
  out[5] = in[7];   // Phi_yy

//  out[0] = in[2];   // u1old
//  out[1] = in[3];   // u2old
//  out[2] = in[4];   // rho_field
//  out[3] = in[5];   // mu_field
//  out[4] = in[6];   // Rx
//  out[5] = in[7];   // Ry
}


void TimeNSRHS_dimensional(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i, N_U;
  double c1, c2, R;
  double Phi_x, Phi_y, Phi_xy, Phi_xx, Phi_yy;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  R = param[0];     // rho_field taken as a param from fe_function in local_assembling
  Phi_x = param[1];
  Phi_y = param[2];
  Phi_xy = param[3];
  Phi_xx = param[4];
  Phi_yy = param[5];

  //Debug code
//  Output::print("p0=R=", R, " p1=Phi_x=", Phi_x, " p2=Phi_y=", Phi_y," p3=Phi_xy=", Phi_xy,
//                " p4=Phi_xx=", Phi_xx, " p5=Phi_yy=",Phi_yy, " p6=", param[6],
//                " p7=",param[7]," p8=" ,param[8]," p9=",param[9]," p10=",
//                param[10]);
//  Output::print("p0=u1=", param[0], " p1=u2=", param[1], " p2=R=", R,
//                " p3=mu=", param[3]," p4=Phi_x=", Phi_x, " p5=Phi_y=",Phi_y,
//                " p6=", param[6]," p7=",param[7]," p8=" ,param[8],
//                " p9=",param[9]," p10=",param[10]);

  /* Curvature kappa of the CSF force */
  double kappa;
  if ((Phi_x*Phi_x+Phi_y*Phi_y) <= 1e-14) //if gradient zero, no interface, try also <= 1e-14
    kappa = 0;
  else
  {
    kappa = -(Phi_xx*Phi_y*Phi_y+Phi_yy*Phi_x*Phi_x-2*Phi_x*Phi_y*Phi_xy)
            /pow(Phi_x*Phi_x+Phi_y*Phi_y,1.5);
  }

  // surface tension coefficients
  double tau = TDatabase::ParamDB->P9;
  double constant = TDatabase::ParamDB->P10;
  if ( constant == 0)
    ErrThrow("ERROR: Parameter P10 is used for surface tension and should not be 0!");
  double surfacetension = tau*kappa/constant;

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*(R*test00*c1 + surfacetension*Phi_x*test00);
    Rhs2[i] += Mult*(R*test00*c2 + surfacetension*Phi_y*test00);
//    Output::print<1>("Kappa ", kappa);// " Tau ", tau, " S ", surfacetension
//                     , " Phi_x ", Phi_x, " Phi_y ", Phi_y,
//                     " Rhs1[i] ", Rhs1[i], " Rhs2[i] ", Rhs2[i]);

  } // endfor i
}



// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// Type 3, Coletti, D(u):D(v)
// Type 3, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType3GalerkinDD_dimensional(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA21 = LocMatrices[2];
  double **MatrixA22 = LocMatrices[3];
  double **MatrixM11 = LocMatrices[4];
  double **MatrixB1  = LocMatrices[5];
  double **MatrixB2  = LocMatrices[6];

  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0];         // u_x
  double *Orig1 = OrigValues[1];         // u_y
  double *Orig2 = OrigValues[2];         // u
  double *Orig3 = OrigValues[3];         // p

  // double c0 = coeff[0];              // nu
  double c1 = coeff[1];                 // f1
  double c2 = coeff[2];                 // f2
  // double c3 = coeff[3];              // rho taken as a coefficient from examples
  // double c4 = coeff[4];              // mu  taken as a coefficient from examples

  double u1 = param[0];                 // u1old
  double u2 = param[1];                 // u2old
  double u3 = param[2];                 // rho_field taken as a param from fe_function in local_assembling
  double u4 = param[3];                 // mu_field taken as a param from fe_function in local_assembling

  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;

  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row;
  double *MatrixRow1, *MatrixRow2;
  double val, val1;

  for(int i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += u3*Mult*test00*c1;
    Rhs2[i] += u3*Mult*test00*c2;

    for(int j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = u3*(u1*ansatz10+u2*ansatz01)*test00;
      val  = u4*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = u4*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = u4*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = u4*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = u3*Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(int i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(int j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -u3*Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -u3*Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// Type 3, Coletti, D(u):D(v)
// Type 3, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD_dimensional(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA21 = LocMatrices[2];
  double **MatrixA22 = LocMatrices[3];
  double **MatrixM11 = LocMatrices[4];
  double **MatrixB1  = LocMatrices[5];
  double **MatrixB2  = LocMatrices[6];
  double **MatrixB1T = LocMatrices[7];
  double **MatrixB2T = LocMatrices[8];

  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0];         // u_x
  double *Orig1 = OrigValues[1];         // u_y
  double *Orig2 = OrigValues[2];         // u
  double *Orig3 = OrigValues[3];         // p

  // double c0 = coeff[0];              // nu
  double c1 = coeff[1];                 // f1
  double c2 = coeff[2];                 // f2
  // double c3 = coeff[3];              // rho taken as a coefficient from examples
  // double c4 = coeff[4];              // mu  taken as a coefficient from examples

  double u1 = param[0];                 // u1old
  double u2 = param[1];                 // u2old
  double u3 = param[2];                 // rho_field taken as a param from fe_function in local_assembling
  double u4 = param[3];                 // mu_field taken as a param from fe_function in local_assembling
  double Phi_x = param[4];              // grad components of phase field phi
  double Phi_y = param[5];              // grad components of phase field phi
  double Phi_xy = param[6];             // grad components of phase field phi
  double Phi_xx = param[7];             // grad components of phase field phi
  double Phi_yy = param[8];             // grad components of phase field phi

  /* Curvature kappa of the CSF force */
  double kappa;
  if ((Phi_x*Phi_x+Phi_y*Phi_y) <= 1e-14) //if gradient zero, no interface, try also <= 1e-14
    kappa = 0;
  else
  {
    kappa = -(Phi_xx*Phi_y*Phi_y+Phi_yy*Phi_x*Phi_x-2*Phi_x*Phi_y*Phi_xy)
            /pow(Phi_x*Phi_x+Phi_y*Phi_y,1.5);
  }
  // surface tension coefficients
  double tau = TDatabase::ParamDB->P9;
  double constant = TDatabase::ParamDB->P10;
  if ( constant == 0)
    ErrThrow("ERROR: Parameter P10 is used for surface tension and should not be 0!");
  double surfacetension = tau*kappa/constant;


  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;

  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row;
  double *MatrixRow1, *MatrixRow2;
  double val, val1;

  for(int i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += u3*Mult*test00*c1+Mult*surfacetension*Phi_x*test00;
    Rhs2[i] += u3*Mult*test00*c2+Mult*surfacetension*Phi_y*test00;

    for(int j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = u3*(u1*ansatz10+u2*ansatz01)*test00;
      val  = u4*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = u4*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = u4*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = u4*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = u3*Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(int j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(int i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(int j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}



// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD_dimensional(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  // double c0 = coeff[0];              // nu
  // double c1 = coeff[1];                 // f1
  // double c2 = coeff[2];                 // f2
  // double c3 = coeff[3];              // rho taken as a coefficient from examples
  // double c4 = coeff[4];              // mu  taken as a coefficient from examples

  double u1 = param[0];                 // u1old
  double u2 = param[1];                 // u2old
  double u3 = param[2];                 // rho_field taken as a param from fe_function in local_assembling
  double u4 = param[3];                 // mu_field taken as a param from fe_function in local_assembling

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = u3*(u1*ansatz10+u2*ansatz01)*test00;
      val  = u4*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = u4*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}

void TimeNSType1_3_4GalerkinDDMass_dimensional(double Mult, double *coeff,
                                         double *param, double hK,
                                         double **OrigValues, int *N_BaseFuncts,
                                         double ***LocMatrices, double **LocRhs)
{
  double **MatrixM;
  double val;
  double *MatrixMRow;
  double ansatz00;
  double test00;
  double *Orig0;
  int i,j,N_U;
  double u3;

  MatrixM = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

//  c0 = coeff[0];                 // nu

//  u1 = param[0];                 // u1old
//  u2 = param[1];                 // u2old
  u3 = param[2];                 // rho_field taken as a param from fe_function in local_assembling
//  u4 = param[3];                 // mu_field taken as a param from fe_function in local_assembling

  /** NOTES: there are 2 ways to consider the property fields in the equations : take it
   * as input from example objects (as a coefficient, written in the example methods or given
   * as user input), use c3 and c4 to use this case and replace them in val, below. The other way
   * is to read them as values from fe_functions taken as a Param in a local Assembling. This is
   * with u3 and u4.
   */

  for(i=0;i<N_U;i++)
  {
    MatrixMRow = MatrixM[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz00 = Orig0[j];

      val  = u3*ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}

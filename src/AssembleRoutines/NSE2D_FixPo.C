// ======================================================================
// @(#)NSE2D_FixPo.C        1.3 06/27/00
//
// common declaration for all Navier-Stokes problems
// fix point iteration
//
// ======================================================================
#include <Convolution.h>
#include <Database.h>
#include <TNSE2D_Routines.h>
// #include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!

// Hotfixglobal_AssembleNSE assemble_nse(Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION);

// ======================================================================
// compute parameter for RFB stabilization
// Brezzi, Marini, Russo, CMAME 194 (2005) 127 - 148
// ======================================================================

double RFB_Parameter(double hK, double eps, double* b)
{
    double peclet, normb,tau;

    normb = sqrt(b[0]*b[0]+b[1]*b[1]);
    if (normb==0)
	return(0.0);
    peclet = normb*hK/(2*eps);
    if (peclet > 1)
        tau = hK/(3*normb);
    else
	tau = hK*hK/eps;

    tau = hK*hK;
    return(tau);
	    
}

// ======================================================================
// compute parameter for SUPG stabilization
// ======================================================================

double SUPG_Parameter(double hK, double eps, double b1, double b2, double c)
{
  double delta, val, norm_b, M, r=2.0, c_inv_2 = 24.0, beta_0 = 1.0, c_pf = 1.0;
  double beta_1 = 1.0, C0, C_h = 1.0, delta_0 = 0.1;
  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  
  switch (TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case 2:
    case 12:
        if (TDatabase::ParamDB->INTERNAL_MESH_CELL_TYPE==3)
        c_inv_2 = 48.0;
  else
        c_inv_2 = 24.0;
      break;
      case 3:
      case 13:
        if (TDatabase::ParamDB->INTERNAL_MESH_CELL_TYPE==3)
        c_inv_2 = (435+sqrt(26025.0))/4.0;
  else
        c_inv_2 = (244+sqrt(9136.0))/3.0;
         break; 
    default: 
         c_inv_2 = (435+sqrt(26025.0))/4.0;
         break;     
  }

  switch (TDatabase::ParamDB->SDFEM_TYPE)
  {
    case 1:
      // standard parameter 
      delta =  hK*hK;
      delta =  delta0*delta;
       //OutPut(delta << endl);
       break;
    case 2:                                      
      // parameter from [BBJL07]
      delta =  hK*hK/(r*r*(eps+c));
      delta =  delta0*delta;
  //OutPut(delta << endl);
       break;
      // paper with Julia Novo
    case 3:
         delta = hK*hK/(3*eps*c_inv_2);
  //OutPut(delta << " ");
  if (c>0)
  {
    val = 1.0/(3.0*c);
    if (val<delta)
      delta = val;
  }
  //OutPut(delta << " ");
  norm_b= fabs(b1);
  if (fabs(b2)>norm_b)
    norm_b = fabs(b2);
  if (norm_b > 0)
  {
  val = beta_0 *hK/(4*norm_b*sqrt(c_inv_2));
  if (val<delta)
          delta = val;
  }
  //OutPut(delta << " ");
  M = 14*eps;
  val = 21 * c * c_pf *  c_pf ;
  if (val > M)
    M = val;
  val = 14 * delta1;
  if (val > M)
    M = val;
  if (c>0)
  {
    val = 21 * norm_b * norm_b/c;
    if (val > M)
      M = val;
  }
  M /= (beta_0*beta_0);
  val = hK * hK /(24 * M * c_inv_2);
  if (val<delta)
          delta = val;

  //OutPut(delta << endl);
     break;
    case 4:
         delta = hK*hK/(2*eps*c_inv_2);
  // OutPut(delta << " ");
  norm_b= fabs(b1);
  if (fabs(b2)>norm_b)
    norm_b = fabs(b2);
  if (norm_b > 0)
  {
          val = beta_1 * hK/(4*norm_b*sqrt(c_inv_2));
    if (val<delta)
            delta = val;
  }
  //OutPut(delta << " ");
  C0 = delta;
  M = 10 * eps * c_inv_2/(C_h*C_h);
  val = 10.0/(C_h*C_h*delta_0);
  if (val > M)
    M = val;
  val = 10 * delta1 * c_inv_2/(C_h*C_h);
  if (val > M)
    M = val;
  val = 10 * C0 * norm_b * norm_b * c_inv_2/(C_h*C_h);
  if (val > M)
    M = val;  
  M /= (beta_1*beta_1);
  val = beta_1 * hK * hK /(16 * M);
  if (val<delta)
          delta = val;

  //OutPut(delta << " " << TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE << endl);
     break;
   case 11:
      // standard parameter 
      if (eps < hK)
        delta = hK;
      else
        delta = hK * hK;
      delta =  delta0*delta;
        //OutPut(delta << endl);
       break;
   case 12:
       delta =  hK;
       delta =  delta0*delta;
         //OutPut(delta << endl);
        break;

    default:
      OutPut("SDFEM_TYPE " << TDatabase::ParamDB->SDFEM_TYPE << " not defined !!!" << endl);
      exit(4711);
  }
  return(delta);
}


// ======================================================================
// Type 1, Standard Galerkin


// ======================================================================
// Type 1, (reduced) SDFEM or (simplified RFB)
// ======================================================================
void NSType1SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double tau, ugrad;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  //delta = delta0*hK;
  tau = RFB_Parameter(hK,c0,&param[0]); // stabilization parameter

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = test00 + tau * (u1*test10+u2*test01);
    Rhs1[i] += Mult*ugrad*c1;
    Rhs2[i] += Mult*ugrad*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*ugrad;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void NSType1Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  //  cout << "c3";

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      //    val += c3 * Orig2[j] *test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 1, Smagorinsky
// ======================================================================
void NSType1Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

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
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 1, VMSProjection, nabla form
// ======================================================================
void NSType1VMSProjection(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11;
  double **MatrixB1, **MatrixB2;
  double **MatrixL, **Matrix_tilde_G11 ;
  double **Matrix_tilde_G22, **Matrix_G11, **Matrix_G22;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P, N_L;
  double c0, c1, c2;
  double u1, u2, mu, delta;
  //cout << "VMS" << endl;
  MatrixA11 = LocMatrices[0];
  MatrixL   = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  Matrix_tilde_G11  = LocMatrices[4];
  Matrix_tilde_G22  = LocMatrices[5];
  Matrix_G11  = LocMatrices[6];
  Matrix_G22  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[2];

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

  // for computational comparisons of Oseen problems
  if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == OSEEN_PROBLEM 
     || TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    u1 = coeff[3];
    u2 = coeff[4];
    param[0] = u1;
    param[1] = u2;
  }

  // compute the size of the mesh lenght for the turbulent viscosity
  delta =  CharacteristicFilterWidth(hK);
  // compute the turbulent viscosity
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    // assemble right hand side
    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      // assemble matrix A (viscous + convective + Smagorinsky term)
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      // assemble matrices B1 and B2, divergence constraint	
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i

  // matrices for VMS
  for(i=0;i<N_U;i++)
  {
    // projections in the momentum balance
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row  = Matrix_tilde_G22[i];
    test10 = Orig0[i];  // velo_x
    test01 = Orig1[i];  // velo_y

     for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];  // large scales 
      Matrix11Row[j] -= Mult * mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * mu * ansatz00 * test01;
    }
  }

  for(i=0;i<N_L;i++)
  {
    // elliptic projection
    Matrix11Row = Matrix_G11[i];
    Matrix22Row = Matrix_G22[i];
    test00 = Orig4[i]; // large scales

    for(j=0;j<N_U;j++)
    {
	ansatz10 = Orig0[j];  // velo_x
	ansatz01 = Orig1[j];  // velo_y

      Matrix11Row[j] -= Mult * ansatz10 * test00;
      Matrix22Row[j] -= Mult * ansatz01 * test00;
    }
  }

  for(i=0;i<N_L;i++)
  {
    // mass matrix of large scale space
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
// Type 1, Standard Galerkin + divergence term 
// ======================================================================
void NSType1GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, du1x, du2y, divu;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u
  Orig1 = OrigValues[1];                          // p
  Orig2 = OrigValues[2];                          // u_x
  Orig3 = OrigValues[3];                          // u_y

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;

      MatrixRow[j] += Mult * val;
    }                                             // endfor j
  }                                               // endfor i
 for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                                             // endfor j

  }                                               // endfor i
  //cout << " end. " << endl;
}


// ======================================================================
// Type 2, Standard Galerkin



// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2, c;
  double u1, u2;
  double delta, ugrad;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  
  // for computational comparisons of Oseen problems
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    u1 = coeff[3];
    u2 = coeff[4];
    c  = coeff[5];
    param[0] = u1;
    param[1] = u2;
  }

  if(c0 < hK)
    delta = delta0*hK*hK;
  else
    delta = delta1*hK*hK;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);
    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];
      // standard terms
      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      // SD term
      val += (-c0*(ansatz20+ansatz02)
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      // for computational comparisons of Oseen problems
      if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
      {
        ansatz00 = Orig0[j];
        val += c * ansatz00 * test00;
        val += c * ansatz00 * ugrad;
      }

      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;

      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void NSType2Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);

      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

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

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

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
void NSType2Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

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
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

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

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 2, Standard Galerkin + divergence term
// ======================================================================
void NSType2GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, du1x, du2y, divu;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u
  Orig1 = OrigValues[1];                          // p
  Orig2 = OrigValues[2];                          // u_x
  Orig3 = OrigValues[3];                          // u_y

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;

      MatrixRow[j] += Mult * val;
    }                                             // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;
      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                                               // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;
      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                                             // endfor j

  }                                               // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================


// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row; // *Matrix21Row, *Matrix12Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

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
void NSType3UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

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
void NSType3Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row; // *Matrix21Row, *Matrix12Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double mu, delta;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

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
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

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
void NSType3SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

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
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = 2*(c0+mu)*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0+mu)*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin + div term, (grad u, grad v)
// ======================================================================
void NSType3GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row; // *Matrix21Row, *Matrix12Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2, du1x, du2y, divu;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1  = LocMatrices[4];
  MatrixB2  = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];                          // u
  Orig1 = OrigValues[1];                          // p
  Orig2 = OrigValues[2];                          // u_x
  Orig3 = OrigValues[3];                          // u_y

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // f1
  c2 = coeff[2];                                  // f2

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
//    Matrix12Row = MatrixA12[i];
//    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;
// OutPut("divu =" << divu << endl);
      Matrix11Row[j] += Mult * val;


      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;
      Matrix22Row[j] += Mult * val;

    }                                             // endfor j
  }                                               // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                                             // endfor j

  }                                               // endfor i
}


// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// ======================================================================



// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// ======================================================================



// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad, r=2.0, maxu, tau;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // for computational comparisons of Oseen problems
   if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == OSEEN_PROBLEM
      || TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    u1 = coeff[3];
    u2 = coeff[4];
  }

  maxu = fabs(u1);
  if (fabs(u2)>maxu)
      maxu = fabs(u2);

  // parameter from [BBJL07]
  delta =  hK*hK/(r*r*(c0+1));
  delta =  delta0*delta;

  tau =  delta1 *(1+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      // viscous term
      val  = c0*(test10*ansatz10+test01*ansatz01);
      // convective term
      val += (u1*ansatz10+u2*ansatz01)*test00;
      // SUPG term
      val +=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * ugrad;
      // term in both diagonal blocks
      Matrix11Row[j] += Mult * (val+tau*test10*ansatz10);
      Matrix22Row[j] += Mult * (val+tau*test01*ansatz01);
      // div-div term
      val = tau * test10*ansatz01; 
      Matrix12Row[j] += Mult * val;

      val = tau * test01*ansatz10; 
      Matrix21Row[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4SDFEMDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  if(c0 < hK)
    delta = delta0*hK*hK;
  else
    delta = delta1*hK*hK;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void NSType4Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

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

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
  
  // including ROSSBY term, ONLY FOR STOKES SO FAR
  if (fabs(TDatabase::ParamDB->ROSSBY_NR)>1e-10)
  {
    for(i=0;i<N_U;i++)
    {
      Matrix12Row = MatrixA12[i];
      Matrix21Row = MatrixA21[i];
      test00 = Orig0[i];
      for(j=0;j<N_U;j++)
      {
        ansatz00 = Orig0[j];
        val  = 2*TDatabase::ParamDB->ROSSBY_NR * ansatz00 * test00;
        Matrix12Row[j] -= Mult * val;
        Matrix21Row[j] += Mult * val;
      }   // endfor j
    }
  }
}


// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void NSType4UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

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

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  
  // including ROSSBY term, ONLY FOR STOKES SO FAR
  if (fabs(TDatabase::ParamDB->ROSSBY_NR)>1e-10)
  {
    for(i=0;i<N_U;i++)
    {
      Matrix12Row = MatrixA12[i];
      Matrix21Row = MatrixA21[i];
      test00 = Orig0[i];
      for(j=0;j<N_U;j++)
      {
        ansatz00 = Orig0[j];
        val  = 2*TDatabase::ParamDB->ROSSBY_NR * ansatz00 * test00;
        Matrix12Row[j] -= Mult * val;
        Matrix21Row[j] += Mult * val;
      }   // endfor j
    }
  }
}


// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType4Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;  // *Matrix21Row, *Matrix12Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double mu, delta;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

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
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

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

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

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
void NSType4SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

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
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = 2*(c0+mu)*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0+mu)*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

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

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================


// ======================================================================
// Type 1,(reduced) SDFEM or (simplified) RFB
// ======================================================================
void NSType1NLSDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow, *Rhs1, *Rhs2;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j,N_U;
  double c0, c1, c2;
  double u1, u2;
  double tau, ugrad;

  MatrixA = LocMatrices[0];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  tau = RFB_Parameter(hK,c0,&param[0]); // stabilization parameter

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = test00+tau * (u1*test10+u2*test01);
  
    Rhs1[i] += Mult*ugrad*c1;
    Rhs2[i] += Mult*ugrad*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*ugrad;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void NSType1_2NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test10, test01;  // test00;
  double *Orig2, *Orig3;  // *Orig0;
  int i,j,N_U;
  double c0;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
//  Orig0 = OrigValues[0];         // u

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
//    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = c0*(test10*ansatz10+test01*ansatz01);

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// ======================================================================
void NSType1_2NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, VMSProjection, nabla form
// Type 2, VMSProjection, nabla form
// ======================================================================
void NSType1_2NLVMSProjection(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double **Matrix_tilde_G11,  **Matrix_tilde_G22;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_L;
  double c0, viscosity, delta;
  double u1, u2, mu;
  //cout << "VMS" << endl;
  MatrixA11 = LocMatrices[0];
  Matrix_tilde_G11  = LocMatrices[1];
  Matrix_tilde_G22  = LocMatrices[2];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // l

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // for computational comparisons of Oseen problems
  if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == OSEEN_PROBLEM
     || TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
      u1 = coeff[3];
      u2 = coeff[4];
      // for the vms coefficient, still the params computed from the code
      // are used
      param[0] = u1;
      param[1] = u2;
      /*param[2] = coeff[5];
      param[3] = coeff[6];
      param[4] = coeff[7];
      param[5] = coeff[8];*/
  }

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  viscosity = mu+c0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;
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
      Matrix11Row[j] -= Mult * mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * mu * ansatz00 * test01;
    }
  }
}

// ======================================================================
// Type 1, Standard Galerkin + div term, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, du1x, du2y, divu;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u
  Orig2 = OrigValues[2];                          // u_x
  Orig3 = OrigValues[3];                          // u_y

  c0 = coeff[0];                                  // nu

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;

      MatrixRow[j] += Mult * val;
    }                                             // endfor j
  }                                               // endfor i
}


// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2NLSDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA = LocMatrices[0];
  MatrixB1T = LocMatrices[1];
  MatrixB2T = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  if(c0 < hK)
    delta = delta0*hK*hK ;
  else
    delta = delta1*hK*hK ;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
	      + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// ======================================================================



// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// ======================================================================



// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void NSType3_4NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig2, *Orig3;
  int i,j,N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void NSType3_4NLUpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig2, *Orig3;
  int i,j,N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void NSType3_4NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double val;
  double *Matrix11Row, *Matrix22Row; // *Matrix21Row, *Matrix12Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j,N_U; // N_P;
  double c0;
  double u1, u2;
  double mu, delta;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];
//  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu

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
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;
      /* do not see why these values are added, 06/03/20
      val  = mu*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = mu*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;
      */
      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void NSType3_4NLSmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2;
  double mu, delta;
  //cout << "hier" << endl;
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu

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
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = 2*(c0+mu)*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*(c0+mu)*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin + div term, (grad u, grad v)
// ======================================================================
void NSType3_4NLGalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, du1x, du2y, divu;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];                          // u
  Orig2 = OrigValues[2];                          // u_x
  Orig3 = OrigValues[3];                          // u_y

  c0 = coeff[0];                                  // nu

  u1 = param[0];                                  // u1old
  u2 = param[1];                                  // u2old
  du1x = param[2];
  du2y = param[5];
  divu = (du1x + du2y)/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;

      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += divu*ansatz00*test00;
      
      Matrix22Row[j] += Mult * val;

    }                                             // endfor j
  }                                               // endfor i
}


// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4NLSDFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad, r=2.0, maxu, tau;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixB1T = LocMatrices[2];
  MatrixB2T = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // for computational comparisons of Oseen problems
   if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == OSEEN_PROBLEM
      || TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
      u1 = coeff[3];
      u2 = coeff[4];
  }

  maxu = fabs(u1);
  if (fabs(u2)>maxu)
      maxu = fabs(u2);

  // parameter from [BBJL07]
  delta =  hK*hK/(r*r*(c0+1));
  delta =  delta0*delta;

  tau =  delta1 *(1+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      // viscous term
      val  = c0*(test10*ansatz10+test01*ansatz01);
      // convective term
      val += (u1*ansatz10+u2*ansatz01)*test00;
      // SUPG term
      val +=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * ugrad;
      // term in both diagonal blocks
      Matrix11Row[j] += Mult * (val+tau*test10*ansatz10);
      Matrix22Row[j] += Mult * (val+tau*test01*ansatz01);
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i
}

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4NLSDFEMDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixB1T = LocMatrices[2];
  MatrixB2T = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  if(c0 < hK)
    delta = delta0*hK*hK;
  else
    delta = delta1*hK*hK;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i
}


// ======================================================================
// auxiliary problem
// ======================================================================
void NSAuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **AuxMatrix;
  double *Rhs1, *Rhs2, val;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j, N_U;
  double mu2, delta;
  double u1, u2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  AuxMatrix = LocMatrices[0];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  //cout << " " << delta;
  mu2 = 0.25*delta*delta/gamma;
  for(i=0;i<N_U;i++)
  {
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*u1;
    Rhs2[i] += Mult*test00*u2;
    //    Rhs1[i] += Mult*test00*coeff[3];
    //    Rhs2[i] += Mult*test00*coeff[3];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// auxiliary problem for differential filter
// ======================================================================
void Filter_Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **AuxMatrix;
  double *Rhs1, *Rhs2, val;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j, N_U;
  double delta;
  double u1, u2;

  AuxMatrix = LocMatrices[0];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);
  delta *= delta;

  for(i=0;i<N_U;i++)
  {
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*u1;
    Rhs2[i] += Mult*test00*u2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = delta*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}

// ======================================================================
// rhs for RFB stabilization
// ======================================================================
void NSRFBRhs(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i,N_U;
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
  }                              // endfor i
}

// ======================================================================
// pressure separation
// ======================================================================
void NSPressSep(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i,N_U;
  double c1, c2;
  double p_x, p_y;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  p_x = param[0];                // p_sep to x
  p_y = param[1];                // p_sep to y

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];
    Rhs1[i] += Mult*test00*(c1-p_x);
    Rhs2[i] += Mult*test00*(c2-p_y);
  }                              // endfor i
}


// ======================================================================
// pressure separation
// ======================================================================
void NSPressSepAuxProb(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, *MatrixRow;
  double *Rhs1;
  double test10,test01,ansatz10,ansatz01;
  double *Orig2, *Orig3;
  double u1, u2, d1u1, d2u1, d1u2, d2u2;
  int i,j, N_U;
  double c1, c2;

  MatrixA = LocMatrices[0];
  Rhs1 = LocRhs[0];

  N_U = N_BaseFuncts[0];

  Orig3 = OrigValues[2];         // u_x
  Orig2 = OrigValues[3];         // u_y

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  if (TDatabase::ParamDB->PRESSURE_SEPARATION==3)
  {
    for(i=0;i<N_U;i++)
    {
      MatrixRow = MatrixA[i];
      test10 = Orig2[i];
      test01 = Orig3[i];
      Rhs1[i] += Mult*(test10*c1+test01*c2);
      for(j=0;j<N_U;j++)
      {
        ansatz10 = Orig2[j];
        ansatz01 = Orig3[j];
        MatrixRow[j] += Mult * (ansatz10*test10+ansatz01*test01);
      }
    }
  }
  else
  {
    u1 = param[0];               // u1old
    u2 = param[1];               // u2old
    d1u1 = param[2];
    d1u2 = param[3];
    d2u1 = param[4];
    d2u2 = param[5];
    for(i=0;i<N_U;i++)
    {
      MatrixRow = MatrixA[i];
      test10 = Orig2[i];
      test01 = Orig3[i];
      Rhs1[i] += Mult*(test10*(c1-u1*d1u1-u2*d2u1)
        +test01*(c2-u1*d1u2-u2*d2u2));
      for(j=0;j<N_U;j++)
      {
        ansatz10 = Orig2[j];
        ansatz01 = Orig3[j];
        MatrixRow[j] += Mult * (ansatz10*test10+ansatz01*test01);
      }
    }
  }
}

// ========================================================================
// routines for FJMT07
// ========================================================================

// ======================================================================
// Type 1, Standard Galerkin
// ======================================================================
void NSType1GalerkinFJMT07(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, tau;
  double sigma = 1/TDatabase::TimeDB->TIMESTEPLENGTH;

  //OutPut("FJMT");

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = coeff[3];                 // u1old
  u2 = coeff[4];                 // u2old

   // SD parameter
  tau = TDatabase::ParamDB->DELTA0*hK*hK;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*(test00*c1 + tau*c1*(u1*test10+u2*test01));
    Rhs2[i] += Mult*(test00*c2 + tau*c2*(u1*test10+u2*test01));

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += tau*(u1*ansatz10+u2*ansatz01)*(u1*test10+u2*test01);
      val += sigma * ansatz00 * test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// ======================================================================
void NSType1_2NLGalerkinFJMT07(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2,tau;
  double sigma = 1/TDatabase::TimeDB->TIMESTEPLENGTH;
  //OutPut("07");
  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu

  u1 = coeff[3];                 // u1old
  u2 = coeff[4];                 // u2old

  // SD parameter
  tau = TDatabase::ParamDB->DELTA0*hK*hK;
  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += sigma*ansatz00*test00;
      val += tau*(u1*ansatz10+u2*ansatz01)*(u1*test10+u2*test01);

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void NSParamsVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
}

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2)
// ========================================================================
void NSParamsVelo_GradVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D10(u1old)
  out[3] = in[5];                // D10(u2old)
  out[4] = in[6];                // D01(u1old)
  out[5] = in[7];                // D01(u2old)
}


void NSParamsVeloExact(double *in, double *out)
{
  double x,y, u1, u2;

  x = in[0];
  y = in[1];
  /*  t1 = x*x;
    t2 = 1.0-x;
    t3 = t2*t2;
    t4 = t1*t3;
    t5 = 1.0-y;
    t6 = t5*t5;
    t7 = y*t6;
    t9 = y*y;
    t10 = t9*t5;
    t13 = x*t3;
    t15 = t1*t2;
    u1 = 2.0*t4*t7-2.0*t4*t10;

    t1 = 1.0-x;
    t2 = t1*t1;
    t3 = x*t2;
    t4 = y*y;
    t5 = 1.0-y;
    t6 = t5*t5;
    t7 = t4*t6;
    t9 = x*x;
    t10 = t9*t1;
    t20 = y*t6;
    t23 = t4*t5;
    u2 = -2.0*t3*t7+2.0*t10*t7;
  */
  u1  =2*Pi*sin(Pi*y)*cos(Pi*y)*sin(Pi*x)*sin(Pi*x)*sin(Pi*x);
  u2 = -3*Pi*sin(Pi*x)*sin(Pi*x)*cos(Pi*x)*sin(Pi*y)*sin(Pi*y);
  out[0] = u1-in[2];             // u1old
  out[1] = u2-in[3];             // u2old
  //out[0] = in[2]; // u1old
  // out[1] = in[3]; // u2old
}


// ========================================================================
// parameters: separated pressure
// ========================================================================
void NSParamsPressSep(double *in, double *out)
{
  out[0] = in[2];                // P_sep to x
  out[1] = in[3];                // P_sep to y
}

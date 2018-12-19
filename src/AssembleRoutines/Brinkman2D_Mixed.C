// ======================================================================
// Brinkman2D_Mixed.C
// Common declaration of Brinkman problems in 2D
// including seperate modules for:
// - equal-order stabilization for P1/P1 
// (non-symmetric GLS, see [Douglas_Wang_An absolutely stabilized finite element method for the Stokes problem_1989]);
// - equal-order stabilization for P2/P2 without P1/P1 terms 
// (non-symmetric GLS, see [Douglas_Wang_An absolutely stabilized finite element method for the Stokes problem_1989]);
// - Grad-Div stabilization;

// Brinkman Problem:
//    -mu_eff \Delta (u1,u2) + \nabla p + sigma (u1,u2) = (f1,f2)
//    \nabla \dcot (u1,u2) = g
//    sigma := mu / K
// ======================================================================

#include "../../include/AssembleRoutines/Brinkman2D_Mixed.h"
#include <Database.h>


// ======================================================================
// Type 4, Standard Galerkin for Brinkman in [p div u] formulation
// Matrix Type:
// A11 0   B1^T
// 0   A22 B2^T
// B1  B2  0
// Type 14, Pressure-Stabilized Galerkin for Brinkman in [p div u] formulation
// Matrix Type:
// A11 0   B1^T
// 0   A22 B2^T
// B1  B2  C
// ======================================================================

// ======================================================================
// Standard Galerkin for Brinkman in [p div v] formulation
// ======================================================================
void BrinkmanType1Galerkin(double Mult, double *coeff,
    double *, double,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs)
{
  double ansatz00, ansatz10, ansatz01;    // ansatz functions
  double test00, test10, test01;          // test functions

  double ** MatrixA11 = LocMatrices[0];
  //double ** MatrixA12 = LocMatrices[1];
  //double ** MatrixA21 = LocMatrices[2];
  double ** MatrixA22 = LocMatrices[3];
 // double ** MatrixC = LocMatrices[4];
   double ** MatrixB1 = LocMatrices[5];
  double ** MatrixB2 = LocMatrices[6];
  double ** MatrixB1T = LocMatrices[7];
  double ** MatrixB2T = LocMatrices[8];

  double * Rhs_u1 = LocRhs[0];            // f_v1
  double * Rhs_u2 = LocRhs[1];            // f_v2
  double * Rhs_p = LocRhs[2];             // f_q

  int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space
  int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space

  // values of FE functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
  // Mult is the quadrature weight (wk)

  double * Orig0 = OrigValues[0];          // u_x (derivative)
  double * Orig1 = OrigValues[1];          // u_y (derivative)
  double * Orig2 = OrigValues[2];          // u
  double * Orig3 = OrigValues[3];          // p

  // values defined in the example
  //double c0 = coeff[0];                 // t^2 = (mueff/mu )* K = mueff/sigma
  double c1 = coeff[1];                   // f1
  double c2 = coeff[2];                   // f2
  double c3 = coeff[3];                   // g (the rhs of incompressibility constraint)
  double mu = coeff[4];                   // fluid viscosity
  double mu_eff = coeff[5];               // effective viscosity
  double K = coeff[6];                    // permeability

// LB Debug New 25.04.18 start:
double approximate_delta_distribution_function = 0;

if (TDatabase::ParamDB->SOURCE_SINK_FUNCTION)
{
	approximate_delta_distribution_function = coeff[9];
}
//cout<< "K::: "<< K << endl;
//cout<< "mueff::: "<< mu_eff << endl;
// LB Debug 25.04.18 end
double val;
  for(int i = 0; i < N_U; i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs_u1[i] += Mult * c1 * test00;      // (f,v) --> (f1, v1)
    Rhs_u2[i] += Mult * c2 * test00;      // (f,v) --> (f2, v2)

    for(int j = 0; j < N_U; j++)
    {
      double ansatz10 = Orig0[j];
      double ansatz01 = Orig1[j];
      double ansatz00 = Orig2[j];

      val  = mu_eff * (test10 * ansatz10 + test01 * ansatz01);     // mueff * (grad u, grad v) --> mueff * (v1_x * u1_x + v1_y * u1_y)
      val += mu/K * (ansatz00 * test00);                           // (u, v) --> (u1, v1)
      MatrixA11[i][j] += Mult * val;

      val  = mu_eff * (test10 * ansatz10 + test01 * ansatz01);     // mueff * (grad u, grad v) --> mueff * (v2_x * u2_x + v2_y * u2_y)
      val += mu/K * (ansatz00 * test00);                           // (u, v) --> (u2, v2)
      MatrixA22[i][j] += Mult * val;
    }

    for(int j = 0; j < N_P; j++)
    {
      ansatz00 = Orig3[j];

      val = -ansatz00 * test10;                      // -(p, div v) --> -p * v1_x
      MatrixB1T[i][j] += Mult * val;

      val = -ansatz00 * test01;                      // -(p, div v) --> -p * v2_y
      MatrixB2T[i][j] += Mult * val;
    }
  }

  for(int i = 0; i < N_P; i++)
  {
    test00 = Orig3[i];      // q

    Rhs_p[i] += Mult * (TDatabase::ParamDB->SIGN_MATRIX_BI) * test00 * c3;        // (g, q) (sign chosen according to stabilization)


    //introduce sources and sinks
    double P = 1;
    Rhs_p[i] += Mult  * (TDatabase::ParamDB->SIGN_MATRIX_BI) * test00 * approximate_delta_distribution_function * P;


    for(int j = 0; j < N_U; j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = TDatabase::ParamDB->SIGN_MATRIX_BI * test00 * ansatz10;                       // +/- (q, div u) --> +/- q * u1_x
      MatrixB1[i][j] += Mult * val;

      val = TDatabase::ParamDB->SIGN_MATRIX_BI * test00 * ansatz01;                       // +/- (q, div u) --> +/- q * u2_y
      MatrixB2[i][j] += Mult * val;
    }
   
   /*
   // consider sources and sinks impacts on the matrix as well 
    for(int j = 0; j < N_P; j++)
    {
      ansatz00 = Orig3[j];   // p

      val = approximate_delta_distribution_function * (ansatz00 * test00);               // source term  delta_fct * ( ph,  qh) 
      MatrixC[i][j] += Mult * val;
    }
  */

  }
}


// ======================================================================
// Standard Galerkin for Brinkman in [u gradq] formulation (macht nicht wirklich Sinn!!!!!!)
// ======================================================================
// To Do: Decide about deleting this assembling
void BrinkmanType2Galerkin(double Mult, double *coeff,
    double *, double,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  //double **MatrixA12 = LocMatrices[1];
  //double **MatrixA21 = LocMatrices[2];
  double **MatrixA22 = LocMatrices[3];
  ////int offset = TDatabase::ParamDB->NSTYPE == 14 ? 1 : 0;
  //double ** MatrixC = LocMatrices[4];
  int offset = 1;
  double **MatrixB1 = LocMatrices[4+offset];
  double **MatrixB2 = LocMatrices[5+offset];
  double **MatrixB1T = LocMatrices[6+offset];
  double **MatrixB2T = LocMatrices[7+offset];

  double *Rhs1 = LocRhs[0];                       // f1
  double *Rhs2 = LocRhs[1];                       // f2
  // double *Rhs_p = LocRhs[2];                   // g_q

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0];         // u_x
  double *Orig1 = OrigValues[1];         // u_y
  double *Orig2 = OrigValues[2];         // u
  //double *Orig3 = OrigValues[3];         // p
  double *Orig4 = OrigValues[4];         // p_x   MUSS NOCH DEFINIERT WERDEN !!!...eventuell
  double *Orig5 = OrigValues[5];         // p_y   MUSS NOCH DEFINIERT WERDEN!!!...eventuell

  //double c0 = coeff[0];                   // nu
  double c1 = coeff[1];                   // f1
  double c2 = coeff[2];                   // f2
  //double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
  double mu = coeff[4];                   // viscosity
  double mu_eff = coeff[5];               // effective viscosity
  double K = coeff[6];                    // permeability

  double val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;

  for(int i = 0; i < N_U; i++)
  {
    double *Matrix11Row = MatrixA11[i];
    // double *Matrix12Row = MatrixA12[i];
    // double *Matrix21Row = MatrixA21[i];
    double *Matrix22Row = MatrixA22[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(int j = 0; j < N_U; j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = mu_eff*(test10*ansatz10+test01*ansatz01);
      val += (mu/K)* (ansatz00*test00);
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = mu_eff*(test10*ansatz10+test01*ansatz01);
      val += (mu/K) *(ansatz00*test00);
      Matrix22Row[j] += Mult * val;
    }

    double *MatrixRow1 = MatrixB1T[i];
    double *MatrixRow2 = MatrixB2T[i];

    for(int j = 0; j < N_P; j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];

      val = Mult*ansatz10*test00;        // p_x*v
      MatrixRow1[j] += val;

      val = Mult*ansatz01*test00;        // p_y*v
      MatrixRow2[j] += val;
    }
  }

  for(int i = 0; i < N_P; i++)
  {
    double *MatrixRow1 = MatrixB1[i];
    double *MatrixRow2 = MatrixB2[i];

    test10 = Orig4[i];
    test01 = Orig5[i];

    for(int j = 0; j < N_U; j++)
    {
      ansatz00 = Orig2[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }
}

// ======================================================================
// Brinkman Problem with GLS Stabilization (only P1/P1)
// assume Delta u = Delta v = 0
// according to [Blank, Caiazzo, Chouly Lozinski, Mura; 2018]
// ======================================================================

void BrinkmanType1GalerkinResidualStabP1(double Mult, double *coeff,
    double *, double hK,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs)
{
  double ansatz00, ansatz10, ansatz01;    // ansatz functions
  double test00, test10, test01;          // test functions

  double ** MatrixA11 = LocMatrices[0];
  // double ** MatrixA12 = LocMatrices[1];
  // double ** MatrixA21 = LocMatrices[2];
  double ** MatrixA22 = LocMatrices[3];
  double ** MatrixC = LocMatrices[4];     // stabilization
  double ** MatrixB1 = LocMatrices[5];
  double ** MatrixB2 = LocMatrices[6];
  double ** MatrixB1T = LocMatrices[7];
  double ** MatrixB2T = LocMatrices[8];

  double * Rhs_u1 = LocRhs[0];            // f_v1
  double * Rhs_u2 = LocRhs[1];            // f_v2
  double * Rhs_p = LocRhs[2];             // f_q resp. g_q

  int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space
  int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space

  // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
  // Mult is the quadrature weight (wk)
  //double * Orig0 = OrigValues[0];          // u_x
  //double * Orig1 = OrigValues[1];          // u_y
  double * Orig2 = OrigValues[2];            // u
  //double * Orig3 = OrigValues[3];          // p
  double * Orig4 = OrigValues[4];            // p_x
  double * Orig5 = OrigValues[5];            // p_y

  // values defined in the example
  //double c0 = coeff[0];                   // t^2
  double c1 = coeff[1];                     // f1
  double c2 = coeff[2];                     // f2
  //double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
  double mu = coeff[4];                     // viscosity
  double mu_eff = coeff[5];                 // effective viscosity
  double K = coeff[6];                      // permeability
  double alpha = coeff[7];                  // equal_order_stab_weight
  //double t = fabs(sqrt((mu_eff/mu) * K));
  double PSPGStab; 
  //double PSPGStab =  alpha * (K/mu) * (hK*hK)/(c0+(hK*hK)); // stabilization = (hK*hK)/(c0*c0+hK*hK) 
  ///double PSPGStab =  alpha * (hK*hK)/(c0+(hK*hK)); // Without Sigma^{-1} as it is in Hannukainen and Sogn   
  ///double PSPGStab =  alpha * (hK*hK)/(c0); // stabilization acccording to Franca and Hughes (basically for stokes) but with t
  ///double PSPGStab =  alpha * (hK*hK)/(mu_eff); // stabilization acccording to Franca and Hughes (basically for stokes) with mueff
  ///double PSPGStab = -alpha * (hK*hK)/(c0+(hK*hK)) ; // stabilization formulation according to Volker ACHTUNG: Hierfür muss SIGN_MATRIX_BI=-1 sein!!!!
  if (TDatabase::ParamDB->l_T == 1)
  {
    PSPGStab = alpha * (hK*hK)/(mu_eff + (mu/K)* (hK*hK)); // Brinkman_P1P1.tex with l_T=h_T 
  }
  else if (TDatabase::ParamDB->l_T == -1)
  {
    //double L_0 = 1.; //0.1 * 1;
    PSPGStab = alpha * (hK*hK)/(mu_eff + (mu/K)* (TDatabase::ParamDB->L_0 * TDatabase::ParamDB->L_0)); // Brinkman_P1P1 with l_T=L_0
  }

  if (TDatabase::ParamDB->SIGN_MATRIX_BI == -1)
  { PSPGStab = alpha * (hK*hK);
    PSPGStab = PSPGStab * (-1);
  }

  double val;
  for(int i = 0; i < N_U; i++)
  {
    test00 = Orig2[i];

    Rhs_u1[i] += Mult * PSPGStab * mu/K * test00 * c1;               // stabilization: (f, sigma v) --> f1 * sigma v1
    Rhs_u2[i] += Mult * PSPGStab * mu/K * test00 * c2;               // stabilization: (f, sigma v) --> f2 * sigma v2

    for(int j = 0; j < N_U; j++)
    {
      ansatz00 = Orig2[j];

      val = PSPGStab * mu/K * mu/K * (ansatz00 * test00);                              // stabilization: (sigma u, sigma v) --> (sigma u1 * sigma v1)
      MatrixA11[i][j] += Mult * val;

      val = PSPGStab * mu/K * mu/K * (ansatz00 * test00);                              // stabilization: (sigma u, sigma v) --> (sigma u2 * sigma v2)
      MatrixA22[i][j] += Mult * val;
    }

    for(int j = 0; j < N_P; j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];

      val = PSPGStab * mu/K *(ansatz10 * test00);                          // stabilization: (grad p, sigma v) --> px * sigma v1
      MatrixB1T[i][j] += Mult * val;

      val = PSPGStab * mu/K * (ansatz01 * test00);                         // stabilization: (grad p, sigma v) --> py * sigma v2
      MatrixB2T[i][j] += Mult * val;
    }
  }

  for(int i = 0; i < N_P; i++)
  {
    test10 = Orig4[i];
    test01 = Orig5[i];

    Rhs_p[i] += Mult * PSPGStab * test10 * c1;                    // stabilization: (f, grad q) --> f1 * qx
    Rhs_p[i] += Mult * PSPGStab * test01 * c2;                    // stabilization: (f, grad q) --> f2 * qy

    for(int j = 0; j < N_U; j++)
    {
      ansatz00 = Orig2[j];

      val = PSPGStab * mu/K * (test10 * ansatz00);                     // stabilization: (sigma u, grad q) --> sigma u1 * qx
      MatrixB1[i][j] += Mult * val;

      val = PSPGStab * mu/K * (test01 * ansatz00);                     // stabilization: (sigma u, grad q) --> sigma u2 * qy
      MatrixB2[i][j] += Mult * val;
    }

    for(int j = 0; j < N_P; j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];

      val = PSPGStab * (ansatz10 * test10 + ansatz01 * test01);               // stabilization: (grad p, grad q) --> px * qx + py * qy
      MatrixC[i][j] += Mult * val;
    }
  }
}


// ======================================================================
// Brinkman Problem with PSPG Stabilization (only the terms for Pk/Pk which are not present for P1/P1)
// takes into account the Laplacian
// see [Blank, Caiazzo, Chouly, Lozinski, Mura; 2018]
// ======================================================================

void BrinkmanType1GalerkinResidualStabP2(double Mult, double *coeff,
    double *, double hK,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs)
{
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;    // ansatz functions
  double test00, test10, test01, test20, test02;              // test functions

  double ** MatrixA11 = LocMatrices[0];
  // double ** MatrixA12 = LocMatrices[1];
  // double ** MatrixA21 = LocMatrices[2];
  double ** MatrixA22 = LocMatrices[3];
  double ** MatrixB1  = LocMatrices[5];
  double ** MatrixB2  = LocMatrices[6];
  double ** MatrixB1T = LocMatrices[7];
  double ** MatrixB2T = LocMatrices[8];

  double * Rhs_u1 = LocRhs[0];            // f_u1           (rhs tested with v1)
  double * Rhs_u2 = LocRhs[1];            // f_u2           (rhs tested with v2)

  int N_U = N_BaseFuncts[0];              // number of basis functions for the velocity space (both for u1 and u2)
  int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space

  // values of FE functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
  // Mult is the quadrature weight (wk)
  double * Orig2 = OrigValues[2];          // u
  double * Orig4 = OrigValues[4];          // p_x
  double * Orig5 = OrigValues[5];          // p_y
  double * Orig6 = OrigValues[6];          // u_xx
  double * Orig7 = OrigValues[7];          // u_yy

  // values defined in the example
  //double c0 = coeff[0];                             // t^2= eps(:=1/Re) ; Re can be set in brinkman2d.dat
  double c1 = coeff[1];                               // f1
  double c2 = coeff[2];                               // f2
  // double c3 = coeff[3];                            // f3 (the rhs of incompressibility constraint)
  double mu = coeff[4];                               // viscosity
  double mu_eff = coeff[5];                           // effective viscosity
  double K = coeff[6];                                // permeability
  // double t = fabs(sqrt((mu_eff/mu)*K));
  double alpha = coeff[7];                            // PSPG Stabilization Parameter
  double PSPGStab; 
  //double PSPGStab =  alpha * (K/mu) * (hK*hK)/(c0+(hK*hK)); // stabilization = (hK*hK)/(c0*c0+hK*hK) ///warum negativ, wenn positiv geht es schief-NOCHMAL TESTEN
  ///double PSPGStab =  alpha * (hK*hK)/(c0+(hK*hK)); // Without Sigma^{-1} as it is in Hannukainen and Sogn   
  ///double PSPGStab =  alpha * (hK*hK)/(c0); // stabilization acccording to Franca and Hughes (basically for stokes) but with t
  ///double PSPGStab =  alpha * (hK*hK)/(mu_eff); // stabilization acccording to Franca and Hughes (basically for stokes) with mueff
  ///double PSPGStab = -alpha * (hK*hK)/(c0+(hK*hK)) ; // stabilization formulation according to Volker ACHTUNG: Hierfür muss SIGN_MATRIX_BI=-1 sein!!!!
  if (TDatabase::ParamDB->l_T == 1)
  {
    PSPGStab = alpha * (hK*hK)/(mu_eff + (mu/K)* (hK*hK)); // Brinkman_P1P1.tex with l_T=h_T 
  }
  else if (TDatabase::ParamDB->l_T == -1)
  {
    //double L_0 = 1.; //0.1 * 1;
    PSPGStab = alpha * (hK*hK)/(mu_eff + (mu/K)* (TDatabase::ParamDB->L_0 * TDatabase::ParamDB->L_0)); // Brinkman_P1P1 with l_T=L_0
  }

  if (TDatabase::ParamDB->SIGN_MATRIX_BI == -1)
  { // TODO LB: Do the stability analysis for the symmetric GLS stabilization and adapt the stabilization parameter  accordingly' 
    PSPGStab = alpha * (hK*hK);
    PSPGStab = PSPGStab * (-1);
  }

  double val;
  for(int i = 0; i < N_U; i++)
  {
    test00 = Orig2[i];          // v
    test20 = Orig6[i];          // v_xx
    test02 = Orig7[i];          // v_yy

    Rhs_u1[i] -= Mult * PSPGStab * mu_eff * c1 * (test20 + test02);  // - Mult * PSPGStab * mu_eff * c1 * (test20 + test02);
    Rhs_u2[i] -= Mult * PSPGStab * mu_eff * c2 * (test20 + test02);  // - Mult * PSPGStab * mu_eff * c2 * (test20 + test02);

    for(int j = 0; j < N_U; j++)
    {
      ansatz00 = Orig2[j];    // u
      ansatz20 = Orig6[j];    // u_xx
      ansatz02 = Orig7[j];    // u_yy

      val = PSPGStab * mu_eff * mu_eff * ((ansatz20 + ansatz02) * (test20 + test02));        // (-mu_eff)^2* (Delta u, Delta v) = (u1_xx + u1_yy)*(v1_xx + v1_yy)
      val += -PSPGStab * mu_eff * (mu/K) * (ansatz20 + ansatz02) * test00;                   // -mu_eff * (mu/K) * (Delta u,v) = (u1_xx + u1_yy) * v1
      val += -PSPGStab * mu_eff * (mu/K)  * (test20 + test02) * ansatz00;                    // -mu_eff * (mu/K) * (Delta v,u) = (v1_xx + v1_yy) * u1
      MatrixA11[i][j] += Mult * val;

      val = PSPGStab * mu_eff * mu_eff * ((ansatz20 + ansatz02) * (test20 + test02));        // (-mu_eff)^2* (Delta u, Delta v) = (u2_xx + u2_yy)*(v2_xx + v2_yy)
      val += -PSPGStab * mu_eff * (mu/K) * (ansatz20 + ansatz02) * test00;                   // -mu_eff * (mu/K) * (Delta u,v) = (u2_xx + u2_yy) * v2
      val += -PSPGStab * mu_eff * (mu/K) * (test20 + test02) * ansatz00;                     // -mu_eff * (mu/K) * (Delta v,u) = (v2_xx + v2_yy) * u2
      MatrixA22[i][j] += Mult * val;
    }

    for(int j = 0; j < N_P; j++)
    {
      ansatz10 = Orig4[j];        // p_x
      ansatz01 = Orig5[j];        // p_y

      val = -PSPGStab * mu_eff * (test20 + test02) * ansatz10;                           // -mu_eff * (Delta v, grad p) = -mu_eff * (v1_xx + v1_yy) * p_x
      MatrixB1T[i][j] += Mult * val;

      val = -PSPGStab * mu_eff * (test20 + test02) * ansatz01;                           // -mu_eff * (Delta v, grad p) = -mu_eff * (v2_xx + v2_yy) * p_y
      MatrixB2T[i][j] += Mult * val;
    }
  }

  for(int i = 0; i < N_P; i++)
  {
    test10 = Orig4[i];          // q_x
    test01 = Orig5[i];          // q_y

    for(int j = 0; j < N_U; j++)
    {
      ansatz20 = Orig6[j];        // u_xx
      ansatz02 = Orig7[j];        // u_yy

      val = -PSPGStab * mu_eff * (ansatz20 + ansatz02) * test10;                    // -mu_eff * (Delta u, grad q) = -mu_eff * (u1_xx + u1_yy) * q_x
      MatrixB1[i][j] += Mult * val;

      val = -PSPGStab * mu_eff * (ansatz20 + ansatz02) * test01;                     // -mu_eff * (Delta u, grad q) = -mu_eff * (u2_xx + u2_yy) * q_y
      MatrixB2[i][j] += Mult * val; 
    }
  }
}

// ======================================================================
// Grad-Div Stabilization
// Additional Terms:  grad_div_stab_weight * (div u - g, div v)
// ======================================================================
void BrinkmanGradDivStab(double Mult, double *coeff,
    double *, double hK,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs)
{
  double ansatz10, ansatz01;    // ansatz functions
  double test10, test01;          // test functions

  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA12 = LocMatrices[1];
  double ** MatrixA21 = LocMatrices[2];
  double ** MatrixA22 = LocMatrices[3];
  //double ** MatrixC = LocMatrices[4];
  //double ** MatrixB1 = LocMatrices[5];
  //double ** MatrixB2 = LocMatrices[6];
  //double ** MatrixB1T = LocMatrices[7];
  //double ** MatrixB2T = LocMatrices[8];

  double * Rhs_u1 = LocRhs[0];            // f_u1
  double * Rhs_u2 = LocRhs[1];            // f_u2
  //double * Rhs_p = LocRhs[2];           // g_q

  int N_U = N_BaseFuncts[0];                // number of basis functions for the velocity space
  //int N_P = N_BaseFuncts[1];              // number of basis functions for the pressure space

  // values of fe functions at quadrature points: Origvalues = f(uk,vk) (Gauß-quadrature: \int f(u,v) = sum wk * f(uk,vk) at Gauß points)
  // Mult is the quadrature weight (wk)
  double * Orig0 = OrigValues[0];          // u_x
  double * Orig1 = OrigValues[1];          // u_y
  //double * Orig2 = OrigValues[2];          // u
  //double * Orig3 = OrigValues[3];          // p

  // values defined in the example
  //double c0 = coeff[0];                   // t^2
  //double c1 = coeff[1];                   // f1
  //double c2 = coeff[2];                   // f2
  double c3 = coeff[3];                   // f3 (the rhs of incompressibility constraint)
  double mu = coeff[4];                   // viscosity
  double mu_eff = coeff[5];               // effective viscosity
  double K = coeff[6];                    // permeability
    
  //double grad_div_stab_weight = coeff[8]; // without considering the units
  //double grad_div_stab_weight = coeff[8] * (mu/K) * hK * hK; // units are fine
  //double grad_div_stab_weight = coeff[8] * mu_eff;  // units are fine
  double grad_div_stab;
  double delta = coeff[8];
  if (TDatabase::ParamDB->l_T == 1)
  {
    grad_div_stab = delta * (mu_eff + (mu/K) * hK * hK); // Brinkman_P1P1.tex with l_T=h_T 
  }
  else if (TDatabase::ParamDB->l_T == -1)
  {
  	grad_div_stab = delta * (mu_eff + (mu/K) * TDatabase::ParamDB->L_0 * TDatabase::ParamDB->L_0); // Brinkman_P1P1.tex with l_T=L_0
  }
  double val;

  if (TDatabase::ParamDB->SIGN_MATRIX_BI == -1)
  { // TODO LB: Do the stability analysis for the symmetric GLS stabilization and adapt the 'grad_div_stab_weight' accordingly.' 
  	grad_div_stab = coeff[8];
  }

  int GradDivStab_propto_h = 0;
  if (GradDivStab_propto_h)
  {
  	int i = 0;
  	if (i == 0)
  	{
  		Output::print("GradDiv-Stab is scaled as O(h).");
  		i += 1;
  	}
  	grad_div_stab = delta * hK;
  }

  // LB Debug New21.08.18 start:
  double approximate_delta_distribution_function = 0;
  if (TDatabase::ParamDB->SOURCE_SINK_FUNCTION)
  {
  	approximate_delta_distribution_function = coeff[9];
  }

  for(int i = 0; i < N_U; i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];

    Rhs_u1[i] += Mult * grad_div_stab * c3 * test10; //(g,div v) = g (dv_1/dx_1 + ...)
    Rhs_u2[i] += Mult * grad_div_stab * c3 * test01; //(g,div v) = g(... + dv_2/dx_2)

    Rhs_u1[i] += Mult * grad_div_stab * approximate_delta_distribution_function * test10; //(g,div v) = g (dv_1/dx_1 + ...)
    Rhs_u2[i] += Mult * grad_div_stab * approximate_delta_distribution_function * test01; //(g,div v) = g(... + dv_2/dx_2)

    for(int j = 0; j < N_U; j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = test10 * ansatz10;        // (div u,div v) = (dv_1/dx_1 * du_1/dx_1 + ... + ...)
      MatrixA11[i][j] += Mult * grad_div_stab * val;

      val = ansatz01 * test10;                                // ... = du_2/dx_2 * dv_1/dx_1
      MatrixA12[i][j] += Mult * grad_div_stab * val;

      val = test01 * ansatz01;        // (div u,div v) = (du_2/dx_2 * dv_2/dx_2 + ... + ...)
      MatrixA22[i][j] += Mult * grad_div_stab * val;

      val = ansatz10 * test01;                                // ... = du_1/dx_1 * dv_2/dx_2
      MatrixA21[i][j] += Mult * grad_div_stab * val;
    }
  }
}



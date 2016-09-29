// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (4*y*(1-y), 0)
// p(x,y) = x-1/2

void ExampleFile()
{
  Output::info<1>("Example", "Poiseuille.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 4*x*(1-x);
  values[1] = 0;
  values[2] = 4-8*x;
  values[3] = -8;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;//x-0.5;
  values[1] = 0;//1;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  if (i==0 || i==1 || i==3)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
//  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=4*Param*(1-Param);
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;
//  double rho = 7800;
//  double mu = 5e-3;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0;//1+8*eps; // f1
    coeff[2] = -9.8*parameters[i][2];// f2
    coeff[3] = parameters[i][2]; //rho; // density
    coeff[4] = parameters[i][3]; //mu; // dynamic viscosity
// cout << coeff[3] << " " << coeff[4] << endl;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void NonLinCoeffs(int n_points, double *x, double *y,
                  double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double c3 = 5.;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = param[0];
    coeff[2] = param[1];
    coeff[3] = c3;
  }
}


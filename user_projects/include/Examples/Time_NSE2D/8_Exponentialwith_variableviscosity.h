// Navier-Stokes problem, with variable viscosity

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;
double USER_parameter3;

void ExampleFile()
{
  Output::print<1>("Example: Variable_viscosity.h");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double a = TDatabase::ParamDB->P1;
  double h = TDatabase::ParamDB->P2;
  double mu0 = TDatabase::ParamDB->P3;
  double val;

    val = h *exp(-a*a)/(2*a*a*mu0);

  values[0] = val*(exp(a*a*(1-(y-h)*(y-h)/(h*h)))-1);
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double a = TDatabase::ParamDB->P1;
  double h = TDatabase::ParamDB->P2;
  double mu0 = TDatabase::ParamDB->P3;
  double val;

  val = h *exp(-a*a)/(2*a*a*mu0);
  
  values[0] = exp(a*a*(1-(y-h)*(y-h)/(h*h)))-1;
  values[1] = 0.0;
  values[2] = -a*a*2*(y-h)/(h*h) *exp(a*a*(1-(y-h)*(y-h)/(h*h)));
  values[3] = -a*a*2/(h*h) *exp(a*a*(1-(y-h)*(y-h)/(h*h)))
              + a*a*2*(y-h)/(h*h) * a*a*2*(y-h)/(h*h) * exp(a*a*(1-(y-h)*(y-h)/(h*h))) ; //Laplacian
              
  values[0] *= val;
  values[1] *= val;
  values[2] *= val;
  values[3] *= val;

}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
//  if (i==1) 
//  {
//      cond = NEUMANN;
//  }
//  else
    cond = DIRICHLET;  
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double a = TDatabase::ParamDB->P1;
  double h = TDatabase::ParamDB->P2;
  double mu0 = TDatabase::ParamDB->P3;
  double val;

  val = h *exp(-a*a)/(2*a*a*mu0);
  
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=val*(exp(a*a*(1-(Param-h)*(Param-h)/(h*h)))-1); // outflow, homogenous NEUMANN or DIRICHLET
            break;
    case 2: value=0;
            break;
    case 3: value=val*(exp(a*a*(1-(Param-h)*(Param-h)/(h*h)))-1); // inflow, DIRICHLET u_x
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0.0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
//  static double eps = 1./TDatabase::ParamDB->RE_NR;
  double eps;
  double derivative;
  int i;
  double *coeff;

  double a = TDatabase::ParamDB->P1;
  double h = TDatabase::ParamDB->P2;
  double mu0 = TDatabase::ParamDB->P3;
  
  double val1[4];   // U1-Exact function and its derivatives
  double val2[4];   // U2-Exact function and its derivatives
  double val3[4];   // P-Exact function and its derivatives

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    
//    eps = 1; //mu0 * exp(a*a*(y[i]-h)*(y[i]-h)/(h*h));
    eps = mu0 * exp(a*a*(y[i]-h)*(y[i]-h)/(h*h));
    // derivative1 = d/dy(eps*val1[2])
    derivative = -1;


    ExactU1(x[i], y[i], val1);
    ExactU2(x[i], y[i], val2);
    ExactP (x[i], y[i], val3);

    coeff[0] = eps;
    coeff[1] = -derivative; //eps*val1[3];// + val3[1]; // f1 non dimensional case
    coeff[2] = 0;//-eps*val2[3] + val3[2]; // f2 non dimensional case

    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
    {
      coeff[1] += val1[0]*val1[1] + val2[0]*val1[2]; // f1 nondimensional
      coeff[2] += val1[0]*val2[1] + val2[0]*val2[2]; // f2 nondimensional
    }
//    coeff[0] = 1;//eps; // THE VARIABLE VISCOSITY SHOULD BE INJECTED IN FRONT OF THE DIFFUSIVE TERM
//    coeff[1] = 2; // f1
//    coeff[2] = 0; // f2
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
//void NonLinCoeffs(int n_points, double *x, double *y,
//                  double **parameters, double **coeffs)
//{
//  static double eps = 1./TDatabase::ParamDB->RE_NR;
//  int i;
//  double *coeff, *param;
//
//  for(i=0;i<n_points;i++)
//  {
//    coeff = coeffs[i];
//    param = parameters[i];
//
//    coeff[0] = eps;
//    coeff[1] = param[0];
//    coeff[2] = param[1];
//  }
//}


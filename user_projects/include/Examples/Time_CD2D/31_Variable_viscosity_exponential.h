
void ExampleFile()
{
  OutPut("Example: 31_Variable_viscosity_exponential.h" << endl);
}

double get_nu()
{
  return 1;
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = false;

// exact solution
void Exact(double x, double y, double *values)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;
  double a = TDatabase::ParamDB->P1;
  double h = TDatabase::ParamDB->P2;
  double mu0 = TDatabase::ParamDB->P3;

  values[0] = mu0*exp(a*a*(y-h)*(y-h)/(h*h));
  values[1] = 0;
  values[2] = mu0*(a*a*2*(y-h)/(h*h))*exp(a*a*(y-h)*(y-h)/(h*h));
  values[3] = mu0*(2*a*a/(h*h) + a*a*2*(y-h)/(h*h)*a*a*2*(y-h)/(h*h))*exp(a*a*(y-h)*(y-h)/(h*h));
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double Param, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;
  double a = TDatabase::ParamDB->P1;
  double h = TDatabase::ParamDB->P2;
  double mu0 = TDatabase::ParamDB->P3;

  switch(BdComp)
  {
    case 0:
      value = mu0*exp(a*a);
    break;
    case 1:
      value = mu0*exp(a*a*(Param-h)*(Param-h)/(h*h));
    break;
    case 2:
      value =  mu0*exp(a*a);
    break;
    case 3:
      value = mu0*exp(a*a*(1-Param-h)*(1-Param-h)/(h*h));
    break;
  } // endswitch
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;
  double a = TDatabase::ParamDB->P1;
  double h = TDatabase::ParamDB->P2;
  double mu0 = TDatabase::ParamDB->P3;

  values[0] = mu0*exp(a*a*(y-h)*(y-h)/(h*h));
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
//  double eps=1/TDatabase::ParamDB->RE_NR;
//  double a=1, b=2, c=1;
  int i;
  double *coeff, y;//x, y;
//  double t=TDatabase::TimeDB->CURRENTTIME;

  double a = TDatabase::ParamDB->P1;
  double h = TDatabase::ParamDB->P2;
  double mu0 = TDatabase::ParamDB->P3;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    x = X[i];
    y = Y[i];

    coeff[0] = 0;
    coeff[1] = 1;
    coeff[2] = 1;
    coeff[3] = 0;

    coeff[4] = mu0*(a*a*2*(y-h)/(h*h))*exp(a*a*(y-h)*(y-h)/(h*h));
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;
  double a = TDatabase::ParamDB->P1;
  double h = TDatabase::ParamDB->P2;
  double mu0 = TDatabase::ParamDB->P3;

  values[0] = mu0*exp(a*a*(y-h)*(y-h)/(h*h));
  values[1] = 0;
  values[2] = mu0*(a*a*2*(y-h)/(h*h))*exp(a*a*(y-h)*(y-h)/(h*h));
  values[3] = mu0*(2*a*a/(h*h) + a*a*2*(y-h)/(h*h)*a*a*2*(y-h)/(h*h))*exp(a*a*(y-h)*(y-h)/(h*h));
}


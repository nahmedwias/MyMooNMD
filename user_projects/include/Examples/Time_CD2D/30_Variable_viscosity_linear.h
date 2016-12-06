
void ExampleFile()
{
  OutPut("Example: 30_Variable_viscosity_linear.h" << endl);
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

  values[0] = 1000*y+700;
  values[1] = 0;
  values[2] = 1000;
  values[3] = 0;
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

  switch(BdComp)
  {
    case 0:
      value = 700;
    break;
    case 1:
      value = 1000*Param+700;
    break;
    case 2:
      value = 1700;
    break;
    case 3:
      value = 1000*(1-Param)+700;
    break;
  } // endswitch
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1000*y+700;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
//  double eps=1/TDatabase::ParamDB->RE_NR;
//  double a=1, b=2, c=1;
  int i;
  double *coeff;// x, y;
//  double t=TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    x = X[i];
//    y = Y[i];

    coeff[0] = 0;
    coeff[1] = 1;
    coeff[2] = 1;
    coeff[3] = 0;

    coeff[4] = 1000;
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1000*y+700;
  values[1] = 1000;
  values[2] = 0;
  values[3] = 0;
}


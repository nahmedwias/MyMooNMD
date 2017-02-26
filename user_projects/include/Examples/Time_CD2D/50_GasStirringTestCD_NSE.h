// This is a test to put gas inflow in a square box filled with liquid
// Domain is UnitSquares
// No Exact solution is known, so it is 0.
// Only Initial solution is known.
// Boundary conditions have to be taken care of.


void ExampleFile()
{
  Output::info<3>("Example: 50_GasStirringTestCD_NSE.h");
}

double get_nu()
{
  return 1;
}

constexpr bool rhs_depends_on_time = false;
constexpr bool coefficients_depend_on_time = true;

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double Param, BoundCond &cond)
{
  if (BdComp == 1 || BdComp==3)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
//  double rho_min = TDatabase::ParamDB->P2;
//  if (BdComp==1 || BdComp==2)
//    value = rho_min;  // this is Dirichlet
//  else

  if (BdComp == 1 || BdComp==3)
  {
    if (t = 0)
      value = 1;
    else
      value = 0;
  }
  else
    value = 0;
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  double x0 = TDatabase::ParamDB->P3; // x position of initial circle
  double y0 = TDatabase::ParamDB->P4; // y position of initial circle
  double radius = TDatabase::ParamDB->P5; // radius of circle
  if ( y <= 0.7 )
    values[0] = 1; // liquid bath
  else
    values[0] = 0; // gas at the top of liquid


//  if ( sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) <= radius )
//  {
//    values[0] = 0; // gas circle in the bottom of the liquid
//  }

//  double A = 1;
//  double a = 20;
//  double b = 0;
//  double c = 20;
//  double x0 = 0.5;
//  double y0 = 0.7;
//  values[0]= A*exp(-a*pow(x-x0,2)+2*b*(x-x0)*(y-y0)-c*pow(y-y0,2));
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
//  double eps=1/TDatabase::ParamDB->RE_NR;
//  double a=1, b=2, c=1;
  double u_x,u_y;
  int i;
  double *coeff;
//  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    x = X[i];
//    y = Y[i];

    u_x = parameters[i][0];
    u_y = parameters[i][1];

    coeff[0] = 0;
    coeff[1] = u_x;
    coeff[2] = u_y;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

// exact solution
//void Initial(double x, double y, double *values)
//{
//  if ( x <= 0.3 && y <= 0.7)
//    values[0] = 1;
//  else
//    values[0] =0;
//}


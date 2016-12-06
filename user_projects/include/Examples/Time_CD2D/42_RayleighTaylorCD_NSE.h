// This is a test to try to define a domain with 2 phases in the square box
// This should reproduce the Rayleigh-Taylor instability!
// No Exact solution is known, so it is 0.
// Only Initial solution is known.
// Boundary conditions have to be taken care of.


void ExampleFile()
{
  Output::info<3>("Example: 42_RayleighTaylorCD_NSE.h");
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
void BoundCondition(int i, double Param, BoundCond &cond)
{
  if (i == 0 || i == 2 )
    cond = DIRICHLET;  // top and bottom
  else
    cond = NEUMANN;    // left and right
}

// value of boundary condition
void BoundValue(int i, double Param, double &value)
{
  if (i == 2)
    value = 1;  // top
  else
    value = 0;
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  if ( y >= 0.5)
    values[0] = 1;
  else
    values[0] = 0;
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
void Initial(double x, double y, double *values)
{
  if ( x <= 0.3 && y <= 0.7)
    values[0] = 1;
  else
    values[0] =0;
}


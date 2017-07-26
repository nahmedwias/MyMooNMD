// Solitary Wave Benchmark Problem - 2-way coupling - Example TCD2D
// See Yue et al 2003, and also Rodriguez thesis 2002.
// Data taken from Yue et al 2003.
// Geometry = SolitaryWave.PRM and SolitaryWave_quad.GEO
// 20*h x 2*h rectangle, where h = 0.1, according to wave
// velocity Cw= 1m/s = sqrt(g*h)
// h is depth of still water. Re= Cw*h/visco_water = 5e4.
// visco ratio air/water=15, density ratio air/water = 1.2e-3
// Initial water surface = Boussinesq profile
// A(x,0) = A0/cosh^2(sqrt(3A0)x/2), where A0 is initial height
// Tested situation = A0/h = 4
// Initial solution is known. The problem is time dependent.
// Only gravity.



void ExampleFile()
{
  Output::info<3>("Example: 44_SolitaryWave_CD_NSE.h");
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
  cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int i, double Param, double &value)
{
  value = 0;
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  double pi = 3.14159265358979;
  double h = 0.1;
  double A0 = TDatabase::ParamDB->P1;

  double phi = y - A0/(cosh(sqrt(3*A0/h)*x/(2*h))*cosh(sqrt(3*A0/h)*x/(2*h))); // initial wave profile

  double eps = 1.2e-2;
  if (phi < -eps)
    values[0] = 1;
  else if ( phi > eps)
    values[0] = 0;
  else
    values[0] =0.5*(1 - (phi/eps) - (1/pi)*sin(phi*pi/eps));
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
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
{//unused function
}


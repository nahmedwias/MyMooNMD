// Example 51 = Gas Stirring in a Ladle using Mazumdar & Guthrie approach,
// where effect of gas plume is modelled through a volume force in the right
// hand side. The gas phase is thus not present in the model, but its effect
// is added through some kind of buyoancy force, depending on gas fraction.

// Note that this can be used as a simple NSE model with extra force
// or with a 2-phase VOF where the second phase is the top slag layer


void ExampleFile()
{
  Output::info<3>("Example: 51_MazumdarGuthrieCD_NSE.h");
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
//  if (BdComp == 2)
    cond = NEUMANN;
//  else
//    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
//  if (BdComp==2)
    value = 0;
//  else
//    value = 1;  // liquid
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  double slag_thickness = 0.05;
  double geometry_height = 1.;

  if (y>=geometry_height - slag_thickness)
    values[0] = 0;
  else
    values[0] = 1;
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
{
}


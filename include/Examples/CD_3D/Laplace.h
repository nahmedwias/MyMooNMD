// ======================================================================
// Sine problem 3D
// ======================================================================
// #include <ConvDiff3D.h>

void ExampleFile()
{
  Output::root_info<1>("Example", "Laplace.h");
}

double PECLET_NUMBER;
// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[1] = Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[2] = Pi*sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[3] = Pi*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[4] = -3*Pi*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  double values[5];
  Exact(x, y, z, values);
  value = values[0];
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
                    double **, double **coeffs)
{
  const double eps = PECLET_NUMBER;
  double values[5];
  for(int i = 0; i < n_points; i++)
  {
    Exact(x[i], y[i], z[i], values);
    coeffs[i][0] = eps;
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;
    coeffs[i][4] = 0;
    coeffs[i][5] = -coeffs[i][0] * values[4] + coeffs[i][1]*values[1]
                   + coeffs[i][2]*values[2] + coeffs[i][3]*values[3]
                   + coeffs[i][4]*values[0];
  }
}


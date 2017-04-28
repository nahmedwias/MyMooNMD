// ==========================================================================
// instationary problem
// ==========================================================================
void ExampleFile()
{
  Output::print<1>("Example: 20_CouplingNSE_CD_Quadratic.h\n");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*t*(x*x+y*z);
  values[1] = 2*x*t*t;
  values[2] = z*t*t;
  values[3] = y*t*t;
  values[4] = 2*t*t;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = t*t*(x*x+y*z);
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*t*(x*x+y*z);
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps = 1;
  double *coeff;
  double t=TDatabase::TimeDB->CURRENTTIME;
  double b1,b2,b3; // convection field (4,3,2)
  
  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    double x = X[i];
    double y = Y[i];
    double z = Z[i];

    b1 = parameters[i][0]; // ux
    b2 = parameters[i][1]; // uy
    b3 = parameters[i][2]; // uz

    coeff[0] = eps; //diffusion coefficient
    coeff[1] = b1;   //ux
    coeff[2] = b2;   //uy
    coeff[3] = b3;   //uy
    coeff[4] = 1;   //reaction coefficient
    coeff[5] = 2*t*(x*x+y*z) -eps * 2*t*t + b1 * 2*x * t*t+ b2 * z *t*t+ b3 * y*t*t + 1 * t*t*(x*x+y*z); //rhs
  }
}


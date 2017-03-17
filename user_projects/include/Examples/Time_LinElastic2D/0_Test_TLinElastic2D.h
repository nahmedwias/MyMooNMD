// ======================================================================
// instationary problem
// ======================================================================


void ExampleFile()
{
  Output::print<1>("Example: 0_Test_TLinElastic2D.h");
//  TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 0;
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = true;

// exact solution
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// initial displacement and velocity
void InitialU1(double x, double y, double *values)
{
  values[0] = x+y;
}
void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}
void InitialV1(double x, double y, double *values)
{
  values[0] = 0;
}
void InitialV2(double x, double y, double *values)
{
  values[0] = 0;
}

// kind of boundary condition (for FE space needed)
// this supposes the same condition for all components u1,u2
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValueU1(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
    break;
    case 1: value = 0;
    break;
    case 2: value = 0;
    break;
    case 3: value = 0;
    break;
  }
}
void BoundValueU2(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
    break;
    case 1: value = 0;
    break;
    case 2: value = 0;
    break;
    case 3: value = 0;
    break;
  }
}

// Coefficients and right hand side
void Coefficients(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;
//  double x, y;// eps, b1, b2, c;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

//    x = X[i];
//    y = Y[i];

    coeff[0] = 0;  // Lame coefficient Lambda, if used
    coeff[1] = 0;  // Lame coefficient mu, if used
    coeff[2] = 0;  // RHS1
    coeff[3] = 0;  // RHS2

  }
}



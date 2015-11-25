// Navier-Stokes problem, sin and cos solution
// 

void ExampleFile()
{
  Output::print<1>("Example: Simple example with sin and cos solution and p=0");
}

// exact solution
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] =          cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[1] =      -Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[2] =       Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[3] =       Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[4] = -3*Pi*Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z); //Laplacien
}
void ExactU2(double x, double y,  double z, double *values)
{
  values[0] =          sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[1] =       Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[2] =      -Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[3] =       Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[4] = -3*Pi*Pi*sin(Pi*x)*cos(Pi*y)*sin(Pi*z); //Laplacien
}
void ExactU3(double x, double y,  double z, double *values)
{
  values[0] =      -2*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[1] =   -2*Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[2] =   -2*Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[3] =    2*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[4] = 6*Pi*Pi*sin(Pi*x)*sin(Pi*y)*cos(Pi*z); //Laplacien
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}
// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  double tol = 1e-10;
  if((fabs(1+x) < tol) || (fabs(1+y) < tol) || (fabs(1+z) < tol)
       || (fabs(1-z) < tol))
    cond = DIRICHLET;
  else
    cond = NEUMANN;

  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;
}
// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  const double eps = 1./TDatabase::ParamDB->RE_NR;
  double tol = 1e-10;
  if((fabs(1+x) < tol) || (fabs(1+y) < tol) || (fabs(1+z) < tol)
       || (fabs(1-z) < tol))
    value = cos(Pi*x)*sin(Pi*y)*sin(Pi*z); //Dirichlet
  else{
    if(fabs(1-x) < tol)
    {
      value = -eps*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z); //Neumann
      if(TDatabase::ParamDB->LAPLACETYPE == 1)
        value += -eps*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
    }
    else
    {
      value = eps*Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z); //Neumann
      if(TDatabase::ParamDB->LAPLACETYPE == 1)
        value += eps*Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
    }
  }
}
void U2BoundValue(double x, double y, double z, double &value)
{
  const double eps = 1./TDatabase::ParamDB->RE_NR;
  double tol = 1e-10;
  if((fabs(1+x) < tol) || (fabs(1+y) < tol) || (fabs(1+z) < tol)
       || (fabs(1-z) < tol))
    value = sin(Pi*x)*cos(Pi*y)*sin(Pi*z); //Dirichlet
  else{
    if(fabs(1-x) < tol)
    {
      value = eps*Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z); //Neumann
      if(TDatabase::ParamDB->LAPLACETYPE == 1)
        value += eps*Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
    }
    else
    {
      value = -eps*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z); //Neumann
      if(TDatabase::ParamDB->LAPLACETYPE == 1)
        value += -eps*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
    }
  }
}
void U3BoundValue(double x, double y, double z, double &value)
{
  const double eps = 1./TDatabase::ParamDB->RE_NR;
  double tol = 1e-10;
  if((fabs(1+x) < tol) || (fabs(1+y) < tol) || (fabs(1+z) < tol)
       || (fabs(1-z) < tol))
    value = -2*sin(Pi*x)*sin(Pi*y)*cos(Pi*z); //Dirichlet
  else{
    if(fabs(1-x) < tol)
    {
      value = -eps*2*Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z); //Neumann
      if(TDatabase::ParamDB->LAPLACETYPE == 1)
        value += eps*Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
    }
    else
      {
      value = -eps*2*Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z); //Neumann
      if(TDatabase::ParamDB->LAPLACETYPE == 1)
        value += eps*Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
      }
  }
}

void LinCoeffs(int n_points, double * X, double * Y, double * Z,
               double **parameters, double **coeffs)
{
  const double eps = 1./TDatabase::ParamDB->RE_NR;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] =  eps;
    coeffs[i][1] =  eps*3*Pi*Pi*cos(X[i]*Pi)*sin(Y[i]*Pi)*sin(Z[i]*Pi); // f1
    coeffs[i][2] =  eps*3*Pi*Pi*sin(X[i]*Pi)*cos(Y[i]*Pi)*sin(Z[i]*Pi); // f2
    coeffs[i][3] = -eps*6*Pi*Pi*sin(X[i]*Pi)*sin(Y[i]*Pi)*cos(Z[i]*Pi); // f3
    coeffs[i][4] = 0; // g
  }
}

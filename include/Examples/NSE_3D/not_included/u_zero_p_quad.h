void ExampleFile()
{
  OutPut("Example: solution not in ansatz space u=(0,0,0), p=x^2 \n");
}

// exact solution
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] =  0;
  values[1] =  0;
  values[2] =  0;
  values[3] =  0;
  values[4] =  0; //Laplacien
}
void ExactU2(double x, double y,  double z, double *values)
{
  values[0] =  0;
  values[1] =  0;
  values[2] =  0;
  values[3] =  0;
  values[4] =  0; //Laplacien
}
void ExactU3(double x, double y,  double z, double *values)
{
  values[0] =  0;
  values[1] =  0;
  values[2] =  0;
  values[3] =  0;
  values[4] =  0; //Laplacien
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = x*x;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}
// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  const double tol = 1e-10;
  if(abs(1+x) < tol)
    cond = DIRICHLET;
  else
    cond = NEUMANN;

  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;

}
// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  const double tol = 1e-10;
    if(abs(1+x) < tol) //x==-1
      value = 0; //Dirichlet
    else
      if (abs(1-x) < tol) //x==1
        value = -x*x; //Neumann
      else
        value = 0; //Neumann
}
void U2BoundValue(double x, double y, double z, double &value)
{
  const double tol = 1e-10;
    if(abs(1+x) < tol) //x==-1
      value = 0; //Dirichlet
    else
      if (abs(1+y) < tol) //y==-1
        value = x*x; //Neumann
      else
        if (abs(1-y) < tol) //y==1
          value = -x*x; //Neumann
        else
          value = 0; //Neumann
}
void U3BoundValue(double x, double y, double z, double &value)
{
  const double tol = 1e-10;
    if(abs(1+x) < tol) //x==-1
      value = 0; //Dirichlet
    else
      if (abs(1+z) < tol) //z==-1
        value = x*x; //Neumann
      else
        if (abs(1-z) < tol) //z==1
          value = -x*x; //Neumann
        else
          value = 0; //Neumann
}

void LinCoeffs(int n_points, double * X, double * Y, double * Z,
               double **parameters, double **coeffs)
{
  const double eps = 1./TDatabase::ParamDB->RE_NR;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = eps;
    coeffs[i][1] = 2*X[i]; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // f3
    coeffs[i][4] = 0; // g
  }
}

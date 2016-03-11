void ExampleFile()
{
  OutPut("Example: solution not in ansatz space u=(x^3+3xy^2,-y^3-3x^2y,0), p=0 \n");
}

// exact solution
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] =  x*x*x+3*x*y*y;
  values[1] =  3*x*x+3*y*y;
  values[2] =  6*x*y;
  values[3] =  0;
  values[4] =  12*x; //Laplacien
}
void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = -y*y*y-3*y*x*x;
  values[1] = -6*x*y;
  values[2] = -3*y*y-3*x*x;
  values[3] =  0;
  values[4] =  -12*y; //Laplacien
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
  values[0] = 0;
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
  const double eps = 1./TDatabase::ParamDB->RE_NR;
  const double tol = 1e-10;
    if(abs(1+x) < tol) //x==-1
      value = x*x*x+3*x*y*y; //Dirichlet
    else
      if(abs(1-x) < tol) //x==1
      {
        value = eps*(3*x*x+3*y*y); //Neumann
        if(TDatabase::ParamDB->LAPLACETYPE == 1)
          value += eps*(3*x*x+3*y*y);
      }
      else
        if(abs(1+y) < tol) //y==-1
        {
          value = -eps*6*x*y; //Neumann
          if(TDatabase::ParamDB->LAPLACETYPE == 1)
            value += eps*6*x*y;
        }
        else
          if(abs(1-y) < tol) //y==1
          {
            value = eps*6*x*y; //Neumann
            if(TDatabase::ParamDB->LAPLACETYPE == 1)
              value += -eps*6*x*y;
          }
          else
            value = 0; //Neumann
}
void U2BoundValue(double x, double y, double z, double &value)
{
  const double eps = 1./TDatabase::ParamDB->RE_NR;
  const double tol = 1e-10;
    if(abs(1+x) < tol) //x==-1
      value = -y*y*y-3*y*x*x; //Dirichlet
    else
      if(abs(1-x) < tol) //x==1
      {
        value = -eps*6*x*y; //Neumann
        if(TDatabase::ParamDB->LAPLACETYPE == 1)
          value += eps*6*x*y;
      }
      else
        if(abs(1+y) < tol) //y==-1
        {
          value = eps*(3*y*y+3*x*x); //Neumann
          if(TDatabase::ParamDB->LAPLACETYPE == 1)
            value += eps*(3*y*y+3*x*x);
        }
        else
          if(abs(1-y) < tol) //y==1
          {
            value = -eps*(3*y*y+3*x*x); //Neumann
            if(TDatabase::ParamDB->LAPLACETYPE == 1)
              value += -eps*(3*y*y+3*x*x);
          }
          else
            value = 0; //Neumann
}
void U3BoundValue(double x, double y, double z, double &value)
{
  const double tol = 1e-10;
    if(abs(1+x) < tol) //x==-1
      value = 0; //Dirichlet
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
    coeffs[i][1] = -eps*12*X[i]; // f1
    coeffs[i][2] = eps*12*Y[i]; // f2
    coeffs[i][3] = 0; // f3
    coeffs[i][4] = 0; // g
  }
}

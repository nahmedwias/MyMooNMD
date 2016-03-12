/**
 * Simple example for Navier-Stokes 3D, used for development, debugging and testing.
 * 
 * Linear velocity solution u = (x,-y,0).
 * Linear pressure solution p = x.
 * Dirichlet and Neumann boundary conditions.
 * 
 * @author ???, Clemens Bartsch imported this from MooNMD
 * @date 2016/03/11 Import to ParMooN.
 * 
 * 
 * @todo Remove global Database dependency.
 */

void ExampleFile()
{
  OutPut("Example: test_u_1_p_1_x. Linear velocity solution u=(x,-y,0), linear pressure solution p=x \n");
}

// exact solution
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] =  x;
  values[1] =  1;
  values[2] =  0;
  values[3] =  0;
  values[4] =  0; //Laplacien
}
void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = -y;
  values[1] =  0;
  values[2] = -1;
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
  values[0] = x;
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
      value = -1; //Dirichlet
    else
      if(abs(1-x) < tol) //x==1
      {
        value = eps-1; //Neumann
        if(TDatabase::ParamDB->LAPLACETYPE == 1)
          value += eps;
      }
      else
       value = 0; //Neumann
}
void U2BoundValue(double x, double y, double z, double &value)
{
  const double eps = 1./TDatabase::ParamDB->RE_NR;
  const double tol = 1e-10;
    if(abs(1+x) < tol) //x==-1
      value = -y; //Dirichlet
    else
      if(abs(1+y) < tol) //y==-1
      {
        value = eps+x; //Neumann
        if(TDatabase::ParamDB->LAPLACETYPE == 1)
          value += eps;
      }
      else
        if(abs(1-y) < tol) //y==1
        {
          value = -eps-x; //Neumann
          if(TDatabase::ParamDB->LAPLACETYPE == 1)
            value += -eps;
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
      if(abs(1+z) < tol) //z==-1
        value = x; //Neumann
      else
        if(abs(1-z) < tol) //z==1
          value = -x; //Neumann
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
    coeffs[i][1] = 1; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // f3
    coeffs[i][4] = 0; // g
  }
}

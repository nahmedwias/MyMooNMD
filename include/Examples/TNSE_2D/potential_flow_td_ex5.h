// u1 = sin(t) sin (PiX) sin(PiY)
// u2 = sin(t) cos(PiX) cos(PiY)
//
double DIMENSIONLESS_VISCOSITY;
void ExampleFile()
{
  OutPut("Example: potential_flow_ex5.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*t*t*(6.*x*y) + (1.-t*t*t)*(-3.*x*x + 3.*y*y);
}

void InitialU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*t*t*(3.*x*x-3.*y*y) + (1.-t*t*t)*(6.*x*y);
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0.;
}
// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*t*t*(6.*x*y) + (1.-t*t*t)*(-3.*x*x + 3.*y*y);
  values[1] = t*t*t*(6.*y) + (1.-t*t*t)*(-6.*x);
  values[2] = t*t*t*(6.*x) + (1.-t*t*t)*(6.*y);
  values[3] = 0.0; // laplacian
  values[4] = 3.*t*t*(6.*x*y) - 3.*t*t*(-3.*x*x + 3.*y*y); // time derivative
}

void ExactU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*t*t*(3.*x*x-3.*y*y) + (1.-t*t*t)*(6.*x*y);
  values[1] = t*t*t*(6.*x) + (1.-t*t*t)*(6.*y);
  values[2] = t*t*t*(-6.*y)+ (1.-t*t*t)*(6.*x);
  values[3] = 0.0; // laplacian
  values[4] = 3.*t*t*(3.*x*x-3.*y*y) - 3.*t*t*(6.*x*y);
}

void ExactP(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = 3.*t*t*(-3.*x*x*y+ y*y*y - x*x*x + 3.*x*y*y);
  values[1] = 3.*t*t*(-6.*x*y - 3.*x*x + 3.*y*y);
  values[2] = 3.*t*t*(-3.*x*x + 3.*y*y + 6.*x*y);
}
// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double x=0, y=0;  
  switch(BdComp)
  {
    case 0: 
      x = Param; y=0;
      break;
    case 1: 
      x = 1; y = Param;
      break;
    case 2: 
      x = 1-Param; y = 1;
      break;
    case 3: 
      x = 0; y = 1-Param;
      break;
    default: cout << "wrong boundary part number" << endl;
      break;
  }
  double val[5];
  ExactU1(x,y,val);
  value=val[0];
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double x=0, y=0;
  
  switch(BdComp)
  {
    case 0: 
      x = Param; y=0;
      break;
    case 1: 
      x = 1; y = Param;
      break;
    case 2: 
      x = 1-Param; y = 1;
      break;
    case 3: 
      x = 0; y = 1-Param;
      break;
    default: cout << "wrong boundary part number" << endl;
      break;
  }
  double val[5];
  ExactU2(x,y,val);
  value=val[0];
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double nu = DIMENSIONLESS_VISCOSITY;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = nu;    
    coeffs[i][1] = 0;
    coeffs[i][2] = 0; 
  }
}



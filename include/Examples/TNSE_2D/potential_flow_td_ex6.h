// u1 = sin(t) sin (PiX) sin(PiY)
// u2 = sin(t) cos(PiX) cos(PiY)
//
double DIMENSIONLESS_VISCOSITY;
void ExampleFile()
{
  Output::print("Example: potential_flow_ex6.h");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*(4.*x*x*x - 12.*x*y*y) + 4.*(1.-t)*(3.*x*x*y-y*y*y);
}

void InitialU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*(-12.*x*x*y + 4.*y*y*y) + 4.*(1.-t)*(x*x*x-3.*x*y*y);
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
  values[0] = t*(4.*x*x*x - 12.*x*y*y) + 4.*(1.-t)*(3.*x*x*y-y*y*y);  
  values[1] = t*(12.*x*x - 12.*y*y) + 4.*(1.-t)*(6.*x*y);
  values[2] = t*(-24.*x*y) + 4.*(1.-t)*(3.*x*x-3.*y*y);
  values[3] = 0.0; // laplacian
  values[4] = (4.*x*x*x - 12.*x*y*y) -4.*(3.*x*x*y-y*y*y); // time derivative
}

void ExactU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*(-12.*x*x*y + 4.*y*y*y) + 4.*(1.-t)*(x*x*x-3.*x*y*y);
  values[1] = t*(-24.*x*y) + 4.*(1.-t)*(3.*x*x-3.*y*y);
  values[2] = t*(-12.*x*x + 12.*y*y) + 4.*(1.-t)*(-6.*x*y);
  values[3] = 0.0; // laplacian
  values[4] = (-12.*x*x*y + 4.*y*y*y) - 4.*(x*x*x-3.*x*y*y);
}

void ExactP(double x, double y, double *values)
{
  values[0] = 6.*x*x*y*y - pow(x,4) - pow(y,4) -4./15. 
              + y*pow(x,3) - x*pow(y,3);
  values[1] = 12.*x*y*y - 4.*pow(x,3) 
              + 3.*y*pow(x,2) - pow(y,3);
  values[2] = 12.*x*x*y - 4.*pow(y,3) 
              + pow(x,3) - 3.*x*pow(y,2);
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
  static double nu = DIMENSIONLESS_VISCOSITY;
  
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = nu;    
    coeffs[i][1] = 0;
    coeffs[i][2] = 0; 
  }
}



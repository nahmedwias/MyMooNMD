// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (t*y, 2*t*y)
// p(x,y) = 0

void ExampleFile()
{
  OutPut("Example: linear_space_cubic_time.h" << endl) ;
}
double DIMENSIONLESS_VISCOSITY;
// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*t*t*y;
  values[1] = 0.;
  values[2] = t*t*t;
  values[3] = 0.;
  values[4] = 3.*t*t*y;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*t*x;
  values[1] = t*t;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 2.*t*x;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
}
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double v[5];
  ExactU1(x,y,v);
  values[0] = v[0];  
}

void InitialU2(double x, double y, double *values)
{
  double v[5];
  ExactU2(x,y,v);
  values[0] = v[0];
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
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
  double x, y;
  switch(BdComp)
  {
    case 0: x = Param; y = 0.;       break;
    case 1: x = 1.;    y = Param;    break;
    case 2: x = 1. - Param; y = 1.;  break;
    case 3: x = 0.;  y = 1. - Param; break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
  double v[5];
  ExactU1(x,y,v);
  value = v[0];
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double x=0, y=0;
  switch(BdComp)
  {
    case 0: x = Param; y = 0.;       break;
    case 1: x = 1.;    y = Param;    break;
    case 2: x = 1. - Param; y = 1.;  break;
    case 3: x = 0.;  y = 1. - Param; break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
  double v[5];
  ExactU2(x,y,v);
  value = v[0];
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = DIMENSIONLESS_VISCOSITY;
  double u1[5], u2[5], p[4];
double t=TDatabase::TimeDB->CURRENTTIME;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = nu;
    ExactU1(X[i], Y[i], u1);
    ExactU2(X[i], Y[i], u2);
    ExactP(X[i], Y[i], p);

    //coeffs[i][1] = u1[4] - nu*u1[3] + u1[0]*u1[1] + u2[0]*u1[2] + p[1];
    //coeffs[i][2] = u2[4] - nu*u2[3] + u1[0]*u2[1] + u2[0]*u2[2] + p[2];
    
coeffs[i][1] =3.*t*t*Y[i];
coeffs[i][2] =2.*t*X[i];
    
    if( (TDatabase::ParamDB->DISCTYPE == RESIDUAL_VMS) && 
        (TDatabase::TimeDB->CURRENTTIME > TDatabase::TimeDB->STARTTIME))
    {
      double ct = TDatabase::TimeDB->CURRENTTIME;
      double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      double pt = ct - dt;
      TDatabase::TimeDB->CURRENTTIME = pt;
      double u1pt[5], u2pt[5], ppt[4];
      ExactU1(X[i], Y[i], u1pt);
      ExactU2(X[i], Y[i], u2pt);
      ExactP(X[i], Y[i], ppt);
      
      coeffs[i][3] = u1pt[4] - nu*u1pt[3] + u1pt[0]*u1pt[1] + u2pt[0]*u1pt[2] + ppt[1];
      coeffs[i][4] = u2pt[4] - nu*u2pt[3] + u1pt[0]*u2pt[1] + u2pt[0]*u2pt[2] + ppt[2];
      TDatabase::TimeDB->CURRENTTIME = ct;
      
    }
    
  }
}
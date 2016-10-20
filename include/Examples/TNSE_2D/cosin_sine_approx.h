#ifndef COSIN_SINE_APPROX_H
#define COSIN_SINE_APPROX_H
#ifndef COSINE_SINE_H
#define COSINE_SINE_H

// u1 = cos(t) * ( sin(\pi * x -0.7) *sin( \pi * y +0.2 ) );
// u2 = cos(t) * ( cos(\pi * x -0.7) *cos( \pi * y +0.2 ) );
// p =  cos(t) * ( sin(x) * cos(y) * (cos(1) * sin(1) - sin(1)) );
//
double DIMENSIONLESS_VISCOSITY;
void ExampleFile()
{
  Output::print("Example: cosine_sin_approx.h") ;
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = (1-t*t)*x*y;
}

void InitialU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0]= (1-t*t)*(1-y*y/2);
}

void InitialP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = (1-t*t) * ( sin(x) * cos(y) + (cos(1) * sin(1) - sin(1)) );
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = (1-t*t)*(x*y);
  values[1] = (1-t*t)*(y);
  values[2] = (1-t*t)*(x);
}

void ExactU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = (1-t*t)*(1-y*y/2.);
  values[1] = 0;
  values[2] = -(1-t*t)*y;
}

void ExactP(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] =  (1-t*t) * ( (x-(x*x*x)/6.0) * (1-(y*y)/2.0));
  values[1] =  (1-t*t) * ( (1.0-3.0*(x*x)/6.0) * (1.-(y*y)/2.0));
  values[2] =  (1-t*t) * ( (x-(x*x*x)/6.0) * (-2.*y/2.0));
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
  double t = TDatabase::TimeDB->CURRENTTIME;
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
  value= (1-t*t)*(x*y);
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
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
 value=(1-t*t)*(1-y*y/2.0);
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = DIMENSIONLESS_VISCOSITY;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double *coeff, x, y;
  double /*u1, u1x, u1y, */u1xx, u1yy;
  double /*u2, u2x, u2y, */u2xx, u2yy;
  double px, py;
  double u1t, u2t;
double valp[4];
  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    ExactP(X[i], Y[i], valp);

    // u1   = cos(t) * ( sin(Pi * x -0.7) *sin( Pi * y +0.2 ) );;
    // u1x  = cos(t) * ( Pi * cos(Pi * x -0.7) *sin( Pi * y +0.2 ) );
    u1xx = -sin(t) * ( Pi *Pi * sin(Pi * x -0.7) *sin( Pi * y +0.2 ) );
    // u1y  = cos(t) * ( sin(Pi * x -0.7) * Pi* cos( Pi * y +0.2 ) );
    u1yy = -sin(t) * ( sin(Pi * x -0.7) * Pi * Pi* sin( Pi * y +0.2 ) );

    u1t = cos(t) * ( sin(Pi * x -0.7) *sin( Pi * y +0.2 ) );;

    // u2   = cos(t) * ( cos(Pi * x -0.7) *cos( Pi * y +0.2 ) );
    // u2x  = -cos(t) * ( Pi * sin(Pi * x -0.7) *cos( Pi * y +0.2 ) );
    u2xx = -sin(t) * ( Pi * Pi* cos(Pi * x -0.7) *cos( Pi * y +0.2 ) );
    // u2y = -cos(t) * ( cos(Pi * x -0.7) *Pi *sin( Pi * y +0.2 ) );
    u2yy = -sin(t) * ( cos(Pi * x -0.7) *Pi *Pi *cos( Pi * y +0.2 ) );

    u2t = cos(t) * ( cos(Pi * x -0.7) *cos( Pi * y +0.2 ) );

    coeff[0] = nu;

    coeff[1] = -2*t*(x*y)      - 0 + valp[1] ;//+px;
    coeff[2] = -2*t*(1-y*y/2.) + nu*(1-t*t)  + valp[2];// py;

  }
}

#endif // COSINE_SINE_H

#endif // COSIN_SINE_APPROX_H

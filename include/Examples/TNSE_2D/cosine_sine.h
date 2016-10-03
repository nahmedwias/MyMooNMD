#ifndef COSINE_SINE_H
#define COSINE_SINE_H

// u1 = cos(t) * ( sin(\pi * x -0.7) *sin ( \pi * y +0.2 ) );
// u2 = cos(t) * ( cos(\pi * x -0.7) *cos ( \pi * y +0.2 ) );
// p =  cos(t) * ( sin(x) * cos(y) * (cos(1) * sin (1) - sin(1)) );
//
double DIMENSIONLESS_VISCOSITY;
void ExampleFile()
{
  Output::print("Example: cosine_sin.h") ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = cos(t) * ( sin(Pi * x -0.7) *sin ( Pi * y +0.2 ) );
}

void InitialU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0]= cos(t) * ( cos(Pi * x -0.7) *cos ( Pi * y +0.2 ) );
}

void InitialP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = cos(t) * ( sin(x) * cos(y) + (cos(1) * sin (1) - sin(1)) );
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = cos(t) * ( sin(Pi * x -0.7) *sin ( Pi * y +0.2 ) );;
  values[1] = cos(t) * ( Pi * cos(Pi * x -0.7) *sin ( Pi * y +0.2 ) );
  values[2] = cos(t) * ( sin(Pi * x -0.7) * Pi* cos ( Pi * y +0.2 ) );
}

void ExactU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = cos(t) * ( cos(Pi * x -0.7) *cos ( Pi * y +0.2 ) );
  values[1] = -cos(t) * ( Pi * sin(Pi * x -0.7) *cos ( Pi * y +0.2 ) );
  values[2] = -cos(t) * ( cos(Pi * x -0.7) *Pi *sin ( Pi * y +0.2 ) );
}

void ExactP(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] =  cos(t) * ( sin(x) * cos(y) + (cos(1) * sin (1) - sin(1)) );
  values[1] =  cos(t) * ( cos(x) * cos(y) );
  values[2] = -cos(t) * ( sin(x) * sin(y) );
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
  value= cos(t) * ( sin(Pi * x -0.7) *sin ( Pi * y +0.2 ) );
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
 value=cos(t) * ( cos(Pi * x -0.7) *cos ( Pi * y +0.2 ) );
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

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    // u1   = cos(t) * ( sin(Pi * x -0.7) *sin ( Pi * y +0.2 ) );;
    // u1x  = cos(t) * ( Pi * cos(Pi * x -0.7) *sin ( Pi * y +0.2 ) );
    u1xx = -cos(t) * ( Pi *Pi * sin(Pi * x -0.7) *sin ( Pi * y +0.2 ) );
    // u1y  = cos(t) * ( sin(Pi * x -0.7) * Pi* cos ( Pi * y +0.2 ) );
    u1yy = -cos(t) * ( sin(Pi * x -0.7) * Pi * Pi* sin ( Pi * y +0.2 ) );

    u1t = -sin(t) * ( sin(Pi * x -0.7) *sin ( Pi * y +0.2 ) );;

    // u2   = cos(t) * ( cos(Pi * x -0.7) *cos ( Pi * y +0.2 ) );
    // u2x  = -cos(t) * ( Pi * sin(Pi * x -0.7) *cos ( Pi * y +0.2 ) );
    u2xx = -cos(t) * ( Pi * Pi* cos(Pi * x -0.7) *cos ( Pi * y +0.2 ) );
    // u2y = -cos(t) * ( cos(Pi * x -0.7) *Pi *sin ( Pi * y +0.2 ) );
    u2yy = -cos(t) * ( cos(Pi * x -0.7) *Pi *Pi *cos ( Pi * y +0.2 ) );

    u2t = -sin(t) * ( cos(Pi * x -0.7) *cos ( Pi * y +0.2 ) );

   //p= cos(t) * ( sin(x) * cos(y) + (cos(1) * sin (1) - sin(1)) );
    px = cos(t) * ( cos(x) * cos(y) );
    py = -cos(t) * ( sin(x) * sin(y) );
    coeff[0] = nu;

    coeff[1] = u1t - nu * (u1xx + u1yy) + px;
    coeff[2] = u2t - nu * (u2xx + u2yy) + py;

  }
}

#endif // COSINE_SINE_H

// Navier-Stokes problem, Driven Cavity
// 
// u(x,y) = ?
// p(x,y) = ?

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;
double pi = 3.14159265358979;

void ExampleFile()
{
  Output::print<1>("Example: DrivenCavity.h");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
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

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double t0 = 0.1;   // parameter for time regularization
  double x1 = 0.1; // distance of BC regularization
  // see Volkers book "Finite Element Methods for
  // Incompressible Flow Problems" - Appendix D

  switch(BdComp)
  {
    case 0:
            value=0;
            break;
    case 1:
            value=0;
            break;
    case 2:  // top moving side velocity
           if ( Param <= x1)
              value = 1-0.25*pow((1-cos((x1-Param)*pi/x1)),2);
            else if (Param >= 1-x1)
              value = 1-0.25*pow((1-cos((Param-(1-x1))*pi/x1)),2);
            else
              value = 1;
            break;
    case 3:
            value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  if (t0 > 0.)
  {
    if (t <= t0)
      value = value*t/t0;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double nu = DIMENSIONLESS_VISCOSITY;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = nu;

    coeff[1] = 0;  // f1
    coeff[2] = 0;  // f2
  }
}



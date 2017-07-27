// Example 51 = Gas Stirring in a Ladle using Mazumdar & Guthrie approach,
// where effect of gas plume is modelled through a volume force in the right
// hand side. The gas phase is thus not present in the model, but its effect
// is added through some kind of buyoancy force, depending on gas fraction.

// Note that this can be used as a simple NSE model with extra force
// or with a 2-phase VOF where the second phase is the top slag layer

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;



void ExampleFile()
{
  Output::info<3>("Example: 51_MazumdarGuthrieNSE_CD.h ") ;
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
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
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  switch(i)
  {
    case 2:
      cond = DIRICHLET;
      break;
    case 0: case 1: case 3:
      cond = DIRICHLET;
//      cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
//      TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
//      TDatabase::ParamDB->FRICTION_TYPE = 1;
//      TDatabase::ParamDB->FRICTION_CONSTANT = 0.0;
//      TDatabase::ParamDB->PENETRATION_CONSTANT = 1.e12;
//      TDatabase::ParamDB->PENETRATION_POWER= -2;
      break;
      /**  END CODE FOR GEOMETRY DAMBREAK_VERYHIGH **/
    default:
      ErrThrow("Unknown BdPart");
  }
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;
}



void U1BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
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
  static double nu = REYNOLDS_number;
//  double t = TDatabase::TimeDB->CURRENTTIME;
  int i;
  double *coeff;
  double  x, y;

  double alpha; // gas fraction depending on x,y
  double Q = 0.004;  double L = 2;  double R = 2;
  double Up = 4.4*pow(Q,0.33)*pow(L,0.25)/pow(R,0.25);
  double rmin; // radius of plume
  double plume_profile; // equation of x and y
  double plume_height = 0.8;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

//    plume_profile = y-5*x+2.5;
    plume_profile = y+5*x-2.5;
    if (y <= plume_height)
    {
      if (y+5*x-2.5 >=0)
      {
        if (y-5*x+2.5 >=0)
          plume_profile = -2*y/5;
        else
          plume_profile = 1.e8;
      }
      else
        plume_profile = 1.e8;
    }
    else
      plume_profile = 1.e8;

//    plume_profile = y-100*(x-0.5)*(x-0.5);
//    plume_profile = x-0.5;
    if ( fabs(plume_profile) <= 1.e-1)
      rmin = 1e-1;
    else
      rmin = plume_profile;

    alpha = Q/(Up*3.14*rmin*rmin);

    coeff[0] = nu;
    coeff[1] = 0;     // f1
    coeff[2] = -10*(1-alpha);     // f2
    coeff[3] = 0;
  }
}

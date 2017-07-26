// Solitary Wave Benchmark Problem - 2-way coupling - Example TNSE2D
// See Yue et al 2003, and also Rodriguez thesis 2002.
// Data taken from Yue et al 2003.
// Geometry = SolitaryWave.PRM and SolitaryWave_quad.GEO
// 20*h x 2*h rectangle, where h = 0.1, according to wave
// velocity Cw= 1m/s = sqrt(g*h)
// h is depth of still water. Re= Cw*h/visco_water = 5e4.
// visco ratio air/water=15, density ratio air/water = 1.2e-3
// Initial water surface = Boussinesq profile
// A(x,0) = A0/cosh^2(sqrt(3A0)x/2), where A0 is initial height
// Tested situation = A0/h = 4
// Initial solution is known. The problem is time dependent.
// Only gravity.

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;



void ExampleFile()
{
  Output::info<3>("Example: 44_SolitaryWave_NSE_CD.h ") ;
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
    case 0: case 1: case 2: case 3:
    case 4: case 5: case 6: case 7:
    case 8: case 9:
      cond = DIRICHLET;  // bottom
      break;
    case 11: case 12: case 13: case 14:
    case 15: case 16: case 17: case 18:
    case 19: case 20:
      cond = NEUMANN;  // top
      break;
    case 10: case 21:  // sides
      cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
      TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
      TDatabase::ParamDB->FRICTION_TYPE = 1;
      TDatabase::ParamDB->FRICTION_CONSTANT = 0.0;
      TDatabase::ParamDB->PENETRATION_CONSTANT = 1.e12;
      TDatabase::ParamDB->PENETRATION_POWER= -2;
      break;
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
  double t = TDatabase::TimeDB->CURRENTTIME;
  int i;
  double *coeff;
//  double  x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    x = X[i];
//    y = Y[i];
//    double rho = parameters[i][2];
//    double mu  = parameters[i][3];

    coeff[0] = nu;

    // Navier-Stokes
    coeff[1] = 0;     // f1
    coeff[2] = -10;   // f2
    coeff[3] = 0;
  }
}


void irgendeinePostProcessFunktion(Time_NSE2D& time_nse2d)
{
}

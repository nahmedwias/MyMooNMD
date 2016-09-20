// u1 = sin(t) sin (PiX) sin(PiY)
// u2 = sin(t) cos(PiX) cos(PiY)
//
void ExampleFile()
{
  OutPut("Example: SinSinSin.h" << endl) ;
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
  values[0]= 0;
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
    double t = TDatabase::TimeDB->CURRENTTIME;
    values[0] = t*(20*x*x*x*y-20*x*y*y*y);
    values[1] = t*(60*x*x*y - 20*y*y*y);
    values[2] = t*(20*x*x*x - 60*x*y*y);
}

void ExactU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;    
  values[0] = t*(5*x*x*x*x + 5*y*y*y*y -30*x*x*y*y);
  values[1] = t*(20*x*x*x - 60*x*y*y);
  values[2] = t*(20*y*y*y - 60*x*x*y);
}

void ExactP(double x, double y, double *values)
{
  values[0] = -(5*x*x*x*x*y  + y*y*y*y*y - 10*x*x*y*y*y);
  values[1] = -(20*x*x*x*y  - 20*x*y*y*y);
  values[2] = -(5*x*x*x*x  + 5*y*y*y*y - 30*x*x*y*y);
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
  value=t*(20*x*x*x*y-20*x*y*y*y);
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
 value=t*(5*x*x*x*x+5*y*y*y*y -30*x*x*y*y);
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = 1/TDatabase::ParamDB->RE_NR;
  
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = nu;    
    coeffs[i][1] = 0;
    coeffs[i][2] = 0; 
  }
}



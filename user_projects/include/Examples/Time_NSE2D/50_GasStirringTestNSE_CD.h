// 

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;



void ExampleFile()
{
  Output::info<3>("Example: 50_GasStirringTestNSE_CD.h ") ;
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
     case 0:  // bottom
       cond = DIRICHLET;
       break;
     case 1: case 2: case 3: // the rest
     case 4: case 5: case 6:
       cond = DIRICHLET;
       break;
     default:
       ErrThrow("Unknown BdComp in example 50.");
   }

//  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0:  // bottom
      value = 0;
      break;
    case 1: case 2: case 3: // the rest
    case 4: case 5: case 6:
      value = 0;
      break;
    default:
      ErrThrow("Unknown BdComp in example 50.");
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double x;
  double first_plug1 = -0.485;
  double first_plug2 = -0.385;
  double second_plug1 = -first_plug2; // plugs are symmetric
  double second_plug2 = -first_plug1;
  switch(BdComp)
  {
    case 0:  // bottom
      x =  -1.329 + Param*(2.658);  // x for bottom part
      if ( ( x >= first_plug1 && x <= first_plug2 ) ||
          ( x >= second_plug1 && x <= second_plug2 ) )
      {
        value = 1;
      }
      else
        value = 0;
      break;
    case 1: case 2: case 3: // the rest
    case 4: case 5: case 6:
      value = 0;
      break;
    default:
      ErrThrow("Unknown BdComp in example 50.");
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = REYNOLDS_number;
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

/*
    // Stokes
    coeff[1] =0;
    coeff[2] =0;
*/

    // Navier-Stokes
    coeff[1] = 0;     // f1
    coeff[2] = -10;     // f2
    coeff[3] = 0;
  }
}

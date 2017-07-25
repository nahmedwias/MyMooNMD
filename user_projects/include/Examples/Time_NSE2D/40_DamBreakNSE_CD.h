// 

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;



void ExampleFile()
{
  Output::info<3>("Example: 40_DamBreakNSE_CD.h ") ;
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
//  switch(i)
//  {
//    case 0: case 1: case 2: case 3:
//    case 4: case 5: case 6: case 7:
//    case 8: case 9:
//      cond = DIRICHLET;
//      break;
//    default:
//      ErrThrow("Unknown BdPart");
//  }

  switch(i)
   {
    /** those lines are for the geometry dambreak_VERYhigh, with
     * Neumann condition at the top
     */
    case 2:
      cond = NEUMANN;
      break;
    case 0: case 1: case 3:
      cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
      TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
      TDatabase::ParamDB->FRICTION_TYPE = 1;
      TDatabase::ParamDB->FRICTION_CONSTANT = 0.0;
      TDatabase::ParamDB->PENETRATION_CONSTANT = 1.e12;
      TDatabase::ParamDB->PENETRATION_POWER= -2;
      break;
      /**  END CODE FOR GEOMETRY DAMBREAK_VERYHIGH **/

//    /** those lines are for the geometry dambreak_high, with
//      * Neumann condition at the top
//      */
//     case 3: case 4:
//       cond = NEUMANN;
//       break;
//     case 0: case 1: case 2:
//     case 5:
//       cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
//       TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
//       TDatabase::ParamDB->FRICTION_TYPE = 1;
//       TDatabase::ParamDB->PENETRATION_CONSTANT = 1.e12;
//       TDatabase::ParamDB->PENETRATION_POWER= -2;
//       break;
//     /**  END CODE FOR GEOMETRY DAMBREAK_HIGH **/

    /** those lines are for the standard geometry dambreak
     * */
//     case 5: case 6: case 7:
//     case 8:
//     cond = NEUMANN;
//       break;
//     case 0: case 1: case 2: case 3:
//     case 9: case 4:
//       cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
//       TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
//       TDatabase::ParamDB->FRICTION_TYPE = 1;
//       TDatabase::ParamDB->FRICTION_CONSTANT = 0.0;
//       TDatabase::ParamDB->PENETRATION_CONSTANT = 1.e12;
//       TDatabase::ParamDB->PENETRATION_POWER= -2;
//       break;
//    /**  END CODE FOR GEOMETRY STANDARD DAMBREAK **/

//    /** those lines are for the geometry dambreak OPENFOAM
//     * */
//     case 6:
//     cond = NEUMANN;
//       break;
//     case 0: case 1: case 2:
//     case 3: case 4: case 5:
//     case 7:
//       cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
//       TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
//       TDatabase::ParamDB->FRICTION_TYPE = 1;
//       TDatabase::ParamDB->PENETRATION_CONSTANT = 1.e12;
//       TDatabase::ParamDB->PENETRATION_POWER= -2;
//       break;
//    /**  END CODE FOR GEOMETRY DAMBREAK OPENFOAM **/

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

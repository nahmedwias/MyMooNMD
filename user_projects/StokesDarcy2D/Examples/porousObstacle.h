/*
An example where the Darcy subdomain is entirely enclosed in the Stokes 
subdomain. No known solution.
*/

void ExampleFile()
{
  Output::print<1>("\nExample: porousObstacle.h");
}

// ============================================================================
// exact solution
// ============================================================================
void ExactP_Darcy(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}
void ExactU1_Darcy(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}
void ExactU2_Darcy(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}
void ExactU1_NSE(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}
void ExactU2_NSE(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}
void ExactP_NSE(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ============================================================================
// boundary conditions
// ============================================================================
void BC_Darcy(int bdComp, double t, BoundCond &cond)  
{
  ErrThrow("BoundCondition_Darcy: should never be called");
}
void BC_Velocity_Darcy(int bdComp, double t, BoundCond &cond)
{
  ErrThrow("BC_Velocity_Darcy: should never be called");
}
void BC_Pressure_Darcy(int bdComp, double t, BoundCond &cond)
{
  ErrThrow("BC_Velocity_Darcy: should never be called");
  cond = NEUMANN;
}
void BoundCondition_NSE(int bdComp, double t, BoundCond &cond)  
{
  cond = DIRICHLET;
  if(bdComp==3)
    cond  = NEUMANN;
}
void BoundConditionPressure_NSE(int bdComp, double t, BoundCond &cond)  
{
  cond = NEUMANN;
}

// ============================================================================
// boundary values
// ============================================================================
void BoundValue_Darcy(int bdComp, double t, double &value)
{
  ErrThrow("BoundValue_Darcy: should never be called");
}
void BoundValueVelocity_Darcy(int bdComp, double t, double &value)
{
  ErrThrow("BoundValueVelocity_Darcy: should never be called");
}
void BoundValuePressure_Darcy(int bdComp, double t, double &value)
{
  ErrThrow("BoundValuePressure_Darcy: should never be called");
}
void U1BoundValue_NSE(int bdComp, double t, double &value)
{
  switch(bdComp)
  {
    case 7:
      value = t*(1-t); // DIRICHLET
      break;
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
      value = 0.;
      break;
    default:
      ErrThrow("U1BoundValue_NSE: wrong boundary part number");
      break;
  }
}
void U2BoundValue_NSE(int bdComp, double t, double &value)
{
  value = 0.;
}
void PressureBoundValue_NSE(int bdComp, double t, double &value)
{
  value = 0;
  //Output::print("PressureBoundValue_NSE: bdComp ", bdComp, ", t ", t);
}

// ============================================================================
// coefficients
// ============================================================================
void LinCoeffs_NSE(int n_points, double * X, double * Y,
               double **parameters, double **coeffs)
{
  double nu= 1/TDatabase::ParamDB->RE_NR;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = nu;
    coeffs[i][1] = 0.; // f1
    coeffs[i][2] = 0.; // f2
    coeffs[i][3] = 0.; // g
  }
}
void LinCoeffs_Darcy(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = K; // diffusion
    coeffs[i][1] = 0; // convection_x
    coeffs[i][2] = 0; // convection_y
    coeffs[i][3] = 0; // reaction
    coeffs[i][4] = 0; // f_D
  }
}
void LinCoeffs_DarcyMixed(int n_points, double *X, double *Y,
               double **parameters, double **coeffs) // for mixed formulation
{
  const double K_1 = 1.0/TDatabase::ParamDB->SIGMA_PERM;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = K_1;
    coeffs[i][1] = 0.; // f1
    coeffs[i][2] = 0.; // f2
    coeffs[i][3] = 0.; // f_D
  }
}

// ============================================================================
// initial solution for iteration
// ============================================================================
void InitialP_Darcy(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;;
}

double prescribedDirichletOutflow()
{
  Output::print("\n\n\n WARNING\n\n Function call to prescribedDirichletOutflow()");
  return 0.0;
}

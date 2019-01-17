// Navier-Stokes problem with sine and cosine functions
// 

double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::print<1>("Example: 2d Stress (Navier-)Stokes simple.h.");
}

// ========================================================================
// exact solution
// ========================================================================
void Stress_XX(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void Stress_XY(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void Stress_YY(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU1(double x, double y, double *values)
{
  values[0] = 1.;
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
void BoundConditionStress(int i, double t, BoundCond &cond)
{
  cond = NEUMANN;
}

void transform(const int BdComp, const double Param, double& x, double& y, 
               double& nx, double& ny)
{
  switch(BdComp)
  {
    case 0:
      x = Param;
      y = 0.;
      nx = 0.;
      ny = -1.;
      break;
    case 1:
      x = 1.;
      y = Param;
      nx = 1.;
      ny = 0.;
      break;
    case 2:
      x = 1. - Param;
      y = 1.;
      nx = 0.;
      ny = 1.;
      break;
    case 3:
      x = 0.;
      y = 1. - Param;
      nx = -1.;
      ny = 0.;
      break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
}

void BoundValue_SXX(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = DIMENSIONLESS_VISCOSITY;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundConditionStress(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double sxx[4];
  double sxy[4];
  double sxz[4];
  Stress_XX(x, y, sxx);
  Stress_XY(x, y, sxy);
  Stress_YY(x, y, sxz);
  
  if(cond == DIRICHLET)
    value = sxx[0];
  else
  {
    ErrThrow("NEUMANN b.c's are not supported");
//     // NEUMANN
//     value = nu * (nx * u1[1] + ny * u1[2]) - sxz[0] * nx;
//     if(lt == 1)
//     {
//       value *= 0.5;
//       value += 0.5 * nu * (nx * u1[1] + ny * sxy[1]);
//     }
  }
  return;
}
void BoundValue_SXY(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = DIMENSIONLESS_VISCOSITY;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundConditionStress(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double sxx[4];
  double sxy[4];
  double sxz[4];
  Stress_XX(x, y, sxx);
  Stress_XY(x, y, sxy);
  Stress_YY(x, y, sxz);
  
  if(cond == DIRICHLET)
    value = sxy[0];
  else
  {
    ErrThrow("NEUMANN b.c's are not supported");
//     // NEUMANN
//     value = nu * (nx * u1[1] + ny * u1[2]) - sxz[0] * nx;
//     if(lt == 1)
//     {
//       value *= 0.5;
//       value += 0.5 * nu * (nx * u1[1] + ny * sxy[1]);
//     }
  }
  return;
}

void BoundValue_SYY(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = DIMENSIONLESS_VISCOSITY;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundConditionStress(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double sxx[4];
  double sxy[4];
  double sxz[4];
  Stress_XX(x, y, sxx);
  Stress_XY(x, y, sxy);
  Stress_YY(x, y, sxz);
  
  if(cond == DIRICHLET)
    value = sxz[0];
  else
  {
    ErrThrow("NEUMANN b.c's are not supported");
//     // NEUMANN
//     value = nu * (nx * u1[1] + ny * u1[2]) - sxz[0] * nx;
//     if(lt == 1)
//     {
//       value *= 0.5;
//       value += 0.5 * nu * (nx * u1[1] + ny * sxy[1]);
//     }
  }
  return;
}

void BoundConditionNS(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = DIMENSIONLESS_VISCOSITY;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundConditionNS(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double u1[4];
  double u2[4];
  double p[4];
  ExactU1(x, y, u1);
  ExactU2(x, y, u2);
  ExactP(x, y, p);
  if(cond == DIRICHLET)
    value = u1[0];
  else
  {
    // NEUMANN
    value = nu * (nx * u1[1] + ny * u1[2]) - p[0] * nx;
    if(lt == 1)
    {
      value *= 0.5;
      value += 0.5 * nu * (nx * u1[1] + ny * u2[1]);
    }
  }
  return;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = DIMENSIONLESS_VISCOSITY;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundConditionNS(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double u1[4];
  double u2[4];
  double p[4];
  ExactU1(x, y, u1);
  ExactU2(x, y, u2);
  ExactP(x, y, p);
  if(cond == DIRICHLET)
    value = u2[0];
  else
  {
    // NEUMANN
    value = nu * (nx * u2[1] + ny * u2[2]) - p[0] * ny;
    if(lt == 1)
    {
      value *= 0.5;
      value += 0.5 * nu * (nx * u1[2] + ny * u2[2]);
    }
  }
  return;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  const double nu = DIMENSIONLESS_VISCOSITY;
  double val1[4];
  double val2[4];
  double val3[4];
  double valsx[4];
  double valsy[4];
  double valsz[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = nu;
    
    ExactU1(X[i], Y[i], val1);
    ExactU2(X[i], Y[i], val2);
    ExactP(X[i], Y[i], val3);
    
    Stress_XX(X[i], Y[i], valsx);
    Stress_XY(X[i], Y[i], valsy);
    Stress_YY(X[i], Y[i], valsz);
    
    
    
    coeffs[i][1] = -nu*val1[3] + val3[1]; // f1
    coeffs[i][2] = -nu*val2[3] + val3[2]; // f2
    
    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
    {
      coeffs[i][1] += val1[0]*val1[1] + val2[0]*val1[2]; // f1
      coeffs[i][2] += val1[0]*val2[1] + val2[0]*val2[2]; // f2
    }
    coeffs[i][3] = val1[1] + val2[2]; // g (divergence)

    // additional coefficient (used only in the Brinkman problem)
    coeffs[i][4] = 0.;
    coeffs[i][5] = 0.;
    coeffs[i][6] = 0.;
  }
  
}

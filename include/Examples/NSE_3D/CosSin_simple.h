/**
 * @file A Navier--Stokes test problem with sine-cosine solution.
 *
 * The boundary data is adapted to the [0.1]^3 unit cube example, which is
 * availabe as default geometry in ParMooN. It will throw an error if you
 * try running it on any other domain - just to make you aware of that fact.

 */

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("EXAMPLE","Simple example with sin and cos solution and p=0");
}

// exact solution
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] =          cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[1] =      -Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[2] =       Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[3] =       Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[4] = -3*Pi*Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z); //Laplacien
}
void ExactU2(double x, double y,  double z, double *values)
{
  values[0] =          sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[1] =       Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[2] =      -Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[3] =       Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[4] = -3*Pi*Pi*sin(Pi*x)*cos(Pi*y)*sin(Pi*z); //Laplacien
}
void ExactU3(double x, double y,  double z, double *values)
{
  values[0] =      -2*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[1] =   -2*Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[2] =   -2*Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[3] =    2*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[4] = 6*Pi*Pi*sin(Pi*x)*sin(Pi*y)*cos(Pi*z); //Laplacien
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}
// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  double diri[5];
  ExactU1(x,y,z,diri);
  value = diri[0]; //Dirichlet value
}
void U2BoundValue(double x, double y, double z, double &value)
{
  double diri[5];
  ExactU2(x,y,z,diri);
  value = diri[0]; //Dirichlet value
}
void U3BoundValue(double x, double y, double z, double &value)
{
  double diri[5];
  ExactU3(x,y,z,diri);
  value = diri[0]; //Dirichlet value
}

void LinCoeffs(int n_points, double * X, double * Y, double * Z,
               double **parameters, double **coeffs)
{
  const double eps = DIMENSIONLESS_VISCOSITY;
  double u1[5], u2[5], u3[5], p[5];
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] =  eps;

    ExactU1(X[i], Y[i], Z[i], u1);
    ExactU2(X[i], Y[i], Z[i], u2);
    ExactU3(X[i], Y[i], Z[i], u3);
    ExactP( X[i], Y[i], Z[i], p);

    coeffs[i][1] = -eps * u1[4] + p[1]; //Stokes: diffusion and pressure gradient
    coeffs[i][2] = -eps * u2[4] + p[2];
    coeffs[i][3] = -eps * u3[4] + p[3];
    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5) //Navier--Stokes:
    {                                              // add convective terms
      coeffs[i][1] += u1[0]*u1[1] + u2[0]*u1[2] + u3[0]*u1[3];
      coeffs[i][2] += u1[0]*u2[1] + u2[0]*u2[2] + u3[0]*u2[3];
      coeffs[i][3] += u1[0]*u3[1] + u2[0]*u3[2] + u3[0]*u3[3];
    }
    coeffs[i][4] = 0; // g
  }
}

// From here it's old stuff: explicitely formulated dirichlet and neumann
// boundary conditions - if you want to reuse them, take care of the code!
// // kind of boundary condition (for FE space needed)
//void BoundCondition(double x, double y, double z, BoundCond &cond)
//{
//  double tol = 1e-10;
//  if((fabs(1+x) < tol) || (fabs(1+y) < tol) || (fabs(1+z) < tol)
//       || (fabs(1-z) < tol))
//    cond = DIRICHLET;
//  else
//    cond = NEUMANN;
//
//}
//
// // former values of boundary condition
//void U1BoundValue(double x, double y, double z, double &value)
//{
//  const double eps = KINEMATIC_VISCOSITY;
//  double tol = 1e-10;
//  if((fabs(1+x) < tol) || (fabs(1+y) < tol) || (fabs(1+z) < tol)
//       || (fabs(1-z) < tol))
//    value = cos(Pi*x)*sin(Pi*y)*sin(Pi*z); //Dirichlet
//  else{
//    if(fabs(1-x) < tol)
//    {
//      value = -eps*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += -eps*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
//    }
//    else
//    {
//      value = eps*Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += eps*Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
//    }
//  }
//}
//void U2BoundValue(double x, double y, double z, double &value)
//{
//  const double eps = KINEMATIC_VISCOSITY;
//  double tol = 1e-10;
//  if((fabs(1+x) < tol) || (fabs(1+y) < tol) || (fabs(1+z) < tol)
//       || (fabs(1-z) < tol))
//    value = sin(Pi*x)*cos(Pi*y)*sin(Pi*z); //Dirichlet
//  else{
//    if(fabs(1-x) < tol)
//    {
//      value = eps*Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += eps*Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
//    }
//    else
//    {
//      value = -eps*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += -eps*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
//    }
//  }
//}
//void U3BoundValue(double x, double y, double z, double &value)
//{
//  const double eps = KINEMATIC_VISCOSITY;
//  double tol = 1e-10;
//  if((fabs(1+x) < tol) || (fabs(1+y) < tol) || (fabs(1+z) < tol)
//       || (fabs(1-z) < tol))
//    value = -2*sin(Pi*x)*sin(Pi*y)*cos(Pi*z); //Dirichlet
//  else{
//    if(fabs(1-x) < tol)
//    {
//      value = -eps*2*Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += eps*Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
//    }
//    else
//      {
//      value = -eps*2*Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z); //Neumann
//      if(TDatabase::ParamDB->LAPLACETYPE == 1)
//        value += eps*Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
//      }
//  }
//}

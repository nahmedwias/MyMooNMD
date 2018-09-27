// Brinkman problem, solution in ansatz space
// velocity pw quadratic, pressure linear
// 
// u(x,y,z) = (x^2+y^2+z^2, x^2+2xy+13, -2xz+5y^2)^T
// p(x,y,z) = 3x-2y+7z-4

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
//double DIMENSIONLESS_VISCOSITY;

double viscosity = -1;
double effective_viscosity = -1;
double permeability = -1;

void ExampleFile()
{
  Output::info<1>("EXAMPLE","AnsatzQuadLin.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = x*x+y*y+z*z;
  values[1] = 2*x;
  values[2] = 2*y;
  values[3] = 2*z;
  values[4] = 6;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = x*x+2*x*z+13;
  values[1] = 2*x+2*z;
  values[2] = 0;
  values[3] = 2*x;
  values[4] = 2;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = -2*x*z+5*y*y;
  values[1] = -2*z;
  values[2] = 10*y;
  values[3] = -2*x;
  values[4] = 10;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 3*x-2*y+7*z-4;
  values[1] = 3;
  values[2] = -2;
  values[3] = 7;
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
  value = x*x+y*y+z*z;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = x*x+2*x*z+13;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  value = -2*x*z+5*y*y;
}
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y,  double z, double *values)
{
  values[0] = x*x+y*y+z*z;
  values[1] = 2*x;
  values[2] = 2*y;
  values[3] = 2*z;
  values[4] = 6;
}

void InitialU2(double x, double y,  double z, double *values)
{
  values[0] = x*x+2*x*z+13;
  values[1] = 2*x+2*z;
  values[2] = 0;
  values[3] = 2*x;
  values[4] = 2;
}

void InitialU3(double x, double y,  double z, double *values)
{
  values[0] = -2*x*z+5*y*y;
  values[1] = -2*z;
  values[2] = 10*y;
  values[3] = -2*x;
  values[4] = 10;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 3*x-2*y+7*z-4;
  values[1] = 3;
  values[2] = -2;
  values[3] = 7;
  values[4] = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
//  static double eps = DIMENSIONLESS_VISCOSITY;

  double *coeff, x, y, z;

    
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 3)// STOKES
  {
    for(int i=0;i<n_points;i++)
    {
        x=X[i];
        y=Y[i];
        z=Z[i];
      coeff = coeffs[i];
        
        //coeff[0] = eps;
        coeff[5] = viscosity;
        coeff[6] = effective_viscosity;
        coeff[7] = permeability;
        coeff[1] = -6*coeff[6] + 3 + (coeff[5]/coeff[7])*(x*x+y*y+z*z); // f1
        coeff[2] = -2*coeff[6] + -2 + (coeff[5]/coeff[7])*(x*x+2*x*z+13); // f2
        coeff[3] = -10*coeff[6] + 7 + (coeff[5]/coeff[7])*(-2*x*z+5*y*y); // f3
        coeff[4] = 0; // g
        coeff[8]=TDatabase::ParamDB->equal_order_stab_weight_PkPk;
        coeff[9]=TDatabase::ParamDB->equal_order_stab_weight_PkPk;
        
    }
  }
  else // Navier-Stokes
  {
    for(int i=0;i<n_points;i++)
    {
        x=X[i];
        y=Y[i];
        z=Z[i];
        coeff = coeffs[i];
        
        //coeff[0] = eps;
        coeff[5] = viscosity;
        coeff[6] = effective_viscosity;
        coeff[7] = permeability;
        coeff[1] = -6*coeff[6] + 3 + (coeff[5]/coeff[7])*(x*x+y*y+z*z); // f1
        coeff[2] = -2*coeff[6] + -2 + (coeff[5]/coeff[7])*(x*x+2*x*z+13); // f2
        coeff[3] = -10*coeff[6] + 7 + (coeff[5]/coeff[7])*(-2*x*z+5*y*y); // f3
        coeff[4] = 0; // g
        coeff[8]=TDatabase::ParamDB->equal_order_stab_weight_PkPk;
        coeff[9]=TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    }
  }
    
}

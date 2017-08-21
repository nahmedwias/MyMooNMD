// Brinkman problem, solution in ansatz space
// velocity pw linear, pressure constant
// 
// u(x,y,z) = (y+z,5x-3z,-x-2y)^T=(u1,u2,u3)^T
// p(x,y,z) = 0

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
//double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::info<1>("EXAMPLE","AnsatzLinConst.h");
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = y+z;      // u1
  values[1] = 0;        // u1_x
  values[2] = 1;        // u1_y
  values[3] = 1;        // u1_z
  values[4] = 0;        // Delta u1=u1_xx+u1_yy+u1_zz
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 5*x-3*z;      // u2
  values[1] = 5;            // u2_x
  values[2] = 0;            // u2_y
  values[3] = -3;           // u2_z
  values[4] = 0;            // Delta u2=u2_xx+u2_yy+u2_zz
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = -x-2*y;       // u3
  values[1] = -1;           // u3_x
  values[2] = -2;           // u3_y
  values[3] = 0;            // u3_z
  values[4] = 0;            // Delta u3=u3_xx+u3_yy+u3_zz
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}
// void ExactNull(double x, double y, double z, double *values)
// {
//   values[0] =0;
//   values[1] =0;
//   values[2] =0;
//   values[3] =0;
//   values[4] =0;
// }

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
    value = y+z;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
    value = 5*x-3*z;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
    value = -x-2*y ;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
//  static double eps = DIMENSIONLESS_VISCOSITY;
    double *coeff;
    
    for(int i=0;i<n_points;i++)
    {
        coeff = coeffs[i];
        
        //coeff[0] = eps;
        coeff[5]=0;//TDatabase::ParamDB->VISCOSITY;//0.;
        coeff[6]= TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
        coeff[7]=TDatabase::ParamDB->PERMEABILITY;
        coeff[1] = (coeff[5]/coeff[7])*(Y[i]+Z[i]); // f1
        coeff[2] = (coeff[5]/coeff[7])*(5*X[i]-3*Z[i]); // f2
        coeff[3] = (coeff[5]/coeff[7])*(-1*X[i]-2*Y[i]); // f3
        coeff[4] = 0; // g
        coeff[8]=TDatabase::ParamDB->equal_order_stab_weight_PkPk;
        coeff[9]=TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    }
 
}

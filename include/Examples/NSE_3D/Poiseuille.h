// Brinkman problem, solution in ansatz space
// Poiseuille (exact solution in P2/P1)
// 
// u(x,y,z) = ( 3*(1-r^2), 0 , 0)^T=(u1,u2,u3)^T
// r^2:= y^2+z^2
// p(x,y,z) = -x+1/2

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("EXAMPLE","Poiseuille.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
    values[0] = 0;           // u1
    values[1] = 0;           // u1_x
    values[2] = 0;           // u1_y
    values[3] = 0;           // u1_z
    values[4] = 0;           // Delta u1
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;            // u2
  values[1] = 0;            // u2_x
  values[2] = 0;            // u2_y
  values[3] = 0;            // u2_z
  values[4] = 0;            // Delta u2=u2_xx+u2_yy+u2_zz
}

void ExactU3(double x, double y,  double z, double *values)
{
    values[0] = 1-(x*x+y*y);//3*(1-(y*y+x*x));      // u3
    values[1] = -2*x; //-6*x*x;        // u3_x
    values[2] = -2*y;        // u3_y
    values[3] = 0;        // u3_z
    values[4] = -4; //-12;        // Delta u3=u3_xx+u3_yy+u3_zz

}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = z-8;
  values[1] = 0;
  values[2] = 0;
  values[3] = 1;
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
    value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
    value = 0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
    value = 1-(x*x+y*y);//3*(1-(y*y+x*x))
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
    double *coeff;
    
    for(int i=0;i<n_points;i++)
    {
        coeff = coeffs[i];
        
        coeff[0] = 1;//eps;
//        coeff[5]=0;//TDatabase::ParamDB->VISCOSITY;//0.;
//        coeff[6]= TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
//        coeff[7]=1e10;//TDatabase::ParamDB->PERMEABILITY;
        coeff[1] = 0; // coeff[6]*(-12)+(-1)+(coeff[5]/coeff[7])*(3*(1-(Y[i]*Y[i]+Z[i]*Z[i]))); // f1
        coeff[2] = 0; // f2
        coeff[3] = 0.;//-coeff[6]*(-4)+(1)+(coeff[5]/coeff[7])*(1-(Y[i]*Y[i]+X[i]*X[i])); // 0; // f3
        coeff[4] = 0; // g
//        coeff[8]=TDatabase::ParamDB->equal_order_stab_weight_P1P1;
//        coeff[9]=TDatabase::ParamDB->equal_order_stab_weight_P2P2;
    }
 
}

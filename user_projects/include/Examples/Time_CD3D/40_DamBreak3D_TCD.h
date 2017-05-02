// ==========================================================================
// instationary problem
// ==========================================================================
void ExampleFile()
{
  Output::print<1>("Example: 40_DamBreak3D_TCD.h\n");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = NEUMANN;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
  if ( (x <= 0.4) && (z <= 0.7) )
    values[0] = 1; // liquid
  else
    values[0] = 0; // gas
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double *coeff;
  double b1,b2,b3;
  
  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    double x = X[i];
//    double y = Y[i];
//    double z = Z[i];

    b1 = parameters[i][0]; // ux
    b2 = parameters[i][1]; // uy
    b3 = parameters[i][2]; // uz

    coeff[0] = 0; //diffusion coefficient
    coeff[1] = b1;   //ux
    coeff[2] = b2;   //uy
    coeff[3] = b3;   //uy
    coeff[4] = 0;   //reaction coefficient

    coeff[5] = 0;
  }
}


// ==========================================================================
// instationary problem
// ==========================================================================

// coordinates of the plugs p1 and p2
//  double p1x = -0.566,  p1y = +0.5;
//  double p2x = -0.468,  p2y = -0.27;
double p1x = 0,  p1y = 0;
double p2x = 0,  p2y = 0;
double p_radius = 0.05;
double height = 0.4;


void ExampleFile()
{
  Output::print<1>("Example: 50_GasStirring_TCD.h\n");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  if (z == 0)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  double in_plug1 = p_radius*p_radius-(x-p1x)*(x-p1x)-(y-p1y)*(y-p1y);
//  double in_plug2 = p_radius*p_radius-(x-p2x)*(x-p2x)-(y-p2y)*(y-p2y);

  if (z == 0)
  {
    if (in_plug1 >= 0)// || in_plug2 >= 0)
      value = 0; // gas injection = DIRICHLET
    else
      value = 1; // Dirichlet, liquid
  }
  else
    value = 0;
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
//  double in_plug1 = p_radius*p_radius-(x-p1x)*(x-p1x)-(y-p1y)*(y-p1y);
//  double in_plug2 = p_radius*p_radius-(x-p2x)*(x-p2x)-(y-p2y)*(y-p2y);

  if ( z <= height )// && z > 0)
  {
    values[0] = 1; // liquid
  }
//  else if (z == 0)
//  {
//    if (in_plug1 >= 0)// || in_plug2 >= 0)
//      values[0] = 0; // gas injection = DIRICHLET
//    else
//      values[0] = 1; // Dirichlet, liquid
//  }
  else
  {
    values[0] = 0; // gas
  }

//  if ( (z == 0) && (in_plug1 >0 || in_plug2 > 0) )
//    values[0] = 0; // gas injection = DIRICHLET

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
    coeff[3] = b3;   //uz
    coeff[4] = 0;   //reaction coefficient

    coeff[5] = 0;
  }
}


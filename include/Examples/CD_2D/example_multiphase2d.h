// Example file for Multiphase 2D

void ExampleFile()
{
  Output::info<3>("Example:", "example_multiphase2d.h");
}

// exact solution
void Exact(double x, double y, double *values)
{
//  const double p = Pi; // 2*Pi;
  values[0] = x+y; //sin(p*x)*sin(p*y); // exact solution
  values[1] = 1; //p*cos(p*x)*sin(p*y); // x-derivative
  values[2] = 1; //p*sin(p*x)*cos(p*y); // y-derivative
  values[3] = 0; //-2*p*p*sin(p*x)*sin(p*y); // laplacian
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  if(BdComp==1)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  if(BdComp==1) // Neumann, x = 1, y = Param
  {
    const double eps = 1. / TDatabase::ParamDB->PE_NR;
    double exact[4];
    Exact(1., Param, exact);
    value = eps * exact[1];
    //value = -eps*Pi*sin(Pi*Param);
  }
  else // Dirichlet
    value = 0;
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
//  const double eps=1/TDatabase::ParamDB->PE_NR;
  double exact[4];
  Output::print<1>("START BILINEAR COEFFICIENTS------------");
  for(int i = 0; i < n_points; i++)
  {
    Output::print<1>("THIS is point number ", i);

    coeffs[i][0] = 0; // diffusion coefficient D
    coeffs[i][1] = 1; // convection coefficient v_x
    coeffs[i][2] = 1; // convection coefficient v_y
    coeffs[i][3] = 0; // reaction coefficient R
    // coeffs[i][4] is the right hand side
    
    Exact(x[i], y[i], exact);
    coeffs[i][4] = -coeffs[i][0]*exact[3]; // diffusion
    coeffs[i][4] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2]; // convection
    coeffs[i][4] += coeffs[i][3]*exact[0]; // reaction
  }
  Output::print<1>("----------------END BILINEAR COEFFICIENTS");
}


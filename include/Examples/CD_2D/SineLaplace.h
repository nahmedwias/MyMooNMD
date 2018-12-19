// ======================================================================
// Sine problem
// ======================================================================


void ExampleFile()
{
  Output::print<1>("Example: SineLaplace.h");
}

// exact solution
void Exact(double x, double y, double *values)
{
  const double p = Pi; // 2*Pi;
  values[0] = sin(p*x)*sin(p*y);
  values[1] = p*cos(p*x)*sin(p*y);
  values[2] = p*sin(p*x)*cos(p*y);
  values[3] = -2*p*p*sin(p*x)*sin(p*y);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double, BoundCond &cond)
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

void BilinearCoeffs(int n_points, double *x, double *y, double **,
                    double **coeffs)
{
  const double eps=1/TDatabase::ParamDB->PE_NR;
  double exact[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;
    coeffs[i][1] = 0;// parameters[i][0]; 
    coeffs[i][2] = 0;//parameters[i][1];
    coeffs[i][3] = 0;
    
    Exact(x[i], y[i], exact);

    coeffs[i][4] = -coeffs[i][0]*exact[3]; // diffusion
    coeffs[i][4] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2]; // convection
    coeffs[i][4] += coeffs[i][3]*exact[0]; // reaction

  // coeffs[i][4] = 0;
  }
}


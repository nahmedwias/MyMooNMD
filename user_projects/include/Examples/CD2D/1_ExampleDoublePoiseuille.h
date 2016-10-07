// Example Convection equation which goes with NSE example 11!!!
// Solves stratified Poiseuille flow, with 2 different fluids.
// The property of the 2 fluids are set in the boundary condition
// in THIS FILE.
double USER1;
double USER2;


void ExampleFile()
{
  Output::info<3>("Example:", "ExampleDoublePoiseuille.h");
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0; // exact solution
  values[1] = 0; // x-derivative
  values[2] = 0; // y-derivative
  values[3] = 0; // laplacian
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  if(BdComp == 1)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  switch (BdComp)
  {
    case 0: value = 0; break;
    case 1: value = 0; break;
    case 2: value = 1; break;
    case 3:
      if (Param <= 0.5)
        value = 1;
      else
        value = 0;
      break;
    default: cout << "wrong boundary part number" << endl;
    break;
  }
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  const double eps=1/TDatabase::ParamDB->PE_NR;
  double exact[4];

  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;              // diffusion coefficient D
    coeffs[i][1] = parameters[i][0]; // convection coefficient v_x
    coeffs[i][2] = parameters[i][1]; // convection coefficient v_y

//    if (fabs(x[i] - 1 ) <= 1e-2 && fabs(y[i] - 0.5) <= 1e-2)
//    {
//      Output::print<1>("x[i]=", x[i], " y[i]=", y[i], " x-velocity=", parameters[i][0], " y-velocity=", parameters[i][1]);
//    }
//    if (fabs(x[i] - 0 ) <= 1e-2 && fabs(y[i] - 0.7) <= 1e-2)
//    {
//      Output::print<1>("x[i]=", x[i], " y[i]=", y[i], " x-velocity=", parameters[i][0], " y-velocity=", parameters[i][1]);
//    }
//    cout << coeffs[i][1] << " ";

    coeffs[i][3] = 0; // reaction coefficient R
    // coeffs[i][4] is the right hand side
    
    Exact(x[i], y[i], exact);
    coeffs[i][4] = -coeffs[i][0]*exact[3]; // diffusion
    coeffs[i][4] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2]; // convection
    coeffs[i][4] += coeffs[i][3]*exact[0]; // reaction
  }
//  Output::print<1>("----------------END BILINEAR COEFFICIENTS");
}


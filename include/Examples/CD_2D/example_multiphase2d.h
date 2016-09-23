// Example file for Multiphase 2D
double USER1;
double USER2;
double USER3;

void ExampleFile()
{
  Output::info<3>("Example:", "example_multiphase2d.h");
}

// exact solution
void Exact(double x, double y, double *values)
{
//  const double p = Pi; // 2*Pi;
  values[0] = 0; //sin(p*x)*sin(p*y); // exact solution
  values[1] = 1; //p*cos(p*x)*sin(p*y); // x-derivative
  values[2] = 1; //p*sin(p*x)*cos(p*y); // y-derivative
  values[3] = 0; //-2*p*p*sin(p*x)*sin(p*y); // laplacian
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  if(BdComp==0 || BdComp ==1 || BdComp == 3)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  if(BdComp==0)
  {
    if ((Param>1.0/3.0)&& (Param<2.0/3.0))
      value = 1;
    else
      value = 0;
  }
  else
    value = 0;
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
//  const double eps=1/TDatabase::ParamDB->PE_NR;
  double exact[4];
//  Output::print<1>("START BILINEAR COEFFICIENTS------------");
  for(int i = 0; i < n_points; i++)
  {
//    Output::print<1>("THIS is point number ", i, " with x and y coordinates "
//                     , x[i]," ", y[i]);

    coeffs[i][0] = 0; // diffusion coefficient D
    coeffs[i][1] = parameters[i][0]; // convection coefficient v_x

//    if (fabs(x[i] - 1 ) <= 1e-2 && fabs(y[i] - 0.5) <= 1e-2)
//    {
//      Output::print<1>("x[i]=", x[i], " y[i]=", y[i], " x-velocity=", parameters[i][0], " y-velocity=", parameters[i][1]);
//    }
//    if (fabs(x[i] - 0 ) <= 1e-2 && fabs(y[i] - 0.7) <= 1e-2)
//    {
//      Output::print<1>("x[i]=", x[i], " y[i]=", y[i], " x-velocity=", parameters[i][0], " y-velocity=", parameters[i][1]);
//    }
//    cout << coeffs[i][1] << " ";
    coeffs[i][2] = parameters[i][1]; // convection coefficient v_y
//    cout << coeffs[i][2] << " ";
    coeffs[i][3] = 0; // reaction coefficient R
    // coeffs[i][4] is the right hand side
    
    Exact(x[i], y[i], exact);
    coeffs[i][4] = -coeffs[i][0]*exact[3]; // diffusion
    coeffs[i][4] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2]; // convection
    coeffs[i][4] += coeffs[i][3]*exact[0]; // reaction
  }
//  Output::print<1>("----------------END BILINEAR COEFFICIENTS");
}


/************************************************
 *  Test for the Darcy limit of the Brinkman problem (mueff = 0)
 *  with sine-cosine solution taken from: DOI: https://doi.org/10.1137/08072632X.
 *
 *  exact solution:
 *  u(x,y) = [u1,u2]^T = [-2 pi cos(2 pi x) * sin(2 pi y),  -2 pi sin(2 pi x) * cos(2 pi y) ]^T
 *  p(x,y) = sigma * sin(2 Pi*x) * sin(2 Pi*y)
 *  u \cdot n = 0 on \partial \Omega
 *
 *  solves the Darcy Problem:
 *  sigma u + p = 0
 *  div u = 0
 *  u \cdot n = ... \partial Omega
 *
 *  The boundary condition is treated as an essential boundary condition.
 **************************************************/

// initialize physical parameters
double effective_viscosity = -1;
double sigma = -1;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_SinCos_DarcyFlow.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = -2 * Pi * cos(2 * Pi * x) * sin(2 * Pi * y);              //u1
  values[1] =  4 * Pi * Pi * sin(2 * Pi * x) * sin(2 * Pi * y);         //u1_x
  values[2] = -4 * Pi * Pi * cos(2 * Pi * x) * cos(2 * Pi * y);         //u1_y
  values[3] = 16 * Pi * Pi * Pi * cos(2 * Pi * x) * sin(2 * Pi * y);    //Delta u1=u1_xx+u1_yy
}

void ExactU2(double x, double y, double *values)
{
  values[0] = -2 * Pi * sin(2 * Pi * x) * cos(2 * Pi * y);            //u2
  values[1] = -4 * Pi * Pi * cos(2 * Pi * x) * cos(2 * Pi * y);       //u2_x
  values[2] =  4 * Pi * Pi * sin(2 * Pi * x) * sin(2 * Pi * y);       //u2_y
  values[3] = 16 * Pi * Pi * Pi * sin(2 * Pi * x) * cos(2 * Pi * y);  //Delta u2=u2_xx + u2_yy
}

void ExactP(double x, double y, double *values)
{
  values[0] = sigma * sin(2 * Pi * x) * sin(2 * Pi * y);                    //p
  values[1] = sigma * 2 * Pi * cos(2 * Pi * x) * sin(2 * Pi * y);           //p_x
  values[2] = sigma * 2 * Pi * sin(2 * Pi * x) * cos(2 * Pi * y);           //p_y
  values[3] = sigma * (-8) * Pi * Pi * sin(2 * Pi * x) * sin(2 * Pi * y);   //Delta p=p_xx+p_yy
}

// ========================================================================
// boundary conditions (Parametrisation of the boundary); Param \in [0,1]
// ========================================================================
void BoundCondition(int i, double Param, BoundCond &cond)
{
  cond = DIRICHLET; // default

  // set Neumann BC
  for (int j = 0; j < TDatabase::ParamDB->n_neumann_boundary; j++)
  {
    if (i == TDatabase::ParamDB->neumann_boundary_id[j])
    {
      cond = NEUMANN;
      return;
    }
  }
  // set Nitsche BC
  for (int j = 0; j < TDatabase::ParamDB->n_nitsche_boundary; j++)
  {
    if (i == TDatabase::ParamDB->nitsche_boundary_id[j])
    {
      cond = DIRICHLET_WEAK;
      return;
    }
  }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
  case 0: value = 0.;
  break;
  case 1: value = -2 * Pi * sin( 2 * Pi * Param );  // u \cdot n = u_1; x = 1 --> cos(2 pi x) = 1
  break;
  case 2: value = 0.;
  break;
  case 3: value = -2 * Pi * sin( 2 * Pi * (1-Param) );  // u \cdot n = - u_1; x = 0 --> cos(2 pi x) = 1
  break;
  default: cout << "No boundary component with this number." << endl;
  break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
  case 0: value = -2 * Pi * sin( 2 * Pi * Param);  // u \cdot n = - u_2; y = 0 --> cos(2 pi y) = 1
  break;
  case 1: value = 0.;
  break;
  case 2: value = -2 * Pi * sin( 2 * Pi * (1-Param));  // u \cdot n = u_2; y = 1 --> cos(2 pi y) = 1
  break;
  case 3: value = 0.;
  break;
  default: cout << "No boundary component with this number." << endl;
  break;
  }
}

// ========================================================================
// coefficients for Brinkman problem:
// mu, f1, f2, g, sigma = mu/permeability
// with:
// -mueff Delta u + grad(p) + sigma u = (f1,f2)
// div(u) = g
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
    double **parameters, double **coeffs)
{
  double val_u1[4];
  double val_u2[4];
  //double val_p[4];

  for(int i = 0; i < n_points; i++)
  {
    ExactU1(x[i], y[i], val_u1);
    ExactU2(x[i], y[i], val_u2);
    //ExactP(x[i], y[i], val_p);

    // physical parameters
    coeffs[i][0] = effective_viscosity;
    coeffs[i][4] = sigma;

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0;   //-coeffs[i][5]*val_u1[3] + val_p[1] + (coeffs[i][4]/coeffs[i][6])*val_u1[0];   // f1
    coeffs[i][2] = 0;   //-coeffs[i][5]*val_u2[3] + val_p[2] + (coeffs[i][4]/coeffs[i][6])*val_u2[0];   // f2

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = val_u1[1] + val_u2[2]; // g

    //coeffs[i][7] = TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    //coeffs[i][8] = TDatabase::ParamDB->grad_div_stab_weight;
  }
}

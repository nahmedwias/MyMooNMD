// Test for the Darcy limit of the Brinkman problem (mueff=0) with sine and cosine functions
/* exact solution:
   u(x,y) = (-2 pi cos(2 pi x) * sin(2 pi y),  -2 pi sin(2 pi x) * cos(2 pi y) )^T = (u1,u2)^T
   p(x,y) = sigma * in(2 Pi*x)*sin(2 Pi*y)
   u \cdot n = 0 on \partial \Omega

   solves the Darcy Problem: 
   sigma u + p = 0
   div u = 0
   u \cdot n = ... \partial Omega

The boundary condition is treated as an essential boundary condition.

   Taken from: 
   @Article{BadiaCodina2009,
   Title       = {{U}NIFIED STABILIZED FINITE ELEMENT FORMULATIONS FOR THE {S}TOKES AND THE {D}ARCY PROBLEMS},
   Author      = {Badia, S. and Codina, R.},
   Journal     = {SIAM Journal on Numerical Analysis},
   Year        = {2009},
   Number      = {3},
   Pages       = {1971-2000},
   Volume      = {47},
   Publisher   = {Society for Industrial and Applied Mathematics}
   }
 */

 double viscosity = -1;
 double effective_viscosity = -1;
 double permeability = -1;

void ExampleFile()
{
  Output::print<1>("Example: SinCos_BadiaCodina_DarcyFlow_Test.h");
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
  //double sigma = TDatabase::ParamDB->VISCOSITY / TDatabase::ParamDB->PERMEABILITY;
  double sigma= viscosity/permeability;

  values[0] = sigma * sin(2 * Pi * x) * sin(2 * Pi * y);                    //p
  values[1] = sigma * 2 * Pi * cos(2 * Pi * x) * sin(2 * Pi * y);           //p_x
  values[2] = sigma * 2 * Pi * sin(2 * Pi * x) * cos(2 * Pi * y);           //p_y
  values[3] = sigma * (-8) * Pi * Pi * sin(2 * Pi * x) * sin(2 * Pi * y);   //Delta p=p_xx+p_yy
}

// ========================================================================
// boundary conditions (Parametrisierung des Randes); Param \in [0,1]
// ========================================================================

void BoundCondition(int i, double Param, BoundCond &cond)
{
  cond = DIRICHLET; // default

  for (int j = 0; j < TDatabase::ParamDB->n_neumann_boundary; j++)
  {
    if (i == TDatabase::ParamDB->neumann_boundary_id[j])
    {
      cond = NEUMANN;
      return;
    }
  }
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
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double val_u1[4];
  double val_u2[4];
  double val_p[4];

  for(int i = 0; i < n_points; i++)
  {
    //coeffs[i][4]= TDatabase::ParamDB->VISCOSITY;
    //coeffs[i][5]= TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
    //coeffs[i][6]= TDatabase::ParamDB->PERMEABILITY;

    coeffs[i][4] = viscosity;
    coeffs[i][5] = effective_viscosity;
    coeffs[i][6] = permeability;


    ExactU1(X[i], Y[i], val_u1);
    ExactU2(X[i], Y[i], val_u2);
    ExactP(X[i], Y[i], val_p);

    coeffs[i][1] = 0; //-coeffs[i][5]*val_u1[3] + val_p[1] + (coeffs[i][4]/coeffs[i][6])*val_u1[0];   // f1
    coeffs[i][2] = 0; //-coeffs[i][5]*val_u2[3] + val_p[2] + (coeffs[i][4]/coeffs[i][6])*val_u2[0];   // f2
    coeffs[i][3] = val_u1[1] + val_u2[2];                                                       // g (divergence)

    coeffs[i][0] = (coeffs[i][5] / coeffs[i][4]) * coeffs[i][6];
    coeffs[i][7] = TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    coeffs[i][8] = TDatabase::ParamDB->grad_div_stab_weight;
    }
}

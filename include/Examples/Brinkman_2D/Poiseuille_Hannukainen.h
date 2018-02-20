// Brinkman 2D problem, 
/* Poiseuille-Problem from the Paper: 
   @Article{Hannukainen2011,
   author="Hannukainen, Antti
   and Juntunen, Mika
   and Stenberg, Rolf",
   title="Computations with finite element methods for the Brinkman problem",
   journal="Computational Geosciences",
   year="2011",
   month="Jan",
   day="01",
   volume="15",
   number="1",
   pages="155--166",
   abstract="Various finite element families for the Brinkman flow (or Stokes--Darcy flow) are tested numerically. Particularly, the effect of small permeability is studied. The tested finite elements are the MINI element, the Taylor--Hood element, and the stabilized equal order methods. The numerical tests include both a priori analysis and adaptive methods.",
   issn="1573-1499",
   doi="10.1007/s10596-010-9204-4",
   url="https://doi.org/10.1007/s10596-010-9204-4"
   }
 */

// authors: Alfonso Caiazzo and Laura Blank

/*
   Domain =  [0,1]x[0,1]
   t^2 = K * (mu_eff/mu)
   u(x,y) = [ K/mu * (1+exp(1/t)-exp((1-y)/t) - exp(y/t)) / (1+exp(1/t)  ,  0 ],   for t>0
   u(x,y) = [ K/mu , 0 ],    for t=0
   p(x,y) = -x + 1/2
 */


void ExampleFile()
{
  Output::print<1>("Example: Poiseuille_Hannukainen.h");
}

// ========================================================================
// exact solution
// ========================================================================

void ExactU1(double x, double y, double *values)
{
  double K = TDatabase::ParamDB->PERMEABILITY;
  double mu = TDatabase::ParamDB->VISCOSITY;
  double mu_eff = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
  double t = fabs(sqrt((mu_eff/mu) * K));

  if (K == 0)
  {
    ErrThrow("The permeability K is zero. Division by zero not allowed. ");
  }
  else if (mu == 0)
  {
    ErrThrow("The fluid viscosity mu is zero. Division by zero not allowed. ");
  }

  if (t == 0)
  {
    values[0] = K/mu;                                                               //u1
    values[1] = 0;                                                                  //u1_x
    values[2] = 0;                                                                  //u1_y
    values[3] = 0;                                                                  //Delta u1
  }
  else
  {
    values[0] = K/mu * (1+exp(1/t)-exp((1-y)/t) - exp(y/t)) / (1+exp(1/t));             //u1
    values[1] = 0;                                                                      //u1_x
    values[2] = K/mu * (exp((1-y)/t)-exp(y/t))/(t * (1+exp(1/t)));                      //u1_y
    values[3] = K/mu * (-exp((1-y)/t)-exp(y/t))/(t * t * (1+exp(1/t)));                 //Delta u1
  }
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;            //u2
  values[1] = 0;            //u2_x
  values[2] = 0;            //u2_y
  values[3] = 0;            //Delta u2
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0.5-x;                      //p
  values[1] = -1;                         //p_x
  values[2] = 0;                          //p_y
  values[3] = 0;                          //Delta p=p_xx+p_yy
}

// ========================================================================
// boundary conditions
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
      // Todo
      cond = DIRICHLET_WEAK;
      return;
    }
  }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double K = TDatabase::ParamDB->PERMEABILITY;
  double mu = TDatabase::ParamDB->VISCOSITY;
  double mu_eff = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
  double t = fabs(sqrt((mu_eff/mu) * K));

  // loop to impose Neumann boundary conditions
  // Since we are using the Neumann boundary condition via boundaryAssembling, the boundvalue here has to be alway zero!!!!
  for (int j = 0; j < TDatabase::ParamDB->n_neumann_boundary; j++)
  {
    if ( BdComp == TDatabase::ParamDB->neumann_boundary_id[j])
    {
      switch(BdComp)
      {
        case 1:
          value = 0.;//TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        case 3:
          value = 0.;//-TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        default:
          Output::print("I cannot impose Neumann boundary condition on component ", BdComp);
          exit(1);
          break;
      }
      return;
    }
  }

  // loop to impose (strong or weak) Dirichlet
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = K/mu * (1+exp(1/t)-exp((1-Param)/t) - exp(Param/t)) / (1+exp(1/t));
            break;
    case 2: value= 0;
            break;
    case 3: value = K/mu * (1+exp(1/t)-exp((Param)/t) - exp((1-Param)/t)) / (1+exp(1/t));
            break;
    default: cout << "No boundary component with this number." << endl;
             break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}


// ========================================================================
// coefficients for Brinkman problem: viscosity, effective viscosity, permeability, f1, f2, g
// (the lhs of the Brinkman problem computed at quadrature points - for the error norms)
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
    double **parameters, double **coeffs)
{
  double val_u1[4];
  double val_u2[4];
  double val_p[4];

  double *coeff;

  for(int i = 0; i < n_points; i++)
  {
    ExactU1(x[i], y[i], val_u1);
    ExactU2(x[i], y[i], val_u2);
    ExactP(x[i], y[i], val_p);

    coeff = coeffs[i];

    coeff[4] = TDatabase::ParamDB->VISCOSITY;
    coeff[5] = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
    coeff[6] = TDatabase::ParamDB->PERMEABILITY;
    coeff[0] = (coeff[5] / coeff[4]) * coeff[6];
    coeff[1] = 0;//-coeff[5] * val_u1[3] - val_p[1] + (coeff[4]/coeff[6]) * val_u1[0];  //(coeff[4]/coeff[6])-1;   // f1 (rhs of Brinkman problem for u1)
    coeff[2] = 0;                                                                        // f2 (rhs of Brinkman problem for u2)
    coeff[3] = val_u1[1] + val_u2[2]; //0;                                               // g (divergence term=u1_x+u2_y)
    coeff[7] = TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    coeff[8] = TDatabase::ParamDB->grad_div_stab_weight;
  }
}



// Brinkman Problem, Poiseuille-Problem

//author: Alfonso Caiazzo and Laura Blank

/*
   Square [0,1]x[0,1]
   effective viscosity = 1
   u(x,y) = (4y(1-y),0),
   p(x,y) = 0.5-x
 */
double viscosity = -1;                                             
double effective_viscosity = -1;                                               
double permeability = -1;

void ExampleFile()
{
  Output::print<1>("Example: Poiseuille.h");
}

// ========================================================================
// exact solution
// ========================================================================

void ExactU1(double, double y, double *values)
{
  values[0] = 4 * y * (1-y);    //u1
  values[1] = 0;                //u1_x
  values[2] = 4 - 8 * y;        //u1_y
  values[3] = -8;               //Delta u1
}

void ExactU2(double, double, double *values)
{
  values[0] = 0;            //u2
  values[1] = 0;            //u2_x
  values[2] = 0;            //u2_y
  values[3] = 0;            //Delta u2
}

void ExactP(double x, double, double *values)
{
  values[0] = 0.5 - x;        //p
  values[1] = -1;   //1;      //p_x
  values[2] = 0;              //p_y
  values[3] = 0;              //Delta p=p_xx+p_yy
}


//=========================================================================
// boundary conditions
//=========================================================================

void BoundCondition(int i, double, BoundCond &cond)
{
  cond = DIRICHLET; // Default

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
    case 1: value = 4 * Param * (1-Param);
            break;
    case 2: value = 0;
            break;
    case 3: value = 4 * (1-Param) * (Param);
            break;
    default: cout << "wrong boundary part number" << endl;
             break;
  }
}


void U2BoundValue(int, double, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// (the lhs of the Brinkman problem computed at quadrature points - for the error norms)
// ========================================================================

void LinCoeffs(int n_points, double *, double *y,
    //std::vector<std::vector<double>> parameters,
    double **,
    //std::vector<std::vector<double>> coeffs,
    double **coeffs)
{
  //std::vector<double> coeff;
  double *coeff;

  for(int i = 0; i < n_points; i++)
  {
    coeff = coeffs[i];

    coeff[4] = viscosity;
    coeff[5] = effective_viscosity;
    coeff[6] = permeability;
    coeff[0] = 0.; // (coeff[5] / coeff[4]) * coeff[6];
    // if f1=coeff[1] is set to zero, we can see the effect of a change in the parameter t^2 in the solution,
    // but convergence rates are no longer meaningful since we do not know the analytical solution

    // if f1=coeff[1] is set to its analytical rhs (8 * coeff[0] - 1 + 4 * y[i] * (1-y[i]);),
    // a cange in the parameter t^2 will not influence the solution
    coeff[1] = 8 * coeff[5] - 1 + (coeff[4]/coeff[6] ) * 4 * y[i] * (1-y[i]);  //0;                    // f1 (rhs of Brinkman problem for u1)
    coeff[2] = 0;                                                           // f2 (rhs of Brinkman problem for u2)
    coeff[3] = 0;                                                           // g (divergence term=u1_x+u2_y)
    coeff[7] = TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    coeff[8] = TDatabase::ParamDB->grad_div_stab_weight;
  }
}



// Brinkman 2D problem,
/* Poiseuille-Problem from the thesis of Discacciati
 */

// authors: Alfonso Caiazzo and Laura Blank

/*
   Domain =  [0,12]x[0,8]
   unknown solution
 */


// physical parameter
// These should be reset when constructing the Example class

double viscosity = -1;
double effective_viscosity = -1;
double permeability = -1;

void ExampleFile()
{
  Output::print<1>("Example: Discacciati_Flow.h");
}

// ========================================================================
// exact solution
// ========================================================================

void ExactU1(double x, double y, double *values)
{
    values[0] = 0;             //u1
    values[1] = 0;             //u1_x
    values[2] = 0;             //u1_y
    values[3] = 0;             //Delta u1
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
  values[0] = 0;                      //p
  values[1] = 0;                      //p_x
  values[2] = 0;                      //p_y
  values[3] = 0;                      //Delta p=p_xx+p_yy
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
//  double K = permeability;
//  double mu = viscosity;
//  double mu_eff = effective_viscosity;


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
        default:
          Output::print("I cannot impose Neumann boundary condition on component ", BdComp);
          exit(1);
          break;
      }
      return;
    }
  }
  double UMAX;
  // loop to impose (strong or weak) Dirichlet
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 0;
               break;
    case 2: value = 0;
            break;
    case 3:
    	      UMAX = 0.1; // m/s
    	     value = 4. * UMAX * (1.-Param) * (Param);
            break;
    default:
    	Output::print("No boundary component with this number.");
    	exit(1);
             break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}


// ========================================================================
// coefficients for Brinkman problem: viscosity, effective viscosity, permeability, f1, f2, g
// (lhs and rhs of the Brinkman problem computed at quadrature points - for the error norms)
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
    double **parameters, double **coeffs)
{

  double val_u1[4];
  double val_u2[4];
  double val_p[4];

  for(int i = 0; i < n_points; i++)
  {

    // physical parameters
    coeffs[i][4] = 1; //viscosity; //TDatabase::ParamDB->VISCOSITY;
    coeffs[i][5] = effective_viscosity; //TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
    coeffs[i][6] = permeability; //TDatabase::ParamDB->PERMEABILITY;

    /*
       if (x[i] < 0.5 && y[i] < 0.5)
       {
       coeffs[i][6] = 0.001;
       }
     */

    /* ****************************************** */
    // Use an analytic_coefficient_function or input sol+mesh for the permeability field
    if (permeability == -1)
    {
      coeffs[i][6] = parameters[i][0];
    }


      // effective_viscosity
  	if ( ( ( (x[i] >=7  ) && (x[i] <= 8  ) )  && ( (y[i] >= 4) && (y[i] <= 6) ) ) )
  	  {
  		coeffs[i][5]= 0;
  	//	coeffs[i][4]= 0.001;
  	  }
  	else if ( ( ( (x[i] >= 9) && (x[i] <= 10) )  && ( (y[i] >= 2) && (y[i] <= 4.2) ) ) )
  	  {
  		coeffs[i][5] = 0;
  	//	coeffs[i][4]= 0.001;
  	  }
    else if ( ( (x[i] >=5  ) && ( (y[i] >= 2) && (y[i] <= 6) ) ) )
    {
    	coeffs[i][5] = 0;
    //	coeffs[i][4]= 0.001;
    }
    else
    {
    	coeffs[i][5] = 0.01;//0.00001;
    //	coeffs[i][4]= 0;
    }


    /* ***************************************** */
    // if sources and sinks via analytic_coefficient_function (delta distr.), then \neq 0
    coeffs[i][9]=0;
    /*
       double use_source_term = 1;
       if (use_source_term == 1)
       {
       coeffs[i][9]= parameters[i][0];
       }
     */


    // Adimensional parameter t^2 = mue/mu*K
    coeffs[i][0] = (coeffs[i][5] / coeffs[i][4]) * coeffs[i][6];


    // (f1,f2)(x,y): RHS for momentum equation
    ExactU1(x[i], y[i], val_u1);
    ExactU2(x[i], y[i], val_u2);
    ExactP(x[i], y[i], val_p);

    coeffs[i][1] = 0;
    //-coeff[5] * val_u1[3] - val_p[1] + (coeff[4]/coeff[6]) * val_u1[0];
    //(coeff[4]/coeff[6])-1;   // f1 (rhs of Brinkman problem for u1)

    coeffs[i][2] = 0;

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = val_u1[1] + val_u2[2]; //0;

    // stabilization parameters
    coeffs[i][7] = TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    coeffs[i][8] = TDatabase::ParamDB->grad_div_stab_weight;
  }
}



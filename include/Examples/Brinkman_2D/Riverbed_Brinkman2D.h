/**
 * Benchmark for Brinkman2d
 * domain Omega = (0, 2) x (0, 2)
 * subdomains pure flow = (0, 2) x (1.5+-eps, 2)
 *            porous = (0, 2) x (0, 1.5+-eps)
 * Interface Gamma = {0, 2} x {1.5+-eps}
 * So the interface starts at (0, 1.5) following some function through 
 * the domain and ending at (2, 1.5)

 Boundary conditions:
      Dirichlet on velocity on (0,2)x{0,2} 
      periodic & Neumann on {0,2}x(0,2)
*/

// physical parameter
// These should be reset when constructing the Example class

double viscosity = -1;
double effective_viscosity = -1;
double permeability = -1;

void ExampleFile()
{
  OutPut("Example: Riverbed_Brinkman2D.h" << endl);
}

// ========================================================================
// exact solution unknown
// ========================================================================

void ExactU1(double, double, double *values)
{
  values[0] = 0; // y/2.*(1-y/2.); 
  values[1] = 0;
  values[2] = 0; // 1./2. - y/2.; 
  values[3] = 0; // -1./2.;         
}
void ExactU2(double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
void ExactP(double, double, double *values)
{
  values[0] = 0; // -x/2+1./2.;   
  values[1] = 0; // -1./2.;      
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================

void BoundCondition(int bdComp, double, BoundCond &cond)  
{
  //cond = DIRICHLET;
  cond = (bdComp==0 || bdComp==2 || bdComp==3 || bdComp==5) ? NEUMANN : DIRICHLET;

for (int j = 0; j < TDatabase::ParamDB->n_nitsche_boundary; j++)
  {
    if (bdComp == TDatabase::ParamDB->nitsche_boundary_id[j])
    {
      // Todo
      cond = DIRICHLET_WEAK;
      return;
    }
  }

 }

// ========================================================================
// boundary values
// ========================================================================

void U1BoundValue(int BdComp, double, double &value)
{
 // loop to impose Neumann boundary conditions
 // Since we are using the Neumann boundary condition via boundaryAssembling2d(), 
 // the boundvalue here has to be alway zero!!!!
 
 for (int j = 0; j < TDatabase::ParamDB->n_neumann_boundary; j++)
  {
    if ( BdComp == TDatabase::ParamDB->neumann_boundary_id[j])
    {
      switch(BdComp)
      {
        case 0:
          value = 0.; //0.5; //TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        case 2:
          value = 0.; //TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        case 3:
          value = 0.; //TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        case 5:
          value = 0.; //0.; //0.5; //TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        default:
          Output::print("I cannot impose Neumann boundary condition on component ", BdComp);
          exit(1);
          break;
      }
      return;
    }
  }
 
  switch(BdComp)
  {
    case 1: 
      value = 0; // 1.; // Dirichlet (reduces to u.n with Nitsche method for mueff=0)
      break;  
    case 4:  
      value = 0; // 1.; // Dirichlet (reduces to u.n with Nitsche method for mueff=0)
      break;
      default: 
      OutPut("BoundValue: wrong boundary part number " << BdComp <<endl);
      break;
  }
}



void U2BoundValue(int BdComp, double, double &value)
{
 // loop to impose Neumann boundary conditions
 // Since we are using the Neumann boundary condition via boundaryAssembling2d(), 
 // the boundvalue here has to be alway zero!!!!
 
 for (int j = 0; j < TDatabase::ParamDB->n_neumann_boundary; j++)
  {
    if ( BdComp == TDatabase::ParamDB->neumann_boundary_id[j])
    {
      switch(BdComp)
      {
        case 0:
          value = 0.; //0.5; //TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        case 2:
          value = 0.; //TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        case 3:
          value = 0.; //TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        case 5:
          value = 0.; //0.5; //TDatabase::ParamDB->neumann_boundary_value[j];
          break;
        default:
          Output::print("I cannot impose Neumann boundary condition on component ", BdComp);
          exit(1);
          break;
      }
      return;
    }
  }

  switch(BdComp)
  {
    case 1: 
      value = 0; // Dirichlet (reduces to u.n with Nitsche method for mueff=0)
      break;  
    case 4:  
      value = 0; // Dirichlet (reduces to u.n with Nitsche method for mueff=0)
      break;
    default: 
      OutPut("BoundValue: wrong boundary part number " << BdComp <<endl);
      break;
  }
}


// ========================================================================
// coefficients
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  //double eps = 1./TDatabase::ParamDB->RE_NR;

  double val_u1[4];
  double val_u2[4];
  double val_p[4];

  for(int i = 0; i < n_points; i++)
  {
   // physical parameters
    coeffs[i][4] = viscosity; 
    coeffs[i][5] = effective_viscosity; 
    coeffs[i][6] = permeability; 

  if (permeability == -2) 
    { 
     // porous region
     if (  ( ( x[i] >= 0. ) && (x[i] <= 0.9) && (y[i] <= 1./9. * x[i] + 1.5) )  ||
           ( ( x[i] > 0.9 ) && (x[i] <= 1.) && (y[i] <= (-1.) * x[i] + 2.5) )  ||
           ( ( x[i] > 1. ) && (x[i] <= 1.9) && (y[i] <= 1./9. * x[i] + 25./18.) ) ||
           ( ( x[i] > 1.9 ) && ( x[i] <= 2.) && (y[i] <= (-1.) * x[i] + 3.5) ) )
       {
       coeffs[i][5] = 0.; //0.0001;
       coeffs[i][4] = 1.;
       coeffs[i][6] = 0.000001;
       }
        else  // flow region
       {
      coeffs[i][4] = 0.;
      coeffs[i][5] = 0.001; //0.000001;
      //coeffs[i][6] = 1000000.0;
      }
 
    }
    /* ****************************************** */
    // Use an analytic_coefficient_function or input sol+mesh for the permeability field
    
   if (permeability == -1) 
    { 
      coeffs[i][6] = parameters[i][0];
    }

    /* ***************************************** */
    // if sources and sinks via analytic_coefficient_function (delta distr.), then \neq 0
    coeffs[i][9]=0;

    // parameter t^2 = mue/mu*K
    coeffs[i][0] = (coeffs[i][5] / coeffs[i][4]) * coeffs[i][6];

    ExactU1(x[i], y[i], val_u1);
    ExactU2(x[i], y[i], val_u2);
    ExactP(x[i], y[i], val_p);

   // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] =  0.; //-coeffs[i][5] * val_u1[3] + val_p[1] + (coeffs[i][4]/coeffs[i][6]) * val_u1[0]; 
   //(coeffs[4]/coeffs[6])-1;   // f1 (rhs of Brinkman problem for u1)  
    
    coeffs[i][2] =  0.; //-coeffs[i][5] * val_u2[3] + val_p[2] + (coeffs[i][4]/coeffs[i][6]) * val_u2[0]; 

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = 0; // val_u1[1] + val_u2[2]; 

    // stabilization parameters
    coeffs[i][7] = TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    coeffs[i][8] = TDatabase::ParamDB->grad_div_stab_weight;
  }
}


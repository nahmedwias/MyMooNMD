/*********************************************************
 * Analytical solution for the Brinkman 3D problem
 * - analogous to a Poiseuille flow
 * Stokes regime (sigma = 0):
   ux = Poiseuille solution (unitary pressure drop)
 * Darcy regime (effective_viscosity = 0):
   ux = 1/sigma
 * Brinkman regime:
   ux = exponential boundary layer
 **********************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1;
double sigma = -1;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman3D_Poiseuille.h");
}

// ========================================================================
// exact solution
// ========================================================================


void ExactU1(double x, double y, double z, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 0.;
}

void ExactU2(double x, double y, double z, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 0.;
}

void ExactU3(double x, double y, double z, double *values)
{
  // ratio between Stokes and Darcy terms
  double t = fabs(sqrt( effective_viscosity/sigma));

  if ( fabs(sigma) < 1e-10 )
  {
    // Stokes regime
    values[0] = 0.; //1./(2.*effective_viscosity) * ((1/4)*x+(1/2) + (1/4)*y+(1/2)) * (1.- ((1/4)*x+(1/2) + (1/4)*y+(1/2)));
    values[1] = 1./(2.*effective_viscosity) * (1. - 2.* (x));
    values[2] = 1./(2.*effective_viscosity) * (1. - 2.* (y));
    values[3] = 0.;
    values[4] = 2*(-1./(effective_viscosity));
  }
  else if  (t == 0)
  {
    // Darcy regime
    values[0] = 1./sigma;  
    values[1] = 0.;
    values[2] = 0.;
    values[3] = 0.;
    values[4] = 0.;
  }
  else
  {
    values[0] = (1./sigma) * (1+exp(2/t)-exp( (2- sqrt(x*x+y*y)  )/t) - exp( sqrt(x*x+y*y)/t)) / (1+exp(2/t));
        //(1./sigma) * (1+exp(1/t)-exp(1- ( (1/4)*x+(1/2)  )/t) - exp( ((1/4)*x+(1/2) )/t)) / (1+exp(1/t))+
                 //(1./sigma) * (1+exp(1/t)-exp(1- (  (1/4)*y+(1/2) )/t) - exp( ( (1/4)*y+(1/2) )/t)) / (1+exp(1/t));
    //(1./sigma) * (1+exp(1/t)-exp(1- ( (1/4)*x+(1/2) + (1/4)*y+(1/2) )/t) - exp( ((1/4)*x+(1/2) + (1/4)*y+(1/2) )/t)) / (1+exp(1/t));
    values[1] = (1./sigma) * (exp(1- ( (1/4)*x+(1/2) + (1/4)*y+(1/2) )/t) - exp( ((1/4)*x+(1/2) + (1/4)*y+(1/2) )/t)) / (4*t*(1+exp(1/t)));
    values[2] = (1./sigma) * (exp(1- ( (1/4)*x+(1/2) + (1/4)*y+(1/2) )/t) - exp( ((1/4)*x+(1/2) + (1/4)*y+(1/2) )/t)) / (4*t*(1+exp(1/t)));
    values[3] = 0.;
    values[4] = 2* (1./sigma) * (exp(1- ( (1/4)*x+(1/2) + (1/4)*y+(1/2) )/t) - exp( ((1/4)*x+(1/2) + (1/4)*y+(1/2) )/t)) / (16*t*t*(1+exp(1/t)));
  }
}

/*
  @attention we assume that pressure has zero average
  this will not be the exact solution (for pressure)
  if pressures with non zero mean are prescribed on 
  the boundaries
 */
void ExactP(double x, double y, double z, double *values)
{
  values[0] = 1-(2/10)*z;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = -(2/10);
  values[4] = 0.;
}

// ========================================================================
// boundary conditions
// ========================================================================


// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{

  cond = DIRICHLET_WEAK;
  //cond = DIRICHLET;

  // set Neumann BC

  // Neumann boundary condition on outlet  (top)
  double _R_CYLINDER = 2.;
  double _HEIGHT = 10.;
  double r = sqrt(x*x+y*y);
  if ((fabs(z-_HEIGHT)<1e-8) || (fabs(z-0)<1e-8) )
  {
    if( fabs(r - _R_CYLINDER)>1e-8  )
    {
      cond = NEUMANN;
    }
  }
//cout << "r: "<< r <<endl;
 /* else if (  ( fabs(r-2.) < 1.e-2 ) && (nitsche_id[0] == 2) )
  {
    cout << "DIICHLET-WEAK"<< endl;
    cout << x <<", "<< y <<", "<< z << endl;
    cond = DIRICHLET_WEAK;
  }
*/
  /*
   for (int j = 0; j < neumann_id.size(); j++)
    {
      if (i == neumann_id[j])
      {
        cond = NEUMANN;
        return;
      }
    }

  // set Nitsche BC
  for (int j = 0; j < nitsche_id.size(); j++)
  {
    if (i == nitsche_id[j])
    {
      cond = DIRICHLET_WEAK;

      return;
    }
  }
  */




  /*
  // Neumann boundary condition on inlet  (bottom)
   if (fabs(z)<1e-8) {
    if ( fabs(r2 - _R_CYLINDER*_R_CYLINDER)>1e-8) {
      cond = NEUMANN;
      //Output::print("NEUMANN on bottom");
    }
  }
  */

}

/*void BoundCondition(double x, double y, double z, BoundCond &cond) // (int i, double Param, BoundCond &cond)
{

  if (z == 0 || z == 10)
    cond = NEUMANN;
  else
    cond = DIRICHLET; // default





  // set Neumann BC
  for (int j = 0; j < neumann_id.size(); j++)
  {
    if (i == neumann_id[j])
    {
      cond = NEUMANN;
      return;
    }
  }

  // set Nitsche BC
  for (int j = 0; j < nitsche_id.size(); j++)
  {
    if (i == nitsche_id[j])
    {
      cond = DIRICHLET_WEAK;
      
      return;
    }
  }

}*/

void U1BoundValue(double x, double y, double z, double &value) // (int BdComp, double Param, double &value)
{
    value = 0.;


  /*
  double t = fabs(sqrt(effective_viscosity/sigma));

  // loop to impose Neumann boundary conditions
  ///@attention we set =0 here, as Neumann BC are imposed using boundary assembling
  for (int j = 0; j < neumann_id.size(); j++)
  {
    if ( BdComp == neumann_id[j])
    {
      switch(BdComp)
      {
      case 1:
        value = 0.; 
        break;
      case 3:
        value = 0.; 
        break;
      default:
        Output::print("I cannot impose Neumann boundary condition on component ",
            BdComp);
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
  case 1:
    if (fabs(sigma)<1e-10) {
      // Stokes solution
      value = 1./(2.*effective_viscosity) * Param * (1.-Param);
    } else {
      value = (1./sigma) *
          (1+exp(1./t)-exp((1.-Param)/t) - exp(Param/t)) / (1.+exp(1./t));
    }
    break;
  case 2: value = 0;
  break;
  case 3:
    if (fabs(sigma)<1e-10) {
      // Stokes solution
      value = 1./(2.*effective_viscosity) * Param * (1.-Param);
    } else {
      value = (1./sigma) *
          (1.+exp(1./t)-exp((Param)/t) - exp((1.-Param)/t)) / (1.+exp(1./t));
    }
    break;
  default: cout << "No boundary component with this number." << endl;
  break;
  }
    */
}

void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// ========================================================================
// coefficients for Brinkman problem:
// mu, f1,f2,g, sigma = mu/permeability
// with:
// -mu Delta u + grad(p) + sigma u = (f1,f2)
// div(u) = g
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
    double **parameters, double **coeffs)
{
  double val_u1[4];
  double val_u2[4];
  double val_u3[4];
  //double val_p[4];

  for(int i = 0; i < n_points; i++)
  {
    ExactU1(x[i], y[i], z[i], val_u1);
    ExactU2(x[i], y[i], z[i], val_u2);
    ExactU3(x[i], y[i], z[i], val_u3);
    //ExactP(x[i], y[i], z[i], val_p);

    // physical parameters
    coeffs[i][0] = effective_viscosity;
    coeffs[i][5] = sigma;

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][4] = 0; //val_u1[1] + val_u2[2] + val_u3[3];
  }
}



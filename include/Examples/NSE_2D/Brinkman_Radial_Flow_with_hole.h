/*********************************************************
 * Analytical solution for the Brinkman 2D problem
 * 
 * Consider a circular domain with radius r_1 and center point (0,0).
 * Inscribed to that circle is another on with smaller radius r_0 and the same center.
 * Both circle have to be considered as boundary.
 * (The appropriate mesh is sources_for_holes.mesh)
 * The radial velocity is fixed at the inner boundary and the radial pressure at the outer boundary as 
 * u(r=r0) = u0
 * p(r=R1) = 0
 * and then the solution is 
 * u(r) = u0 * r_0/r_1
 * p(r) = -sigma * u0 * r_0 * log(r/r_1) 
 *
 * r := sqrt(x^2 + y^2)
 * Hence,
 * u_1 = u0 * r_0 * x/(x^2 + y^2);
 * u_2 = u0 * r_0 * y/(x^2 + y^2);
 * p = -sigma * u0 * r_0 * 0.5 * log( (x^2 + y^2)/r_1^2 );
 * du_1/dx = u0 * r_0 * ( (-x^2 +y^2)/(x^2 + y^2)^2 );
 * d^2u_1/dx^2 = u0 * r_0 * ( 2 * x * (x^2-3y^2)/((x^2+y^2)^3))
 * du_1/dy = u0 * r_0 * ( -2xy/((x^2+y^2)^2) ); 
 * d^2u_1/dy^2 = u0 * r_0 * ( 2*x*(3y^2-x^2)/((x^2+y^2)^3) )
 * 
 * Delta (u_1,u_2) = 0
 * (f_1,f_2) = (0,0);
 **********************************************************/

// initialize physical parameters
// These should be reset when constructing the Example class
double effective_viscosity = -1.;
double sigma = -1.;
std::vector<size_t> neumann_id;
std::vector<size_t> nitsche_id;

double u0 = 1.;
double r_1 = 1.;
double r_0 = 0.05;

void ExampleFile()
{
  Output::print<1>("Example: Brinkman_Radial_Flow_with_hole.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    values[0] = u0 * r_0 * x/(x*x+y*y);        
    values[1] = u0 * r_0 * ( (-x*x+y*y)/((x*x+y*y)*(x*x+y*y)) );
    values[2] = u0 * r_0 * ( -2*x*y/((x*x+y*y)*(x*x+y*y)) );    
    values[3] = 0.; 
// values[3] = u0 * r_0 * (  2*x*(x*x-3*y*y)/((x*x + y*y)*(x*x + y*y)*(x*x + y*y)) + 2*x*(3*y*y-x*x)/((x*x + y*y)*(x*x + y*y)*(x*x + y*y))  );  
}

void ExactU2(double x, double y, double *values)
{
  values[0] = u0 * r_0 * y/(x*x + y*y);
  values[1] = u0 * r_0 * (-2*y*x/((x*x+y*y)*(x*x+y*y)));
  values[2] = u0 * r_0 * (y*y - x*x)/((x*x+y*y)*(x*x+y*y));
  values[3] = 0.; 
// values[3] = u0 * r_0 * (  2*y*(3*x*x-y*y)/((x*x+y*y)*(x*x+y*y)*(x*x+y*y)) + 2*y*(y*y-3*x*x)/((x*x+y*y)*(x*x+y*y)*(x*x+y*y))  );
}

/*
  @attention we assume that pressure has zero average
  this will not be the exact solution (for pressure)
  if pressures with non zero mean are prescribed on 
  the boundaries
 */

void ExactP(double x, double y, double *values)
{
  values[0] = -sigma * u0 * r_0 * 0.5 * log( (x*x + y*y)/(r_1*r_1) );     
  values[1] = -sigma * u0 * r_0 * 0.5 * 2*x/(x*x+y*y);
  values[2] = -sigma * u0 * r_0 * 0.5 * 2*y/(x*x+y*y);
  values[3] = 0.;
// values[3] = -sigma * u0 * r_0 * 0.5 * (  2*(-x*x+y*y)/((x*x+y*y)*(x*x+y*y)) + 2*(-y*y+x*x)/((x*x+y*y)*(x*x+y*y)) );
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double Param, BoundCond &cond)
{
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
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t = fabs(sqrt(effective_viscosity/sigma));

  // loop to impose Neumann boundary conditions
  ///@attention we set =0 here, as Neumann BC are imposed using boundary assembling
  for (int j = 0; j < neumann_id.size(); j++)
  {
    if ( BdComp == neumann_id[j])
    {
      switch(BdComp)
      {
      case 0:
        value = 0.;  //p(r=R1) = 0
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
  case 1: value = u0;  //u(r=r0) = u0
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
// coefficients for Brinkman problem:
// mu, f1,f2,g, sigma = mu/permeability
// with:
// -mu Delta u + grad(p) + sigma u = (f1,f2)
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

    // (f1,f2)(x,y): RHS for momentum equation
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;

    //g(x,y):  RHS for mass conservation equation
    coeffs[i][3] = val_u1[1] + val_u2[2];
    coeffs[i][4] = sigma;
  }
}



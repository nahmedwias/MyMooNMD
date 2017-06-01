// This is a test to try to define a domain with 2 phases in the square box
// This should be submitted to a velocity field and see how it behaves
// No Exact solution is known, so it is 0.
// Only Initial solution is known.
// Boundary conditions have to be taken care of.


void ExampleFile()
{
  Output::info<3>("Example: 40_DamBreakCD_NSE.h");
}

double get_nu()
{
  return 1;
}

constexpr bool rhs_depends_on_time = false;
constexpr bool coefficients_depend_on_time = true;

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double Param, BoundCond &cond)
{
//  if (BdComp == 1 || BdComp==2)
//    cond = DIRICHLET;
//  else
    cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
//  double rho_min = TDatabase::ParamDB->P7;
//  if (BdComp==1 || BdComp==2)
//    value = rho_min;  // this is Dirichlet
//  else
    value =0;
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
//  double vmin=0;
//  double vmax=1;
//  values[0] = vmin+(vmax-vmin)*exp(-10e13*(pow((x-0.3),10)+5*pow((y-0.7),10)));
//  double columnsize = 0.05715;
//  if ( x <= columnsize && y <= columnsize)
//  double density_ratio = TDatabase::ParamDB->P7;
//  double rho_min = TDatabase::ParamDB->P7; // this is the density of the bottom fluid

//  double phi = (1/0.05715)*(0.05715*0.05715-x*x) - y;
//  if (phi >=0)

  double dam_height = 0.05715;
  double dam_width  = 0.05715;

  /* Code for the sharp dam */
  if ( x <= dam_width && y <= dam_height)
    values[0] = 1;
  else
    values[0] = 0;

  /* Code for a smoother dam column */
//  // note that when height=width, the corner is a circle
//  // otherwise, one would need an ellipse
//  double corner_radius = 0.03;
//  double x0 = dam_width - corner_radius;
//  double y0 = dam_height - corner_radius;
//
//  if (x <= dam_width)
//  {
//    if (y <= y0)
//    {
//      values[0] = 1;
//    }
//    else if (y <= dam_height) // y is between y0 and height
//    {
//      if (x <= x0)
//      {
//        values[0] = 1;
//      }
//      else // x should now be between x0 and width
//      {
//        if (corner_radius*corner_radius - (x-x0)*(x-x0) - (y-y0)*(y-y0) >= 0)
//          values[0] = 1;
//        else
//          values[0] = 0;
//      }
//    }
//    else  // y should now be > height
//      values[0] = 0;
//  }
//  else
//    values[0] = 0;


  /* Code for a quarter of circle dam */
//  if ( x*x + y*y - dam_height*dam_height <= 0)
//    values[0]=1;
//  else
//    values[0]=0;

//  double A = 1;
//  double a = 20;
//  double b = 0;
//  double c = 20;
//  double x0 = 0.05;
//  double y0 = 0.07;
//  values[0]= A*exp(-a*pow(x-x0,2)+2*b*(x-x0)*(y-y0)-c*pow(y-y0,2));
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
//  double eps=1/TDatabase::ParamDB->RE_NR;
//  double a=1, b=2, c=1;
  double u_x,u_y;
  int i;
  double *coeff;
//  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    x = X[i];
//    y = Y[i];

    u_x = parameters[i][0];
    u_y = parameters[i][1];

    coeff[0] = 0;
    coeff[1] = u_x;
    coeff[2] = u_y;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
}


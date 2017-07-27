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
    cond = NEUMANN;
//  else
//    cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  if (BdComp==1 || BdComp==2)
    value = 0;
  else
    value = 0;
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  double dam_height = 0.05715;
  double dam_width  = 0.05715;
//  double dam_height_australia = 0.6;
//  double dam_width_australia  = 1.2;
//  double dam_height_openfoam = 0.292;
//  double dam_width_openfoam  = 0.1461;

//  /* Code for the sharp dam OPENFOAM example */
//  if ( x <= dam_width_openfoam && y <= dam_height_openfoam)
//    values[0] = 1;
//  else
//    values[0] = 0;

//  /* Code for the sharp dam */
//  if ( x <= dam_width && y <= dam_height)
//    values[0] = 1;
//  else
//    values[0] = 0;

//  /* Code for the sharp dam - variant */
//  double pi = 3.14159265358979;
//  double max = fabs(x) > fabs(y) ? fabs(x) : fabs(y);
//  double phi = max - dam_height;
//  double eps = 0.02;
//  double regularized_phi = 0.5*(1 - (phi/eps) - (1/pi)*sin(phi*pi/eps));
//  if ( phi < - eps)
//    values[0] = 1;
//  else if ( phi > eps)
//    values[0] = 0;
//  else
//    values[0] = regularized_phi;


//  /* Code for the sharp dam AUSTRALIA PAPER - FLUENT */
//  if ( x <= dam_width_australia && y <= dam_height_australia)
//    values[0] = 1;
//  else
//    values[0] = 0;

//  /* Code for a smoother dam column */
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
  double pi = 3.14159265358979;
  double phi = (x)*(x) + y*y - dam_height*dam_height;
  double eps = 0.001;
  double regularized_phi = 0.5*(1 - (phi/eps) - (1/pi)*sin(phi*pi/eps));
  if ( phi < - eps)
    values[0] = 1;
  else if ( phi > eps)
    values[0] = 0;
  else
    values[0] = regularized_phi;
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

void dambreak_postprocess(Time_CD2D& time_cd2d)
{
  const TFESpace2D *Space;
  TFEFunction2D& fct(time_cd2d.get_function());
  Space = fct.GetFESpace2D();
//  TBaseCell *cell;
//  TCollection *Coll;
//  Coll = Space->GetCollection();

  int N_DOFS = fct.GetLength();
  double x, y;

  int n_dofs_interface = 0;
  double epsilon = 1.e-1;

  double val[4];

  for(int i=0; i<N_DOFS;i++)
  {
    Space->GetDOFPosition(i,x,y);
    fct.FindGradient(x,y,val);

    if (val[0] > epsilon && val[0] < 1-epsilon)
    {
      Output::print<1>("DOF ", i, " c = ", val[0], " (x,y) = (", x, ", ", y, ")" );
      n_dofs_interface++;
    }
  }
  Output::print<1>("Number of dofs on interface = ", n_dofs_interface);

//    if ( (x-x0)*(x-x0)+(y-y0)*(y-y0)-R*R <= - epsilon )
//    {
//      pfct.FindGradient(x,y,val);
//      pmax += val[0];
////      Output::print<1>("dof i = ", i, " x = ", x, " y = ", y, " pmax = ", val[0]);
//      n_dofs_in++;
//    }
//    else if ( (x-x0)*(x-x0)+(y-y0)*(y-y0)-R*R >= epsilon )
//    {
//      pfct.FindGradient(x,y,val);
//      pmin += val[0];
////      Output::print<1>("dof i = ", i, " x = ", x, " y = ", y, " pmin = ", val[0]);
//      n_dofs_out++;
//    }
//  pmax /= n_dofs_in;
//  pmin /= n_dofs_out;
//
//  Output::print<1>("######## Pin = ", pmax);
//  Output::print<1>("######## Pout = ", pmin);
//  Output::print<1>("######## Delta P = ", pmax - pmin);
//  Output::print<1>("######## Relative error = ", ((pmax-pmin)-1.1805)*100/1.1805, "%");
}




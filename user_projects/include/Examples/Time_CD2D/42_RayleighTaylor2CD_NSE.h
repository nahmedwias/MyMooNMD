// This is a test to try to define a domain with 2 phases in the square box
// This should reproduce the Rayleigh-Taylor instability!
// No Exact solution is known, so it is 0.
// Only Initial solution is known.
// Boundary conditions have to be taken care of.
// IT IS A DIFFERENT SET UP THAN EXAMPLE 41.
// HERE, WE ARE BASED ON POCHET ET AL, 2013


void ExampleFile()
{
  Output::info<3>("Example: 42_RayleighTaylor2CD_NSE.h");
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
void BoundCondition(int i, double Param, BoundCond &cond)
{
  switch(i)
  {
    case 0:  // bottom
    case 9:  // top
      cond = DIRICHLET;
      break;
    case  1: case  2: case  3:
    case  4: case  5: case  6:
    case  7: case  8: case 10:
    case 11: case 12: case 13:
    case 14: case 15: case 16:
    case 17: case 18:
      cond = NEUMANN;     // right and left
     break;
    default:
      ErrThrow("Problem with Boundary number.");
  }
}

// value of boundary condition
void BoundValue(int i, double Param, double &value)
{
  // NOTE: THIS PART DEPENDS HIGHLY ON THE GEOMETRY
  // CHECK THE BOUNDARY IF IT IS DIVIDED INTO SEVERAL
  // PARTS!!
//  double rho_min = TDatabase::ParamDB->P7; // this is the density of the bottom fluid
  switch(i)
  {
    case 0:  // bottom
      value = 0;
      break;
    case 9:  // top
      value = 1;
      break;
    case  1: case  2: case  3:
    case  4: case  5: case  6:
    case  7: case  8: case 10:
    case 11: case 12: case 13:
    case 14: case 15: case 16:
    case 17: case 18:
       value = 0;     // right and left
     break;
    default:
      ErrThrow("Problem with Boundary number.");
  }
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
//  double rho_min = TDatabase::ParamDB->P7; // this is the density of the bottom fluid
  double pi = 3.14159265358979;
  double d = 1.0; // this is the width of the rectangular domain
  double interface0 = 2-(-0.05*d*(cos(2*pi*(x)/d)) ); // this is the initial interface, as given in (Pochet,13)

  // Initial fluids density field
  values[0] = 1*(0.5+0.5*tanh((y-interface0)/(0.01*d)));
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
    coeff[1] = u_x;//0;
    coeff[2] = u_y;//0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
  double rho_min = TDatabase::ParamDB->P7; // this is the density of the bottom fluid
  double pi = 3.14159265358979;
  double d = 1.; // this is the width of the rectangular domain
  double interface0 = -0.1*d*cos(2*pi*x/d); // this is the initial interface, as given in (Fraigneau et al,01)

  // Initial fluids density field
  values[0] = rho_min*(2+tanh((y-interface0)/(0.01*d)));
}


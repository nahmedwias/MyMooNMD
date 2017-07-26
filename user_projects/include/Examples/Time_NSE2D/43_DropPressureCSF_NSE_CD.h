// Navier-Stokes problem for the drop pressure 43_TCD2D
// Should calculate the drop pressure in a fluid at rest.
// The drop is in middle and subject to surface tension only.
// initial velocity =0 everyhwere. No gravity
// See paper Brackbille et al. 1996

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;



void ExampleFile()
{
  Output::info<3>("Example: 43_DropPressureCSF_NSE_CD.h ") ;
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = NEUMANN;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = REYNOLDS_number;
  double t = TDatabase::TimeDB->CURRENTTIME;
  int i;
  double *coeff;
//  double  x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    x = X[i];
//    y = Y[i];
//    double rho = parameters[i][2];
//    double mu  = parameters[i][3];

    coeff[0] = nu;

/*
    // Stokes
    coeff[1] =0;
    coeff[2] =0;
*/

    // Navier-Stokes
    coeff[1] = 0;     // f1
    if (t < 0.01)
      coeff[2] = 0;     // f2
    else
      coeff[2] = -10;     // f2
    coeff[3] = 0;
  }
}


void compute_pressure_drop(Time_NSE2D& time_nse2d)
{
//  const TFESpace2D *PSpace;
//  TFEFunction2D& pfct(time_nse2d.get_pressure());
//  PSpace = pfct.GetFESpace2D();
////  TBaseCell *cell;
////  TCollection *Coll;
////  Coll = PSpace->GetCollection();
//
//  int N_DOFS = pfct.GetLength();
//  double x, y;
//  double pmin = 0., pmax = 0.;
//  int n_dofs_in = 0, n_dofs_out = 0;
//  double x0 = TDatabase::ParamDB->P4;
//  double y0 = TDatabase::ParamDB->P5; // center of unit square=center of drop
//  double R = TDatabase::ParamDB->P6; // radius of drop, 2cm
//  double epsilon = 0.0001;
//  double val[4];
//
//  for(int i=0; i<N_DOFS;i++)
//  {
//    PSpace->GetDOFPosition(i,x,y);
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
//  }
//
//  pmax /= n_dofs_in;
//  pmin /= n_dofs_out;
//
//  Output::print<1>("######## Pin = ", pmax);
//  Output::print<1>("######## Pout = ", pmin);
//  Output::print<1>("######## Delta P = ", pmax - pmin);
//  Output::print<1>("######## Relative error = ", ((pmax-pmin)-1.1805)*100/1.1805, "%");
}

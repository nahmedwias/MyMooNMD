// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (4*y*(1-y), 0):=(u1,u2)
// p(x,y) = x-1/2

void ExampleFile()
{
  Output::print<1>("Example: Poiseuille.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 4*y*(1-y);    // u1
  values[1] = 0;            //u1_x
  values[2] = 4-8*y;        //u1_y
  values[3] = -8;           //Delta u1
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
  values[0] = x-0.5;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp) //counterclockwise numbering
  {
    case 0: value=0;
            break;
    case 1: value=4*Param*(1-Param); //Param \b [0,1]
            break;
    case 2: value=0;
            break;
    case 3: value=4*Param*(1-Param);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 1+8*eps + 4*y[i]*(1-y[i]); // f1 // 1 wegen p_x // 8eps= -epsDelta u
    coeff[2] = 0; // f2
    coeff[3] = 0; //f3  rechte Seite von Massenerhaltungsgl.
  }
}


// Brinkman problem, Poiseuille-Problem
// 
// u(x,y) = (4*y*(1-y), 0) = (u1,u2)
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
  values[0] = 4*y*(1-y);    //u1
  values[1] = 0;            //u1_x
  values[2] = 4-8*y;        //u1_y
  values[3] = -8;           //Delta u1
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
  values[0] = x-0.5;        //p
  values[1] = 1;            //p_x
  values[2] = 0;            //p_y
  values[3] = 0;            //Delta p=p_xx+p_yy
}

// ========================================================================
// boundary conditions (The boundary is Parametrized using Param \in [0,1])
// ========================================================================
void BoundCondition(int i, double Param, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=4*Param*(1-Param);
            break;
    case 2: value=0;
            break;
    case 3: value=4*(1-Param)*(Param);
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
// (the lhs of the Brinkman problem computed at quadrature points - for the error norms)
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

    coeff[0] = eps;                             //WAS IST DAS?
    //coeff[1] = 1 + 8*eps + 4*y*(1-y);
    coeff[1] = 1 + 8*eps + 4*y[i]*(1-y[i]);     // f1 (lhs of Brinkman problem for u1)
    coeff[2] = 0;                               // f2 (lhs of Brinkman problem for u2)
    coeff[3] = 0;                               // g (divergence term=u1_x+u2_y)

  }

}



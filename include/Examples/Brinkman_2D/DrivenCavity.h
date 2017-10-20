// Navier-Stokes problem, Driven cavity
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  Output::print<1>("Example: DrivenCavity.h");
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
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 0;
            break;
    case 2: if(Param<0.00001 || Param>0.99999) 
              value = 0;
            else
              value = 1;
            break;
    case 3: value = 0;
            break;
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
  if(BdComp>3) cout << "wrong boundary part number" << endl;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

 //   coeff[0] = eps;
 //   coeff[1] = 0; // f1
 //   coeff[2] = 0; // f2
        coeff[4] = TDatabase::ParamDB->VISCOSITY;
        coeff[5] = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
        coeff[6] = TDatabase::ParamDB->PERMEABILITY;
        coeff[0] = (coeff[5] / coeff[4]) * coeff[6];
        coeff[1] = 0;//-coeff[5] * val_u1[3] - val_p[1] + (coeff[4]/coeff[6]) * val_u1[0];  //(coeff[4]/coeff[6])-1;                       // f1 (rhs of Brinkman problem for u1)
        coeff[2] = 0;                                               // f2 (rhs of Brinkman problem for u2)
        coeff[3] = 0;                                               // g (divergence term=u1_x+u2_y)
        coeff[7] = TDatabase::ParamDB->equal_order_stab_weight_PkPk;
        coeff[8] = TDatabase::ParamDB->grad_div_stab_weight;
     }
}


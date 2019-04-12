// Brinkman problem, Benchmark problem
//
// u(x,y) = unknown
// p(x,y) = unknown

double viscosity = -1;                                                        
double effective_viscosity = -1;                                              
double permeability = -1;

void ExampleFile()
{
  OutPut("Example: Benchmark_Neum.h" << endl) ;
}

#define __BENCH__

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
  if (i == 1)
  {
    cond = NEUMANN;
  }
  else
    cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value= 0;
            break;
    case 2: value = 0;
            break;
    case 3: value = 1.2 * Param * (1 - Param); // 4*0.3
            break;
    case 4: value = 0;
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }  
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
  if(BdComp > 4) cout << "wrong boundary part number: " << BdComp << endl;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
 // double eps = 1/TDatabase::ParamDB->RE_NR;
  double *coeff;

  for(int i = 0; i < n_points; i++)
  {
    coeff = coeffs[i];

    coeff[4] = TDatabase::ParamDB->VISCOSITY;
    coeff[5] = TDatabase::ParamDB->EFFECTIVE_VISCOSITY;
    coeff[6] = TDatabase::ParamDB->PERMEABILITY;
    coeff[0] = (coeff[5] / coeff[4]) * coeff[6];
    coeff[1] = 0;//-coeff[5] * val_u1[3] - val_p[1] + (coeff[4]/coeff[6]) * val_u1[0];  //(coeff[4]/coeff[6])-1; // f1 (rhs of Brinkman problem for u1)
    coeff[2] = 0;                                               // f2 (rhs of Brinkman problem for u2)
    coeff[3] = 0;                                               // g (divergence term=u1_x+u2_y)
    coeff[7] = TDatabase::ParamDB->equal_order_stab_weight_PkPk;
    coeff[8] = TDatabase::ParamDB->grad_div_stab_weight;
  }
}

// Brinkman problem, First Problem from Hannukainen:
// Computations with Finite Element Methods for the Brinkman
// 
// u(x,y) = (cos(x)*sinh(y) , sin(x)*cosh(y)) = (u1,u2)
// p(x,y) = -sin(x)*sinh(y)+C
// C is chosen such that p \in L_0^2
// div u=0
// Au=0
// ========================================================================


void ExampleFile()
{
  Output::print<1>("Example: hannukainen.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    values[0] = cos(x)*sinh(y);                                                //u1
    values[1] = -sin(x)*sinh(y);                                               //u1_x
    values[2] = cos(x)*cosh(y);                                                //u1_y
    values[3] = 0;                 //-cos(x)*sinh(y) + sinh(y)*cos(x)          //Delta u1
}

void ExactU2(double x, double y, double *values)
{
    values[0] = sin(x)*cosh(y);                                                //u2
    values[1] = cos(x)*cosh(y);                                                //u2_x
    values[2] = sin(x)*sinh(y);                                                //u2_y
    values[3] = 0;                 //-sin(x)*cosh(y) + cosh(y)*sin(x)          //Delta u2
}

void ExactP(double x, double y, double *values)
{
    double C = cosh(1)+cos(1)-cos(1)*cosh(1)-1; //chosen s.t. p \in L_0^2
    values[0] = -sin(x)*sinh(y)+C;                                              //p
    values[1] = -cos(x)*sinh(y);                                                //p_x
    values[2] = -sin(x)*cosh(y);                                                //p_y
    values[3] = 0;                 // sin(x)*sinh(y) -sinh(y)*sin(x)            //Delta p=p_xx+p_yy
}

// ========================================================================
// boundary conditions (The boundary is parametrized using Param \in [0,1].)
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
    case 0: value=0;                        //cos(Param)*sinh(0)=cos(Param)*0=0
            break;
    case 1: value=cos(1)*sinh(Param);
            break;
    case 2: value=cos(1-Param)*sinh(1);
            break;
    case 3: value=cos(0)*sinh(1-Param);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
    {
    case 0: value=sin(Param);               //sin(Param)*cosh(0)=sin(Param)*1
            break;
    case 1: value=sin(1)*cosh(Param);
            break;
    case 2: value=sin(1-Param)*cosh(1);
            break;
    case 3: value=sin(0)*cosh(1-Param);
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
    static double eps = 1./TDatabase::ParamDB->RE_NR;                       // fixed by Reynolds number in ".dat"
  static double K= TDatabase::ParamDB->SIGMA_PERM;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;                                                         // viscosity
    coeff[1] = -cos(x[i])*sinh(y[i]) + K * cos(x[i])*sinh(y[i]);            // f1 (lhs of Brinkman problem for u1)
    coeff[2] = -sin(x[i])*cosh(y[i]) + K * sin(x[i])*cosh(y[i]);            // f2 (lhs of Brinkman problem for u2)
    coeff[3] = 0;                                                           // g (divergence term=u1_x+u2_y)
    coeff[4] = K;                                                           // permeability

  }

}



// Brinkman problem, Poiseuille-Problem from Hannukainen: Computations...
// 
// u(x,y) = (cos(x)*sinh(y) , sin(x)*cosh(y))           // u(x,y) = (4*y*(1-y), 0) = (u1,u2)
// p(x,y) = -sin(x)*sinh(y)+C           // p(x,y) = x-1/2
// C is chosen such that p \in L_0^2
// div u=0
// Au=0




void ExampleFile()
{
  Output::print<1>("Example: Poiseuille_Hannukainen.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    values[0] = cos(x)*sinh(y);         //4*y*(1-y);    //u1
    values[1] = -sin(x)*sinh(y);        //0;            //u1_x
    values[2] = cos(x)*cosh(y);         //4-8*y;        //u1_y
    values[3] = 0;            //-cos(x)*sinh(y) + sinh(y)*cos(x)          //-8;           //Delta u1
}

void ExactU2(double x, double y, double *values)
{
    values[0] = sin(x)*cosh(y);         //0;            //u2
    values[1] = cos(x)*cosh(y);         //0;            //u2_x
    values[2] = sin(x)*sinh(y);         //0;            //u2_y
    values[3] = 0;        // -sin(x)*cosh(y) + cosh(y)*sin(x)    //Delta u2
}

void ExactP(double x, double y, double *values)
{
    double C = cosh(1)+cos(1)-cos(1)*cosh(1)-1; //chosen s.t. p \in L_0^2
    values[0] = -sin(x)*sinh(y)+C;          //x-0.5;        //p
    values[1] = -cos(x)*sinh(y);          //1;            //p_x
    values[2] = -sin(x)*cosh(y);          //0;            //p_y
    values[3] = 0;         // sin(x)*sinh(y) -sinh(y)*sin(x)   //Delta p=p_xx+p_yy
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
    case 0: value=0;            //sin(Param)*cosh(0)=sin(Param)*0=0
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
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  static double K= TDatabase::ParamDB->SIGMA_PERM;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;                             //WAS IST DAS?
      coeff[1] = -cos(x[i])*sinh(y[i]) + K * cos(x[i])*sinh(y[i]);
   // coeff[1] = 1 + 8*eps + 4*y[i]*(1-y[i]);     // f1 (lhs of Brinkman problem for u1)
      coeff[2] = -sin(x[i])*cosh(y[i]) + K * sin(x[i])*cosh(y[i]);
   // coeff[2] = 0;                               // f2 (lhs of Brinkman problem for u2)
    coeff[3] = 0;                                 // g (divergence term=u1_x+u2_y)
    coeff[4] = K;

  }

}



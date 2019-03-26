// ======================================================================
// convection skew to the domain
// Hughes, Mallet, Mizukami 1986
// ======================================================================
#define __HMM_1986__

void ExampleFile()
{
  OutPut("Example: HMM1986.h" << endl) ;
  TDatabase::ParamDB->INTERNAL_OUTFLOW_BOUNDARY[0]= 1;
  TDatabase::ParamDB->INTERNAL_OUTFLOW_BOUNDARY[1]= 1;
  TDatabase::ParamDB->INTERNAL_OUTFLOW_BOUNDARY[2]= 0;
  TDatabase::ParamDB->INTERNAL_OUTFLOW_BOUNDARY[3]= 0;
}
// exact solution (this is the solution for eps = 0)
void Exact(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

// kind of boundary condition
void BoundCondition(int i, double t, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
    switch(BdComp)
    {
	case 0:
	case 1:
	    value = 0;
	    break;
	case 2:
	    if (Param < 1e-6)
		value = 0;
	    else
		value = 1;
	    break;
	case 3:
	    if (Param<0.3)
		value = 1;
	    else
		value = 0;
	    break;
    }
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, arg;

  arg = -Pi/3;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = cos(arg);
    coeff[2] = sin(arg);
    coeff[3] = 0;
    coeff[4] = 0;
  }
}
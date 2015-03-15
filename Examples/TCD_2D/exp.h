// ======================================================================
// instationary problem
// ======================================================================

/// ========================================================================
// example file
// ========================================================================

#define __SIN3__

void ExampleFile()
{
  OutPut("Example: exp.h" << endl); 
}

// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
  values[1] = exp(t)*2*Pi*cos(2*Pi*x)*sin(2*Pi*y);
  values[2] = exp(t)*2*Pi*sin(2*Pi*x)*cos(2*Pi*y);
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
 cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
}


void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->PE_NR;
  double b1=1, b2=-1, c=1;
  int i;
  double *coeff;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;

    coeff[4] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y))
	- eps * exp(t)*4*Pi*Pi*(-sin(2*Pi*x)*sin(2*Pi*y)-sin(2*Pi*x)*sin(2*Pi*y))
       + b1 * exp(t)*2*Pi*cos(2*Pi*x)*sin(2*Pi*y)
       + b2 * exp(t)*2*Pi*sin(2*Pi*x)*cos(2*Pi*y)
	+ c *  exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
  }
}

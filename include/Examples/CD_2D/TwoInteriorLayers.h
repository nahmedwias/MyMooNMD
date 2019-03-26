// ======================================================================
// sharp characteristic interior layer
// Knopp, Lube, Rapin, CMAME 2002
// ======================================================================
#define __TWO_INTERIOR_LAYERS__

void ExampleFile()
{
  Output::print<1>("Example: TwoInteriorLayers.h");
}
// exact solution
void Exact(double, double, double *values)
{

  values[0]= 0;
  values[1]= 0;
  values[2]= 0;
  values[3]= 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double, BoundCond &cond)
{
    if (i==3)
	cond = NEUMANN;
    else
	cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
   if (BdComp==0)
   {
      if ((Param>1.0/3.0)&& (Param<2.0/3.0))
         value = 1;
      else
         value = 0;
   }
   else
      value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y, double **,
                    double **coeffs)
{
  double eps = 1./TDatabase::ParamDB->PE_NR;
  for(int i = 0; i < n_points; i++)
  {
    double x = X[i];
    double y = Y[i];
    coeffs[i][0] = eps;
    coeffs[i][1] = -y;
    coeffs[i][2] = x;
    coeffs[i][3] = 0;
    coeffs[i][4] = 0;
    coeffs[i][5] = std::sqrt(x*x + y*y);
  }
}


 

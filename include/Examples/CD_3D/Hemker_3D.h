/**
 * This is the originally 2D Hemker example adopted to the 3D problem
 *
 */
double PECLET_NUMBER;

void ExampleFile()
{
  Output::print("Example: Hemker_3D.h");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


// Helper function - finds out if (x,y,z) is on the interior of the outflow face
bool on_outflow_interior(double x, double y, double z)
{
  double tol = 1e-10;
  return fabs(9-x) < tol
      && !(fabs(y -3) < tol)  && !(fabs(y+3) < tol)
      && !(fabs(z -3) < tol)  && !(fabs(z+3) < tol);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  if(on_outflow_interior(x,y,z)) // Neumann only on interior of the outflow face 
    cond = NEUMANN;
  else
  {
    cond = DIRICHLET;
  }
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  value = 0;
  if(fabs(x*x + y*y - 1) < 1e-10)
  {
    value = 1;
  }
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  double eps= PECLET_NUMBER;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps; //diffusion coefficient
    coeff[1] = 1;   //ux
    coeff[2] = 0;   //uy
    coeff[3] = 0;   //uz
    coeff[4] = 0;   //reaction coefficient
    coeff[5] = 0;   //rhs
  }
}


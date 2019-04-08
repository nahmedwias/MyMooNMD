//Example 6 from Allendes, Barrenechea and Rankin, SIAM-SCI17

double PECLET_NUMBER;

void ExampleFile()
{
  Output::root_info("Example", "boundary_layer_known_3d.h");
}

// exact solution is not known
void Exact(double x, double y, double z, double *values)
{
  double fx, fe,gx;
  fx=exp((x-1.0)/PECLET_NUMBER);
  fe=1.0/(exp(-1.0/PECLET_NUMBER)-1.0);
  gx=(x+(fx-exp(-1./PECLET_NUMBER)*fe));
  
  values[0] = y*z*(1-y)*(1-z)*gx; //original solution
  values[1] = y*z*(1-x)*(1-y)*(1-z)*(1+(fe*fx)/PECLET_NUMBER);  //x-derivative
  values[2] = z*(1-y)*(1-z)*gx-y*z*(1-z)*gx;  //y-derivative
  values[3] = y*(1-y)*(1-z)*gx-y*z*(1-y)*gx;  //z-derivative
  values[4] = -2*z*(1-z)*gx-2*y*(1-y)*gx+(y*z*(1-y)*(1-z)*fe*fx)/(PECLET_NUMBER*PECLET_NUMBER);    //laplacian
}

// kind of boundary condition (needed for FE space)
void BoundCondition(double, double, double, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double, double, double, double &value)
{
  value=0;
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
                    double **, double **coeffs)
{
  double exact[5];
  for(int i = 0; i < n_points; ++i)
  {
    coeffs[i][0] = PECLET_NUMBER; //diffusion coefficient
    coeffs[i][1] = 1;   //x coordinate of convection b
    coeffs[i][2] = 0;   //y coordinate of convection b
    coeffs[i][3] = 0;   //z coordinate of convection b
    coeffs[i][4] = 0;   //reaction coefficient c
    
    
    Exact(x[i],y[i],z[i],exact);
    coeffs[i][5] = -coeffs[i][0]*exact[4]; // diffusion
    coeffs[i][5] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2] + coeffs[i][3]*exact[3]; // convection
    coeffs[i][5] += coeffs[i][4]*exact[0]; // reaction
  }
}


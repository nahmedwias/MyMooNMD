//Simple Example on the cube [0,1]^3
//with Dirichlet boundary condition

#include <cmath>

void ExampleFile()
{
  Output::print("A simple SinCos Darcy program on [0,1]^3 with Dirichlet ",
                "boundary condition using mixed finite elements") ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[1] = Pi*Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[2] = Pi*Pi*sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[3] = Pi*Pi*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[4] = -3*Pi*Pi*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
}

void ExactU2(double x, double y,double z, double *values)
{
  values[0] = -Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[1] = Pi*Pi*sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[2] = Pi*Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[3] = -Pi*Pi*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[4] = 3*Pi*Pi*Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = -Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[1] = Pi*Pi*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[2] = -Pi*Pi*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[3] = Pi*Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[4] = 3*Pi*Pi*Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
}


void ExactP(double x, double y, double z, double *values)
{
  values[0] = cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[1] = -Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[2] = Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[3] = Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[4] = -3*Pi*Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void FluxBoundValue(double x, double y, double z, double &value)
{
  double tol = 1e-10;
   bool pointOnBoundary = true;
   if(std::abs(x) <= tol) // (x == 0.0)
   {
     if( (y>=0 && y<=1) && ( z>=0 && z<=1 ))
       value = 0; // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(x-1.0) <= tol) //(x == 1.0)
   {
     if( (y>=0 && y<=1) && ( z>=0 && z<=1 ))
       value = 0; // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(y) <= tol) // (y == 0.0)
   {
     if( (x>=0 && x<=1) && ( z>=0 && z<=1 ))
       value = Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z); // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(y-1.0) <= tol) // (y == 1.0)
   {
     if( (x>=0 && x<=1) && ( z>=0 && z<=1 ))
       value = -Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z); // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(z) <= tol) // (z == 0.0)
   {
     if( (x>=0 && x<=1) && ( y>=0 && y<=1 ))
       value = Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z); // Dirichlet
     else
       pointOnBoundary = false;
   }
   else if(std::abs(z-1.0) <= tol) // (z == 1.0)
   {
     if( (x>=0.0 && x<=1.0) && ( y>=0.0 && y<=1.0 ))
       value = -Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z); // Dirichlet
     else
       pointOnBoundary = false;
   }
   else
     pointOnBoundary = false;

   if(!pointOnBoundary)
   {
     ErrMsg("trying to evaluate boundary data at a point not belonging to the "
            << "boundary, (" << x << "," << y << "," << z << ")\n");
     exit(0);
   }
}

void BoundConditionPressure(double, double, double, BoundCond &cond)
{
  cond = NEUMANN; 
}

void PressureBoundValue(double, double, double, double &value)
{
  value = 0;
}

// coefficients in the pde
void LinCoeffs(int n_points, double *X, double *Y,double *Z, double **,
               double **coeffs)
{
  const double eps = 1./TDatabase::ParamDB->SIGMA_PERM;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = eps;      // permeability
    coeffs[i][1] = 0;        // f1, forces
    coeffs[i][2] = 0;        // f2, forces
    coeffs[i][3] = 0;        // f3, forces
    //source terms:
    coeffs[i][4] = 3 * Pi * Pi * cos(Pi*X[i]) * sin(Pi*Y[i]) * sin(Pi*Z[i]);
  }
}



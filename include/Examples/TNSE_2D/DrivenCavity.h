#ifndef DRIVEN_CAVITY_H
#define DRIVEN_CAVITY_H


// Navier-Stokes problem, Driven Cavity
// 
// u(x,y) = ?
// p(x,y) = ?
double DIMENSIONLESS_VISCOSITY;


void ExampleFile()
{
  Output::print("Example: DrivenCavity.h");  
}
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
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
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double x_1=0.1;
  double x = 1- Param;
  
  if(BdComp>3)
  {
    ErrThrow( "ERROR in file " , __FILE__ , ", line: ",  __LINE__ ,
              ": wrong boundary part number: " , BdComp);
  }
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2:
    {
       if (x <= x_1)
       {
          value = cos(Pi*(x_1-x)/x_1);
          value= 1 -(1 - value)*(1- value)/4;
       }
       else
       {
         if (x >= 1-x_1)
         {
          value = cos(Pi*(x-(1-x_1))/x_1);
          value= 1 -(1 - value)*(1- value)/4;
         }
         else
           value = 1.0;
       }
            break;
    }
    case 3: value=0;
            break;
    default: Output::print( "wrong boundary part number" );
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  for(int i=0;i<n_points;i++)
  {
    double *coeff = coeffs[i];

    coeff[0] = DIMENSIONLESS_VISCOSITY;

    coeff[1] = 0;  // f1
    coeff[2] = 0;  // f2
  }
}


#endif

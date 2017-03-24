#ifndef CHANNELSTEP_H
#define CHANNELSTEP_H

double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::print("Example: ChannelStep.h (Dirichlet boundary conditions)");
}
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0 ;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}

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

void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double y;

  y = 2*(1-Param);

  switch(BdComp)
  {
    case 0: 
    case 1: 
    case 2: 
    case 3: 
    case 4: 
    case 6: 
      value=0;
      break;
     case 5:
     case 7:  
       if (t <= 4)
         value = sin(Pi*t/8) * y *(2-y)*1.5;
       else
         value = y *(2-y)*1.5;
      break;
    default: cout << "wrong boundary part number" << endl;
    exit(4711);
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value=0;
}

void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  for(int i=0; i<n_points; i++)
  {
    double *coeff = coeffs[i];

    coeff[0] = DIMENSIONLESS_VISCOSITY;
    coeff[1] = 0;
    coeff[2] = 0;
  }
}
#endif // CHANNELSTEP_H

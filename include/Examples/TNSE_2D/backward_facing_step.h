#include "NSE_2D/backward_facing_step.h"

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = y >= 0. ? 4*y*(1-y) : 0.;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}

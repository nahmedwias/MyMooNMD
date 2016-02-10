/**
 * An example for two reacting species in a flow, quickly cobbled together
 * for to test the coupled solver module.
 * Draws "inspiration" and some of its functions from the example file in
 * Example/CD_2D/SharpBoundaryLayer.
 */

#ifndef USER_PROJECTS_EXAMPLES_TEST_TIME_H_
#define USER_PROJECTS_EXAMPLES_TEST_TIME_H_

void ExampleFile()
{
  Output::print<1>(
      "Example: Test_Time.h. Is supposed to be used in time dependent case. "
      "No exact solutions known, put to 0. "
      "Interpret error values as norm of solution.");
}

// terms responsible for the coupled reaction
double couplingTerm1(const double* const c)
{
  double rate = 3;
  return rate * c[0] * c[1]; //consumption
}

double couplingTerm2(const double* const c)
{
  double rate = 3;
  return - rate * c[0] * c[1]; //production
}

// Unknown exact solution - put to zero.
void ExactC1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// Unknown exact solution - put to zero.
void ExactC2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition - same for both species
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  // numbering of boundary components starts with 0 at the y==0 bdry
  // and proceeds counterclockwise for UnitSquare.PRM
  if(BdComp==2)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

// value of boundary condition for species 1
void BoundValueC1(int BdComp, double Param, double &value)
{
  //on one Dirichlet boundary: jump function
  if(BdComp == 3){
    if(Param < 0.5)
      value = 1;
    else
      value = 0;
  }
  //on the other boundaries zero
  else
    value = 0;
}

// value of boundary condition for species 2
void BoundValueC2(int BdComp, double Param, double &value)
{
  //on one Dirichlet boundary: jump function
  if(BdComp == 3){
    if(Param < 0.5)
      value = 0;
    else
      value = 1;
  }
  //on the other boundaries zero
  else
    value = 0;
}

void BilinearCoeffsC1(int n_points, double *x, double *y,
                      double **parameters, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    double diffusion_1 = 0.1;
    coeffs[i][0] = diffusion_1; //diffusion coefficient
    coeffs[i][1] = cos(Pi/18);//convection in x direction: cos(10 deg)
    coeffs[i][2] = sin(Pi/18);//convection in y direction: sin(10 deg)
    coeffs[i][3] = 0; //reaction coefficient

    coeffs[i][4] = 0; //rhs, outer body force coefficients
  }
}

void BilinearCoeffsC2(int n_points, double *x, double *y,
                      double **parameters, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    double diffusion_2 = 0.2;
    coeffs[i][0] = diffusion_2; //diffusion coefficient
    coeffs[i][1] = cos(Pi/18);//convection in x direction: cos(10 deg)
    coeffs[i][2] = sin(Pi/18);//convection in y direction: sin(10 deg)
    coeffs[i][3] = 0; //reaction coefficient

    coeffs[i][4] = 0; //rhs, outer body force coefficients
  }
}

void InitialConditionC1(double x,  double y, double *values)
{
  values[0] = 0;
}

void InitialConditionC2(double x,  double y, double *values)
{
  values[0] = 0;
}

/****************************************************************************************
 * ParameterFunction and AssemblingFunctions used in the "Linearized Decoupled" strategy.
 ****************************************************************************************/

void ParameterFunction(double* in, double* out){
  // Just skip the first two entries - this is where TAuxParam2D->GetParameters() places the x and y value,
  // and pass on as many values as nCoupled_ - for this example that would be 2.
  for (size_t i = 0; i<2 ; ++i){
    out[i] = in [i+2];
  }
}

void AssemblingFunctionC1(
    double Mult, double *coeff, double *param,
    double hK, double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs)
{
  // For convenience: rename the place to write to (Rhs) and the place to read from (Orig,
  // which suposedly contains values of test functions)
  double* Rhs = LocRhs[0];
  double* Orig = OrigValues[0];

  // Calculate the value of the coupling term at the quad point.
  // This relies on the correct interaction with the ParameterFunction and the
  // input order of fe functions to the aux object.

  // coupled term evaluated at the quad point this
  // AssembleFctParam2D is called upon
  double coupledTerm= couplingTerm1(param);

  //Loop over all local base functions.
  for(int i=0;i<N_BaseFuncts[0];i++)
  {
    Rhs[i] -= Mult*coupledTerm*Orig[i]; //(Mult contains quad weigth, and maybe even the determinant due to trafo)
  }
}

void AssemblingFunctionC2(
    double Mult, double *coeff, double *param,
    double hK, double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs)
{
  // For convenience: rename the place to write to (Rhs) and the place to read from (Orig,
  // which suposedly contains values of test functions)
  double* Rhs = LocRhs[0];
  double* Orig = OrigValues[0];

  // Calculate the value of the coupling term at the quad point.
  // This relies on the correct interaction with the ParameterFunction and the
  // input order of fe functions to the aux object.

  // coupled term evaluated at the quad point this
  // AssembleFctParam2D is called upon
  double coupledTerm= couplingTerm2(param);

  //Loop over all local base functions.
  for(int i=0;i<N_BaseFuncts[0];i++)
  {
    Rhs[i] -= Mult*coupledTerm*Orig[i]; //(Mult contains quad weigth, and maybe even the determinant due to trafo)
  }
}



#endif /* USER_PROJECTS_EXAMPLES_TEST_TIME_H_ */

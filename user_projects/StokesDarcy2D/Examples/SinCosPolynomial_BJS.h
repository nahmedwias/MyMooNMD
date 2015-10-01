/*
 Benchmark for Stokes-Darcy
 domain Omega = (0,pi)x(-1,1)
 subdomains Omega_Darcy = (0,Pi)x(-1,0)
            Omega_Stokes= (0,Pi)x(0,1)
 Interface Gamma = (0,Pi)x{0}
 Boundary conditions:
      Dirichlet on Stokes-velocity
      Dirichlet on Darcy-pressure on (0,Pi)x{-1} 
      Neumann on {0,Pi}x(-1,0)


 This Example has been taken from
 "A parallel Robin-Robin domain decomposition method for the Stokes-Darcy 
 system" by Wenbin Chen, Max Gunzburger, Fei Hua, Xiaoming Wang
 
 In the definition of v there is a K^2 instead of a K in the paper. This must 
 be a mistake, I guess, because only with a K the Beavers-Joseph-Saffmann (BJS) 
 condition is fulfilled
*/

void ExampleFile()
{
  OutPut("Example: SinCosPolynomial_BJS.h\n");
  
  // set the boundary condition u.t + alpha t.T.n = 0 on the interface (BJS)
  TDatabase::ParamDB->StoDa_interfaceType = 0;
  if(TDatabase::ParamDB->SIGMA_PERM!=1.0 || TDatabase::ParamDB->RE_NR!=1.0)
  {
    OutPut("\n\nWARNING: This example requires reynolds number = 1 and "
         << "Permeability = 1.\n\n\n");
  }
  TDatabase::ParamDB->SIGMA_PERM = 1;
  TDatabase::ParamDB->RE_NR = 1;
}

// ============================================================================
// exact solution
// ============================================================================
void ExactP_Darcy(double x, double y, double *values)
{
  //const double K = TDatabase::ParamDB->SIGMA_PERM;
  //const double nu= 1/TDatabase::ParamDB->RE_NR;
  values[0] = exp(y)*sin(x);
  values[1] = exp(y)*cos(x);
  values[2] = exp(y)*sin(x);
  values[3] = 0;
}
void ExactU1_Darcy(double x, double y, double *values)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  values[0] = -K*exp(y)*cos(x);
  values[1] = K*exp(y)*sin(x);
  values[2] = -K*exp(y)*cos(x);
  values[3] = 0; // Laplace not needed
}
void ExactU2_Darcy(double x, double y, double *values)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  values[0] = -K*exp(y)*sin(x);
  values[1] = -K*exp(y)*cos(x);
  values[2] = -K*exp(y)*sin(x);
  values[3] = 0; // Laplace not needed
}
void ExactU1_NSE(double x, double y, double *values)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nu= 1/TDatabase::ParamDB->RE_NR;
  const double alpha = TDatabase::ParamDB->StoDa_alpha;
  values[0] = cos(x)*( - 1/(2*nu) + (-alpha/(2*nu*nu) + K)*y);
  values[1] = -sin(x)*( - 1/(2*nu) + (-alpha/(2*nu*nu) + K)*y);
  values[2] = cos(x)*( (-alpha/(2*nu*nu) + K));
  values[3] = -cos(x)*( - 1/(2*nu) + (-alpha/(2*nu*nu) + K)*y); //Laplacien
}
void ExactU2_NSE(double x, double y, double *values)
{
  double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu= 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  values[0] = sin(x)*(-K - y/(2*nu) + (-alpha/(4*nu*nu) + K/2)*y*y);
  values[1] = cos(x)*(-K - y/(2*nu) + (-alpha/(4*nu*nu) + K/2)*y*y);;
  values[2] = sin(x)*( - 1/(2*nu) + (-alpha/(2*nu*nu) + K)*y);;
  values[3] = 0; //Laplacien
}
void ExactP_NSE(double x, double y, double *values)
{
  //static double K = TDatabase::ParamDB->SIGMA_PERM;
  //static double nu= 1/TDatabase::ParamDB->RE_NR;
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ============================================================================
// boundary conditions
// ============================================================================
void BC_Darcy(int bdComp, double t, BoundCond &cond)  
{
  switch(bdComp)
  {
    case 0: 
      //cond = NEUMANN;
      cond = DIRICHLET; 
      break;      
    case 1:  
      cond = DIRICHLET;  
      break;
    case 2: 
      //cond = NEUMANN;
      cond = DIRICHLET;  
      break;
    default: cout << "BoundCondition_Darcy: wrong boundary part number" << endl;
      break;
  }
}
void BC_Velocity_Darcy(int bdComp, double t, BoundCond &cond)
{
  switch(bdComp)
    {
      case 0: 
        cond = DIRICHLET; 
        break;      
      case 1:  
        cond = NEUMANN;  
        break;
      case 2: 
        cond = DIRICHLET;  
        break;
      default: cout << "BC_Velocity_Darcy: wrong boundary part number" << endl;
        break;
    }
}
void BC_Pressure_Darcy(int bdComp, double t, BoundCond &cond)
{ cond = NEUMANN; }
void BoundCondition_NSE(int bdComp, double t, BoundCond &cond)  
{
  cond = DIRICHLET;
  if(bdComp==0 || bdComp==1 || bdComp==2)
    cout << "BoundCondition_NSE: wrong boundary part number";
}
void BoundConditionPressure_NSE(int bdComp, double t, BoundCond &cond)  
{
  cond = NEUMANN;
  if(bdComp==0 || bdComp==1 || bdComp==2)
    cout << "BoundConditionPressure_NSE: wrong boundary part number";
}

// ============================================================================
// boundary values
// ============================================================================
void BoundValue_Darcy(int bdComp, double t, double &value)
{
  switch(bdComp)
  {
    case 0: // x==0, y==-t
      //value = -exp(-t); // NEUMANN  
      value = 0; // DIRICHLET
      break;      
    case 1: //y==-1 x==t*pi
      // value = -exp(-1)*sin(Pi*t); // NEUMANN
      value = exp(-1)*sin(Pi*t); // DIRICHLET
      break;
    case 2: // x==Pi, y==-1+t
      //value = -exp(-1+t); // NEUMANN
      value = 0; //DIRICHLET 
      break;
    default: cout << "BoundValue_Darcy: wrong boundary part number" << endl;
      break;
  }
}
void BoundValueVelocity_Darcy(int bdComp, double t, double &value)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  switch(bdComp)
  {
    case 0: // DIRICHLET
      value = K*exp(-t);
      break;      
    case 1: // NEUMANN
      value = exp(-1.0)*sin(Pi*t);
      break;
    case 2: // DIRICHLET
      value = K*exp(-1+t);
      break;
    default: cout << "BoundValue_Darcy: wrong boundary part number" << endl;
      break;
  }
}
void BoundValuePressure_Darcy(int bdComp, double t, double &value)
{
  value = 0;
}
void U1BoundValue_NSE(int bdComp, double t, double &value)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nu= 1/TDatabase::ParamDB->RE_NR;
  const double alpha = TDatabase::ParamDB->StoDa_alpha;
  double x,y;
  switch(bdComp)
  {
    case 3: //Dirichlet
      x=Pi; y=t; // x==Pi, y==t
      break;      
    case 4: //Dirichlet
      x=Pi*(1-t); y=1; // x==Pi*(1-t), y==1
      break;
    case 5: //Dirichlet
      x=0; y=1-t; // x==0, y==1-t 
      break;
    default: cout << "U1BoundValue_NSE: wrong boundary part number" << endl;
      break;
  }
  value = (-(1/(2*nu)) + 2*(K/2 - alpha/(4*nu*nu))*y)*cos(x);
}
void U2BoundValue_NSE(int bdComp, double t, double &value)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nu= 1/TDatabase::ParamDB->RE_NR;
  const double alpha = TDatabase::ParamDB->StoDa_alpha;
  double x,y;
  switch(bdComp)
  {
  case 3: //Dirichlet
    x=Pi; y=t; // x==Pi, y==t
    break;      
  case 4: //Dirichlet
	  x=Pi*(1-t); y=1; // x==Pi*(1-t), y==1
    break;
  case 5: //Dirichlet
    x=0; y=1-t; // x==0, y==1-t 
    break;
  default: cout << "U1BoundValue_NSE: wrong boundary part number" << endl;
    break;
  }
  value = (-K - y/(2*nu) + (K/2 - alpha/(4*nu*nu))*y*y)*sin(x);
}
void PressureBoundValue_NSE(int bdComp, double t, double &value)
{
  value = 0;
  //cout << "PressureBoundValue_NSE: bdComp " << bdComp << ", t " << t << endl;
}

// ============================================================================
// coefficients
// ============================================================================
void LinCoeffs_NSE(int n_points, double * X, double * Y,
               double **parameters, double **coeffs)
{
  double *coeff, x, y;
  double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu= 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];

    coeff[0] = nu;
    // RHS for exact solution
    coeff[1] = (-0.5 - alpha*y/(2*nu) + K*nu*y)*cos(x); // f1
    coeff[2] = -(alpha*(-2+y*y)+2*nu*(y-K*nu*(-4 + y*y)))*sin(x)/(4*nu); // f2
    coeff[3] = 0; // g
  }
}
void LinCoeffs_Darcy(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = K; // diffusion
    coeffs[i][1] = 0; // convection_x
    coeffs[i][2] = 0; // convection_y
    coeffs[i][3] = 0; // reaction
    coeffs[i][4] = 0; // f_D
  }
}
void LinCoeffs_DarcyMixed(int n_points, double *X, double *Y,
               double **parameters, double **coeffs) // for mixed formulation
{
  const double K_1 = 1.0/TDatabase::ParamDB->SIGMA_PERM;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = K_1;
    // for RHS
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // f_D
  }
}

// ============================================================================
// initial solution for iteration
// ============================================================================
void InitialP_Darcy(double x, double y, double *values)
{
  values[0] = y*(2/Pi-0.5)*6*x*(1-x)+0.5;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  //values[0] = (2/Pi)*cos(x*Pi/2)*cos(y*Pi/2) - y*(x-0.5);
  //values[1] = -sin(x*Pi/2)*cos(y*Pi/2) - y;
  //values[2] = -cos(x*Pi/2)*sin(y*Pi/2) - (x-0.5);
  //values[3] = 0;
}

double prescribedDirichletOutflow()
{
  OutPut("\n\n\n WARNING\n\n Function call to prescribedDirichletOutflow()\n");
  return 0.0;
}

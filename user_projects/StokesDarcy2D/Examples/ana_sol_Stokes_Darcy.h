// Benchmark for Stokes-Darcy
// domain Omega = (0,1)x(0,2)
// subdomains Omega_Darcy = (0,1)x(0,1)
//            Omega_Stokes= (0,1)x(1,2)
// Interface Gamma = (0,1)x{1}
// Boundary conditions:
//      Dirichlet on Stokes-velocity on \partial Omega_Stokes \setminus Gamma
//      Dirichlet on Darcy-pressure on (0,1)x{0} 
//      Neumann on {0,1}x(0,1)

void ExampleFile()
{
  OutPut("Example: ana_sol_Stokes_Darcy.3.h" << endl);
  //TDatabase::ParamDB->StoDa_alpha = 1; //alpha
  //TDatabase::ParamDB->P1 = 1; //a1
  //TDatabase::ParamDB->P2 = 0; //a2
  //TDatabase::ParamDB->P3 = 1; //a3
  //TDatabase::ParamDB->P4 = 0; //a4
  //TDatabase::ParamDB->P5 = 0; //a5
  
  if(TDatabase::ParamDB->StoDa_interfaceType != 0)
  {
    OutPut("WARNING: This example requires Beaver--Joseph--Saffman interface ");
    OutPut("conditions.\n         The parameter 'StoDa_interfaceType' has ");
    OutPut("been changed accordingly (to be 0)\n");
  }  
  TDatabase::ParamDB->StoDa_interfaceType = 0; // ( Beavers-Joseph-Saffman)
}

// ========================================================================
// exact solution
// ========================================================================
void ExactP_Darcy(double x, double y, double *values)
{
  double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  //double a2 = TDatabase::ParamDB->P2;
  double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  double a5 = TDatabase::ParamDB->P5;
  
  values[0] = 2*nu*a1 + nu*a3*(2*x-1) - 0.5*a1/alpha
             +(1/K)*(1-y)*(a4*x + a1*x*x/(2*alpha*nu) + a5 - a1);
  values[1] = 2*nu*a3 + (1/K)*(1-y)*(a4 + a1*x/(alpha*nu));
  values[2] = -(1/K)*(a4*x + a1*x*x/(2*alpha*nu) + a5 - a1);
  values[3] = (1/K)*(1-y)*a1/(alpha*nu);
}
void ExactU1_Darcy(double x, double y, double *values)
{
  double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  //double a2 = TDatabase::ParamDB->P2;
  double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  //double a5 = TDatabase::ParamDB->P5;
  
  values[0] = -2*K*nu*a3 - (1-y)*(a4 + a1*x/(alpha*nu));
  values[1] = -(1-y)*(a1/(alpha*nu));
  values[2] = a4 + a1*x/(alpha*nu);
  values[3] = 0;
}
void ExactU2_Darcy(double x, double y, double *values)
{
  //double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  //double a2 = TDatabase::ParamDB->P2;
  //double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  double a5 = TDatabase::ParamDB->P5;
  
  values[0] = (a4*x + a1*x*x/(2*alpha*nu) + a5 - a1);
  values[1] = a4 + a1*x/(alpha*nu);
  values[2] = 0;
  values[3] = a1/(alpha*nu);
}
void ExactU1_NSE(double x, double y, double *values)
{
  //double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  double a2 = TDatabase::ParamDB->P2;
  double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  //double a5 = TDatabase::ParamDB->P5;
  
  values[0] = a1*x + (y-1+alpha*nu)*a2 + (y*y-1+2*alpha*nu)*a3 + alpha*nu*a4;
  values[1] = a1;
  values[2] = a2 + 2*y*a3;
  values[3] = 2*a3; //Laplacien
}
void ExactU2_NSE(double x, double y, double *values)
{
  //double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  //double a2 = TDatabase::ParamDB->P2;
  //double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  double a5 = TDatabase::ParamDB->P5;

  values[0] = a4*x + a1*x*x/(2*alpha*nu) + a5 - a1*y;
  values[1] = a4 + a1*x/(alpha*nu);
  values[2] = -a1;
  values[3] = a1/(alpha*nu); //Laplacien
}
void ExactP_NSE(double x, double y, double *values)
{
  //double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  
  double a1 = TDatabase::ParamDB->P1;
  //double a2 = TDatabase::ParamDB->P2;
  double a3 = TDatabase::ParamDB->P3;
  //double a4 = TDatabase::ParamDB->P4;
  //double a5 = TDatabase::ParamDB->P5;
  
  values[0] = nu*a3*(2*x-1) + (a1/alpha)*(y-1.5);
  values[1] = 2*nu*a3;
  values[2] = a1/alpha;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BC_Darcy(int bdComp, double t, BoundCond &cond)  
{
  cond = DIRICHLET;
  if(bdComp == 0) cond = NEUMANN;
  if(bdComp == 3 || bdComp == 4 || bdComp == 5)
    cout << "BC_Darcy: wrong boundary part number" << endl;
}
void BC_Velocity_Darcy(int bdComp, double t, BoundCond &cond)
{
  cond = (bdComp==1) ? NEUMANN : DIRICHLET;
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
  //cond = DIRICHLET;
}
void BC_Pressure_Darcy(int bdComp, double t, BoundCond &cond)
{
  cond = NEUMANN;
}
void BoundCondition_NSE(int bdComp, double t, BoundCond &cond)  
{
  cond = DIRICHLET;
  if(bdComp == 3)  cond = NEUMANN;
  if(bdComp==0 || bdComp==1 || bdComp==2)
    cout << "BoundCondition_NSE: wrong boundary part number";
}
void BoundConditionPressure_NSE(int bdComp, double t, BoundCond &cond)  
{
  cond = NEUMANN;
  if(bdComp==0 || bdComp==1 || bdComp==2)
    cout << "BoundConditionPressure_NSE: wrong boundary part number";
}

// ========================================================================
// boundary values
// ========================================================================
void BoundValue_Darcy(int bdComp, double t, double &value)
{
  double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  //double a2 = TDatabase::ParamDB->P2;
  double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  double a5 = TDatabase::ParamDB->P5;

  switch(bdComp)
  {
    case 0: // x==0, y==1-t
      //value = 2*nu*a1-nu*a3 - 0.5*a1/alpha + (1/K)*t*(a5-a1);// DIRICHLET
      value = -2*K*nu*a3 - t*a4;// NEUMANN
      break;      
    case 1: //y==0 x==t
      value = 2*nu*a1 + nu*a3*(2*t-1) - 0.5*a1/alpha
              +(1/K)*(a4*t + a1*t*t/(2*alpha*nu) + a5 - a1); // DIRICHLET
      //value = (1/K)*(a4*t + a1*t*t/(2*alpha*nu) + a5 - a1); // Neumann
      break;
    case 2: // x==1, y==t
      value = 2*nu*a1 + nu*a3 - 0.5*a1/alpha
             +(1/K)*(1-t)*(a4 + a1/(2*alpha*nu) + a5 - a1); // DIRICHLET
      //value = 2*K*nu*a3 + (1-t)*(a4+a1/(alpha*nu)); // NEUMANN
      break;
    default: cout << "BoundValue_Darcy: wrong boundary part number" << endl;
      break;
  }
}
void BoundValueVelocity_Darcy(int bdComp, double t, double &value)
{
  double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  //double a2 = TDatabase::ParamDB->P2;
  double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  double a5 = TDatabase::ParamDB->P5;
  switch(bdComp)
  {
    case 0: // DIRICHLET
      value = 2*nu*K*a3 + t*a4;
      break;      
    case 1: // NEUMANN
      value = 2*nu*a1 + nu*a3*(2*t-1) - 0.5*a1/alpha
              +(1/K)*(a4*t + a1*t*t/(2*alpha*nu) + a5 - a1);
      //value = -(a4*t + a1*t*t/(2*alpha*nu) + a5 - a1); // DIRICHLET
      break;
    case 2: // DIRICHLET
      value = -2*nu*K*a3 - (1-t)*(a4 + a1/(alpha*nu));
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
  //double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1.0/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  double a2 = TDatabase::ParamDB->P2;
  double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  //double a5 = TDatabase::ParamDB->P5;
  
  switch(bdComp)
  {
    case 3: // x==1, y==1+t
      // DIRICHLET
      //value = a1 + (t+alpha*nu)*a2 + ((1+t)*(1+t)-1+2*alpha*nu)*a3 + alpha*nu*a4;
      value = 2*a1*nu - a3*nu - a1*(-0.5 + t)/alpha; // NEUMANN
      break;
    case 4: // x==1-t, y==2
      // DIRICHLET
      value = a1*(1-t) + (1+alpha*nu)*a2 + (3+2*alpha*nu)*a3 + alpha*nu*a4;
      //value = nu*(a2 + 4*a3 + a4 + (a1 * (1-t)) / (alpha * nu)); // NEUMANN
      break;
    case 5: // x==0, y==2-t
      // DIRICHLET
      value = (1-t+alpha*nu)*a2 + ((2-t)*(2-t)-1+2*alpha*nu)*a3 + alpha*nu*a4;
      //value = -2*a1*nu - a3*nu + a1*(0.5 - t)/alpha; // NEUMANN
      break;
    default: cout << "U1BoundValue_NSE: wrong boundary part number" << endl;
      break;
  }
}
void U2BoundValue_NSE(int bdComp, double t, double &value)
{
  //double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  double a2 = TDatabase::ParamDB->P2;
  double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  double a5 = TDatabase::ParamDB->P5;

  switch(bdComp)
  {
    case 3: // x==1, y==1+t
      //value = a4 + a1/(2*alpha*nu) + a5 - a1*(1+t); // DIRICHLET
      value = nu*(a2 + a4 + a1/(alpha*nu) + 2*a3*(1+t)); // NEUMANN
      break;      
    case 4: // x==1-t, y==2
      value = a4*(1-t) + a1*(1-t)*(1-t)/(2*alpha*nu) + a5 - 2*a1; // DIRICHLET
      //value = -0.5 * a1 / alpha - 2*a1*nu - a3 * nu * (2*(1-t) - 1); // NEUMANN
      break;
    case 5: // x==0, y==2-t
      value = a5 - a1*(2-t); // DIRICHLET
      //value = -nu*(a2 + a4 + 2*a3*(2-t)); // NEUMANN
      break;
    default: cout << "U2BoundValue_NSE: wrong boundary part number" << endl;
      break;
  }
}
void PressureBoundValue_NSE(int bdComp, double t, double &value)
{
  value = 0;
  //cout << "PressureBoundValue_NSE: bdComp " << bdComp << ", t " << t << endl;
}

// ========================================================================
// coefficients
// ========================================================================
void LinCoeffs_NSE(int n_points, double * X, double * Y,
               double **parameters, double **coeffs)
{
  double nu = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = nu;
    // RHS for exact solution
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // g
  }
}
void LinCoeffs_Darcy(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double K = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1    = TDatabase::ParamDB->P1;
  
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = K;// diffusion
    coeffs[i][1] = 0;    // convection_x
    coeffs[i][2] = 0;    // convection_y
    coeffs[i][3] = 0;    // reaction
    coeffs[i][4] = -(1-Y[i])*a1/(alpha*nu);    // f_D
  }
}
void LinCoeffs_DarcyMixed(int n_points, double *X, double *Y,
               double **parameters, double **coeffs) // for mixed formulation
{
  const double K_1 = 1./TDatabase::ParamDB->SIGMA_PERM;
  const double nu    = 1/TDatabase::ParamDB->RE_NR;
  const double alpha = TDatabase::ParamDB->StoDa_alpha;
  const double a1    = TDatabase::ParamDB->P1;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = K_1;
    // for RHS
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = -(1-Y[i])*(a1/(alpha*nu)); // g = div u
  }
}


double prescribedDirichletOutflow()
{
  double sigma = TDatabase::ParamDB->SIGMA_PERM;
  double nu    = 1/TDatabase::ParamDB->RE_NR;
  double alpha = TDatabase::ParamDB->StoDa_alpha;
  double a1 = TDatabase::ParamDB->P1;
  //double a2 = TDatabase::ParamDB->P2;
  //double a3 = TDatabase::ParamDB->P3;
  double a4 = TDatabase::ParamDB->P4;
  double a5 = TDatabase::ParamDB->P5;
  return -sigma*(0.5*a4 + a1/(6*alpha*nu) + a5 - a1);
}

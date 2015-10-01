// Benchmark for Stokes-Darcy
// domain Omega = (0,1)x(0,2)
// subdomains Omega_Darcy = (0,1)x(0,1)
//            Omega_Stokes= (0,1)x(1,2)
// Interface Gamma = (0,1)x{1}
// Boundary conditions:
//      Dirichlet on Stokes-velocity on {0}x(1,2) and (0,1)x{2}
//      Neumann for Stokes on {1}x(1,2)
//      Dirichlet on Darcy-pressure on (0,1)x{0} 
//      Neumann on {0,1}x(0,1)

void ExampleFile()
{
  OutPut("Example: Polynomial_ut0.h\n");
  
  // set the boundary condition u_s.t = 0 on the interface
  TDatabase::ParamDB->StoDa_interfaceType = 1;
  // pressure projection into L^2_0 is set for Dirichlet problems only
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0; 
}

// ============================================================================
// exact solution
// ============================================================================
void ExactP_Darcy(double x, double y, double *values)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nu= 1./TDatabase::ParamDB->RE_NR;
  values[0] = (x*(1 - x)*(y - 1) + y*y*y/3 - y*y + y)/K + 2*nu*x;
  values[1] = ((1 - 2*x)*(y - 1))/K + 2*nu;
  values[2] = (x*(1 - x) + (y-1)*(y-1))/K;
  values[3] = 0;
}
void ExactU1_Darcy(double x, double y, double *values)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nu= 1./TDatabase::ParamDB->RE_NR;
  values[0] = -(1 - 2*x)*(y - 1) - 2*nu*K;
  values[1] = 2*(y-1);
  values[2] = 2*x-1;
  values[3] = 0;
}
void ExactU2_Darcy(double x, double y, double *values)
{
  //const double K = TDatabase::ParamDB->SIGMA_PERM;
  //const double nu= 1/TDatabase::ParamDB->RE_NR;
  values[0] = -(x*(1 - x) + (y-1)*(y-1));
  values[1] = 2*x-1;
  values[2] = -2*(y-1);
  values[3] = 0;
}
void ExactU1_NSE(double x, double y, double *values)
{
  values[0] = y*y-2*y+1;
  values[1] = 0;
  values[2] = 2*y-2;
  values[3] = 2; //Laplacien
}
void ExactU2_NSE(double x, double y, double *values)
{
  values[0] = x*x-x;
  values[1] = 2*x-1;
  values[2] = 0;
  values[3] = 2; //Laplacien
}
void ExactP_NSE(double x, double y, double *values)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nu= 1/TDatabase::ParamDB->RE_NR;
  values[0] = 2*nu*(x+y-1) + 1/(3*K);
  values[1] = 2*nu;
  values[2] = 2*nu;
  values[3] = 0;
}

// ============================================================================
// boundary conditions
// ============================================================================
void BC_Darcy(int bdComp, double t, BoundCond &cond)  
{
  cond = (bdComp == 1) ? DIRICHLET : NEUMANN;
  if(bdComp != 0 && bdComp != 1 && bdComp != 2)
    cout << "BC_Darcy: wrong boundary part number " << bdComp << endl;
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
{
  cond = NEUMANN;
}
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
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nu= 1/TDatabase::ParamDB->RE_NR;
  switch(bdComp)
  {
    case 0: // x==0, y==1-t
      value = t-2*nu*K; // NEUMANN
      //value = ((1-t) * (1-t) * (1-t)/3 - (1-t) * (1-t) + (1-t))/K; // DIRICHLET
      break;      
    case 1: //y==0 x==t
      value = -t*(1-t)/K + 2*t*nu; // DIRICHLET
      break;
    case 2: // x==1, y==t
      value = (1-t)+2*nu*K; // NEUMANN
      //value = (t*t*t/3 - t*t + t)/K + 2*nu; // DIRICHLET
      break;
    default: cout << "BoundValue_Darcy: wrong boundary part number" << endl;
      break;
  }
}
void BoundValueVelocity_Darcy(int bdComp, double t, double &value)
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const double nu= 1.0/TDatabase::ParamDB->RE_NR;
  switch(bdComp)
  {
    case 0: // x==0, y==1-t
      value = -t + 2*nu*K; // Dirichlet
      break;      
    case 1: //y==0 x==t
      value = -t*(1-t)/K + 2*nu*t; // Neumann
      break;
    case 2: // x==1, y==t
      value = -1+t - 2*nu*K; // Dirichlet
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
  //const double nu = 1.0/TDatabase::ParamDB->RE_NR;
  //const double K = TDatabase::ParamDB->SIGMA_PERM;
  switch(bdComp)
  {
    case 3: // x==1, y==1+t
      value = (1+t)*(1+t) - 2*(1+t) + 1; // Dirichlet
      //value = -1/(3*K) - 2*nu*(1+t); // Neumann
      break;      
    case 4: // x==1-t, y==2
      value = 1; // Dirichlet
      //value = nu*(1+2*(1-t)); // Neumann
      break;
    case 5: // x==0, y==2-t
      value = (2-t)*(2-t) - 2*(2-t) + 1; // Dirichlet
      //value = 1/(3*K) + 2*nu*(1-t); // Neumann
      break;
    default: cout << "U1BoundValue_NSE: wrong boundary part number" << endl;
      break;
  }
}
void U2BoundValue_NSE(int bdComp, double t, double &value)
{
  //const double nu = 1.0/TDatabase::ParamDB->RE_NR;
  //const double K  = TDatabase::ParamDB->SIGMA_PERM;
  switch(bdComp)
  {
    case 3: // x==1, y==1+t
      value = 0; // Dirichlet
      //value = nu*(-1 + 2*(1+t)); // Neumann
      break;
    case 4: // x==1-t, y==2
      value = (1-t)*(1-t) - (1-t); // Dirichlet
      //value = -1/(3*K) - 2*nu*(1+(1-t)); // Neuamann
      break;
    case 5: // x==0, y==2-t
      value = 0; // Dirichlet 
      //value = -nu*(-3 + 2*(2-t)); // Neumann
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

// ============================================================================
// coefficients
// ============================================================================
void LinCoeffs_NSE(int n_points, double * X, double * Y,
               double **parameters, double **coeffs)
{
  const double nu= 1.0/TDatabase::ParamDB->RE_NR;
  bool nonlin = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE!=0;
  for(int i=0;i<n_points;i++)
  {
    const double ugradu_x = nonlin ? (2*Y[i]-2)*(X[i]*X[i]-X[i]) : 0.0;
    const double ugradu_y = nonlin ? (2*X[i]-1)*(Y[i]*Y[i]-2*Y[i]+1) : 0.0;
    coeffs[i][0] = nu;
    // RHS for exact solution
    coeffs[i][1] = ugradu_x; // f1
    coeffs[i][2] = ugradu_y; // f2
    coeffs[i][3] = 0; // g
  }
}
void LinCoeffs_Darcy(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double K = TDatabase::ParamDB->SIGMA_PERM;
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
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = 1.0/K;
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

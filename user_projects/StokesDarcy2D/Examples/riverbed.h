/*
 Benchmark for Stokes-Darcy
 domain Omega = (0,10)x(0,5)
 subdomains Omega_Darcy = (0,10)x(0,phi)
            Omega_Stokes= (0,10)x(phi,5)
 Interface Gamma = {0,10}x{phi(t),t in [0,10]}
 phi is some function with phi(0)=2=phi(10) and 0<phi<5
 So the interface starts at (0,2) following some function through the domain
 and ending at (10,2)
 Boundary conditions:
      Dirichlet on Stokes-velocity on {0}x(2,5) and (0,10)x{5}
      Neumann (homogeneous on Stokes-stress on {10}x(2,5)
      Dirichlet on Darcy-pressure on (0,10)x{0} 
      Neumann on {0,20}x(0,2)
*/

void ExampleFile()
{
  OutPut("Example: riverbed_flow_Stokes_Darcy.h" << endl) ;
  TDatabase::ParamDB->StoDa_periodicBoundary = 1;
  //TDatabase::ParamDB->StoDa_periodicBoundaryPressureDrop = 0.001;
}

// ========================================================================
// exact solution unknown
// ========================================================================
void ExactP_Darcy(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
void ExactU1_Darcy(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
void ExactU2_Darcy(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
void ExactU1_NSE(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
void ExactU2_NSE(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
void ExactP_NSE(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BC_Darcy(int bdComp, double t, BoundCond &cond)  
{
  cond = (bdComp==0 || bdComp==2)? NEUMANN : NEUMANN;
  if (bdComp==3 || bdComp==4 || bdComp==5)
    OutPut("BC_Darcy: wrong boundary part number " << bdComp << endl);
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
  //cond = DIRICHLET;
  cond = (bdComp==3 || bdComp==5) ? NEUMANN : DIRICHLET;
  if(bdComp==0 || bdComp==1 || bdComp==2)
    OutPut("BoundCondition_NSE: wrong boundary part number " << bdComp <<endl);
}
void BoundConditionPressure_NSE(int bdComp, double t, BoundCond &cond)  
{
  
  cond = (bdComp==3 || bdComp==5) ? ROBIN : NEUMANN;
  if(bdComp==0 || bdComp==1 || bdComp==2)
    OutPut("BoundConditionPressure_NSE: wrong boundary part number " << bdComp 
           << endl);
}

// ========================================================================
// boundary values
// ========================================================================
void BoundValue_Darcy(int bdComp, double t, double &value)
{
  //static const double K = TDatabase::ParamDB->SIGMA_PERM;
  switch(bdComp)
  {
    case 0: // x==0, y==1-t
      value = 0; // NEUMANN
      //value = 0; // Dirichlet
      break;
    case 1: //y==0 x==t 
      value = 0; // NEUMANN
      //value = (1-t)*TDatabase::ParamDB->StoDa_periodicBoundaryPressureDrop;
      break;
    case 2: // x==1, y==t
      value = 0; // Neumann
      //value = 0; // DIRICHLET
      break;
    default: 
      OutPut("BoundValue_Darcy: wrong boundary part number " << bdComp <<endl);
      break;
  }
}
void BoundValueVelocity_Darcy(int bdComp, double t, double &value)
{
  switch(bdComp)
    {
      case 0: // DIRICHLET
        value = 0;
        break;      
      case 1: // NEUMANN
        value = 0; 
        break;
      case 2: // DIRICHLET
        value = 0;
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
  switch(bdComp)
  {
    case 3:
      //value = 4*t*(1-t); // Dirichlet
      value = 0.0; // Neumann
      break;      
    case 4:  
      value = 0;
      break;
    case 5: 
      //value = 4*t*(1-t); // Dirichlet
      value = TDatabase::ParamDB->StoDa_periodicBoundaryPressureDrop; // Neumann
      break;
    default: 
      OutPut("BoundValue_Darcy: wrong boundary part number " << bdComp <<endl);
      break;
  }
}
void U2BoundValue_NSE(int bdComp, double t, double &value)
{
  switch(bdComp)
  {
    case 3: 
      value = 0;
      break;      
    case 4:  
      value = 0;
      break;
    case 5: 
      value = 0;
      break;
    default: 
      OutPut("BoundValue_Darcy: wrong boundary part number " << bdComp <<endl);
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
  double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    // RHS for exact solution
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // g
  }
}
void LinCoeffs_Darcy(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double sigma = TDatabase::ParamDB->SIGMA_PERM;

  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = sigma;// diffusion
    coeffs[i][1] = 0;    // convection_x
    coeffs[i][2] = 0;    // convection_y
    coeffs[i][3] = 0;    // reaction
    coeffs[i][4] = 0; // f_D
  }
}
void LinCoeffs_DarcyMixed(int n_points, double *X, double *Y,
               double **parameters, double **coeffs) // for mixed formulation
{
  double sigma = TDatabase::ParamDB->SIGMA_PERM;
  //double nu    = 1/TDatabase::ParamDB->RE_NR;
  //double alpha = TDatabase::ParamDB->StoDa_alpha;
  
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = 1/sigma;
    // for RHS
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // f_D
  }
}


double prescribedDirichletOutflow()
{
  OutPut("\n\n\n WARNING\n\n Function call to prescribedDirichletOutflow()\n");
  return 0.0;
}



#include <math.h> //pow

#include <array>

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

//side effect: sets the global parameter
void ExampleFile()
{
  Output::root_info<1>("Example: CylinderSquare22000.h");
  Output::print("Example: CylinderSquare22000.h ", TDatabase::ParamDB->P7, 
         " % noise (only for U1 !!!, see [JK05])");
  if (!TDatabase::ParamDB->FRICTION_TYPE)
  {
    Output::print("no slip b.c.");
  }
  else
     OutPut("slip+friction b.c." << endl);
   TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 22000; 
   /*if (TDatabase::ParamDB->CELL_MEASURE!=2)
   {
       TDatabase::ParamDB->CELL_MEASURE = 2;
       OutPut("CELL_MEASURE changed to " << 
       TDatabase::ParamDB->CELL_MEASURE << endl);
   }*/
   if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 100)
   {
      TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE = 101;
      Output::print("TURBULENT_VISCOSITY_TYPE changed to ", 
		    TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE);
   }
   Output::print("CYLINDER_22000_YPLUS_SIDES: ", TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES); 
   Output::print("CYLINDER_22000_YPLUS_FRONT: ", TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT); 
   Output::print("CYLINDER_22000_YPLUS_BACK:  ", TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK); 
}

// ========================================================================
// initial condition
// ========================================================================
void InitialU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
     cond  =  SLIP_FRICTION_PENETRATION_RESISTANCE;
     TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;

    // outflow boundary condition
    if (fabs(x-2.5)<1e-6)
    {
	cond = NEUMANN;
    } 
    // inflow boundary condition
    if (fabs(x)<1e-6)
      cond  = DIRICHLET;
    
    // boundary conditions at the column
    if ((fabs(x-0.45)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
    {
      cond  = DIRICHLET;
      //cout << "left ";
    }
    if ((fabs(x-0.55)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
    {
      cond  = DIRICHLET;
      //cout << "right ";
    }
    if ((fabs(y-0.65)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
    {
      cond  = DIRICHLET;
      //cout << "lower ";
    }
    if ((fabs(y-0.75)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
    {
      cond  = DIRICHLET;
      //cout << "upper ";
    }
  
  /*  // no slip b.c. on the bottom, top, column
    if (TDatabase::ParamDB->FRICTION_TYPE==0)
    {
	TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY = 3;
	// bottom and top 
	   if (fabs(z)<1e-6)
	   cond  = DIRICHLET;
	   if (fabs(z-0.4)<1e-6)
	   cond  = DIRICHLET;
	// boundary conditions at the column
	if ((fabs(x-0.45)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
	{
	    cond  = DIRICHLET;
	    //cout << "left ";
	}
	if ((fabs(x-0.55)<1e-6) && (y >= 0.65-1e-6) && (y<=0.75+1e-6))
	{
	    cond  = DIRICHLET;
	    //cout << "right ";
	}
	if ((fabs(y-0.65)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
	{
	    cond  = DIRICHLET;
	    //cout << "lower ";
	}
	if ((fabs(y-0.75)<1e-6) && (x >= 0.45-1e-6) && (x<=0.55+1e-6))
	{
	    cond  = DIRICHLET;
	    //cout << "upper ";
	}
        // this is for the circular cylinder
	if (fabs((x-0.5)*(x-0.5) +(y-0.7)*(y-0.7)-2.5e-3)<=1e-6)
	    cond  = DIRICHLET;
    }
    else
    {
	// slip with friction boundary conditions on all walls
	TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY = 3;
    }*/
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  double noise = TDatabase::ParamDB->P7/100.0;
  // double eps = 1e-4;
  //if( (fabs(x)<1e-6) || (fabs(x-2.5)<1e-6) )
  if((fabs(x)<1e-6))
  {
    //if ((fabs(z)>eps)&&(fabs(z-0.4)>eps)&&(fabs(y)>eps)&&(fabs(y-1.4)>eps))
    value = 1 + noise * ((double)rand()/RAND_MAX-0.5);
    //else
    //value = 1;
  }
  else
  {
    value = 0;
  }
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  // double noise = TDatabase::ParamDB->P7/1000.0;
  if((fabs(x)<1e-6))
  {
    //value = noise* ((double)rand()/RAND_MAX-0.5);
    value = 0;
  }
  else
  {
    value = 0;
  }  
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  //double noise = TDatabase::ParamDB->P7/1000.0;
  if((fabs(x)<1e-6))
  {
    //value = noise * ((double)rand()/RAND_MAX-0.5);
    value = 0;
  }
  else
  {
    value = 0;
  }    
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{

  double eps = DIMENSIONLESS_VISCOSITY; // the kinematic viscosity (1e-3 in the paper cited above)

  for(int i=0;i<n_points;i++)
  {
    double* coeff = coeffs[i];
      
    coeff[0] = eps; // kinematic viscosity
    coeff[1] = 0;   // f1
    coeff[2] = 0;   // f2
    coeff[3] = 0;   // f3
  }
}
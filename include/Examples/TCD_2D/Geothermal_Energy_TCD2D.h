// TCD2D problem, 
/* Geothermal Energy
 * author: Laura Blank
 * date: 10/09/2018
 */



// physical parameters
// These should be reset when constructing the Example class
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double viscosity = -1;
double effective_viscosity = -1;
double permeability = -1;

void ExampleFile()
{
  Output::print<1>("Example: Geothermal_Energy_TCD2D.h");
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constexpr bool rhs_depends_on_time = false;
constexpr bool coefficients_depend_on_time = false;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0; 
  values[2] = 0;
  values[3] = 0; 
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
 // if(BdComp==1)
 //   cond = NEUMANN;
 // else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
   
value = 150;

/*
switch(BdComp)
  {
    case 0: value = 420; // u.n=0
            break;
    case 1: value = 120*(1-Param) + 300; //K/mu * (1+exp(1/t)-exp((1-Param)/t) - exp(Param/t)) / (1+exp(1/t)); //p=...
            //value = 0; // TEST:sources/sinks,
            break;
    case 2: value = 300; // u.n=0
            break;
    case 3: value = 120*(Param) + 300; //K/mu * (1+exp(1/t)-exp((Param)/t) - exp((1-Param)/t)) / (1+exp(1/t)); // p=...
            //value = 0; // TEST: sources/sinks,
            break;
    default: cout << "No boundary component with this number." << endl;
             break;
  }
*/

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// initial condition
void InitialCondition(double x, double y, double *values)
{
  //double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = 150;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  const double eps=1/TDatabase::ParamDB->PE_NR;
  //double t = TDatabase::TimeDB->CURRENTTIME;

  double T_in = 50.; 
  double a = 0.05;
  std::vector<double> center_source = {4.5, 3}; 

  double exact[4];

  for(int i = 0; i < n_points; i++)
  {
    Exact(x[i], y[i], exact);

    coeffs[i][0] = eps;  //diffusivity
    coeffs[i][1] = parameters[i][0]; //0;  //b1
    coeffs[i][2] = parameters[i][1]; //0;  //b2
    coeffs[i][3] = 0;  //reaction coefficient

   coeffs[i][4] = -coeffs[i][0]*exact[3]; // diffusion
   coeffs[i][4] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2]; // convection
   coeffs[i][4] += coeffs[i][3]*exact[0]; // reaction

   // coeffs[i][4] = 0;  //f

if ( (((x[i]-center_source[0])/a) < 1 && ((x[i]-center_source[0])/a) > -1) && 
     (((y[i]-center_source[1])/a) < 1 && ((y[i]-center_source[1])/a) > -1)    )
{
       coeffs[i][3]  +=   1/(a*a)* ( (cos(Pi*(x[i]-center_source[0])/a)+1)/2 )* ( (cos(Pi*(y[i]-center_source[1])/a)+1)/2 ); // source
       coeffs[i][4]  +=   T_in * 1/(a*a)* ( (cos(Pi*(x[i]-center_source[0])/a)+1)/2 )* ( (cos(Pi*(y[i]-center_source[1])/a)+1)/2 ); // source
}

coeffs[i][5] = 0;

  }
}


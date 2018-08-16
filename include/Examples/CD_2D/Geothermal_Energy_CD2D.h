// CD2D problem, 
/* Geothermal Energy
 */

// author: Laura Blank

/*
   Domain =  [0,1]x[0,1]
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
  Output::print<1>("Example: Geothermal_Energy_CD2D.h");
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// exact solution
void Exact(double x, double y, double *values)
{
  const double p = Pi; // 2*Pi;
  values[0] = 0; //sin(p*x)*sin(p*y);
  values[1] = 0; //p*cos(p*x)*sin(p*y);
  values[2] = 0; //p*sin(p*x)*cos(p*y);
  values[3] = 0; //-2*p*p*sin(p*x)*sin(p*y);
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
  //if(BdComp==1) // Neumann, x = 1, y = Param
  //{
  //  const double eps = 1. / TDatabase::ParamDB->PE_NR;
  //  double exact[4];
  //  Exact(1., Param, exact);
  //  value = eps * exact[1];
  //  //value = -eps*Pi*sin(Pi*Param);
  //}
  //else // Dirichlet
   
// value = 40;

switch(BdComp)
  {
    case 0: value = 200; // u.n=0
            break;
    case 1: value = 160*(1-Param) + 40; //K/mu * (1+exp(1/t)-exp((1-Param)/t) - exp(Param/t)) / (1+exp(1/t)); //p=...
            //value = 0; // TEST:sources/sinks,
            break;
    case 2: value = 40; // u.n=0
            break;
    case 3: value = 160*Param + 40; //K/mu * (1+exp(1/t)-exp((Param)/t) - exp((1-Param)/t)) / (1+exp(1/t)); // p=...
            //value = 0; // TEST: sources/sinks,
            break;
    default: cout << "No boundary component with this number." << endl;
             break;
  }


}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  const double eps=1/TDatabase::ParamDB->PE_NR;
  double exact[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;  //diffusivity
    coeffs[i][1] = parameters[i][0]; //0;  //b1
    coeffs[i][2] = parameters[i][1]; //0;  //b2
    coeffs[i][3] = 0;  //reaction coefficient
    
    Exact(x[i], y[i], exact);

 //   coeffs[i][4] = -coeffs[i][0]*exact[3]; // diffusion
 //   coeffs[i][4] += coeffs[i][1]*exact[1] + coeffs[i][2]*exact[2]; // convection
 //   coeffs[i][4] += coeffs[i][3]*exact[0]; // reaction

 //  coeffs[i][4] = 0;  //f
 coeffs[i][4] = -( 1/sqrt((Pi*0.01*0.01) )* exp(-(((x[i]-0.4)*(x[i]-0.4))+((y[i]-0.4)*(y[i]-0.4)))/(0.01*0.01)))*10000;
// (1/sqrt((Pi*0.01*0.01) )* exp(-(((x[i]-0.6)*(x[i]-0.6))+((y[i]-0.6)*(y[i]-0.6)))/(0.01*0.01)) )*10000;


  }
}


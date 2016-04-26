// Navier-Stokes problem, Driven cavity
//
// u(x,y) = unknown
// p(x,y) = unknown

#define __UREA__
#define __UREA__PIPE__

#include <Urea_3d4d.h>

void ExampleFile()
{
  OutPut("Example: Urea_pipe_nse_mult_twist.h with " << TDatabase::ParamDB->N_CELL_LAYERS  << " cell layers in flow direction"<< endl);

 /* if ((TDatabase::ParamDB->N_CELL_LAYERS != 21) && (TDatabase::ParamDB->N_CELL_LAYERS != 84))
  {
    OutPut("no initial data for N_CELL_LAYERS: " << TDatabase::ParamDB->N_CELL_LAYERS);
    exit(4711);
  }*/

  TDatabase::ParamDB->DRIFT_Z = 3380;
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1356;
  TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 1;
  TDatabase::ParamDB->UREA_PIPE_RADIUS = 0.3;

  OutPut("UREA_PIPE_RADIUS: " << TDatabase::ParamDB->UREA_PIPE_RADIUS << endl);
  OutPut("UREA_REACTION_DISC: " << TDatabase::ParamDB->UREA_REACTION_DISC << endl);
  OutPut("UREA_PB_DISC: " << TDatabase::ParamDB->UREA_PB_DISC << endl);
  OutPut("UREA_PB_DISC_STAB: " << TDatabase::ParamDB->UREA_PB_DISC_STAB<<endl);
  OutPut("UREA_SOLD_PARAMETER_TYPE: "<< TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE <<endl);
  OutPut("UREA_MODEL: " << TDatabase::ParamDB->UREA_MODEL << endl);
  OutPut("UREA_CONC_TOL: " << TDatabase::ParamDB->UREA_CONC_TOL << endl);
  OutPut("UREA_CONC_MAXIT: " << TDatabase::ParamDB->UREA_CONC_MAXIT << endl);

  OutPut("UREA_AGGR_BROWNIAN: " << TDatabase::ParamDB->UREA_AGGR_BROWNIAN<< endl);
  OutPut("UREA_AGGR_POL_ORDER: " << TDatabase::ParamDB->UREA_AGGR_POL_ORDER << endl);
  OutPut("UREA_AGGR_BROWNIAN_TEMP: " << TDatabase::ParamDB->UREA_AGGR_BROWNIAN<< endl);
  OutPut("UREA_AGGR_BROWNIAN_SCAL: " << TDatabase::ParamDB->UREA_AGGR_BROWNIAN_SCAL<< endl);
  OutPut("UREA_AGGR_SHEAR_FACTOR_TYPE: " << TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR_TYPE << endl);
  OutPut("UREA_AGGR_SHEAR_FACTOR: " << TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR << endl);

  OutPut("UREA_l_infty: " << TDatabase::ParamDB->UREA_l_infty <<endl);
  OutPut("UREA_u_infty: " << TDatabase::ParamDB->UREA_u_infty <<endl);
  OutPut("UREA_c_infty: " << TDatabase::ParamDB->UREA_c_infty <<endl);
  OutPut("UREA_temp_infty: " << TDatabase::ParamDB->UREA_temp_infty <<endl);
  OutPut("UREA_f_infty: " << TDatabase::ParamDB->UREA_f_infty<<endl);
  OutPut("UREA_nu: " << TDatabase::ParamDB->UREA_nu<<endl);
  OutPut("UREA_rho: " << TDatabase::ParamDB->UREA_rho<<endl);
  OutPut("UREA_c_p: " << TDatabase::ParamDB->UREA_c_p<<endl);
  OutPut("UREA_lambda: " << TDatabase::ParamDB->UREA_lambda<<endl);
  OutPut("UREA_D_P_0: " << TDatabase::ParamDB->UREA_D_P_0<<endl);
  OutPut("UREA_D_P_MAX: " << TDatabase::ParamDB->UREA_D_P_MAX <<endl);
  OutPut("UREA_k_v: " << TDatabase::ParamDB->UREA_k_v<<endl);
  OutPut("UREA_m_mol: " << TDatabase::ParamDB->UREA_m_mol<<endl);
  OutPut("UREA_D_J: " << TDatabase::ParamDB->UREA_D_J<<endl);
  OutPut("UREA_rho_d: " << TDatabase::ParamDB->UREA_rho_d <<endl);
  OutPut("UREA_k_g: " << TDatabase::ParamDB->UREA_k_g<<endl);
  OutPut("UREA_g: " << TDatabase::ParamDB->UREA_g <<endl);
  OutPut("UREA_rho_sat_1: " << TDatabase::ParamDB->UREA_rho_sat_1 <<endl);
  OutPut("UREA_rho_sat_2: " << TDatabase::ParamDB->UREA_rho_sat_2<<endl);
  OutPut("UREA_beta_nuc: " << TDatabase::ParamDB->UREA_beta_nuc<<endl);
  OutPut("UREA_alfa_nuc: " << TDatabase::ParamDB->UREA_alfa_nuc<<endl);
  OutPut("UREA_INFLOW_SCALE: " << TDatabase::ParamDB->UREA_INFLOW_SCALE <<endl);
  OutPut("UREA_inflow_time: " << TDatabase::ParamDB->UREA_inflow_time <<endl);
  OutPut("PB_DISC_TYPE: " << TDatabase::ParamDB->PB_DISC_TYPE <<endl);
  OutPut("PB_TIME_DISC: " << TDatabase::ParamDB->PB_TIME_DISC <<endl);
  OutPut("INTERNAL_STEADY_STATE_MATRICES_OR_RHS: " << TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS <<endl);

  if(!TDatabase::ParamDB->UREA_PIPE)
  {
    TDatabase::ParamDB->UREA_PIPE = 1;
    OutPut("Setting UREA_PIPE to " << TDatabase::ParamDB->UREA_PIPE  << endl);
  }
  else
    OutPut("UREA_PIPE: " << TDatabase::ParamDB->UREA_PIPE << endl);

  // set some parameters
  //TDatabase::ParamDB->GRID_TYPE = 3;
  //OutPut("GRID_TYPE set to " << TDatabase::ParamDB->GRID_TYPE << endl);
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  // Poiseuille flow with R = 1, dp/dz = -0.1
  // v = -1/(4 eta) (R^2-r^2) * nabla p
  // ParamDB->UREA_INFLOW_SCALE  in ml/min !!!
  // based on derivation 13/12/06
  double R, R2, r2yz;
  double l_infty = TDatabase::ParamDB->UREA_l_infty;
 
  R = 0.3;
  R2=R*R;

  //if (y > 0)
  //  values[0] =
  //    (R2-(y-R_bent)*(y-R_bent)-z*z) * V_in /(3000 * Pi * R2*R2)/TDatabase::ParamDB->UREA_u_infty;
  //else
  //  values[0] =
  //    -(R2-(y+R_bent)*(y+R_bent)-z*z) * V_in /(3000 * Pi * R2*R2)/TDatabase::ParamDB->UREA_u_infty;
  double vFlux = 7.2;
  r2yz = (y-5.65)*(y-5.65)+z*z; 

  // factor covering inflow from boundary, dimensionless representation, etc.
  double uFactor = 2. * vFlux * TDatabase::ParamDB->UREA_INFLOW_SCALE
    / (Pi * R2 * TDatabase::ParamDB->UREA_u_infty);
  values[0] = uFactor * (1-r2yz/R2);
  //if (TDatabase::TimeDB->CURRENTTIME < 1.0)
  //  values[0] *= TDatabase::TimeDB->CURRENTTIME; 
  
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  
  //OutPut(values[0] << " ");
 // if (abs(y-5.65)>0.3)
  //  exit(1);
  //OutPut(" in " << x << " "  << y << " " << z << " " << values[0] << " : ");
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


// ========================================================================
// initial solution
// ========================================================================

void InitialU1(double x, double y, double z, double *values)
{
  // Poiseuille flow with R = 1, dp/dz = -0.1
  // v = -1/(4 eta) (R^2-r^2) * nabla p
  double val[5];
  
  ExactU1(x,y,z,val);
  //if ((x<50)&& (y > 0))
  if (x<50)
  {
    values[0] = val[0];
    //OutPut(y << " " << val[0] << "::");
  }
  else
    values[0] = 0;

 values[0] = 0.0;
 //OutPut(values[0] << " ");
}


void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}


void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
}


void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}


// ========================================================================
// boundary conditions
// ========================================================================

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  double eps = 1e-6, r2yz;

  cond = DIRICHLET;

  r2yz =  (y-5.65)*(y-5.65)+(z-5.3333333333333325)*(z-5.3333333333333325);

  if (x>53.75-25)
  {
     cond = NEUMANN;
     TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
     OutPut("neum " << x <<" " << y << " " << z << endl);
  }  
  //if ((fabs(x-43.43137254901961)< eps))//&& ((0.09-r2yz) < 1e-3))
  //if (x>53.75)
  //  OutPut("bdry " <<  setprecision(16) << x <<" " << y << " " <<  setprecision(16) << z << " "  << r2yz << endl);
  //if ((fabs(60-x)<eps))
  //if (z>5.03333)
 // {
    // outflow
   // cond = NEUMANN;
    //OutPut("neum " << x <<" " << y << " " << z << endl);
   // TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
 // }
}


// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  double eps = 1e-6, val[5];

  value = 0.0;

  // inflow
  //if ((fabs(x)<eps)&&(y>0))
  if ((fabs(x)<eps))
  {
    ExactU1(x,y,z,val);
    value = val[0];
  }
  if (x>53.75-25)
  {
    ExactU1(0,y,z-5.3333333333333325,val);
    value = val[0];
    value = 0;
    OutPut("out " << x << " "  << y << " " << z << " " << value << endl);
  }
  //OutPut(value << " ");
}


// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
double **parameters, double **coeffs)
{
  double eps;
  int i;
  double *coeff;
  // reference values
  double L_infty = TDatabase::ParamDB->UREA_l_infty;
  double U_infty = TDatabase::ParamDB->UREA_u_infty;
  double eta = 1.19e-3;
  double rho = 1100;
  double nu;
  
  nu = eta/rho;

  eps =  nu/(L_infty*U_infty);
  //OutPut(nu <<  " " <<  L_infty << " "  << U_infty << " " << 1.0/eps<< endl);
  //eps = 1/TDatabase::ParamDB->RE_NR;
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0;                                 // f1
    coeff[2] = 0;                                 // f2
    coeff[3] = 0;                                 // f3
  }
}


// ========================================================================
// switch of x- and z- coordinate
// ========================================================================

void SwitchCoords_xz(TCollection *coll)
{
  int N_Cells, N_V, i, j;
  double x, y, z, nx, ny, nz;
  TVertex *vertex;
  TBaseCell *cell;

  N_Cells = coll->GetN_Cells();

  // initialise ClipBoard to see which vertices are still untouched
  for(i=0 ; i<N_Cells ; i++)
  {
    cell = coll->GetCell(i);
    N_V = cell->GetN_Vertices();

    for (j=0 ; j<N_V ; j++)
    {
      vertex = cell->GetVertex(j);
      vertex->SetClipBoard(0);
    }
  }

  // switch coordinates and 'normal' coordinates if ClipBoard=0
  for(i=0 ; i<N_Cells ; i++)
  {
    cell = coll->GetCell(i);
    N_V = cell->GetN_Vertices();

    for (j=0 ; j<N_V ; j++)
    {
      vertex = cell->GetVertex(j);
      if (!vertex->GetClipBoard())
      {
        vertex->GetCoords(x, y, z);
        vertex->GetNormalCoords(nx, ny, nz);
        vertex->SetCoords(z, y, x);
        vertex->SetNormal(nz, ny, nx);
        vertex->SetClipBoard(1);
      }
    }
  }
}


void CoordsTrafo_bentPipe(double x, double y, double &x_bent_curve, double &y_bent_curve)
{
  double l_infty = TDatabase::ParamDB->UREA_l_infty;
  double  pi = 3.1415926535897;
  double pipe_radius = TDatabase::ParamDB->UREA_PIPE_RADIUS;
  double bent_radius = 1;
  double  x_start_curve = 50;
  double  x_end_curve = 160;

  if(x < x_start_curve)
  {
    x_bent_curve=x;
    y_bent_curve=bent_radius+y;
  }
  else
  {
    if(x<=x_start_curve+(x_end_curve-x_start_curve)/2)
    {
      x_bent_curve=x_start_curve+(bent_radius+y) *
        sin(pi/(x_end_curve-x_start_curve)*(x-x_start_curve));
      y_bent_curve=(bent_radius+y) *
        cos(pi/(x_end_curve-x_start_curve)*(x-x_start_curve));
    }
    else
    {
      if(x<x_start_curve+(x_end_curve-x_start_curve))
      {
        x_bent_curve=x_start_curve+(bent_radius+y) *
          cos(pi/(x_end_curve-x_start_curve)*(x-x_start_curve-(x_end_curve-x_start_curve)/2));
        y_bent_curve=-(bent_radius+y) *
          sin(pi/(x_end_curve-x_start_curve)
          *(x-x_start_curve-(x_end_curve-x_start_curve)/2));
      }
      else
      {
        x_bent_curve=x_start_curve-(x-x_end_curve);
        y_bent_curve=-(bent_radius+y);
      }

    }
  }
}


void CoordsTrafo_bentPipe_3d_old(double x, double y,double z, double &x_bent_curve, double &y_bent_curve, double &z_bent_curve)
{
  double l_infty = TDatabase::ParamDB->UREA_l_infty;
  double  pi = 3.1415926535897;
  double pipe_radius = TDatabase::ParamDB->UREA_PIPE_RADIUS;
  double bent_radius = 5.6;
  double  x_start_curve = 30;
  double  x_end_curve = 100;
  double lead = 1.5;

  if(x < x_start_curve)
  {
    x_bent_curve=x;
    y_bent_curve=bent_radius+y;
    z_bent_curve=z;
  }
  else
  {
    if(x<=x_start_curve+(x_end_curve-x_start_curve)/2||(x>=170&&x<=205)||(x>=310&&x<=345))
    {
      x_bent_curve=x_start_curve+(bent_radius+y) *
        sin(pi/(x_end_curve-x_start_curve)*(x-x_start_curve));
      y_bent_curve=(bent_radius+y) *
        cos(pi/(x_end_curve-x_start_curve)*(x-x_start_curve));
      z_bent_curve=z+lead/140*(x-x_start_curve);
    }
    else
    {
      if(x<x_start_curve+(x_end_curve-x_start_curve)||(x>=205&&x<=240)||(x>=345&&x<=380))
      {
        x_bent_curve=x_start_curve+(bent_radius+y) *
          cos(pi/(x_end_curve-x_start_curve)*(x-x_start_curve-(x_end_curve-x_start_curve)/2));
        y_bent_curve=-(bent_radius+y) *
          sin(pi/(x_end_curve-x_start_curve)
          *(x-x_start_curve-(x_end_curve-x_start_curve)/2));
        z_bent_curve=z+lead/140*(x-x_start_curve);
      }
      else
      {
        if(x<170||(x>=240&&x<=310)||(x>=380&&x<=450))
        {
          x_bent_curve=x_start_curve-(bent_radius+y) *sin(pi/70*(x-x_end_curve));
          y_bent_curve=-(bent_radius+y) *
            cos(pi/(70)
            *(x-x_end_curve));
          z_bent_curve=z+lead/140*(x-x_start_curve);

        }
        else
	{
          x_bent_curve=x-420;
          y_bent_curve=(bent_radius+y);
          z_bent_curve=z+lead/140*(x-x_start_curve);
        }
      }
    }
  }
}

void CoordsTrafo_bentPipe_3d(double x, double y,double z, double &x_bent_curve, double &y_bent_curve, double &z_bent_curve)
{
  double l_infty = TDatabase::ParamDB->UREA_l_infty;
  double pipe_radius = TDatabase::ParamDB->UREA_PIPE_RADIUS;
  double length_pipe = TDatabase::ParamDB->DRIFT_Z;

  // determined by experiment
  // not yet used
  int number_of_circles = 4;

  double bent_radius = 5.65;
  int number_cells_per_circle = 20; // was 32

  // rounding to nearest positive integer
  int number_inlet_cells = 1; // was 3

  /*
   * This somehow configures the distance between the circles
   * I assume it is the ratio of 2x Wall_Thickness to pipe radius.
   * (Here 0.4/0.3)
   */
  double lead = 4./3.;

  // distance in cm from one sandwich cell to another
  double delta_x = length_pipe/TDatabase::ParamDB->N_CELL_LAYERS;
  // shorten inlet
  double x_prime;

  // outer radius
  double radius = bent_radius+pipe_radius;
 
  // the inlet part is straight
  if(x <= number_inlet_cells*delta_x)
  {
    x_bent_curve = x;
    y_bent_curve = bent_radius+y;
    z_bent_curve = z;
  }
  else
  {
	// on which circle is this x lying (counting from bottom to top)
    int number_circle = (int)(x/delta_x - number_inlet_cells)/(number_cells_per_circle) ;

    // where within that cell is x located ?
    x_prime = (x/delta_x-number_inlet_cells) - number_circle*number_cells_per_circle;

    //outlet
    if( TDatabase::ParamDB->N_CELL_LAYERS - number_inlet_cells
       -(number_circle+1) * number_cells_per_circle < 0 )
    {
      x_bent_curve = number_inlet_cells*delta_x+x_prime+10;
      y_bent_curve = bent_radius+y;
      z_bent_curve = z+(number_circle)*lead;
    }
    // helically coiled middle part
    else
    {
      if( x_prime <= number_cells_per_circle )
      {
	    x_bent_curve = number_inlet_cells*delta_x+(bent_radius+y)
	    		     * sin(x_prime*2*Pi/number_cells_per_circle);
        y_bent_curve = (bent_radius+y)
        		     * cos(x_prime*2*Pi/number_cells_per_circle);
        z_bent_curve = z+(number_circle)*lead
        		     + lead*x_prime/number_cells_per_circle;
      }
    }
  }
  if (x_bent_curve>1)
    x_bent_curve -= 25;
  if ((x_bent_curve>25)&&(x_bent_curve<26))
    x_bent_curve = 5./4 * 14.4;
  if ((x_bent_curve>26)&&(x_bent_curve<27))
    x_bent_curve = 6./4 * 14.4;
  if ((x_bent_curve>27)&&(x_bent_curve<28))
    x_bent_curve = 7./4 * 14.4;
  
}

void SwitchCoords_bent_pipe(TCollection *coll)
{
  int N_Cells, N_V, i, j;
  double x, y, z, x_bent_curve, y_bent_curve, z_bent_curve;
  TVertex *vertex;
  TBaseCell *cell;

  N_Cells = coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0 ; i<N_Cells ; i++)
  {
    cell = coll->GetCell(i);
    N_V = cell->GetN_Vertices();

    for (j=0 ; j<N_V ; j++)
    {
      vertex = cell->GetVertex(j);
      vertex->SetClipBoard(0);
    }
  }

  // switch coordinates and 'normal' coordinates if ClipBoard=0
  for(i=0 ; i<N_Cells ; i++)
  {
    cell = coll->GetCell(i);
    N_V = cell->GetN_Vertices();

    for (j=0 ; j<N_V ; j++)
    {
      vertex = cell->GetVertex(j);
      if (!vertex->GetClipBoard())
      {
        vertex->GetCoords(x, y, z);
        CoordsTrafo_bentPipe_3d(x, y, z, x_bent_curve, y_bent_curve,z_bent_curve );
        vertex->SetCoords(x_bent_curve, y_bent_curve, z_bent_curve);
        vertex->SetNormal(x_bent_curve, y_bent_curve, z_bent_curve);
        vertex->SetClipBoard(1);
      }
    }
  }
}


// ========================================================================
// calculate square coordinates for given circle coordinates
// ========================================================================

void CoordsTrafo_CircleToSquare(double x, double y, double& nx, double& ny)
{
  double phi, r, eps=1e-8;

  // calculate polar coordinates
  r = sqrt(x*x + y*y)*0.7;
  // case x = y = 0
  if (fabs(r)<eps)
  {
    nx = 0;
    ny = 0;
    return;
  }
  else                                            //r>0
  {
    if (y>=0)
      phi = acos(x/r);
    else
      phi = -acos(x/r);

    // transformation onto square coordinates
    if ((fabs(y) - fabs(x))<=eps && x>0 )
    {
      ny = phi*r*4/Pi;
      nx = r;
      return;
    }

    if ((fabs(y)-fabs(x))>eps && y>0)
    {
      nx = (Pi/2 - phi)*r*4/Pi;
      ny = r;
      return;
    }

    if((fabs(y)-fabs(x))<=eps && x<0)
    {
      if (phi<0)
        phi = phi + 2*Pi;
      ny = (Pi - phi)*r*4/Pi;
      nx = -r;
      return;
    }

    if((fabs(y)-fabs(x))>eps && y<0)
    {
      if (phi<0)
        phi = phi + 2*Pi;
      nx = (-3*Pi/2 + phi)*r*4/Pi;
      ny = -r;
      return;
    }
  }
}


// ========================================================================
// save square coordinates for given circle coordinates
// ========================================================================

void SquareCoords(TCollection *coll)
{
  int N_Cells, N_V, i, j;
  double x = 0, y = 0, z = 0, nx, ny;
  TVertex *vertex;
  TBaseCell *cell;

  N_Cells = coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0 ; i<N_Cells ; i++)
  {
    cell = coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    for (j=0 ; j<N_V ; j++)
    {
      vertex = cell->GetVertex(j);
      vertex->SetClipBoard(0);
    }
  }

  // get corresponding coordinates in the square for each vertex
  // save them in 'normal' coordinates of vertex
  for(i=0 ; i<N_Cells ; i++)
  {
    cell = coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    for (j=0 ; j<N_V ; j++)
    {
      vertex = cell->GetVertex(j);
      if (!vertex->GetClipBoard())
      {
        vertex->GetCoords(x, y, z);

        CoordsTrafo_CircleToSquare(x, y, nx, ny);
        vertex->SetNormal(nx, ny, z);
        vertex->SetClipBoard(1);
      }
    }
  }
}

void FindValue(TFEFunction3D *u1)
{
  double val[4];
  u1->FindGradient(60,6.05,4.92,val);
  OutPut("center of outflow " << val[0] << endl);
  //exit(1);
}

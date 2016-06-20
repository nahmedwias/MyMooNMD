#include <Constants.h>


void ExampleFile()
{
  Output::print("Example :  Tunnel1.h");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  values[0]=0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  auto boundBox = [] (const double x, const double y, const double z,
                      const double xmin, const double xmax, 
                      const double ymin, const double ymax, 
                      const double zmin, const double zmax )
  {
    if( (x<xmax && x > xmin )  &&
        (y<ymax && y > ymin)  &&
        (z<zmax && z > zmin))
    {
      return true;
    }
    else
      return false;
  };
  
  auto boundcomp = [] (const double x, const double y, const double z, 
                      const double nx, const double ny, const double nz,
                      const double px, const double py, const double pz)
  {
    if(fabs((x-px)*nx + (y-py)*ny + (z-pz)*nz) < 1e-6)
    {
      return true;
    }
    else 
      return false;
  };
  
  // set: NEUMAN//  Component: 2
  double px2, py2, pz2, nx2, ny2, nz2;
  double xmin2, ymin2, xmax2, ymax2, zmin2, zmax2;
  px2= 810.5529785;  py2=-2.605742216;  pz2=132.7259827;
  nx2=-0.06652197599633372; ny2=-5.165893689914315e-06; nz2=0.9977849601406391;
  
  xmax2= 814.3061157000001; xmin2= 809.2967896;
  ymax2= 0.689115262;       ymin2= -4.330946064;
  zmax2= 133.0695435;       zmin2= 132.5489105;
  
  //NEUMANN Component: 4
  double px4, py4, pz4, nx4, ny4, nz4;
  double xmin4, ymin4, xmax4, ymax4, zmin4, zmax4;
  px4 = 846.5855713; py4=-2.184515238; pz4=163.4323883;
  nx4 =0.08391559697625746; ny4=2.066068724861857e-06; nz4=-0.9964728659526307;
  
  xmax4 = 847.2120716000001; xmin4 = 842.1840576;
  ymax4 = 0.7077504158; ymin4 = -4.366158581;
  zmax4 = 163.5767266;  zmin4 = 162.9701141;
  
  // NEUMANN Component: 5
  double px5, py5, pz5, nx5, ny5, nz5;
  double xmin5, ymin5, xmax5, ymax5, zmin5, zmax5;
  px5 = 788.1939697; py5=-3.196940899; pz5=131.7538757;
  nx5 =-0.06651716722641562; ny5=5.332620354083918e-06; nz5=0.9977852807271393;
  
  xmax5= 789.9502808;  xmin5= 784.9191499999999;
  ymax5= 0.6843026042; ymin5= -4.349986157;
  zmax5= 131.9642883; zmin5= 131.4422146;
  
  // NEUMANN Component: 6
  double px6, py6, pz6, nx6, ny6, nz6;
  double xmin6, ymin6, xmax6, ymax6, zmin6, zmax6;
  px6 = 870.7096558;  py6=-2.310198307; pz6=165.7664032;
  nx6 = -0.09038615888166045; ny6= 0.0002674819171914294; nz6= 0.9959067580532041;
  
  xmax6= 871.4936768; xmin6= 866.5301724;
  ymax6= 0.6559827089;ymin6= -4.366750336;
  zmax6= 165.9283539; zmin6= 165.3228221;  
 
  int bbox2=boundBox(x,y,z,xmin2,xmax2,ymin2,ymax2,zmin2,zmax2);
  int bdcomp2=boundcomp(x,y,z,nx2,ny2,nz2,px2,py2,pz2);
  // 3 is the inflow 
  int bbox4=boundBox(x,y,z,xmin4,xmax4,ymin4,ymax4,zmin4,zmax4);
  int bdcomp4=boundcomp(x,y,z,nx4,ny4,nz4,px4,py4,pz4);
  
  int bbox5=boundBox(x,y,z,xmin5,xmax5,ymin5,ymax5,zmin5,zmax5);
  int bdcomp5=boundcomp(x,y,z,nx5,ny5,nz5,px5,py5,pz5);
  
  int bbox6=boundBox(x,y,z,xmin6,xmax6,ymin6,ymax6,zmin6,zmax6);
  int bdcomp6=boundcomp(x,y,z,nx6,ny6,nz6,px6,py6,pz6);
  
  if( ( bbox2 && bdcomp2) || (bbox4 && bdcomp4) ||
      ( bbox5 && bdcomp5) || (bbox6 && bdcomp6) )
  { 
    cond = NEUMANN;
  }
  else
  {
    cond = DIRICHLET;
  }
 
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  auto boundBox = [] (const double x, const double y, const double z,
                      const double xmin, const double xmax, 
                      const double ymin, const double ymax, 
                      const double zmin, const double zmax )
  {
    if( (x<xmax && x > xmin )  &&
        (y<ymax && y > ymin)  &&
        (z<zmax && z > zmin))
    {
      return true;
    }
    else
      return false;
  };
  
  auto boundcomp = [] (const double x, const double y, const double z, 
                      const double nx, const double ny, const double nz,
                      const double px, const double py, const double pz)
  {
    if(fabs((x-px)*nx + (y-py)*ny + (z-pz)*nz) < 1e-6)
    {
      return true;
    }
    else 
      return false;
  };

  /// 3 is the inflow part 
  // Inflow part of bc Component: 3
  double px3, py3, pz3, nx3, ny3, nz3;
  double xmin3, ymin3, xmax3, ymax3, zmin3, zmax3;
  px3 =5.466983376; py3=34.4424; pz3=-0.4399964431;
  nx3 =-3.883987658032814e-06; ny3=-0.9999999999924574; nz3=-1.939282693995977e-09; 
  xmax3= 5.566983376; xmin3= -5.583245613;
  ymax3= 34.54240189; ymin3= 34.3424;
  zmax3= 5.583052254; zmin3= -5.576516723999999;
  
  int bbox3=boundBox(x,y,z,xmin3,xmax3,ymin3,ymax3,zmin3,zmax3);
  int bdcomp3=boundcomp(x,y,z,nx3,ny3,nz3,px3,py3,pz3);
  
  // TODO: ask Abid about the inflow b.c
  if( bbox3 && bdcomp3)
  {
    value = 1; // what is the inflow here ??? 
  }
  else
    value =0; 
}

void U2BoundValue(double x, double y, double z, double &value)
{
  value =0;
}

void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}
// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  // TODO: what will be the reynold's number
  const double eps = 1/TDatabase::ParamDB->RE_NR;
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;
    // f1
    coeffs[i][1] = 0;
    // f2
    coeffs[i][2] = 0;
    // f3
    coeffs[i][3] = 0;
  }
}

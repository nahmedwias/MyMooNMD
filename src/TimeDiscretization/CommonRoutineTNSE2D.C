#include <CommonRoutineTNSE2D.h>
#include <Database.h>
//=============================================================================
/******************************************************************************/
// Stabilization parameters for residual-based VMS methods
// of Bazilevs et al. (2007)
// references go to Ahmed et. al, Arch. Computat. Methods Eng. 24, 115 - 164 (2017)
/******************************************************************************/
void stabilization_parameters_equal_order(double Mult, double* u, double* coeff,
                                          double* params)
{
  double x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
  double y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
  double x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
  double y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
  double x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
  double y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
 
  double d11, d12, d21, d22;
  double rec_detjk = 1./coeff[19];
  
  
  // triangle
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[3]== -4711)
  {
    d11 = (y2-y0) * rec_detjk;  //dxi/dx
    d12 = (x0-x2) * rec_detjk;  //dxi/dy
    d21 = (y0-y1) * rec_detjk;  //deta/dx
    d22 = (x1-x0) * rec_detjk;  //deta/dy
  }
  else
  {
    double x3 = TDatabase::ParamDB->INTERNAL_VERTEX_X[3];
    double y3 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[3];
    double xc1=(-x0 + x1 + x2 - x3) * 0.25;
    double xc2=(-x0 - x1 + x2 + x3) * 0.25;
    double yc1=(-y0 + y1 + y2 - y3) * 0.25;
    double yc2=(-y0 - y1 + y2 + y3) * 0.25;
    
    double rec_detjk = 1./(xc1*yc2 - xc2*yc1);
    // quadrilateral
    d11 = yc2 * rec_detjk;  //dxi/dx
    d12 = -xc2 * rec_detjk;  //dxi/dy
    d21 = -yc1 * rec_detjk;  //deta/dx
    d22 = xc1 * rec_detjk;  //deta/dy
  }

  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double cinv = TDatabase::ParamDB->DELTA0;
  double tau_m, tau_c;

  double g11 = d11*d11 + d21*d21;
  double g12 = d11*d12 + d21*d22;
  double g22 = d12*d12 + d22*d22;

  tau_m = g11*g11 + 2.*g12*g12 + g22*g22; // G:G
  tau_m *= cinv*coeff[0]*coeff[0];
  tau_m +=  4./(dt*dt);
  tau_m += u[0] * (g11*u[0]+g12*u[1]) + u[1]*(g12*u[0]+g22*u[1]);
  tau_m = 1./sqrt(tau_m);

  tau_c = (d11+d21)*(d11+d21)+(d12+d22)*(d12+d22);
  tau_c *= tau_m;
  tau_c = 1./tau_c;

  params[0] = tau_m;
  params[1] = tau_c;
  
  // for output 
  TDatabase::ParamDB->P14 = tau_m;
  TDatabase::ParamDB->P15 = tau_c;
}

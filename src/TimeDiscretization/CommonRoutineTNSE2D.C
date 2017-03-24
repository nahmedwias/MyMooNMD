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
    d11 = (y2-y1) * rec_detjk;  //dxi/dx
    d12 = (x1-x2) * rec_detjk;  //dxi/dy
    d21 = (y0-y1) * rec_detjk;  //deta/dx
    d22 = (x1-x0) * rec_detjk;  //deta/dy
  }
  else
  {
    // quadrilateral
    d11 = (y2-y1) * 0.5 * rec_detjk;  //dxi/dx
    d12 = (x1-x2) * 0.5 * rec_detjk;  //dxi/dy
    d21 = (y0-y1) * 0.5 * rec_detjk;  //deta/dx
    d22 = (x1-x0) * 0.5 * rec_detjk;  //deta/dy
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
}

// compute square of the Frobenius norm of the tensor
// use deformation tensor
double frobeniusNormTensor(double *gradu)
{
int viscosity_tensor = TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR;
double frobenius_norm_tensor = 0;
if (viscosity_tensor==0)
{
  // compute (grad(u)+grad(u)^T)/2
  double a = gradu[0]+gradu[0];
  double b = gradu[1]+gradu[2];
  double c = gradu[3]+gradu[3];
  frobenius_norm_tensor = (a*a+ 2*b*b + c*c)/4.0;
}
// use grad u
else
{
  frobenius_norm_tensor
     =  gradu[0]*gradu[0] + gradu[1]*gradu[1]
        + gradu[2]*gradu[2] + gradu[3]*gradu[3];
}
return frobenius_norm_tensor;
}

double turbulentviscosity(double hK, double* gradU, double* u, double* uConv)
{
/// compute Characteristic Filter Width
double filter_constant = TDatabase::ParamDB->FILTER_WIDTH_CONSTANT;
double filter_power = TDatabase::ParamDB->FILTER_WIDTH_POWER;

double delta;
if (filter_power==1)
  delta = filter_constant * hK;
else
  delta = filter_constant*pow(hK,filter_power);

int viscosity_type = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
double frobenius_norm_tensor = frobeniusNormTensor(gradU);

double viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
double viscosity_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
double viscosity_sigma= TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
double viscosity;

switch(viscosity_type){
  case 0: //
    viscosity = 0;
    break;
  case 1:
    viscosity = viscosity_constant * delta * delta * sqrt(frobenius_norm_tensor);
    break;
  case 2:
    viscosity = viscosity_constant * delta * delta *
                pow(frobenius_norm_tensor, viscosity_power/2.0);
    viscosity = viscosity_constant * delta * delta *
                pow(frobenius_norm_tensor, viscosity_power-2.0);
    break;
  case 3:
    viscosity = viscosity_constant * pow(delta,viscosity_sigma) *
                  pow(frobenius_norm_tensor,(viscosity_power-2)/2.0)
                  / pow(fabs(log(delta)),0.5*(viscosity_power-1));
    break;
  case 4:
    // \|u-g_\delta \ast u\|_2
    viscosity = viscosity_constant * delta * sqrt((u[0]-uConv[0])*(u[0]-uConv[0])+
               (u[1]-uConv[1])*(u[1]-uConv[1]));
    break;
  case 5:
    viscosity = viscosity_constant;
    break;
  case 6:
    if(TDatabase::ParamDB->DEFECT_CORRECTION_TYPE)
      viscosity = delta;
    else
      viscosity = 0;
    break;
  default:
    ErrThrow("This type of turbulent viscosity is not implemented !!!" );
    break;
}

return viscosity;
}


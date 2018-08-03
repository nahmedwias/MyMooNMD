#include <CommonRoutineTNSE2D.h>
#include <Database.h>

// ======================================================================
// compute turbulent viscosity for LES
// ======================================================================
double turbulentviscosity(double hK, double* gradU, double* u, double* uConv)
{
  double delta;
  double filter_constant = TDatabase::ParamDB->FILTER_WIDTH_CONSTANT;
  double filter_power = TDatabase::ParamDB->FILTER_WIDTH_POWER;

  if (filter_power==1)
    delta = filter_constant *hK;
  else
    delta =  filter_constant*pow(hK,filter_power);
  
  int nu_type = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
  double nu_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  int nu_tensor =  TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR;
  double nu_power, nu_sigma;
  double frobenius_norm_tensor,nu,a,b,c;

  // compute turbulent viscosity
  switch(nu_type)
  {
      // no turbulent viscosity
      case 0:
	  nu = 0;
	  return(nu);
	  break;
      case 5:
	  nu = nu_constant;
	  return(nu);
	  break;
	  // defect correction  
      case 6:
	  if (TDatabase::ParamDB->DEFECT_CORRECTION_TYPE)	  
	      nu = delta;
	  else
	      nu = 0;  // bwe
	  return(nu);
	  break;
  }

  // compute square of the Frobenius norm of the tensor
  // use deformation tensor
  if (nu_tensor==0)
  {
    // compute (grad(u)+grad(u)^T)/2
    a = gradU[0]+gradU[0];
    b = gradU[1]+gradU[2];
    c = gradU[3]+gradU[3];
    frobenius_norm_tensor = (a*a+ 2*b*b + c*c)/4.0;
  }
  // use grad u
  else
    frobenius_norm_tensor =  gradU[0]*gradU[0] + gradU[1]*gradU[1] +
      gradU[2]*gradU[2] + gradU[3]*gradU[3];

  // compute turbulent viscosity
  switch(nu_type)
  {
    case 1:                      // Smagorinsky
    case 17:                     // Smagorinsky
      nu = nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
      break;
    case 2:                      // p laplacian
      nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
      nu = nu_constant * delta * delta * pow(frobenius_norm_tensor, nu_power/2.0);
      nu = nu_constant * delta * delta * pow(frobenius_norm_tensor, nu_power-2.0);
      break;
    case 3:                      // Layton SIAM J. Sci. Comput. 1996
      // nu = nu_0 |ln(h)|^(-0.5*(p-1)) h^sigma \|\nabla u\|^(p-2)
      nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
      nu_sigma = TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
      nu = nu_constant * pow(delta,nu_sigma) *
        pow(frobenius_norm_tensor,(nu_power-2)/2.0)
        / pow(fabs(log(delta)),0.5*(nu_power-1));
      break;
    case 4:                      // \|u-g_\delta \ast u\|_2
      nu = nu_constant * delta * sqrt((u[0]-uConv[0])*(u[0]-uConv[0])+
        (u[1]-uConv[1])*(u[1]-uConv[1]));
      //   OutPut(u[0] << " " << uConv[0] << "        "<< endl);
      //    OutPut(u[1] << " " << uConv[1] << endl);
      break;
    default:
      OutPut("This type of turbulent viscosity is not implemented !!!" << endl);
      exit(4711);
      break;
  }
  return(nu);
}

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

// ========================================================================
// parameters: u1old, u2old
// used for : GALERKIN
//            COLETTI (non linear steps)
// ========================================================================
void TimeNSParams2(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
}

// ========================================================================
// u1old, u2old, and the previous time solutions for BDF2 time 
// stepping scheme
// ========================================================================
void TimeNSType14ParamsSUPG(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // u1-previous time
  out[3] = in[5];                // u2-previous time
}

// ========================================================================
// parameters: u1old, u2old,
// used for : COLETTI, Smagorinsky
// ========================================================================
void TimeNSParamsVelo_GradVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2
}

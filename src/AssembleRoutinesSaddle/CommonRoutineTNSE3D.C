#include <CommonRoutineTNSE3D.h>

#include <Database.h>
#include <Convolution.h>
#include <ConvDiff.h>

double standard(double delta, double frobeniusNorm)
{
  double viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * delta * sqrt(frobeniusNorm);
  return viscosity;
}

double laplacian(double delta, double frobeniusNorm)
{
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double  viscosity_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;

  double viscosity = viscosity_constant * delta * delta * 
                      pow(frobeniusNorm, viscosity_power/2.0);
  return viscosity;
}

// Layton SIAM J. Sci. Comput. 1996 
double getViscosityLayton96(double delta, double frobeniusNorm)
{
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double  viscosity_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
  double  viscosity_sigma = TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
  
  double viscosity = viscosity_constant * pow(delta,viscosity_sigma) * 
                 pow(frobeniusNorm,(viscosity_power-2)/2.0)
                 /pow(fabs(log(delta)),2*(viscosity_power-1)/3);
  return viscosity;
}
// \|u-g_\delta \ast u\|_2 
// give better name
double normof_velocityminus_gdeltaTimesuStar(double delta, double* u, double *uConv)
{
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * 
              sqrt((u[0]-uConv[0])*(u[0]-uConv[0]) +
                   (u[1]-uConv[1])*(u[1]-uConv[1]) +
                   (u[2]-uConv[2])*(u[2]-uConv[2]));
  return viscosity;
}
// \|Du^h-G^H\|_2
// the values of the off diagonal entries of G^H has to be divided by 2
// since the basis functions are l^H/2 !!!
double normof_Defvelocityminus_gdeltaTimesuStar(double delta, double* gradu,
                                                double *uConv)
{
  double a11, a12, a13, a22, a23, a33;
  a11 = gradu[0] - uConv[0];
  a12 = (gradu[1]+gradu[3]-uConv[1])/2.0;
  a13 = (gradu[2]+gradu[6]-uConv[2])/2.0;
  a22 = gradu[4]-uConv[3];
  a23 = (gradu[5]+gradu[7]-uConv[4])/2.0;
  a33 = gradu[8]-uConv[5];
  
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * delta * 
              sqrt(a11*a11+2*a12*a12+2*a13*a13 
                   +a22*a22+2*a23*a23+a33*a33);
  return viscosity;
}
// all parameters free
double parametersFree(double delta, double frobeniusNorm)
{
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double  viscosity_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
  double viscosity = viscosity_constant * pow(delta,viscosity_power) * 
      pow(frobeniusNorm,viscosity_power/2.0);
  return viscosity;
}
// walls at z=0, and z = 2
double vanDriestDampingChannel(double delta, double frobeniusNorm, double* z)
{
  // walls at z=0 and z=2
  double reynolds_number = TDatabase::ParamDB->RE_NR;
  double viscosity;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double van_driest_damping = 26.0;
  double zplus = reynolds_number*(1-fabs(1-z[0]));
  if(zplus < 5)
  {
    viscosity = viscosity_constant * delta * delta * (1-exp(-zplus/van_driest_damping))
                * (1-exp(-zplus/van_driest_damping)) * sqrt(frobeniusNorm);
  }
  else
  {
    viscosity = viscosity_constant * delta * delta * sqrt(frobeniusNorm);
  }
  return viscosity;
}

// van Driest damping for cylinder with squared cross--section
// left and right wall at the cylinder
double vanDriestDampingCylinder(double delta, double frobeniusNorm, double* x, 
                                double* y)
{
  double x0, y0;
  double van_driest_damping = 26.0;
  double eps = 1e-06;
  double zplus = 1000;
  
  if ((x[0] > 0.45 - eps) && (x[0] < 0.55 + eps))
  {
    // distance to the wall
    if (y[0] > 0.7)
      y0 = y[0] - 0.75;
    else
      y0 = 0.65 - y[0];
    // wall units 
    zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES * y0;
  }
  if ((y[0] > 0.65 - eps) && (y[0] < 0.75 + eps))
  {
    // distance to the wall
    if (x[0] < 0.5)
    {
      x0 = 0.45 - x[0];
      // wall units 
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT * x0;
    }
    else
    {
      x0 = x[0] - 0.55;
      // wall units 
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK * x0;
    }
  }
  double viscosity;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  if(zplus < 5)
  {
    viscosity = viscosity_constant * delta * delta * (1-exp(-zplus/van_driest_damping)) *
                (1-exp(-zplus/van_driest_damping)) * sqrt(frobeniusNorm);
  }
  else
  {
    viscosity = viscosity_constant * delta * delta * sqrt(frobeniusNorm);
  }
  return viscosity;
}
// walls at z=0 and z=2
double vanDriestDampingChannelContinuous(double delta, double frobeniusNorm, double* z)
{
  double reynolds_number = TDatabase::ParamDB->RE_NR;
  double zplus = reynolds_number*(1-fabs(1-z[0]));
  double van_driest_damping = 26.0;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * delta * (1-exp(-zplus/van_driest_damping)) *
              (1-exp(-zplus/van_driest_damping)) * sqrt(frobeniusNorm);

  return viscosity;
}

// van Driest damping for channel flow (paper: Rudman, Blackburn'99)
  // walls at z=0 and z=2
double vanDriestDampingChannelRB99(double delta, double frobeniusNorm, double* z)
{
  double reynolds_number = TDatabase::ParamDB->RE_NR;
  double zplus = reynolds_number*(1-fabs(1-z[0]));
  double van_driest_damping = 26.0;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  
  double viscosity = viscosity_constant * delta * delta * 
              (1-exp(-(zplus/van_driest_damping)*(zplus/van_driest_damping) *
              (zplus/van_driest_damping))) * sqrt(frobeniusNorm);

  return viscosity;
}

// van Driest damping for channel flow (paper: Rudman, Blackburn'99) with diff A+
double vanDriestDampingChannelRB99APlus(double delta, double frobeniusNorm, double* z)
{
  double reynolds_number = TDatabase::ParamDB->RE_NR;
  double van_driest_damping = 17.0;
  // walls at z=0 and z=2
  double zplus = reynolds_number*(1-fabs(1-z[0]));
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * delta * (1-exp(-(zplus/van_driest_damping) *
              (zplus/van_driest_damping)*(zplus/van_driest_damping))) * 
              sqrt(frobeniusNorm);
  return viscosity;
}

// van Driest damping (continuous, classical) for cylinder with squared cross--section
double vanDriestDampingCylinderContinuous(double delta, double frobeniusNorm, double* x, 
                                double* y)
{
  // left and right wall at the cylinder
  double x0, y0;
  double zplus = 1000;
  double eps = 1e-06;
  
  if ((x[0] > 0.45 - eps) && (x[0] < 0.55 + eps))
  {
    // distance to the wall
    if (y[0] > 0.7)
      y0 = y[0] - 0.75;
    else
      y0 = 0.65 - y[0];
    // wall units
    zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES * y0;
  }
  if ((y[0] > 0.65 - eps) && (y[0] < 0.75 + eps))
  {
    // distance to the wall
    if (x[0] < 0.5)
    {
      x0 = 0.45 - x[0];
      // wall units
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT * x0;
    }
    else
    {
      x0 = x[0] - 0.55;
      // wall units
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK * x0;
    }
  }
  double van_driest_damping = 26.0;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * delta * (1-exp(-zplus/van_driest_damping)) *
              (1-exp(-zplus/van_driest_damping)) * sqrt(frobeniusNorm);
  
  return viscosity;
}
// van Driest damping (paper: Rudman, Blackburn'99) for cylinder with 
// squared cross--section
double vanDriestDampingCylinderRB99(double delta, double frobeniusNorm, double* x,
                                    double *y)
{
  // left and right wall at the cylinder
  double x0, y0, zplus = 1000;
  double eps = 1e-06;
  
  if ((x[0] > 0.45 - eps) && (x[0] < 0.55 + eps))
  {
    // distance to the wall
    if (y[0] > 0.7)
      y0 = y[0] - 0.75;
    else
      y0 = 0.65 - y[0];
    // wall units
    zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES * y0;
  }
  if ((y[0] > 0.65 - eps) && (y[0] < 0.75 + eps))
  {
    // distance to the wall
    if (x[0] < 0.5)
    {
      x0 = 0.45 - x[0];
      // wall units
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT * x0;
    }
    else
    {
      x0 = x[0] - 0.55;
      // wall units
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK * x0;
    }
  }
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double van_driest_damping = 26.0;
  
  double viscosity = viscosity_constant* delta * delta *
             (1-exp(-(zplus/van_driest_damping)*(zplus/van_driest_damping) *
             (zplus/van_driest_damping))) *  sqrt(frobeniusNorm);
  return viscosity;
}
/// Verstappen model (J Sci Comput'11) 
/// C = 1/mu_max as on page 107
double VerstappenViscositymodelCmuMax(double delta, double* gradu, double frobeniusNorm)
{
  double a11, a12, a13, a22, a23, a33;
  // compute filter width in coordinate directions
  // hk is just the value if the filter width would be zero (this should not happen)
  double hk = delta/2.0;
  double delta_x = Mesh_size_in_convection_direction(hk,1,0,0);
  delta_x *= delta_x;
  double delta_y = Mesh_size_in_convection_direction(hk,0,1,0);
  delta_y *= delta_y;
  double delta_z = Mesh_size_in_convection_direction(hk,0,0,1);
  delta_z *= delta_z;
  
  double mu_max = 4 * ( 1./delta_x + 1./delta_y + 1./delta_z );
  if(TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR == 0)
  {
    // 0.5*(grad(u) + grad(u)^T)
    a11 = gradu[0]+gradu[0];
    a11 /= 2.;
    a12 = gradu[1]+gradu[3];
    a12 /= 2.;
    a13 = gradu[2]+gradu[6];
    a13 /= 2.;
    a22 = gradu[4]+gradu[4];
    a22 /= 2.;
    a23 = gradu[5]+gradu[7];
    a23 /= 2.;
    a33 = gradu[8]+gradu[8];
    a33 /= 2.;
  }
  else
    ErrThrow("ERROR: Verstappen model needs a symmetric stress tensor!");
  
  double invariant_3 = - (a11*a22*a33 + 2.*a12*a23*a13 - a12*a12*a33 - 
                          a23*a23*a11 - a13*a13*a22);
  double invariant_2 = 0.5 * frobeniusNorm;
      
  double viscosity = (1.5 * fabs(invariant_3) ) / (mu_max * invariant_2);
  return viscosity;
}
// Verstappen model (J Sci Comput'11) 
// C = (h/pi)^2 as on page as on page 97, where Delta = (h_x*h_y*h_z)^(1/3) 
double VerstappenViscositymodelChPhi(double delta, double* gradu, double frobeniusNorm)
{
  double a11, a12, a13, a22, a23, a33;
  // compute filter width in coordinate directions
  // hk is just the value if the filter width would be zero (this should not happen)
  double hk = delta/2.0;
  double delta_x = Mesh_size_in_convection_direction(hk,1,0,0);
  double delta_y = Mesh_size_in_convection_direction(hk,0,1,0);
  double delta_z = Mesh_size_in_convection_direction(hk,0,0,1);
  
  // TODO: change cell width hk using CELL_MEASURE (more elegant), now too slow! 
  
  hk = delta_x*delta_y*delta_z;
  hk = pow(hk,1.0/3.0);
  
  if(TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR==0)
  {
    // 0.5*(grad(u) + grad(u)^T)
    a11 = gradu[0]+gradu[0];
    a11 /= 2.;
    a12 = gradu[1]+gradu[3];
    a12 /= 2.;
    a13 = gradu[2]+gradu[6];
    a13 /= 2.;
    a22 = gradu[4]+gradu[4];
    a22 /= 2.;
    a23 = gradu[5]+gradu[7];
    a23 /= 2.;
    a33 = gradu[8]+gradu[8];
    a33 /= 2.;
  }
  else
    ErrThrow("ERROR: Verstappen model needs a symmetric stress tensor!");
  
  double invariant_3 = - (a11*a22*a33 + 2.*a12*a23*a13 - a12*a12*a33 -
                          a23*a23*a11 - a13*a13*a22);
  double invariant_2 = 0.5 * frobeniusNorm;
      
  double viscosity = ( 1.5 * hk * hk * fabs(invariant_3) )  / ( Pi * Pi * invariant_2 );  
  return viscosity;
}
// eddy viscosity model: Vreman, Phys. Fluids 16 (10), 3670 -3681, 2004
// frobenius norm of gradient of velocity
// use same notations as in paper
double vermanViscosityModel(double delta, double *gradu)
{
  double hk, delta_x, delta_y, delta_z;
  double b11, b12, b13, b22, b23, b33;
  double Alpha, Bbeta;
  
  Alpha = gradu[0]*gradu[0]+gradu[1]*gradu[1] +
                 gradu[2]*gradu[2]+gradu[3]*gradu[3] +
                 gradu[4]*gradu[4]+gradu[5]*gradu[5] +
                 gradu[6]*gradu[6]+gradu[7]*gradu[7] +
                 gradu[8]*gradu[8];
  double viscosity;
  if(fabs(Alpha)<1e-12)
  {
    viscosity = 0;
  }
  else
  {
    // compute filter width in coordinate directions
    // hk is just the value if the filter width would be zero 
    // (this should not happen)
    hk = delta/2.0;
    delta_x = Mesh_size_in_convection_direction(hk,1,0,0);
    delta_x *=delta_x;
    delta_y = Mesh_size_in_convection_direction(hk,0,1,0);
    delta_y *=delta_y;
    delta_z = Mesh_size_in_convection_direction(hk,0,0,1);
    delta_z *=delta_z;
    
    // compute second invariant of gradient of velocity, scaled with 
    // filter widht in coordinate directions
    b11 = delta_x*gradu[0]*gradu[0]+delta_y*gradu[3]*gradu[3]+delta_z*gradu[6]*gradu[6];
    b22 = delta_x*gradu[1]*gradu[1]+delta_y*gradu[4]*gradu[4]+delta_z*gradu[7]*gradu[7];
    b33 = delta_x*gradu[2]*gradu[2]+delta_y*gradu[5]*gradu[5]+delta_z*gradu[8]*gradu[8];
    b12 = delta_x*gradu[0]*gradu[1]+delta_y*gradu[3]*gradu[4]+delta_z*gradu[6]*gradu[7];
    b13 = delta_x*gradu[0]*gradu[2]+delta_y*gradu[3]*gradu[5]+delta_z*gradu[6]*gradu[8];
    b23 = delta_x*gradu[1]*gradu[2]+delta_y*gradu[4]*gradu[5]+delta_z*gradu[7]*gradu[8];
    Bbeta = b11*b22 - b12*b12 +b11*b33 - b13*b13 + b22*b33 - b23*b23;
      // check for round-off errors
    if (Bbeta<0)
    {
      Bbeta = 0;
    }
    double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
    viscosity = viscosity_constant*sqrt(Bbeta/Alpha);
  }
  return viscosity;
}

double VerstappenModelSimple(double delta, double* gradu, double frobeniusNorm)
{
  double a11, a12, a13, a22, a23, a33;
  if(TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR==0)
  {
    a11 = gradu[0]+gradu[0];
    a12 = gradu[1]+gradu[3];
    a13 = gradu[2]+gradu[6];
    a22 = gradu[4]+gradu[4];
    a23 = gradu[5]+gradu[7];
    a33 = gradu[8]+gradu[8];
  }
  else
    ErrThrow("ERROR: Verstappen model needs a symmetric stress tensor!");
  double invariant_3 = - (a11 * a22 * a33 + 2.* a12 * a23 * a13 - 
                          a12 * a12 * a33 - a23 * a23 * a11 - a13 * a13 * a22);
  invariant_3 /= 8.0;
  double invariant_2 = frobeniusNorm;
  double viscosity = (6.0 * delta * delta * fabs(invariant_3) );
  viscosity /= (Pi * Pi * invariant_2);
}

double frobeniusNormTensor(double* u, double* gradu, double* uConv, int proj_space)
{
  int viscosityTensor = TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR;
  double a11, a12, a13,a22, a23, a33;
  double frobenius_norm_tensor;
  
  if(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 0)
  { 
    // no turbulent viscosity
    return (0);
  }
  
  // compute square of the Frobenius norm of the tensor using 
  // deformation tensor
  switch(viscosityTensor)
  {
    case 0: // 0.5*(grad(u) + grad(u)^T)
      a11 = gradu[0]+gradu[0];
      a12 = gradu[1]+gradu[3];
      a13 = gradu[2]+gradu[6];
      a22 = gradu[4]+gradu[4];
      a23 = gradu[5]+gradu[7];
      a33 = gradu[8]+gradu[8];
      
      frobenius_norm_tensor = 2*(a12*a12 + a13*a13 + a23*a23);
      frobenius_norm_tensor += a11*a11 + a22*a22 + a33*a33;
      frobenius_norm_tensor /= 4;
      break;
    case 1: // grad u form
      frobenius_norm_tensor = gradu[0]*gradu[0] + gradu[1]*gradu[1] + 
                              gradu[2]*gradu[2] + gradu[3]*gradu[3] + 
                              gradu[4]*gradu[4] + gradu[5]*gradu[5] + 
                              gradu[6]*gradu[6] + gradu[7]*gradu[7] +
                              gradu[8]*gradu[8];
      break;
    case 2:
      // deformation tensor of small scales
      // compute (grad(u)+grad(u)^T)/2 - G^H works only with VMS methods
      if (uConv == nullptr)
      {
        OutPut("TURBULENT_VISCOSITY_TENSOR 2 works only with VMS methods !!!" << endl);
        exit(4711);
      }
      a11 = (gradu[0]+gradu[0])/2.0 - uConv[0];
      a12 = (gradu[1]+gradu[3])/2.0 - uConv[1];
      a13 = (gradu[2]+gradu[6])/2.0 - uConv[2];
      a22 = (gradu[4]+gradu[4])/2.0 - uConv[3];
      a23 = (gradu[5]+gradu[7])/2.0 - uConv[4];
      a33 = (gradu[8]+gradu[8])/2.0 - uConv[5];
      frobenius_norm_tensor = 2*(a12*a12 + a13*a13 + a23*a23);
      frobenius_norm_tensor += a11*a11 + a22*a22 + a33*a33;
      break;
    default:
      ErrThrow("Turbulent viscosity tnsor: ", viscosityTensor, " is not implemented so far");
  }
  return frobenius_norm_tensor;
}


double turbulentViscosity3D(double hK, double* u, double* gradu, 
                            double* uConv, double* x, double* y, double* z, double proj_space)
{
  /// compute Characteristic Filter Width
  double filter_constant = TDatabase::ParamDB->FILTER_WIDTH_CONSTANT;
  double filter_power = TDatabase::ParamDB->FILTER_WIDTH_POWER;
  
  double delta, viscosity;

  if (filter_power==1)
    delta = filter_constant * hK;
  else
    delta = filter_constant*pow(hK,filter_power);
  
  // compute the frobenius norm tensor
  double frobenius_norm_tensor=frobeniusNormTensor(u, gradu, nullptr);
  
  int viscosityType = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
  switch(viscosityType)
  {
    case 1:// Smagorinsky
      viscosity = standard(delta, frobenius_norm_tensor);
      break;
    case 2: // p laplacian
      viscosity = laplacian(delta, frobenius_norm_tensor);
      break;
    case 3: // Layton SIAM J. Sci. Comput. 1996 
      viscosity = getViscosityLayton96(delta, frobenius_norm_tensor);
      break;
    case 4: // \|u-g_\delta \ast u\|_2
      viscosity = normof_velocityminus_gdeltaTimesuStar(delta, u, uConv);
      break;
    case 5: // \|Du^h-G^H\|_2
      viscosity = normof_Defvelocityminus_gdeltaTimesuStar(delta, gradu, uConv);
      break;
    case 6: // all parameters free
      viscosity = parametersFree(delta, frobenius_norm_tensor);
      break;
    case 7: // van Driest damping for channel
      viscosity = vanDriestDampingChannel(delta, frobenius_norm_tensor, z);
      break;
    case 8: // van Driest damping for cylinder with squared cross--section
            // left and right wall at the cylinder
      viscosity = vanDriestDampingCylinder(delta, frobenius_norm_tensor, x, y);
      break;
    case 9:  // van Driest damping for channel flow (continuous)
      viscosity = vanDriestDampingChannelContinuous(delta, frobenius_norm_tensor, z);
      break;
    case 10: // van Driest damping for channel flow (paper: Rudman, Blackburn'99)
      viscosity = vanDriestDampingChannelRB99(delta, frobenius_norm_tensor, z);
      break;
    case 11:// van Driest damping for channel flow 
            // (paper: Rudman, Blackburn'99) with diff A+
      viscosity = vanDriestDampingChannelRB99APlus(delta, frobenius_norm_tensor, z);
      break;
    case 12: // van Driest damping (continuous, classical) for cylinder 
             // with squared cross--section
      viscosity = vanDriestDampingCylinderContinuous(delta, frobenius_norm_tensor, x, y);
      break;
    case 13:  // van Driest damping (paper: Rudman, Blackburn'99) 
              // for cylinder with squared cross--section
      viscosity = vanDriestDampingCylinderRB99(delta, frobenius_norm_tensor, x, y);
      break;
    case 14: // Verstappen model (J Sci Comput'11)
      viscosity = VerstappenViscositymodelCmuMax(delta, gradu, frobenius_norm_tensor);
      break;
    case 15: // Verstappen model (J Sci Comput'11)
      viscosity = VerstappenViscositymodelChPhi(delta, gradu, frobenius_norm_tensor);
      break;
    case 16: // eddy viscosity model: Vreman, Phys. Fluids 16 (10), 3670 -3681, 2004
             // frobenius norm of gradient of velocity
             // use same notations as in paper
      viscosity = vermanViscosityModel(delta, gradu);
      break;
    case 17: // simple Verstappen model (J Sci Comput'11, p. 97) 
      viscosity = VerstappenModelSimple(delta, gradu, frobenius_norm_tensor);
      break;
    default:
      ErrThrow("Turbulent viscosity for viscosity type : ", viscosityType, " is not implemented");
  }
  return (viscosity);
}

double DivDivStab3D(double u1, double u2, double u3, double hK, double eps) 
{
  int divdiv_type = TDatabase::ParamDB->DIV_DIV_STAB_TYPE;
  double c_1 = TDatabase::ParamDB->DIV_DIV_STAB_C1;
  double c_2 = TDatabase::ParamDB->DIV_DIV_STAB_C2;
  double tau;
    
  switch(divdiv_type)
  {
    // constant
    case 0:
      tau = c_2;
      break;
      // for non inf-sup stable fe, Codina, Gravemeier
    case 1:
      tau = c_2*sqrt(u1*u1+u2*u2+u3*u3)*hK/c_1;
      tau = tau*tau + eps*eps ;
      tau = sqrt(tau);
      break;
    case 2:
      tau = sqrt(u1*u1+u2*u2+u3*u3)*hK/c_1;
      tau = tau*tau + eps*eps ;
      tau = c_2*sqrt(tau);
      break;
    default:
      OutPut("div-div stabilization " << divdiv_type << " not implemented !!!");
      exit(4711);
  }
  return(tau);
}


void stabilization_parameters_equal_order3D(double Mult, double* u, double* coeff, 
                                          double* params)
{
  double eps  = 1e-12;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double C_I = TDatabase::ParamDB->DELTA0;
  
  double x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
  double y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
  double z0 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[0];
  
  double x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
  double y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
  double z1 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[1];
  
  double x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
  double y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
  double z2 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[2];
  
  
  double nu = coeff[0];
  double rec_detjk = coeff[19];
  rec_detjk = 1/rec_detjk;
  double d11, d12, d13, d21, d22, d23, d31, d32, d33;
  // tetrahedron
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[4] == -4711)
  {
    double x3 = TDatabase::ParamDB->INTERNAL_VERTEX_X[3];
    double y3 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[3];
    double z3 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[3];
    
    d11 = ((y2-y0)*(z3-z0)+(y3-y0)*(z0-z2)) * rec_detjk;  //dxi/dx
    d12 = ((x3-x0)*(z2-z0)+(x2-x0)*(z0-z3)) * rec_detjk;  //dxi/dy
    d13 = ((x2-x0)*(y3-y0)+(x3-x0)*(y0-y2)) * rec_detjk;  //dxi/dz
    
    d21 = ((y3-y0)*(z1-z0)+(y1-y0)*(z0-z3)) * rec_detjk;  //deta/dx
    d22 = ((x1-x0)*(z3-z0)+(x3-x0)*(z0-z1)) * rec_detjk;  //deta/dy
    d23 = ((x3-x0)*(y1-y0)+(x1-x0)*(y0-y3)) * rec_detjk;  //deta/dz
    
    d31 = ((y1-y0)*(z2-z0)+(y2-y0)*(z0-z1)) * rec_detjk;  //dzeta/dx
    d32 = ((x2-x0)*(z1-z0)+(x1-x0)*(z0-z2)) * rec_detjk;  //dzeta/dy
    d33 = ((x1-x0)*(y2-y0)+(x2-x0)*(y0-y1)) * rec_detjk;  //dzeta/dz	
  }
  else // hexahedron
  {
    double x4 = TDatabase::ParamDB->INTERNAL_VERTEX_X[4];
    double y4 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[4];
    double z4 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[4];
    
    d11 = ((y1-y0)*(z4-z0)+(y0-y4)*(z1-z0)) * 0.5 * rec_detjk;  //dxi/dx
    d12 = ((x4-x0)*(z1-z0)+(x1-x0)*(z0-z4)) * 0.5 * rec_detjk;  //dxi/dy
    d13 = ((x1-x0)*(y4-y0)+(x0-x4)*(y1-y0)) * 0.5 * rec_detjk;  //dxi/dz
    
    d21 = ((y4-y0)*(z1-z2)+(y1-y2)*(z0-z4)) * 0.5 * rec_detjk;  //deta/dx
    d22 = ((x1-x2)*(z4-z0)+(x0-x4)*(z1-z2)) * 0.5 * rec_detjk;  //deta/dy
    d23 = ((x4-x0)*(y1-y2)+(x1-x2)*(y0-y4)) * 0.5 * rec_detjk;  //deta/dz
    
    d31 = ((y1-y2)*(z1-z0)+(y1-y0)*(z2-z1)) * 0.5 * rec_detjk;  //dzeta/dx
    d32 = ((x1-x0)*(z1-z2)+(x1-x2)*(z0-z1)) * 0.5 * rec_detjk;  //dzeta/dy
    d33 = ((x1-x2)*(y1-y0)+(x1-x0)*(y2-y1)) * 0.5 * rec_detjk;  //dzeta/dz
  }
  
  double g11 = d11*d11 + d21*d21 + d31*d31;
  double g12 = d11*d12 + d21*d22 + d31*d32;
  double g13 = d11*d13 + d21*d23 + d31*d33;
  double g22 = d12*d12 + d22*d22 + d32*d32;
  double g23 = d12*d13 + d22*d23 + d32*d33;
  double g33 = d13*d13 + d23*d23 + d33*d33;
  
  // G : G
  double tau_m = g11*g11 + 2*g12*g12 + 2*g13*g13 + g22*g22 + 2*g23*g23 + g33*g33; 
  
  tau_m *= C_I*nu*nu;
  tau_m +=  4/(time_step*time_step); 
  tau_m += u[0] * (g11*u[0]+g12*u[1]+g13*u[2]) + u[1]*(g12*u[0]+g22*u[1]+g23*u[2])
           + u[2]*(g13*u[0]+g23*u[1]+g33*u[2]);
  if (tau_m < eps)
  {
    params[0] = 0;
    params[1] = 0;
    return;
  }
  
  // parameter for the mementum equation
  tau_m = 1./sqrt(tau_m);
  double tau_c;
  tau_c  = (d11+d21+d31)*(d11+d21+d31);
  tau_c += (d12+d22+d32)*(d12+d22+d32);
  tau_c += (d13+d23+d33)*(d13+d23+d33);
  
  tau_c *= tau_m;
  tau_c = 1./tau_c;
  
  params[0] = tau_m;
  params[1] = tau_c;
}

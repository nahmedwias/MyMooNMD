// ======================================================================
// @(#)TNSE3D_FixPo.C        1.3 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#include <Database.h>
#include <Convolution.h>
#include <MooNMD_Io.h>
#include <ConvDiff.h>
#include <stdlib.h>
#include <TNSE3D_Routines.h>
#include <MainUtilities.h>
// ======================================================================
// compute turbulent viscosity for LES 
// ======================================================================
double TurbulentViscosity3D(double delta, double* gradU, double* u, 
			    double* uConv, double* x, double* y, double* z, double proj_space)
{
  int nu_type = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
  double nu_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  int nu_tensor =  TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR;
  double Re, zplus, eps = 1e-6, x0, y0;
  double Alpha, b11, b22, b33, b12, b13, b23, Bbeta;
  double nu_power, nu_sigma;
  double delta_x, delta_y, delta_z, hk;
  double mu_max, invariant_2, invariant_3;
  double frobenius_norm_tensor,nu,a11=0,a12=0,a13=0,a22=0,a23=0,a33=0;
  
  // van Driest damping, Do 09.02.06
  double A =26.0;
  //OutPut(proj_space << " ");
  // no turbulent viscosity
  if ((proj_space==2)||(nu_type==0))
  {
    nu = 0;
    //OutPut("ps "  << proj_space << " ");
    return(nu);
  }

  // compute square of the Frobenius norm of the tensor
  // use deformation tensor
  switch(nu_tensor)
  {
    case 0:
      // compute (grad(u)+grad(u)^T)/2
      
      a11 = gradU[0]+gradU[0];
      a12 = gradU[1]+gradU[3];
      a13 = gradU[2]+gradU[6];
      a22 = gradU[4]+gradU[4];
      a23 = gradU[5]+gradU[7];
      a33 = gradU[8]+gradU[8];
      frobenius_norm_tensor = 2*(a12*a12 + a13*a13 + a23*a23);
      frobenius_norm_tensor += a11*a11 + a22*a22 + a33*a33;
      frobenius_norm_tensor /= 4;
    
      break;

    case 1:
      // use grad u
      frobenius_norm_tensor =  gradU[0]*gradU[0] + gradU[1]*gradU[1] +
      gradU[2]*gradU[2] + gradU[3]*gradU[3]  + gradU[4]*gradU[4] + 
      gradU[5]*gradU[5] + gradU[6]*gradU[6] + gradU[7]*gradU[7] +
      gradU[8]*gradU[8];
  
      break;

    case 2:
    // deformation tensor of small scales
    // compute (grad(u)+grad(u)^T)/2 - G^H
    // works only with VMS methods
    if (uConv == NULL)
    {
       OutPut("TURBULENT_VISCOSITY_TENSOR 2 works only with VMS methods !!!" << endl);
       exit(4711);
    }
    a11 = (gradU[0]+gradU[0])/2.0 - uConv[0];
    a12 = (gradU[1]+gradU[3])/2.0 - uConv[1];
    a13 = (gradU[2]+gradU[6])/2.0 - uConv[2];
    a22 = (gradU[4]+gradU[4])/2.0 - uConv[3];
    a23 = (gradU[5]+gradU[7])/2.0 - uConv[4];
    a33 = (gradU[8]+gradU[8])/2.0 - uConv[5];
    frobenius_norm_tensor = 2*(a12*a12 + a13*a13 + a23*a23);
    frobenius_norm_tensor += a11*a11 + a22*a22 + a33*a33;
    
    break;
    
    default:
      OutPut("TURBULENT_VISCOSITY_TENSOR " << TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR  <<
             " not implemented !!!" << endl);
      exit(4711);
  }

  // compute turbulent viscosity
  switch(nu_type)
  {
  case 1: // Smagorinsky
  case 17: // Smagorinsky
    nu = nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
    break;
  case 2: // p laplacian
    nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
    nu = nu_constant * delta * delta * pow(frobenius_norm_tensor, nu_power/2.0);
    break;
  case 3: // Layton SIAM J. Sci. Comput. 1996 
    // nu = nu_0 |ln(h)|^(-2/3*(p-1)) h^sigma \|\nabla u\|^(p-2)
    nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
    nu_sigma = TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
    nu = nu_constant * pow(delta,nu_sigma) * 
      pow(frobenius_norm_tensor,(nu_power-2)/2.0) 
      /pow(fabs(log(delta)),2*(nu_power-1)/3);
    break;
  case 4: // \|u-g_\delta \ast u\|_2
    nu = nu_constant * delta * sqrt((u[0]-uConv[0])*(u[0]-uConv[0])+
                                    (u[1]-uConv[1])*(u[1]-uConv[1])
                                    +(u[2]-uConv[2])*(u[2]-uConv[2]));
//    OutPut(" nu " << nu << " " << fabs((u[0]-uConv[0])) << " " <<  fabs((u[1]-uConv[1]))
//           << " " << fabs((u[2]-uConv[2])) << endl);
    break;
  case 5: // \|Du^h-G^H\|_2
      // the values of the off diagonal entries of G^H has to be divided by 2
      // since the basis functions are l^H/2 !!!
      a11 = gradU[0] - uConv[0];
      a12 = (gradU[1]+gradU[3]-uConv[1])/2.0;
      a13 = (gradU[2]+gradU[6]-uConv[2])/2.0;
      a22 = gradU[4]-uConv[3];
      a23 = (gradU[5]+gradU[7]-uConv[4])/2.0;
      a33 = gradU[8]-uConv[5];
      nu = nu_constant * delta * delta * sqrt(a11*a11+2*a12*a12+2*a13*a13+a22*a22+2*a23*a23+a33*a33);
      //   OutPut(" nu " << nu);// << " " << fabs((u[0]-uConv[0])) << " " <<  fabs((u[1]-uConv[1]))
//           << " " << fabs((u[2]-uConv[2])) << endl);
    break;
  case 6: // all parameters free
    nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
    nu_sigma = TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
    nu = nu_constant * pow(delta,nu_sigma) * 
      pow(frobenius_norm_tensor,nu_power/2.0) ;
    break;
  case 100: // van Driest damping for channel
 	Re = TDatabase::ParamDB->RE_NR;
	// walls at z=0 and z=2
	zplus = Re*(1-fabs(1-z[0]));
	if(zplus < 5)
	  {
	    nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
	     (1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);
	  }
	else
	  {
	    nu  =  nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
	  }
	break;
  case 101: // van Driest damping for cylinder with squared cross--section
	  // left and right wall at the cylinder
	  zplus = 1000;
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
	  
	  if(zplus < 5)
	    {
	      nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
		(1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);
	      //OutPut("nu " << x[0] << " " << y[0] << " " << x0 << " " << y0 << " " << zplus << endl);
	    }
	  else
	    {
	      nu  =  nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
	    }
	  /*    
	      
	  r = sqrt(x0*x0+y0*y0);
	  // find the distance of the center to the boundary of the cylinder
	  // in the direction of (x[0],y[0])
	  // compute the parameter of the line from the center to (x[0],y[0])
	  // which determines the intersection with the boundary of the cylinder
	  int found = 0;
	  if (fabs(x0)>eps)
	  {
	      lambda = 0.05/x0;
	      val = 0.5 + lambda * y0;
	      if ((val>0.45-eps)&&(val<0.55+eps)&&(lambda>0))
		  found = 1;
	  }
	  if ((fabs(x0)>eps)&&(!found))
	  {
	      lambda = -0.05/x0;
	      val = 0.5 + lambda * y0;
	      if ((val>0.45-eps)&&(val<0.55+eps)&&(lambda>0))
		  found = 1;
	  }
	  if ((fabs(y0)>eps)&&(!found))
	  {
	      lambda = 0.05/y0;
	      val = 0.5 + lambda * x0;
	      if ((val>0.45-eps)&&(val<0.55+eps)&&(lambda>0))
		  found = 1;
	  }
	  if ((fabs(y0)>eps)&&(!found))
	  {
	      lambda = -0.05/y0;
	      val = 0.5 + lambda * x0;
	      if ((val>0.45-eps)&&(val<0.55+eps)&&(!found))
		  found = 1;
	  }
	  if (!found)
	  {
	      OutPut("intersection not found !" << endl);
	      exit(4711);
	  }
	  // compute distance
	  r = (1-lambda)*r;
	  //OutPut("aft " << x[0] << " " << y[0] << " " << r <<endl);
 	  Re = TDatabase::ParamDB->RE_NR;
	  // THIS HAS TO BE CORRECTED
	  zplus = Re*r;
	  if(zplus < 5)
	    {
	      nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
		(1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);
	    }
	  else
	    {
	      nu  =  nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
	      }*/
	  break;
  case 102:                                     // van Driest damping for channel flow (continuous)
      Re = TDatabase::ParamDB->RE_NR;
      // walls at z=0 and z=2
      zplus = Re*(1-fabs(1-z[0]));
      nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
        (1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);

      /*OutPut("Van Driest: " << (1-exp(-zplus/A)) *
      (1-exp(-zplus/A)) << endl);*/
    break;
  case 103:                                     // van Driest damping for channel flow (paper: Rudman, Blackburn'99)
      Re = TDatabase::ParamDB->RE_NR;
      // walls at z=0 and z=2
      zplus = Re*(1-fabs(1-z[0]));
      nu = nu_constant * delta * delta * (1-exp(-(zplus/A)*(zplus/A)*(zplus/A))) *  sqrt(frobenius_norm_tensor);

      /*OutPut("Van Driest: " << (1-exp(-zplus/A)) *
      (1-exp(-zplus/A)) << endl);*/
    break;
  case 104:                                     // van Driest damping for channel flow (paper: Rudman, Blackburn'99) with diff A+
      Re = TDatabase::ParamDB->RE_NR;
      A =17.0;
      // walls at z=0 and z=2
      zplus = Re*(1-fabs(1-z[0]));
      nu = nu_constant * delta * delta * (1-exp(-(zplus/A)*(zplus/A)*(zplus/A))) *  sqrt(frobenius_norm_tensor);

      /*OutPut("Van Driest: " << (1-exp(-zplus/A)) *
      (1-exp(-zplus/A)) << endl);*/
    break;

  case 105:
      // eddy viscosity model: Vreman, Phys. Fluids 16 (10), 3670 -3681, 2004
      // frobenius norm of gradient of velocity
      // use same notations as in paper
      
      Alpha = gradU[0]*gradU[0]+gradU[1]*gradU[1]+gradU[2]*gradU[2]+gradU[3]*gradU[3]
              +gradU[4]*gradU[4]+gradU[5]*gradU[5]
              +gradU[6]*gradU[6]+gradU[7]*gradU[7]+gradU[8]*gradU[8];
      if (fabs(Alpha)<1e-12)
      {
        nu = 0;
        break;
      }

      // compute filter width in coordinate directions
      // hk is just the value if the filter width would be zero (this should not happen)
      
      hk = delta/2.0;
      delta_x = Mesh_size_in_convection_direction(hk,1,0,0);
      delta_x *=delta_x;
      delta_y = Mesh_size_in_convection_direction(hk,0,1,0);
      delta_y *=delta_y;
      delta_z = Mesh_size_in_convection_direction(hk,0,0,1);
      delta_z *=delta_z;

      // compute second invariant of gradient of velocity, scaled with filter widht in coordinate directions
      b11 = delta_x*gradU[0]*gradU[0]+delta_y*gradU[3]*gradU[3]+delta_z*gradU[6]*gradU[6];
      b22 = delta_x*gradU[1]*gradU[1]+delta_y*gradU[4]*gradU[4]+delta_z*gradU[7]*gradU[7];
      b33 = delta_x*gradU[2]*gradU[2]+delta_y*gradU[5]*gradU[5]+delta_z*gradU[8]*gradU[8];
      b12 = delta_x*gradU[0]*gradU[1]+delta_y*gradU[3]*gradU[4]+delta_z*gradU[6]*gradU[7];
      b13 = delta_x*gradU[0]*gradU[2]+delta_y*gradU[3]*gradU[5]+delta_z*gradU[6]*gradU[8];
      b23 = delta_x*gradU[1]*gradU[2]+delta_y*gradU[4]*gradU[5]+delta_z*gradU[7]*gradU[8];
      Bbeta = b11*b22 - b12*b12 +b11*b33 - b13*b13 + b22*b33 - b23*b23;
      // check for round-off errors
      if (Bbeta<0)
      {
        Bbeta = 0;
      }

      // scale, in Vreman (2004) it is recommended to use 2.5 times Smagorinsky constant
      nu = nu_constant*sqrt(Bbeta/Alpha);

    break;
  case 106:                                     // van Driest damping (continuous, classical) for cylinder with squared cross--section
      // left and right wall at the cylinder
      zplus = 1000;
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
      nu = nu_constant * delta * delta * (1-exp(-zplus/A)) *
        (1-exp(-zplus/A)) * sqrt(frobenius_norm_tensor);
      //OutPut("nu " << x[0] << " " << y[0] << " " << x0 << " " << y0 << " " << zplus << endl);
    break;
  case 107:                                     // van Driest damping (paper: Rudman, Blackburn'99) for cylinder with squared cross--section
      // left and right wall at the cylinder
      zplus = 1000;
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
      nu = nu_constant * delta * delta * (1-exp(-(zplus/A)*(zplus/A)*(zplus/A))) *  sqrt(frobenius_norm_tensor);
      //OutPut("nu " << x[0] << " " << y[0] << " " << x0 << " " << y0 << " " << zplus << endl);
    break;
      
  case 108: /** Verstappen model (J Sci Comput'11) */
      
      /* C = 1/mu_max as on page 107 */
      
      // compute filter width in coordinate directions
      // hk is just the value if the filter width would be zero (this should not happen)
      hk = delta/2.0;
      delta_x = Mesh_size_in_convection_direction(hk,1,0,0);
      delta_x *= delta_x;
      delta_y = Mesh_size_in_convection_direction(hk,0,1,0);
      delta_y *= delta_y;
      delta_z = Mesh_size_in_convection_direction(hk,0,0,1);
      delta_z *= delta_z;
      
      mu_max = 4 * ( 1./delta_x + 1./delta_y + 1./delta_z );
      
      switch(nu_tensor)
      {
      case 0:
    	  a11 = a11/2.;
    	  a12 = a12/2.;
    	  a13 = a13/2.;
    	  a22 = a22/2.;
    	  a23 = a23/2.;
    	  a33 = a33/2.;
    	break;
    
      case 1:
    	  OutPut("ERROR: Verstappen model needs a symmetric stress tensor!" << endl);
        exit(0);
      }
      
      //invariant_3 = - det(D(u))
      invariant_3 = - (a11*a22*a33 + 2.*a12*a23*a13 - a12*a12*a33 - a23*a23*a11 - a13*a13*a22);
      invariant_2 = 0.5 * frobenius_norm_tensor;
      
      nu = (1.5 * fabs(invariant_3) ) / (mu_max * invariant_2);      
    break;
      
  case 109:  /** Verstappen model (J Sci Comput'11) */
      
      /* C = (h/pi)^2 as on page as on page 97, where Delta = (h_x*h_y*h_z)^(1/3) */
      
      // compute filter width in coordinate directions
      // hk is just the value if the filter width would be zero (this should not happen)
      hk = delta/2.0;
      delta_x = Mesh_size_in_convection_direction(hk,1,0,0);
      delta_y = Mesh_size_in_convection_direction(hk,0,1,0);
      delta_z = Mesh_size_in_convection_direction(hk,0,0,1);
      
      /* TODO: change cell width hk using CELL_MEASURE (more elegant), now too slow! */
      
      hk = delta_x*delta_y*delta_z;
      hk = pow(hk,1.0/3.0);
      
      switch(nu_tensor)
      {
  case 0:
    a11 = a11/2.;
    a12 = a12/2.;
    a13 = a13/2.;
    a22 = a22/2.;
    a23 = a23/2.;   
    a33 = a33/2.; 
    break;
    
  case 1:
    OutPut("ERROR: Verstappen model needs a symmetric stress tensor!" << endl);
    exit(0);
      }
      
      /* invariant_3 = - det(D(u)) */
      invariant_3 = - (a11*a22*a33 + 2.*a12*a23*a13 - a12*a12*a33 - a23*a23*a11 - a13*a13*a22);
      invariant_2 = 0.5 * frobenius_norm_tensor;
      
      nu = ( 1.5 * hk * hk * fabs(invariant_3) )  / ( Pi * Pi * invariant_2 );      
      
      break;
  default:
    OutPut("This type of turbulent viscosity is not implemented !!!" << endl);
    exit(4711);
  break;
  }
  
  return(nu);
}

/******************************************************************************/
//
// computation of SUPG parameter following 
// Bazilevs, Calo, Cottrell, Hughes, Reali, Scovazzi
//
/******************************************************************************/

void SUPG_Param3D(double u1, double u2, double u3, double* coeff, double* params)
{
  ErrThrow("not tested and/or adjusted yet: ");
    double x1, x2, x3, x0, y1, y2, y3, y0, z1, z2, z3, z0, x4, y4, z4;
    double d11, d12, d13, d21, d22, d23, d31, d32, d33, nu;
    double g11, g12, g13, g22, g23, g33;
    double rec_detjk, tau_c, tau_m;
    double eps  = 1e-12;
    double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;    
    double C_I = TDatabase::ParamDB->DELTA0;
    
    nu = coeff[0];             
    rec_detjk = coeff[19];
    rec_detjk = 1/rec_detjk;

    // tetrahedron
    if (TDatabase::ParamDB->INTERNAL_VERTEX_X[4] == -4711)
    {
	x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
	y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
	z0 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[0];
	x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
	y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
	z1 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[1];
	x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
	y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
	z2 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[2];
	x3 = TDatabase::ParamDB->INTERNAL_VERTEX_X[3];
	y3 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[3];
	z3 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[3];
	
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
    else
    {
	// hexahedron
	x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
	y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
	z0 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[0];
	x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
	y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
	z1 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[1];
	x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
	y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
	z2 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[2];
	x4 = TDatabase::ParamDB->INTERNAL_VERTEX_X[4];
	y4 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[4];
	z4 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[4];
	
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

    g11 = d11*d11 + d21*d21 + d31*d31;
    g12 = d11*d12 + d21*d22 + d31*d32;
    g13 = d11*d13 + d21*d23 + d31*d33;
    g22 = d12*d12 + d22*d22 + d32*d32;
    g23 = d12*d13 + d22*d23 + d32*d33;
    g33 = d13*d13 + d23*d23 + d33*d33;
    
    tau_m = g11*g11 + 2*g12*g12 + 2*g13*g13 + g22*g22 + 2*g23*g23 + g33*g33; // G:G
    //OutPut("det "  << rec_detjk << " " << tau_m << endl);
   
    tau_m *= C_I*nu*nu;
    tau_m +=  4/(time_step*time_step); 
    tau_m += u1 * (g11*u1+g12*u2+g13*u3) + u2*(g12*u1+g22*u2+g23*u3)
	+ u3*(g13*u1+g23*u2+g33*u3);
    if (tau_m < eps)
    {
	params[0] = 0;
	params[1] = 0;
	return;
    }
    tau_m = 1/sqrt(tau_m); // this is the parameter for the momentum equation
    
    tau_c = (d11+d21+d31)*(d11+d21+d31)+(d12+d22+d32)*(d12+d22+d32)
	+(d13+d23+d33)*(d13+d23+d33);
//    OutPut(" taucbef " << tau_c << " "  << d11+d21+d31 << " " << d12+d22+d32 << " " << d13+d23+d33 <<
//	   " " << rec_detjk << ":");
    tau_c *= tau_m;
    if (tau_c < eps)
    {
	tau_c = 0;
    }
    else
	tau_c = 1/tau_c;

    //  OutPut(" tauc " << tau_c << " "  << tau_m);
    params[0] = tau_m;
    params[1] = tau_c;

/*
  delta = (d11*d11+d21*d21+d31*d31)*(d11*d11+d21*d21+d31*d31) + 2*(d11*d12+d21*d22+d31*d32)*(d11*d12+d21*d22+d31*d32) +  // G:G
          (d12*d12+d22*d22+d32*d32)*(d12*d12+d22*d22+d32*d32) + 2*(d11*d13+d21*d23+d31*d33)*(d11*d13+d21*d23+d31*d33) +
          2*(d12*d13+d22*d23+d32*d33)*(d12*d13+d22*d23+d32*d33) + (d13*d13+d23*d23+d33*d33)*(d13*d13+d23*d23+d33*d33);
  delta *= C_I*nu*nu;         
  delta += 4/(time_step*time_step);  
  delta += u1*u1*(d11*d11+d21*d21+d31*d31)+2*u1*u2*(d11*d12+d21*d22+d31*d32)+2*u1*u3*(d11*d13+d21*d23+d31*d33)+
           u2*u2*(d12*d12+d22*d22+d32d32)+2*u2*u3*(d12*d13+d22*d23+d32*d33)+u3*u3*(d13*d13+d23*d23+d33d33);            // uGu
  delta = sqrt(delta);
  
  delta /= (d11+d21+d31)*(d11+d21+d31)+(d12+d22+d32)*(d12+d22+d32)+(d13+d23+d33)*;   // gg
  
  return delta;
*/
}

// ======================================================================
// compute parameter for Leray-alpha model
// ======================================================================
double LerayAlpha_Param3D(double hK) 
{
  double c = TDatabase::ParamDB->DELTA0;
    
  return(c*hK);
}


// ======================================================================
// Type 1, Standard Galerkin
// Type 1, ClassicalLES
// Type 1, GL00Convolution
// ======================================================================
void TimeNSType1Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;
    
    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
       
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;
    
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// for Type 1 is no SDFEM available
// ======================================================================

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
        test001*ansatz001);

      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, Smagorinsky
// the nonlinear viscosity is treated implicitly
// ======================================================================
void TimeNSType1Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, Galdi-Layton 98 Model auxiliary problem
// ======================================================================
void TimeNSType1GL00AuxProblem3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM, **AuxMatrix;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double *AuxMatrixRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, mu2, delta;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[3];
  MatrixB2 = LocMatrices[4];
  MatrixB3 = LocMatrices[5];
  AuxMatrix = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);
  Output::print("Please adjust 'TurbulentViscosity3D to the turbulentViscosity3D'" 
   " and check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    AuxMatrixRow = AuxMatrix[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;

      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
            
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, Standard Galerkin
// Type 2, ClassicalLES
// Type 2, GL00Convolution
// ======================================================================
void TimeNSType2Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];
  MatrixB1T = LocMatrices[5];
  MatrixB2T = LocMatrices[6];
  MatrixB3T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// for Type 2 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void TimeNSType2Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];
  MatrixB1T = LocMatrices[5];
  MatrixB2T = LocMatrices[6];
  MatrixB3T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void TimeNSType2Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB3 = LocMatrices[4];
  MatrixB1T = LocMatrices[5];
  MatrixB2T = LocMatrices[6];
  MatrixB3T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i
  
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 2, GL00AuxProblem
// ======================================================================
void TimeNSType2GL00AuxProblem3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixB3, **MatrixM, **AuxMatrix;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixRow3, *MatrixMRow;
  double *AuxMatrixRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, mu2, delta;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[3];
  MatrixB2 = LocMatrices[4];
  MatrixB3 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];
  MatrixB3T = LocMatrices[8];
  AuxMatrix = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f2

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    AuxMatrixRow = AuxMatrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;

      val = ansatz000*test000;
      MatrixMRow[j] += Mult * val;

      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
            
    } // endfor j
    
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// Type 3, ClassicalLES, (grad u, grad v)
// Type 3, GL00Convolution, (grad u, grad v)`
// ======================================================================
void TimeNSType3Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11; // **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row; // *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
//  MatrixM22 = LocMatrices[10];
//  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[10];
  MatrixB2  = LocMatrices[11];
  MatrixB3  = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
//    MatrixM22Row  = MatrixM22[i];
//    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
//      MatrixM22Row[j] += val;
//      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// Type 3, ClassicalLES, D(u):D(v)
// Type 3, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType3GalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = Mult*c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;

    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3UpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, mu;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM = LocMatrices[9];
  MatrixB1  = LocMatrices[10];
  MatrixB2  = LocMatrices[11];
  MatrixB3  = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM[i];
    // MatrixM22Row  = MatrixM22[i];
    // MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      // MatrixM22Row[j] += val;
      // MatrixM33Row[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixMRow;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, viscosity;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM = LocMatrices[9];
  MatrixB1  = LocMatrices[10];
  MatrixB2  = LocMatrices[11];
  MatrixB3  = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
  
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixMRow  = MatrixM[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixMRow[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i      
}

// ======================================================================
// Type 3, GL00AuxProblem (grad u, grad v)
// ======================================================================
void TimeNSType3GL00AuxProblem3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *AuxMatrixRow;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;  

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;
 
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);
  Output::print("Please adjust 'TurbulentViscosity3D to the turbulentViscosity3D'" 
   " and check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;

      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType3GL00AuxProblemDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *AuxMatrixRow;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;
  double delta, mu2, viscosity;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);
  Output::print("Please adjust 'TurbulentViscosity3D to the turbulentViscosity3D'" 
   " and check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");
  
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
 
      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
            
   } // endfor j
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// Type 4, ClassicalLES, (grad u, grad v)
// Type 4, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType4Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11; // **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2, **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row; // *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
//  MatrixM22 = LocMatrices[10];
//  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[10];
  MatrixB2  = LocMatrices[11];
  MatrixB3  = LocMatrices[12];
  MatrixB1T  = LocMatrices[13];
  MatrixB2T  = LocMatrices[14];
  MatrixB3T  = LocMatrices[15];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
//    MatrixM22Row  = MatrixM22[i];
//    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
//      MatrixM22Row[j] += val;
//      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// Type 4, ClassicalLES, D(u):D(v)
// Type 4, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixB1  = LocMatrices[10];
  MatrixB2  = LocMatrices[11];
  MatrixB3  = LocMatrices[12];
  MatrixB1T = LocMatrices[13];
  MatrixB2T = LocMatrices[14];
  MatrixB3T = LocMatrices[15];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// for Type 4 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void TimeNSType4Upwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T  = LocMatrices[15];
  MatrixB2T  = LocMatrices[16];
  MatrixB3T  = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void TimeNSType4UpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[12];
  MatrixB2  = LocMatrices[13];
  MatrixB3  = LocMatrices[14];
  MatrixB1T = LocMatrices[15];
  MatrixB2T = LocMatrices[16];
  MatrixB3T = LocMatrices[17];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4Smagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double **MatrixM;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T, **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixMRow;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM = LocMatrices[9];
  MatrixB1  = LocMatrices[10];
  MatrixB2  = LocMatrices[11];
  MatrixB3  = LocMatrices[12];
  MatrixB1T  = LocMatrices[13];
  MatrixB2T  = LocMatrices[14];
  MatrixB3T  = LocMatrices[15];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
  
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixMRow  = MatrixM[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixMRow[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType4SmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixMRow;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, viscosity;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM = LocMatrices[9];
  MatrixB1  = LocMatrices[10];
  MatrixB2  = LocMatrices[11];
  MatrixB3  = LocMatrices[12];
  MatrixB1T = LocMatrices[13];
  MatrixB2T = LocMatrices[14];
  MatrixB3T = LocMatrices[15];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
  
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixMRow  = MatrixM[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixMRow[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, GL00AuxProblem, (grad u, grad v)
// ======================================================================
void TimeNSType4GL00AuxProblem3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double *AuxMatrixRow;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3;
  double delta, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;  

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[4];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;
 
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);
  Output::print("Please adjust 'TurbulentViscosity3D to the turbulentViscosity3D'" 
   " and check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;

      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType4GL00AuxProblemDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *AuxMatrixRow;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu;
  double delta, mu2, viscosity;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           -4711);
  mu = mu/2.0;
  viscosity = c0+mu;
  
  Output::print("Please adjust 'TurbulentViscosity3D to the turbulentViscosity3D'" 
   " and check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
 
      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;            
   } // endfor j 

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 4, VMS_Projection, D(u):D(v)
// ======================================================================

void TimeNSType4VMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double **MatrixL, **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33;
  double **Matrix_G11, **Matrix_G22,   **Matrix_G33;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_L;
  double c0, c1, c2, c3;
  double u1, u2, u3, mu, viscosity;
  
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixL   = LocMatrices[10];
  MatrixB1  = LocMatrices[11];
  MatrixB2  = LocMatrices[12];
  MatrixB3  = LocMatrices[13];
  MatrixB1T = LocMatrices[14];
  MatrixB2T = LocMatrices[15];
  MatrixB3T = LocMatrices[16];
  Matrix_tilde_G11  = LocMatrices[17];
  Matrix_tilde_G22  = LocMatrices[18];
  Matrix_tilde_G33  = LocMatrices[19];
  Matrix_G11  = LocMatrices[20];
  Matrix_G22  = LocMatrices[21];
  Matrix_G33  = LocMatrices[22];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p
  Orig5 = OrigValues[5]; // l
  
  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old


  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  double *projection_space_label = &param[15];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, projection_space_label[0]);
  
  mu = mu/2.0;
  viscosity = c0+mu;
  
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];
      val1 = Mult*ansatz000;
      
      val = -val1*test100;
      MatrixRow1[j] += val;
      val = -val1*test010;
      MatrixRow2[j] += val;
      val = -val1*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];
     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig5[j];
        val =  Mult * 2*mu * ansatz000;
        Matrix11Row[j] -= val * test100;
        Matrix22Row[j] -= val * test010;
        Matrix33Row[j] -= val * test001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     Matrix11Row = Matrix_G11[i];
     Matrix22Row = Matrix_G22[i];
     Matrix33Row = Matrix_G33[i];
     test000 = Orig5[i];     
     val =  Mult * test000;

     for(j=0;j<N_U;j++)
     {        
        ansatz100 = Orig0[j];
        ansatz010 = Orig1[j];
        ansatz001 = Orig2[j];

        Matrix11Row[j] -= val * ansatz100;
        Matrix22Row[j] -= val * ansatz010;
        Matrix33Row[j] -= val * ansatz001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     test000 = Orig5[i];     
     MatrixRow1 = MatrixL[i];
     for(j=0;j<N_L;j++)
     {
        ansatz000 = Orig5[j];
        MatrixRow1[j] += Mult * test000 * ansatz000;
     }
  }
}
// ======================================================================
// Type 4, VMS_Projection with streamline formulation D(u):D(v)
// ======================================================================

void TimeNSType4VMS_ProjectionStreamlineDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double **MatrixL, **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33;
  double **Matrix_G11, **Matrix_G22,   **Matrix_G33;
  double *Rhs1, *Rhs2, *Rhs3, val, val1, test_stream, ansatz_stream, conv;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_L;
  double c0, c1, c2, c3, delta;
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixL   = LocMatrices[12];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  Matrix_tilde_G11  = LocMatrices[19];
  Matrix_tilde_G22  = LocMatrices[20];
  Matrix_tilde_G33  = LocMatrices[21];
  Matrix_G11  = LocMatrices[22];
  Matrix_G22  = LocMatrices[23];
  Matrix_G33  = LocMatrices[24];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p
  Orig5 = OrigValues[5]; // l
  
  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           param[21]);
Output::print("Please adjust 'TurbulentViscosity3D to the turbulentViscosity3D'" 
   " and check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    val1 = Mult*test000;

    Rhs1[i] += val1*c1;
    Rhs2[i] += val1*c2;
    Rhs3[i] += val1*c3;
    test_stream = (u1*test100+u2*test010+u3*test001);

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      ansatz_stream = (u1*ansatz100+u2*ansatz010+u3*ansatz001);
      conv = ansatz_stream * test000;
      ansatz_stream *= mu*test_stream;     // sd projection term 

      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += conv + ansatz_stream;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += conv + ansatz_stream;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += conv + ansatz_stream;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];
      val1 = Mult*ansatz000;
      
      val = -val1*test100;
      MatrixRow1[j] += val;
      val = -val1*test010;
      MatrixRow2[j] += val;
      val = -val1*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];
    val1 = Mult*test000;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixRow1[j] += val;

      val = -val1*ansatz010;
      MatrixRow2[j] += val;

      val = -val1*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];
     val1 =u1*test100+u2*test010+u3*test001;

     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig5[j];

        val =  Mult * mu * ansatz000 *val1;
        Matrix11Row[j] -= val * u1;
        Matrix22Row[j] -= val * u2;
        Matrix33Row[j] -= val * u3;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     Matrix11Row = Matrix_G11[i];
     Matrix22Row = Matrix_G22[i];
     Matrix33Row = Matrix_G33[i];
     test000 = Orig5[i];     
     val =  Mult * test000;

     for(j=0;j<N_U;j++)
     {        
        ansatz100 = Orig0[j];
        ansatz010 = Orig1[j];
        ansatz001 = Orig2[j];

        Matrix11Row[j] -= val * ansatz100;
        Matrix22Row[j] -= val * ansatz010;
        Matrix33Row[j] -= val * ansatz001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     test000 = Orig5[i];     
     MatrixRow1 = MatrixL[i];
     for(j=0;j<N_L;j++)
     {
        ansatz000 = Orig5[j];
        MatrixRow1[j] += Mult * test000 * ansatz000;
     }
  }
}

// ======================================================================
// Type 4, LerayAlpha, D(u):D(v)
// ======================================================================
void TimeNSType4LerayAlphaDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33, **AuxMatrix;
  double **MatrixM11, **MatrixM22, **MatrixM33;
  double **MatrixB1, **MatrixB2,  **MatrixB3;
  double **MatrixB1T, **MatrixB2T,  **MatrixB3T;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row, *AuxMatrixRow;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double *MatrixRow1, *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3;
  double u1, u2, u3, alpha;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixM11 = LocMatrices[9];
  MatrixM22 = LocMatrices[10];
  MatrixM33 = LocMatrices[11];
  MatrixB1  = LocMatrices[13];
  MatrixB2  = LocMatrices[14];
  MatrixB3  = LocMatrices[15];
  MatrixB1T = LocMatrices[16];
  MatrixB2T = LocMatrices[17];
  MatrixB3T = LocMatrices[18];
  AuxMatrix = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3
 
  u1 = param[0]; // u1 filtered
  u2 = param[1]; // u2 filtered
  u3 = param[2]; // u3 filtered
  
  alpha = LerayAlpha_Param3D(hK);
 Output::print("check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    MatrixM33Row  = MatrixM33[i];
    AuxMatrixRow = AuxMatrix[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = c0*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix22Row[j] += Mult * val;

      val  = c0*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = c0*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = c0*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      Matrix33Row[j] += Mult * val;

      val = Mult*(ansatz000*test000);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
      MatrixM33Row[j] += val;
 
      val  = alpha*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;            
   } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    MatrixRow3 = MatrixB3T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j];

      val = -Mult*ansatz000*test100;
      MatrixRow1[j] += val;
      val = -Mult*ansatz000*test010;
      MatrixRow2[j] += val;
      val = -Mult*ansatz000*test001;
      MatrixRow3[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];
    MatrixRow3 = MatrixB3[i];

    test000 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -Mult*test000*ansatz100;
      MatrixRow1[j] += val;

      val = -Mult*test000*ansatz010;
      MatrixRow2[j] += val;

      val = -Mult*test000*ansatz001;
      MatrixRow3[j] += val;
    } // endfor j
  } // endfor i
}


// ======================================================================
// assemble matrix for auxiliary problem
// ======================================================================

void MatrixAuxiliaryProblem(double Mult, double *coeff, 
                            double *param, double hK, 
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double *AuxMatrixRow, **AuxMatrix;;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U; 
  double delta, mu2, val;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;
 
  // solution does not need to be convolved
  AuxMatrix = LocMatrices[0];
  N_U = N_BaseFuncts[0];
     
  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
    
  // filter width
  delta =  CharacteristicFilterWidth(hK);
  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;
  ErrThrow("not tested and/or adjusted yet: ");
  for(i=0;i<N_U;i++)
  {
    AuxMatrixRow = AuxMatrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
       
      val  = mu2*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += ansatz000*test000;
      AuxMatrixRow[j] += Mult * val;
      
    } // endfor j
  } // endfor i
}

// ======================================================================
// Assembling routine for all nonlinear matrices
// ======================================================================

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// Type 2, Standard Galerkin, only nonlinear part
// ======================================================================
void TimeNSType1_2NLGalerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu
 
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// Type 2, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1_2NLUpwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  c0 = coeff[0]; // nu

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// Type 2, Smagorinsky, only nonlinear part
// Type 1, ClassicalLES, only nonlinear part
// Type 2, ClassicalLES, only nonlinear part
// Type 1, GL00Convolution, only nonlinear part
// Type 2, GL00Convolution, only nonlinear part
// Type 1, GL00AuxProblem, only nonlinear part
// Type 2, GL00AuxProblem, only nonlinear part
// ======================================================================
void TimeNSType1_2NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0, mu;
  double u1, u2, u3;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu
 
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
  
  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
  

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val1;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = Mult*coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = c0*Orig0[i];
    test010 = c0*Orig1[i];
    test001 = c0*Orig2[i];
    test000 = Mult*Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 += (test100*ansatz100+test010*ansatz010 +test001*ansatz001);
      Matrix11Row[j] += test100*ansatz100+val1;
      Matrix22Row[j] += test010*ansatz010+val1;
      Matrix33Row[j] += test001*ansatz001+val1;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// Type 4, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3_4NLUpwind3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  c0 = coeff[0]; // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = Mult*c0*(test100*ansatz100+test010*ansatz010+
                      test001*ansatz001);
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// Type 4, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3_4NLUpwindDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  c0 = coeff[0]; // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      Matrix11Row[j] += Mult * val;
      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      Matrix22Row[j] += Mult * val;
      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v), only nonlinear part
// Type 4, Smagorinsky, (grad u, grad v), only nonlinear part
// Type 3, ClassicalLES, (grad u, grad v), only nonlinear part
// Type 4, ClassicalLES, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 3, GL00AuxProblem, (grad u, grad v), only nonlinear part
// Type 4, GL00AuxProblem, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLSmagorinsky3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix22Row,  *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;
  double c0, mu;
  double u1, u2, u3;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixA33 = LocMatrices[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
 
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = (c0+mu)*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
      Matrix33Row[j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 3, ClassicalLES, D(u):D(v), only nonlinear diagonal blocks
// Type 4, ClassicalLES, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1,val2, val3, val4;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;
  double c0, viscosity, dummy_param[6] = {0, 0, 0, 0, 0, 0};
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &dummy_param[0];

  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
  
  viscosity = Mult*(mu/2.0+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = viscosity*Orig0[i];
    test010 = viscosity*Orig1[i];
    test001 = viscosity*Orig2[i];
    test000 = Mult*Orig3[i];

    for(j=0;j<N_U;j++)
    {      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      //t100a100 = test100*ansatz100;
      //t010a010 = test010*ansatz010;
      //t001a001 = test001*ansatz001;

      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      /* val1 += (test100*ansatz100+test010*ansatz010+test001*ansatz001);
      Matrix11Row[j] += test100*ansatz100+val1;
      Matrix12Row[j] += test100*ansatz010; 
      Matrix13Row[j] += test100*ansatz001;
      Matrix21Row[j] += test010*ansatz100;
      Matrix22Row[j] += test010*ansatz010+val1;
      Matrix23Row[j] += test010*ansatz001;
      Matrix31Row[j] += test001*ansatz100;
      Matrix32Row[j] += test001*ansatz010;
      Matrix33Row[j] += test001*ansatz001+val1;*/
      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;
      val1 += val2+val3+val4;
      Matrix11Row[j] += val2+val1;
      Matrix12Row[j] += test010*ansatz100; 
      Matrix13Row[j] += test001*ansatz100;
      Matrix21Row[j] += test100*ansatz010;
      Matrix22Row[j] += val3+val1;
      Matrix23Row[j] += test001*ansatz010;
      Matrix31Row[j] += test100*ansatz001;
      Matrix32Row[j] += test010*ansatz001;
      Matrix33Row[j] += val4+val1;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLVMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1,val2, val3, val4, valu1, valu2, valu3;
  double **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33; 
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;  // double ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U, N_L;  // int N_P;
  double c0, viscosity;
  double u1, u2, u3, mu;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  Matrix_tilde_G11  = LocMatrices[9];
  Matrix_tilde_G22  = LocMatrices[10];
  Matrix_tilde_G33  = LocMatrices[11];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // l

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  double *x = &param[12];
  double *y = &param[13];
  double *z = &param[14];
  double *u = &param[0];
  double *gradu = &param[3];
  double *uConv = &param[0];
  double *projection_space_label = &param[15];
  
  mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, projection_space_label[0]);
  viscosity = Mult*(mu/2.0+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = viscosity*Orig0[i];
    test010 = viscosity*Orig1[i];
    test001 = viscosity*Orig2[i];
    test000 = Mult*Orig3[i];
    valu1 = u1 * test000;
    valu2 = u2 * test000;
    valu3 = u3 * test000;

    for(j=0;j<N_U;j++)
    {      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      //val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val1 = valu1*ansatz100+valu2*ansatz010+valu3*ansatz001;

      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;
      val1 += val2+val3+val4;
      Matrix11Row[j] += val2+val1;
      Matrix12Row[j] += test010*ansatz100; 
      Matrix13Row[j] += test001*ansatz100;
      Matrix21Row[j] += test100*ansatz010;
      Matrix22Row[j] += val3+val1;
      Matrix23Row[j] += test001*ansatz010;
      Matrix31Row[j] += test100*ansatz001;
      Matrix32Row[j] += test010*ansatz001;
      Matrix33Row[j] += val4+val1;
    } // endfor j
  } // endfor i

  val2 = Mult * mu;
  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];

     for(j=0;j<N_L;j++)
     {       
	 //ansatz000 = Orig4[j];
        val1 = val2 * Orig4[j];
        Matrix11Row[j] -= val1 * test100;
        Matrix22Row[j] -= val1 * test010;
        Matrix33Row[j] -= val1 * test001;
     }
  }   
}

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), adaptive coarse space
// Type 4, VMS_Projection, D(u):D(v), adaptive coarse space
// ======================================================================
void TimeNSType3_4NL_Adap_VMS_ProjectionDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double **MatrixL, **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33;
  double **Matrix_G11, **Matrix_G22,   **Matrix_G33;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double *MatrixRow1; // double *MatrixRow2, *MatrixRow3;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;  // double *Orig5;
  int i,j,N_U, N_L;  // int N_P;
  double c0, delta, val, val1;  // double c1,c2,c3;
  double u1, u2, u3, mu, viscosity;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  MatrixL   = LocMatrices[9];
  Matrix_tilde_G11  = LocMatrices[10];
  Matrix_tilde_G22  = LocMatrices[11];
  Matrix_tilde_G33  = LocMatrices[12];
  Matrix_G11  = LocMatrices[13];
  Matrix_G22  = LocMatrices[14];
  Matrix_G33  = LocMatrices[15];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // l
  
  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[0],&param[12],&param[13],&param[14],
                           param[21]);
  mu = mu/2.0;
  viscosity = c0+mu;
Output::print("Please adjust 'TurbulentViscosity3D to the turbulentViscosity3D'" 
   " and check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];

    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
 
    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      val  = viscosity*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz100);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz100);
      Matrix13Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz010);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = viscosity*(test001*ansatz010);
      Matrix23Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz001);
      Matrix31Row[j] += Mult * val;

      val  = viscosity*(test010*ansatz001);
      Matrix32Row[j] += Mult * val;

      val  = viscosity*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001);
      val += val1;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];

     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig4[j];
        val =  Mult * 2*mu * ansatz000;
        Matrix11Row[j] -= val * test100;
        Matrix22Row[j] -= val * test010;
        Matrix33Row[j] -= val * test001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     Matrix11Row = Matrix_G11[i];
     Matrix22Row = Matrix_G22[i];
     Matrix33Row = Matrix_G33[i];
     test000 = Orig4[i];     
     val =  Mult * test000;

     for(j=0;j<N_U;j++)
     {        
        ansatz100 = Orig0[j];
        ansatz010 = Orig1[j];
        ansatz001 = Orig2[j];

        Matrix11Row[j] -= val * ansatz100;
        Matrix22Row[j] -= val * ansatz010;
        Matrix33Row[j] -= val * ansatz001;
     }
  }   

  for(i=0;i<N_L;i++)
  {
     test000 = Orig4[i];     
     MatrixRow1 = MatrixL[i];
     for(j=0;j<N_L;j++)
     {
        ansatz000 = Orig4[j];
        MatrixRow1[j] += Mult * test000 * ansatz000;
     }
  }
}

// ======================================================================
// Type 3, VMS_Projection explicit, only Matrix_tilde_G??
// Type 4, VMS_Projection explicit, only Matrix_tilde_G??
// ======================================================================
void TimeNSType3_4VMS_ProjectionExpl3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double val1;
  double **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33; 
  double *Matrix11Row, *Matrix22Row, *Matrix33Row;
  double ansatz000;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig4;
  int i,j,N_U, N_L;
  double delta, mu;

  Matrix_tilde_G11  = LocMatrices[0];
  Matrix_tilde_G22  = LocMatrices[1];
  Matrix_tilde_G33  = LocMatrices[2];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig4 = OrigValues[3]; // l
 //OutPut(param[12] << " "<< param[13] << " " <<param[14] << " ");
 //OutPut(param[0] << " "<< param[1] << " " <<param[2] << " "<<endl);
  OutPut("CHECK IF PARAM[21] IS SET !!!"<<endl);
  exit(4711);
  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[15],&param[12],&param[13],&param[14],
                           param[21]);
  Output::print("Please adjust 'TurbulentViscosity3D to the turbulentViscosity3D'" 
   " and check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");
  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];

     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig4[j];
        val1 = Mult * mu * ansatz000;
        Matrix11Row[j] -= val1 * test100;
        Matrix22Row[j] -= val1 * test010;
        Matrix33Row[j] -= val1 * test001;
     }
  }   
}

// ======================================================================
// Type 3, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMS_Projection, D(u):D(v), only nonlinear diagonal blocks
// streamline projection
// ======================================================================
void TimeNSType3_4NLVMS_ProjectionStreamlineDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1,val2, val3, val4;
  double **Matrix_tilde_G11, **Matrix_tilde_G22, **Matrix_tilde_G33; 
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j,N_U,N_L;  // int N_P;
  double c0, delta, test_stream, ansatz_stream;
  double u1, u2, u3, mu, val;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  Matrix_tilde_G11  = LocMatrices[9];
  Matrix_tilde_G22  = LocMatrices[10];
  Matrix_tilde_G33  = LocMatrices[11];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[2];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  Orig4 = OrigValues[4]; // l

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity3D(delta,&param[3],&param[0],&param[12],&param[12],&param[13],&param[14],
                           param[21]);
Output::print("Please adjust 'TurbulentViscosity3D to the turbulentViscosity3D'" 
   " and check the matrices as well espcially in NSTYPE 3 and 4");
  ErrThrow("not tested and/or adjusted yet: ");
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test_stream = Mult*(u1*test100+u2*test010+u3*test001);
    test100 *= Mult*c0;
    test010 *= Mult*c0;
    test001 *= Mult*c0;
    test000 = Mult*Orig3[i];

    for(j=0;j<N_U;j++)
    {      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz_stream = (u1*ansatz100+u2*ansatz010+u3*ansatz001); // sd direction
      val1 = ansatz_stream*test000;           // convective term
      ansatz_stream *= mu * test_stream;      // sd projection term 

      val2 = test100*ansatz100;
      val3 = test010*ansatz010;
      val4 = test001*ansatz001;
      val1 += val2+val3+val4;
      Matrix11Row[j] += val2+val1+ansatz_stream;
      Matrix12Row[j] += test010*ansatz100; 
      Matrix13Row[j] += test001*ansatz100;
      Matrix21Row[j] += test100*ansatz010;
      Matrix22Row[j] += val3+val1+ansatz_stream;
      Matrix23Row[j] += test001*ansatz010;
      Matrix31Row[j] += test100*ansatz001;
      Matrix32Row[j] += test010*ansatz001;
      Matrix33Row[j] += val4+val1+ansatz_stream;
    } // endfor j
  } // endfor i

  for(i=0;i<N_U;i++)
  {
     Matrix11Row = Matrix_tilde_G11[i];
     Matrix22Row = Matrix_tilde_G22[i];
     Matrix33Row = Matrix_tilde_G33[i];
     test100 = Orig0[i];
     test010 = Orig1[i];
     test001 = Orig2[i];
     val1 =u1*test100+u2*test010+u3*test001;

     for(j=0;j<N_L;j++)
     {       
        ansatz000 = Orig4[j];
        val = Mult * mu * ansatz000 * val1;
        Matrix11Row[j] -= val * u1;
        Matrix22Row[j] -= val * u2;
        Matrix33Row[j] -= val * u3;
     }
  }   
}

// ======================================================================
// Type 3, div-div stabilization
// Type 4, div-div stabilization
// ======================================================================
void TimeNSType3_4NLDivDivDD3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21;
  double **MatrixA22, **MatrixA23, **MatrixA31, **MatrixA32;
  double **MatrixA33;
  double val1;   // double val2, val3, val4;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double tautest100, tautest010, tautest001;
  double c0test100, c0test010, c0test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;  // double *Orig4;
  int i,j,N_U; // int N_P;
  double c0, tau;
  double u1, u2, u3; // double mu;
  double val;
  double theta1 = TDatabase::TimeDB->THETA1;
  
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  tau = DivDivStab3D(u1,u2,u3,hK,c0);
  //OutPut(" "<<tau);
  // for time-dependent problems, the div-div term is incorporated
  // into the matrix A which will be multiplied with theta1
  if (theta1 > 0)
      tau /= theta1;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    tautest100 = tau*test100;
    tautest010 = tau*test010;
    tautest001 = tau*test001;
    c0test100 = c0*test100;
    c0test010 = c0*test010;
    c0test001 = c0*test001;


    for(j=0;j<N_U;j++)
    {      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      // convection term
      val1 = (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;

      val  = c0*(2*test100*ansatz100+test010*ansatz010
                   +test001*ansatz001)+tautest100*ansatz100;
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0test010*ansatz100+tautest100*ansatz010;
      Matrix12Row[j] += Mult * val;

      val  = c0test001*ansatz100+tautest100*ansatz001;
      Matrix13Row[j] += Mult * val;

      val  = c0test100*ansatz010+tautest010*ansatz100;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+2*test010*ansatz010
                   +test001*ansatz001)+tautest010*ansatz010;
      val += val1;
      Matrix22Row[j] += Mult * val;

      val  = c0test001*ansatz010+tautest010*ansatz001;
      Matrix23Row[j] += Mult * val;

      val  = c0test100*ansatz001+tautest001*ansatz100;
      Matrix31Row[j] += Mult * val;

      val  = c0test010*ansatz001+tautest001*ansatz010;
      Matrix32Row[j] += Mult * val;

      val  = c0*(test100*ansatz100+test010*ansatz010
                   +2*test001*ansatz001)+tautest001*ansatz001;
      val += val1;
      Matrix33Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}

// ======================================================================
// Type 4, div-div, SUPG
// ======================================================================
void TimeNSType4Params_SUPG(double *in, double *out)
{
  out[0] = in[3];
  out[1] = in[4];
  out[2] = in[5];
  // u1old, u2old, u3old previous time 
  out[3] = in[6]; 
  out[4] = in[7]; 
  out[5] = in[8];
}
// ======================================================================
void TimeNSType4_SUPGDD3D(double Mult, double *coeff, double *param, double hK, 
     double **OrigValues, int *N_BaseFuncts,double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = LocMatrices[2];
  double **MatrixA21 = LocMatrices[3];
  double **MatrixA22 = LocMatrices[4];
  double **MatrixA23 = LocMatrices[5];
  double **MatrixA31 = LocMatrices[6];
  double **MatrixA32 = LocMatrices[7];
  double **MatrixA33 = LocMatrices[8];
  double **MassMatrix = LocMatrices[9];
  double **MatrixB1  = LocMatrices[10];
  double **MatrixB2  = LocMatrices[11];
  double **MatrixB3  = LocMatrices[12];
  double **MatrixB1T = LocMatrices[13];
  double **MatrixB2T = LocMatrices[14];
  double **MatrixB3T = LocMatrices[15];

  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_y
  double *Orig3 = OrigValues[3]; // u
  
  double *Orig4 = OrigValues[4]; // p
  double *Orig5 = OrigValues[5]; // p_x
  double *Orig6 = OrigValues[6]; // p_y
  double *Orig7 = OrigValues[7]; // p_z
  
  //double *Orig8 = OrigValues[8]; // u_xx
  //double *Orig9 = OrigValues[9]; // u_yy
  //double *Orig10 = OrigValues[10]; // u_yy

  double c0 = coeff[0]; // nu
  double c1 = coeff[1]; // f1
  double c2 = coeff[2]; // f2
  double c3 = coeff[3]; // f3

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old
  
  double u1_pt = param[3];
  double u2_pt = param[4];
  double u3_pt = param[5];

  double val;
  double test000, test100, test010, test001;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  //double ansatz200, ansatz020, ansatz002;
  //TODO: specify the parameter accordingly
  double tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  double tau_c = TDatabase::ParamDB->DELTA1;
  double tau_m_ugradv;
  
  for(int i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    tau_m_ugradv = tau_m*(u1_pt*test100 + u2_pt*test010 + u3_pt*test001);
    
    Rhs1[i] += Mult*(test000+tau_m_ugradv)*c1;
    Rhs2[i] += Mult*(test000+tau_m_ugradv)*c2;
    Rhs3[i] += Mult*(test000+tau_m_ugradv)*c3;

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      //ansatz200 = Orig8[j];
      //ansatz020 = Orig9[j];
      //ansatz002 = Orig10[j];
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*tau_m_ugradv; 
      val += tau_c * test100 * ansatz100;
      MatrixA11[i][j] += Mult * val;

      val  = c0*(test010*ansatz100) + tau_c * test100 * ansatz010;
      MatrixA12[i][j] += Mult * val;

      val  = c0*(test001*ansatz100) + tau_c * test100 * ansatz001;
      MatrixA13[i][j] += Mult * val;

      val  = c0*(test100*ansatz010) + tau_c * test010 * ansatz100;
      MatrixA21[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*tau_m_ugradv; 
      val += tau_c * test010 * ansatz010;
      MatrixA22[i][j] += Mult * val;

      val  = c0*(test001*ansatz010) + tau_c * test010 * ansatz001;
      MatrixA23[i][j] += Mult * val;

      val  = c0*(test100*ansatz001) + tau_c * test001 * ansatz100;
      MatrixA31[i][j] += Mult * val;

      val  = c0*(test010*ansatz001) + tau_c * test001 * ansatz010;
      MatrixA32[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*tau_m_ugradv;
      val += tau_c * test001 * ansatz001;
      MatrixA33[i][j] += Mult * val;
      // weighted mass matrix
      val = Mult*ansatz000*(test000 + tau_m_ugradv);
      MassMatrix[i][j] += val;
    } 
    for(int j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j]; // p
      ansatz100 = Orig5[j];
      ansatz010 = Orig6[j];
      ansatz001 = Orig7[j];
      // B1T
      val  = -ansatz000 * test100;
      val +=  ansatz100 * tau_m_ugradv;
      MatrixB1T[i][j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * tau_m_ugradv;
      MatrixB2T[i][j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * tau_m_ugradv;
      MatrixB3T[i][j] += Mult*val;
    }
  } 

  for(int i=0;i<N_P;i++)
  {
    test000 = Orig4[i];
    double val1 = Mult*test000;

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixB1[i][j] += val;
      
      val = -val1*ansatz010;
      MatrixB2[i][j] += val;
      
      val = -val1*ansatz001;
      MatrixB3[i][j] += val;
    } // endfor j
  } // endfor i
}

// ======================================================================
void TimeNSType4NL_SUPGDD3D(double Mult, double *coeff, double *param, double hK, 
         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA22 = LocMatrices[1];
  double **MatrixA33 = LocMatrices[2];
  
  int N_U = N_BaseFuncts[0];

  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_y
  double *Orig3 = OrigValues[3]; // u
  
  // double *Orig4 = OrigValues[4]; // u_xx
  // double *Orig5 = OrigValues[5]; // u_yy
  // double *Orig6 = OrigValues[6]; // u_zz

  double c0 = coeff[0]; // nu

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old

  double u1_pt = param[3];
  double u2_pt = param[4];
  double u3_pt = param[5];

  double val;
  double test000, test100, test010, test001;
  double ansatz100, ansatz010, ansatz001;
  // double ansatz200, ansatz020, ansatz002;
  //TODO: specify the parameter accordingly
  double tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  double tau_c = TDatabase::ParamDB->DELTA1;
  double tau_m_ugradv;
    
  for(int i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    tau_m_ugradv = tau_m*(u1_pt*test100 + u2_pt*test010 + u3_pt*test001);   

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      //ansatz000 = Orig3[j];
      
      // ansatz200 = Orig4[j];
      // ansatz020 = Orig5[j];
      // ansatz002 = Orig6[j];
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*tau_m_ugradv; 
      val += tau_c * test100 * ansatz100;
      MatrixA11[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*tau_m_ugradv;
      val += tau_c * test010 * ansatz010;
      MatrixA22[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*tau_m_ugradv; 
      val += tau_c * test001 * ansatz001;
      MatrixA33[i][j] += Mult * val;
    } 
  }
}

void TimeNSType4_SUPGExtraDD3D(double Mult, double* coeff, double* param, double hK, 
     double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA22 = LocMatrices[1];
  double **MatrixA33 = LocMatrices[2];
  double **MassMatrix = LocMatrices[3];
  double **MatrixB1T = LocMatrices[4];
  double **MatrixB2T = LocMatrices[5];
  double **MatrixB3T = LocMatrices[6];

  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_y
  double *Orig3 = OrigValues[3]; // u
  
  double *Orig4 = OrigValues[4]; // 
  double *Orig5 = OrigValues[5]; // 
  double *Orig6 = OrigValues[6]; //
  double *Orig7 = OrigValues[7]; // 
  
  // double *Orig8 = OrigValues[8]; // 
  // double *Orig9 = OrigValues[9]; // 
  // double *Orig10 = OrigValues[10]; // 
  
  

  double c0 = coeff[0]; // nu
  double c1 = coeff[1]; // f1
  double c2 = coeff[2]; // f2
  double c3 = coeff[3]; // f3

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old
  
  // double u1_pt = param[3];
  // double u2_pt = param[4];
  // double u3_pt = param[5];

  double val;
  double test000, test100, test010, test001;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  //double ansatz200, ansatz020, ansatz002;
  //TODO: specify the parameter accordingly
  double tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  double tau_c = TDatabase::ParamDB->DELTA1;
  double tau_m_ugradv;
  
  for(int i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    tau_m_ugradv = tau_m*(u1*test100 + u2*test010 + u3*test001);
    Rhs1[i] += Mult*(test000+tau_m_ugradv)*c1;
    Rhs2[i] += Mult*(test000+tau_m_ugradv)*c2;
    Rhs3[i] += Mult*(test000+tau_m_ugradv)*c3;
    
    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      //ansatz200 = Orig8[j];
      //ansatz020 = Orig9[j];
      //ansatz002 = Orig10[j];
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*tau_m_ugradv; 
      val += tau_c * test100 * ansatz100;
      MatrixA11[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*tau_m_ugradv; 
      val += tau_c * test010 * ansatz010;
      MatrixA22[i][j] += Mult * val;
     
      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*tau_m_ugradv; 
      val += tau_c * test001 * ansatz001;
      MatrixA33[i][j] += Mult * val;
      // weighted mass matrix
      val = Mult*ansatz000*(test000 + tau_m_ugradv);
      MassMatrix[i][j] += val;
    } 

    for(int j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j]; 
      ansatz100 = Orig5[j];
      ansatz010 = Orig6[j];
      ansatz001 = Orig7[j];
      // B1T
      double val  = -ansatz000 * test100;
      val +=  ansatz100 * tau_m_ugradv;
      MatrixB1T[i][j] += Mult*val;
      
      val  = -ansatz000 * test010;
      val +=  ansatz010 * tau_m_ugradv;
      MatrixB2T[i][j] += Mult*val;

      val  = -ansatz000 * test001;
      val +=  ansatz001 * tau_m_ugradv;
      MatrixB3T[i][j] += Mult*val;
    }
  } 
}

// ======================================================================
void TimeNSRhs_SUPGDD3D(double Mult, double *coeff, double *param, double hK, 
     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];
  
  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_z
  double *Orig3 = OrigValues[3]; // u
  
  double c1 = coeff[1];
  double c2 = coeff[2];
  double c3 = coeff[3];

  double u1=param[0]; // u1old
  double u2=param[1]; // u2old
  double u3=param[2]; // u3old
  
  int N_U = N_BaseFuncts[0];
  
  //TODO: specify the parameter accordingly
  double tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  
  for(int i=0; i<N_U; ++i)
  {
    double test100 = Orig0[i];
    double test010 = Orig1[i];
    double test001 = Orig2[i];
    double test000 = Orig3[i];
    
    double tau_m_ugradv = tau_m*(u1*test100 + u2*test010 + u3*test001);
    Rhs1[i] += Mult*(test000 + tau_m_ugradv)*c1;
    Rhs2[i] += Mult*(test000 + tau_m_ugradv)*c2;
    Rhs3[i] += Mult*(test000 + tau_m_ugradv)*c3;
  }
}
// ======================================================================
// Type 4 : Residual-based VMS
// assembling routine for residual-based VMS method of Bazilevs et al. (2007)
// references go to Ahmed et. al, Arch. Computat. Methods Eng. 24, 115 - 164 (2017)
// all matrices and right-hand side
// ======================================================================
void TimeNSType4Residual_VMSDD3D(double Mult, double* coeff, double* param, double hK, 
     double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs)
{ // A-blocks
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = LocMatrices[2];
  double **MatrixA21 = LocMatrices[3];
  double **MatrixA22 = LocMatrices[4];
  double **MatrixA23 = LocMatrices[5];
  double **MatrixA31 = LocMatrices[6];
  double **MatrixA32 = LocMatrices[7];
  double **MatrixA33 = LocMatrices[8];
  // M-blocks
  double **MatrixM11 = LocMatrices[9];
  double **MatrixM12 = LocMatrices[10];
  double **MatrixM13 = LocMatrices[11];
  double **MatrixM21 = LocMatrices[12];
  double **MatrixM22 = LocMatrices[13];
  double **MatrixM23 = LocMatrices[14];
  double **MatrixM31 = LocMatrices[15];
  double **MatrixM32 = LocMatrices[16];
  double **MatrixM33 = LocMatrices[17];
  // B, BT-blocks
  double **MatrixB1  = LocMatrices[18];
  double **MatrixB2  = LocMatrices[19];
  double **MatrixB3  = LocMatrices[20];
  double **MatrixB1T = LocMatrices[21];
  double **MatrixB2T = LocMatrices[22];
  double **MatrixB3T = LocMatrices[23];

  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_y
  double *Orig3 = OrigValues[3]; // u
  
  double *Orig4 = OrigValues[4]; // p
  double *Orig5 = OrigValues[5]; // p_x
  double *Orig6 = OrigValues[6]; // p_y
  double *Orig7 = OrigValues[7]; // p_z
  
  double *Orig8 = OrigValues[8]; // u_xx
  double *Orig9 = OrigValues[9]; // u_yy
  double *Orig10 = OrigValues[10]; // u_yy

  double c0 = coeff[0]; // nu
  double c1 = coeff[1]; // f1
  double c2 = coeff[2]; // f2
  double c3 = coeff[3]; // f3
  
  double c1_old = coeff[4]; // f1_previous time 
  double c2_old = coeff[5]; // f2_previous time 
  double c3_old = coeff[6]; // f3_previous time 

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old
  // additional terms for the residual computations
  double u1_min1 = param[3];
  double u2_min1 = param[4];
  double u3_min1 = param[5];
  double u1x = param[6];
  double u2x = param[7];
  double u3x = param[8];
  double u1y = param[9];
  double u2y = param[10];
  double u3y = param[11];
  double u1z = param[12];
  double u2z = param[13];
  double u3z = param[14];
  double u1xx = param[15];
  double u2xx = param[16];
  double u3xx = param[17];
  double u1yy = param[18];
  double u2yy = param[19];
  double u3yy = param[20];
  double u1zz = param[21];
  double u2zz = param[22];
  double u3zz = param[23];
  double px = param[24];
  double py = param[25];
  double pz = param[26];
  
  double u1_min2 = param[27]; // previous time solution u1
  double u2_min2 = param[28]; // previous time solution u2
  double u3_min2 = param[29]; // previous time solution u3
  
  double val;
  double test000, test100, test010, test001;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double ansatz200, ansatz020, ansatz002;

  double tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  double tau_c = TDatabase::ParamDB->DELTA1;

  double tau_m_ugradv;
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
   
  // old residual to be used Eq(51)
  double res1 = tau_m *(c1_old - 1./dt*(u1_min1-u1_min2) + c0*(u1xx + u1yy + u1zz)
                          -(u1*u1x + u2*u1y + u3*u1z) -px );
  double res2 = tau_m *(c2_old - 1./dt*(u2_min1-u2_min2) + c0*(u2xx + u2yy + u2zz)
                          -(u1*u2x + u2*u2y + u3*u2z) -py );
  double res3 = tau_m *(c3_old - 1./dt*(u3_min1-u3_min2) + c0*(u3xx + u3yy + u3zz)
                          -(u1*u3x + u2*u3y + u3*u3z) -pz );

  
  for(int i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    // streamline derivative with extrapolated velocity from previous time
    tau_m_ugradv = tau_m*(u1_min1*test100 + u2_min1*test010 + u3_min1*test001);
    // contribution to rhs from Galerkin discretization and SUPG term 
    Rhs1[i] += Mult*(test000+tau_m_ugradv)*c1;
    Rhs2[i] += Mult*(test000+tau_m_ugradv)*c2;
    Rhs3[i] += Mult*(test000+tau_m_ugradv)*c3;
    
    // contribution from second nonlinear term
    Rhs1[i] += Mult*tau_m*u1_min1*(c1*test100 + c2*test010 + c3*test001);
    Rhs2[i] += Mult*tau_m*u2_min1*(c1*test100 + c2*test010 + c3*test001);
    Rhs3[i] += Mult*tau_m*u3_min1*(c1*test100 + c2*test010 + c3*test001);
    
    // contribution from third nonlinear term 
    Rhs1[i] += Mult*tau_m*res1*(c1*test100 + c2*test010 + c3*test001);
    Rhs2[i] += Mult*tau_m*res2*(c1*test100 + c2*test010 + c3*test001);
    Rhs3[i] += Mult*tau_m*res3*(c1*test100 + c2*test010 + c3*test001);

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      ansatz200 = Orig8[j];
      ansatz020 = Orig9[j];
      ansatz002 = Orig10[j];
      
      double laplacian = -c0*(ansatz200 + ansatz020 + ansatz002);
      double ugradu = (u1*ansatz100 + u2*ansatz010 + u3*ansatz001);
      
      // stiffness matrix blocks
      // viscous term (deformation tensor)
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      // convective term 
      val += ugradu*test000;
      // velocity contribution of SUPG term 
      val += (laplacian + ugradu)*tau_m_ugradv;
      // grad-div term 
      val += tau_c * test100 * ansatz100;
      // FOR EFFICIENCY: MERGE SEOCND CROSS AND SUBGRID TERM 
      // second cross term 
      val += tau_m * (laplacian + ugradu) * u1_min1 * test100;
      // subgrid term 
      val += tau_m * (laplacian + ugradu) * res1 * test100;
      MatrixA11[i][j] += Mult * val;

      // viscous term and grad-div term 
      val  = c0*(test010*ansatz100) + tau_c * test100 * ansatz010;
      // second cross term
      val += tau_m * ugradu * u1_min1 * test010;
      // subgrid term 
      val += tau_m * ugradu * res1 * test010;
      MatrixA12[i][j] += Mult * val;

      // viscous term and grad-div term 
      val  = c0*(test001*ansatz100) + tau_c * test100 * ansatz001;
      // second cross term
      val += tau_m * ugradu * u1_min1 * test001;
      // subgrid term 
      val += tau_m * ugradu * res1 * test001;
      MatrixA13[i][j] += Mult * val;

      val  = c0*(test100*ansatz010) + tau_c * test010 * ansatz100;
      val += tau_m * ugradu * u2_min1 * test100;
      val += tau_m * ugradu * res2 * test100;
      MatrixA21[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += ugradu*test000;
      val += ugradu*tau_m_ugradv; 
      val += tau_c * test010 * ansatz010;
      val += tau_m * (laplacian + ugradu) * u2_min1 * test010;
      val += tau_m * (laplacian + ugradu) * res2 * test010;
      MatrixA22[i][j] += Mult * val;

      val  = c0*(test001*ansatz010) + tau_c * test010 * ansatz001;
      val += tau_m * ugradu * u2_min1 * test001;
      val += tau_m * ugradu * res2 * test001;
      MatrixA23[i][j] += Mult * val;

      val  = c0*(test100*ansatz001) + tau_c * test001 * ansatz100;
      val += tau_m * ugradu * u3_min1 * test100;
      val += tau_m * ugradu * res3 * test100;
      MatrixA31[i][j] += Mult * val;

      val  = c0*(test010*ansatz001) + tau_c * test001 * ansatz010;
      val += tau_m * ugradu * u3_min1 * test010;
      val += tau_m * ugradu * res3 * test010;
      MatrixA32[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += ugradu*test000;
      val += (laplacian +ugradu)*tau_m_ugradv;
      val += tau_c * test001 * ansatz001;
      val += tau_m * (laplacian + ugradu) * u3_min1 * test001;
      val += tau_m * (laplacian + ugradu) * res3 * test001;
      MatrixA33[i][j] += Mult * val;
      
      // mass matrix blocks
      // Galerkin + SUPG term 
      val = ansatz000*(test000 + tau_m_ugradv);
      // second cross term 
      val += tau_m * u1_min1   * ansatz000 * test100;
      // subgrid term 
      val += tau_m * res1 * ansatz000 * test100;
      MatrixM11[i][j] += Mult * val;

      val = tau_m * u1_min1   * ansatz000 * test010;
      val += tau_m * res1 * ansatz000 * test010;
      MatrixM12[i][j] += Mult*val;
      
      val = tau_m * u1_min1   * ansatz000 * test001;
      val += tau_m * res1 * ansatz000 * test001;
      MatrixM13[i][j] += Mult*val;
      
      val = tau_m * u2_min1 * ansatz000 * test100;
      val += tau_m * res2 * ansatz000 * test100;
      MatrixM21[i][j] += Mult * val;
      
      val = ansatz000*(test000 + tau_m_ugradv);
      val += tau_m * u2_min1 * ansatz000 * test010;
      val += tau_m * res2 * ansatz000 * test010;
      MatrixM22[i][j] += Mult * val;
      
      val = tau_m * u2_min1   * ansatz000 * test001;
      val += tau_m * res2 * ansatz000 * test001;
      MatrixM23[i][j] += Mult * val ;
      
      val = tau_m * u3_min1   * ansatz000 * test100;
      val += tau_m * res3 * ansatz000 * test100;
      MatrixM31[i][j] += Mult * val;
      
      val = tau_m * u3_min1   * ansatz000 * test010;
      val += tau_m * res3 * ansatz000 * test010;
      MatrixM32[i][j] += Mult * val;
      
      val = ansatz000*(test000 + tau_m_ugradv);
      val += tau_m * u3_min1 * ansatz000 * test001;
      val += tau_m * res3 * ansatz000 * test001;
      MatrixM33[i][j] += Mult * val;
    }

    // coupling pressure (ansatz) - velocity (test)
    for(int j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j]; // p
      ansatz100 = Orig5[j];
      ansatz010 = Orig6[j];
      ansatz001 = Orig7[j];
      
      // B1T
      // Galerkin
      val  = -ansatz000 * test100;
      // SUPG
      val +=  ansatz100 * tau_m_ugradv;
      // second cross term 
      val += tau_m * u1_min1   * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      // subgrid term 
      val += tau_m * res1 * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixB1T[i][j] += Mult*val;
      
      // B2T 
      val  = -ansatz000 * test010;
      val +=  ansatz010 * tau_m_ugradv;
      val += tau_m * u2_min1 * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      val += tau_m * res2 * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixB2T[i][j] += Mult*val;

      // B3T 
      val  = -ansatz000 * test001;
      val +=  ansatz001 * tau_m_ugradv;
      val += tau_m * u3_min1 * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      val += tau_m * res3 * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
      MatrixB3T[i][j] += Mult*val;
    }
  } 

  // coupling velocity (ansatz) - pressure (test)
  for(int i=0;i<N_P;i++)
  {
    test000 = Orig4[i];
    double val1 = Mult*test000;

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];

      val = -val1*ansatz100;
      MatrixB1[i][j] += val;
      
      val = -val1*ansatz010;
      MatrixB2[i][j] += val;
      
      val = -val1*ansatz001;
      MatrixB3[i][j] += val;
    } // endfor j
  } // endfor i
}
// ======================================================================
void TimeNSType4NLResidual_VMSDD3D(double Mult, double *coeff, double *param, 
 double hK, double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
 double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = LocMatrices[2];
  double **MatrixA21 = LocMatrices[3];
  double **MatrixA22 = LocMatrices[4];
  double **MatrixA23 = LocMatrices[5];
  double **MatrixA31 = LocMatrices[6];
  double **MatrixA32 = LocMatrices[7];
  double **MatrixA33 = LocMatrices[8];

  int N_U = N_BaseFuncts[0];

  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_y
  double *Orig3 = OrigValues[3]; // u
  
  double *Orig4 = OrigValues[4]; // u_xx
  double *Orig5 = OrigValues[5]; // u_yy
  double *Orig6 = OrigValues[6]; // u_yy

  double c0 = coeff[0]; // nu
  
  double c1_old = coeff[4]; // f1_previous time 
  double c2_old = coeff[5]; // f2_previous time 
  double c3_old = coeff[6]; // f3_previous time 

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old
  // additional terms for the residual computations
  double u1_min1 = param[3];
  double u2_min1 = param[4];
  double u3_min1 = param[5];
  double u1x = param[6];
  double u2x = param[7];
  double u3x = param[8];
  double u1y = param[9];
  double u2y = param[10];
  double u3y = param[11];
  double u1z = param[12];
  double u2z = param[13];
  double u3z = param[14];
  double u1xx = param[15];
  double u2xx = param[16];
  double u3xx = param[17];
  double u1yy = param[18];
  double u2yy = param[19];
  double u3yy = param[20];
  double u1zz = param[21];
  double u2zz = param[22];
  double u3zz = param[23];
  double px = param[24];
  double py = param[25];
  double pz = param[26];
  
  double u1_min2 = param[27]; // previous time solution u1
  double u2_min2 = param[28]; // previous time solution u2
  double u3_min2 = param[29]; // previous time solution u3

//   for(int i=0; i<30; i++)
//     cout<< param[i] << "  ";
//   cout << endl;
  double val;
  double test000, test100, test010, test001;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double ansatz200, ansatz020, ansatz002;
  //TODO: specify the parameter accordingly
  double tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  double tau_c = TDatabase::ParamDB->DELTA1;
  double tau_m_ugradv;
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  
  for(int i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    tau_m_ugradv = tau_m*(u1_min1*test100 + u2_min1*test010 + u3_min1*test001);
    
    // old residual in Eq(51)
    double res1 = tau_m *(c1_old - 1./dt*(u1_min1-u1_min2) + c0*(u1xx + u1yy + u1zz)
                          -(u1*u1x + u2*u1y + u3*u1z) -px );
    double res2 = tau_m *(c2_old - 1./dt*(u2_min1-u2_min2) + c0*(u2xx + u2yy + u2zz)
                          -(u1*u2x + u2*u2y + u3*u2z) -py );
    double res3 = tau_m *(c3_old - 1./dt*(u3_min1-u3_min2) + c0*(u3xx + u3yy + u3zz)
                          -(u1*u3x + u2*u3y + u3*u3z) -pz );
  
    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      ansatz200 = Orig4[j];
      ansatz020 = Orig5[j];
      ansatz002 = Orig6[j];
      
      double laplacian = -c0*(ansatz200 + ansatz020 + ansatz002);
      double ugradu = (u1*ansatz100 + u2*ansatz010 + u3*ansatz001);
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += ugradu*test000;
      val += (laplacian + ugradu)*tau_m_ugradv;
      val += tau_c * test100 * ansatz100;
      val += tau_m * (laplacian + ugradu) * u1_min1 * test100;
      val += tau_m * (laplacian + ugradu) * res1 * test100;
      MatrixA11[i][j] += Mult * val;

      val  = c0*(test010*ansatz100) + tau_c * test100 * ansatz010;
      val += tau_m * ugradu * u1_min1 * test010;
      val += tau_m * ugradu * res1 * test010;
      MatrixA12[i][j] += Mult * val;

      val  = c0*(test001*ansatz100) + tau_c * test100 * ansatz001;
      val += tau_m * ugradu * u1_min1 * test001;
      val += tau_m * ugradu * res1 * test001;
      MatrixA13[i][j] += Mult * val;

      val  = c0*(test100*ansatz010) + tau_c * test010 * ansatz100;
      val += tau_m * ugradu * u2_min1 * test100;
      val += tau_m * ugradu * res2 * test100;
      MatrixA21[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += ugradu*test000;
      val += (laplacian + ugradu)*tau_m_ugradv; 
      val += tau_c * test010 * ansatz010;
      val += tau_m * (laplacian + ugradu) * u2_min1 * test010;
      val += tau_m * (laplacian + ugradu) * res2 * test010;
      MatrixA22[i][j] += Mult * val;

      val  = c0*(test001*ansatz010) + tau_c * test010 * ansatz001;
      val += tau_m * ugradu * u2_min1 * test001;
      val += tau_m * ugradu * res2 * test001;
      MatrixA23[i][j] += Mult * val;

      val  = c0*(test100*ansatz001) + tau_c * test001 * ansatz100;
      val += tau_m * ugradu * u3_min1 * test100;
      val += tau_m * ugradu * res3 * test100;
      MatrixA31[i][j] += Mult * val;

      val  = c0*(test010*ansatz001) + tau_c * test001 * ansatz010;
      val += tau_m * ugradu * u3_min1 * test010;
      val += tau_m * ugradu * res3 * test010;
      MatrixA32[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += ugradu*test000;
      val += (laplacian +ugradu)*tau_m_ugradv;
      val += tau_c * test001 * ansatz001;
      val += tau_m * (laplacian + ugradu) * u3_min1 * test001;
      val += tau_m * (laplacian + ugradu) * res3 * test001;
      MatrixA33[i][j] += Mult * val;
    }
  }
}
// ======================================================================
void TimeNSType4Residual_VMS_RhsDD3D(double Mult, double* coeff, double* param, 
 double hK, double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, 
 double** LocRhs)
{
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];

  int N_U = N_BaseFuncts[0];

  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_y
  double *Orig3 = OrigValues[3]; // u
  
  double c0 = coeff[0]; // nu
  double c1 = coeff[1]; // f1
  double c2 = coeff[2]; // f2
  double c3 = coeff[3]; // f3
  
  double c1_old = coeff[4]; // f1_previous time 
  double c2_old = coeff[5]; // f2_previous time 
  double c3_old = coeff[6]; // f3_previous time 

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old
  // additional terms for the residual computations
  double u1m1 = param[3];
  double u2m1 = param[4];
  double u3m1 = param[5];
  double u1x = param[6];
  double u2x = param[7];
  double u3x = param[8];
  double u1y = param[9];
  double u2y = param[10];
  double u3y = param[11];
  double u1z = param[12];
  double u2z = param[13];
  double u3z = param[14];
  double u1xx = param[15];
  double u2xx = param[16];
  double u3xx = param[17];
  double u1yy = param[18];
  double u2yy = param[19];
  double u3yy = param[20];
  double u1zz = param[21];
  double u2zz = param[22];
  double u3zz = param[23];
  double px = param[24];
  double py = param[25];
  double pz = param[26];
  
  double u1m2 = param[27]; // previous time solution u1
  double u2m2 = param[28]; // previous time solution u2
  double u3m2 = param[29]; // previous time solution u3

  double test000, test100, test010, test001;
  // double ansatz200, ansatz020, ansatz002;
  //TODO: specify the parameter accordingly
  double tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  double tau_m_ugradv;
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  
  for(int i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    tau_m_ugradv = tau_m*(u1m1*test100 + u2m1*test010 + u3m1*test001);
    Rhs1[i] += Mult*(test000+tau_m_ugradv)*c1;
    Rhs2[i] += Mult*(test000+tau_m_ugradv)*c2;
    Rhs3[i] += Mult*(test000+tau_m_ugradv)*c3;
    
    // old residual in Eq(51)
    double res1=tau_m*(c1_old-1./dt*(u1m1-u1m2)+c0*(u1xx+u1yy+u1zz)-(u1m1*u1x+u2m1*u1y+u3m1*u1z)-px);
    double res2=tau_m*(c2_old-1./dt*(u2m1-u2m2)+c0*(u2xx+u2yy+u2zz)-(u1m1*u2x+u2m1*u2y+u3m1*u2z)-py);
    double res3=tau_m*(c3_old-1./dt*(u3m1-u3m2)+c0*(u3xx+u3yy+u3zz)-(u1m1*u3x+u2m1*u3y+u3m1*u3z)-pz);
    //cout <<"rhs: " <<  res1 << "  " << res2 << "  " << res3 << endl;
    // contribution from second nonlinear term
    Rhs1[i] += Mult*tau_m*u1m1*(c1*test100 + c2*test010 + c3*test001);
    Rhs2[i] += Mult*tau_m*u2m1*(c1*test100 + c2*test010 + c3*test001);
    Rhs3[i] += Mult*tau_m*u3m1*(c1*test100 + c2*test010 + c3*test001);
    // contribution from third nonlinear term 
    Rhs1[i] += Mult*tau_m*res1*(c1*test100 + c2*test010 + c3*test001);
    Rhs2[i] += Mult*tau_m*res2*(c1*test100 + c2*test010 + c3*test001);
    Rhs3[i] += Mult*tau_m*res3*(c1*test100 + c2*test010 + c3*test001);
  }
}
// ======================================================================
void TimeNSType4_Residual_VMS_ExtraDD3D(double Mult, double* coeff, 
 double* param, double hK, double** OrigValues, int* N_BaseFuncts, 
 double*** LocMatrices, double** LocRhs)
{
  // A-blocks
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = LocMatrices[2];
  double **MatrixA21 = LocMatrices[3];
  double **MatrixA22 = LocMatrices[4];
  double **MatrixA23 = LocMatrices[5];
  double **MatrixA31 = LocMatrices[6];
  double **MatrixA32 = LocMatrices[7];
  double **MatrixA33 = LocMatrices[8];
  // M-blocks
  double **MatrixM11 = LocMatrices[9];
  double **MatrixM12 = LocMatrices[10];
  double **MatrixM13 = LocMatrices[11];
  double **MatrixM21 = LocMatrices[12];
  double **MatrixM22 = LocMatrices[13];
  double **MatrixM23 = LocMatrices[14];
  double **MatrixM31 = LocMatrices[15];
  double **MatrixM32 = LocMatrices[16];
  double **MatrixM33 = LocMatrices[17];
  // BT-blocks
  double **MatrixB1T = LocMatrices[18];
  double **MatrixB2T = LocMatrices[19];
  double **MatrixB3T = LocMatrices[20];

  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = LocRhs[2];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0]; // u_x
  double *Orig1 = OrigValues[1]; // u_y
  double *Orig2 = OrigValues[2]; // u_y
  double *Orig3 = OrigValues[3]; // u
  
  double *Orig4 = OrigValues[4]; // p
  double *Orig5 = OrigValues[5]; // p_x
  double *Orig6 = OrigValues[6]; // p_y
  double *Orig7 = OrigValues[7]; // p_z
  
  double *Orig8 = OrigValues[8]; // u_xx
  double *Orig9 = OrigValues[9]; // u_yy
  double *Orig10 = OrigValues[10]; // u_yy

  double c0 = coeff[0]; // nu
  double c1 = coeff[1]; // f1
  double c2 = coeff[2]; // f2
  double c3 = coeff[3]; // f3

  double u1 = param[0]; // u1old
  double u2 = param[1]; // u2old
  double u3 = param[2]; // u3old
  // additional terms for the residual computations
  double u1x = param[3];
  double u2x = param[4];
  double u3x = param[5];
  double u1y = param[6];
  double u2y = param[7];
  double u3y = param[8];
  double u1z = param[9];
  double u2z = param[10];
  double u3z = param[11];
  double u1xx = param[12];
  double u2xx = param[13];
  double u3xx = param[14];
  double u1yy = param[15];
  double u2yy = param[16];
  double u3yy = param[17];
  double u1zz = param[18];
  double u2zz = param[19];
  double u3zz = param[20];
  double px = param[21];
  double py = param[22];
  double pz = param[23];
  double u1m1 = param[24]; // previous time solution u1
  double u2m1 = param[25]; // previous time solution u2
  double u3m1 = param[26]; // previous time solution u3

  double val;
  double test000, test100, test010, test001;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double ansatz200, ansatz020, ansatz002;
  //TODO: specify the parameter accordingly
  double tau_m = TDatabase::ParamDB->DELTA0*hK*hK;
  double tau_c = TDatabase::ParamDB->DELTA1;
  double tau_m_ugradv;
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  for(int i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    tau_m_ugradv = tau_m*(u1*test100 + u2*test010 + u3*test001);
    Rhs1[i] += Mult*(test000+tau_m_ugradv)*c1;
    Rhs2[i] += Mult*(test000+tau_m_ugradv)*c2;
    Rhs3[i] += Mult*(test000+tau_m_ugradv)*c3;
    
    // old residual in Eq(51)
    double res1=tau_m*(c1-1./dt*(u1-u1m1)+c0*(u1xx+u1yy+u1zz)-(u1*u1x+u2*u1y+u3*u1z)-px);
    double res2=tau_m*(c2-1./dt*(u2-u2m1)+c0*(u2xx+u2yy+u2zz)-(u1*u2x+u2*u2y+u3*u2z)-py);
    double res3=tau_m*(c3-1./dt*(u3-u3m1)+c0*(u3xx+u3yy+u3zz)-(u1*u3x+u2*u3y+u3*u3z)-pz);
    
    // contribution from second nonlinear term
    Rhs1[i] += Mult*tau_m*u1*(c1*test100 + c2*test010 + c3*test001);
    Rhs2[i] += Mult*tau_m*u2*(c1*test100 + c2*test010 + c3*test001);
    Rhs3[i] += Mult*tau_m*u3*(c1*test100 + c2*test010 + c3*test001);
    // contribution from third nonlinear term 
    Rhs1[i] += Mult*tau_m*res1*(c1*test100 + c2*test010 + c3*test001);
    Rhs2[i] += Mult*tau_m*res2*(c1*test100 + c2*test010 + c3*test001);
    Rhs3[i] += Mult*tau_m*res3*(c1*test100 + c2*test010 + c3*test001);

    for(int j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      ansatz200 = Orig8[j];
      ansatz020 = Orig9[j];
      ansatz002 = Orig10[j];
      
      double laplacian = -c0*(ansatz200 + ansatz020 + ansatz002);
      double ugradu = (u1*ansatz100 + u2*ansatz010 + u3*ansatz001);
      
      val  = 2*c0*(test100*ansatz100+0.5*test010*ansatz010
                   +0.5*test001*ansatz001);
      val += ugradu*test000;
      val += (laplacian + ugradu)*tau_m_ugradv;
      val += tau_c * test100 * ansatz100;
      val += tau_m * (laplacian + ugradu) * u1 * test100;
      val += tau_m * (laplacian + ugradu) * res1 * test100;
      MatrixA11[i][j] += Mult * val;

      val  = c0*(test010*ansatz100) + tau_c * test100 * ansatz010;
      val += tau_m * ugradu * u1 * test010;
      val += tau_m * ugradu * res1 * test010;
      MatrixA12[i][j] += Mult * val;

      val  = c0*(test001*ansatz100) + tau_c * test100 * ansatz001;
      val += tau_m * ugradu * u1 * test001;
      val += tau_m * ugradu * res1 * test001;
      MatrixA13[i][j] += Mult * val;

      val  = c0*(test100*ansatz010) + tau_c * test010 * ansatz100;
      val += tau_m * ugradu * u2 * test100;
      val += tau_m * ugradu * res2 * test100;
      MatrixA21[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+test010*ansatz010
                   +0.5*test001*ansatz001);
      val += ugradu*test000;
      val += ugradu*tau_m_ugradv; 
      val += tau_c * test010 * ansatz010;
      val += tau_m * (laplacian + ugradu) * u2 * test010;
      val += tau_m * (laplacian + ugradu) * res2 * test010;
      MatrixA22[i][j] += Mult * val;

      val  = c0*(test001*ansatz010) + tau_c * test010 * ansatz001;
      val += tau_m * ugradu * u2 * test001;
      val += tau_m * ugradu * res2 * test001;
      MatrixA23[i][j] += Mult * val;

      val  = c0*(test100*ansatz001) + tau_c * test001 * ansatz100;
      val += tau_m * ugradu * u3 * test100;
      val += tau_m * ugradu * res3 * test100;
      MatrixA31[i][j] += Mult * val;

      val  = c0*(test010*ansatz001) + tau_c * test001 * ansatz010;
      val += tau_m * ugradu * u3 * test010;
      val += tau_m * ugradu * res3 * test010;
      MatrixA32[i][j] += Mult * val;

      val  = 2*c0*(0.5*test100*ansatz100+0.5*test010*ansatz010
                   +test001*ansatz001);
      val += ugradu*test000;
      val += (laplacian +ugradu)*tau_m_ugradv;
      val += tau_c * test001 * ansatz001;
      val += tau_m * (laplacian + ugradu) * u3 * test001;
      val += tau_m * (laplacian + ugradu) * res3 * test001;
      MatrixA33[i][j] += Mult * val;
      
      // weighted mass matrix (galerkin + supg ) terms
      val = ansatz000*(test000 + tau_m_ugradv);
      val += tau_m * u1   * ansatz000 * test100;
      val += tau_m * res1 * ansatz000 * test100;
      MatrixM11[i][j] += Mult * val;

      val = tau_m * u1   * ansatz000 * test010;
      val += tau_m * res1 * ansatz000 * test010;
      MatrixM12[i][j] += Mult*val;
      val = tau_m * u1   * ansatz000 * test001;
      val += tau_m * res1 * ansatz000 * test001;
      MatrixM13[i][j] += Mult*val;
      
      val = tau_m * u2 * ansatz000 * test100;
      val += tau_m * res2 * ansatz000 * test100;
      MatrixM21[i][j] += Mult * val;
      
      val = ansatz000*(test000 + tau_m_ugradv);
      val += tau_m * u2 * ansatz000 * test010;
      val += tau_m * res2 * ansatz000 * test010;
      MatrixM22[i][j] += Mult * val;
      
      val = tau_m * u2   * ansatz000 * test001;
      val += tau_m * res2 * ansatz000 * test001;
      MatrixM23[i][j] += Mult * val ;
      
      val = tau_m * u3   * ansatz000 * test100;
      val += tau_m * res3 * ansatz000 * test100;
      MatrixM31[i][j] += Mult * val;
      
      val = tau_m * u3   * ansatz000 * test010;
      val += tau_m * res3 * ansatz000 * test010;
      MatrixM32[i][j] += Mult * val;
      
      val = ansatz000*(test000 + tau_m_ugradv);
      val += tau_m * u3 * ansatz000 * test001;
      val += tau_m * res3 * ansatz000 * test001;
      MatrixM33[i][j] += Mult * val;
    } 

    for(int j=0;j<N_P;j++)
    {
      ansatz000 = Orig4[j]; // p
      ansatz100 = Orig5[j];
      ansatz010 = Orig6[j];
      ansatz001 = Orig7[j];
      // B1T
      val  = -ansatz000*test100;
      val += ansatz100 *tau_m_ugradv;
      val += tau_m*u1  *(ansatz100*test100+ansatz010*test010+ansatz001*test001);
      val += tau_m*res1*(ansatz100*test100+ansatz010*test010+ansatz001*test001);
      MatrixB1T[i][j] += Mult*val;
      
      // B2T 
      val  = -ansatz000* test010;
      val +=  ansatz010* tau_m_ugradv;
      val += tau_m*u2  *(ansatz100*test100+ansatz010*test010+ansatz001*test001);
      val += tau_m*res2*(ansatz100*test100+ansatz010*test010+ansatz001*test001);
      MatrixB2T[i][j] += Mult*val;

      // B3T 
      val  = -ansatz000* test001;
      val +=  ansatz001* tau_m_ugradv;
      val += tau_m*u3  *(ansatz100*test100+ansatz010*test010+ansatz001*test001);
      val += tau_m*res3*(ansatz100*test100+ansatz010*test010+ansatz001*test001);
      MatrixB3T[i][j] += Mult*val;
    }
  } 
}

// ======================================================================
void TimeNSType4Params_Residual_VMS(double *in, double *out)
{
  out[0] = in[3];
  out[1] = in[4];
  out[2] = in[5];
  // u1x, u2x, u3x old
  out[3] = in[6]; 
  out[4] = in[7]; 
  out[5] = in[8]; 
  // u1y. u2y, u3y old
  out[6] = in[9]; 
  out[7] = in[10]; 
  out[8] = in[11]; 
  // u1z, u2z, u3z, old
  out[9] = in[12]; 
  out[10] = in[13]; 
  out[11] = in[14]; 
  // u1xx, u2xx, u3xx old
  out[12] = in[15]; 
  out[13] = in[16]; 
  out[14] = in[17]; 
  // u1yy, u2yy, u3yy, old
  out[15] = in[18]; 
  out[16] = in[19]; 
  out[17] = in[20]; 
  // u1zz, u2zz, u3zz, old
  out[18] = in[21];  
  out[19] = in[22];  
  out[20] = in[23];  
  // p_x, p_y, p_z
  out[21] = in[24]; 
  out[22] = in[25]; 
  out[23] = in[26]; 
  // u1, u2, u3, previous time sol's
  out[24] = in[27]; 
  out[25] = in[28]; 
  out[26] = in[29]; 
}

void TimeNSType4Params_Residual_VMS_Extrapolate(double* in, double* out)
{
  out[0] = in[3];
  out[1] = in[4];
  out[2] = in[5];
  // previous time solution and derivatives
  out[3] = in[6]; 
  out[4] = in[7]; 
  out[5] = in[8]; 

  out[6] = in[9]; 
  out[7] = in[10]; 
  out[8] = in[11]; 

  out[9] = in[12]; 
  out[10] = in[13]; 
  out[11] = in[14]; 

  out[12] = in[15]; 
  out[13] = in[16]; 
  out[14] = in[17]; 

  out[15] = in[18]; 
  out[16] = in[19]; 
  out[17] = in[20]; 

  out[18] = in[21];  
  out[19] = in[22];  
  out[20] = in[23];  

  out[21] = in[24]; 
  out[22] = in[25]; 
  out[23] = in[26]; 
  
  // p_x, p_y, p_z
  out[24] = in[27]; 
  out[25] = in[28]; 
  out[26] = in[29]; 
  // u1, u2, u3, previous time sol's
  out[27] = in[30]; 
  out[28] = in[31]; 
  out[29] = in[32]; 
}


// ======================================================================
// ROSENBROCK
// ======================================================================
void TimeNSType1GalerkinJ3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  Output::print("Nothing has been tested yer: !! :( ");
  ErrThrow("not tested and/or adjusted yet: ");
  
  double **MatrixA;
  double val;
  double *MatrixRow;  // double *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz100, ansatz010, ansatz001;  // douible ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U;  // int N_P;
  double c0;
  double u1, u2, u3;
  
  MatrixA = LocMatrices[0];
     
  N_U = N_BaseFuncts[0];
  //N_P = N_BaseFuncts[1];
  
  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u
     
  c0 = coeff[0]; // nu
  
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+
                 test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      MatrixRow[j] += Mult * val;
    } // endfor j
  } // endfor i
}

void TimeNSGalerkinC3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  Output::print("Nothing has been tested yer: !! :( ");
  ErrThrow("not tested and/or adjusted yet: ");
  
  double *Rhs1, *Rhs2, *Rhs3;
  double test000;
  double *Orig3;
  int i, N_U;
  double c4, c5, c6;
  
  //cout << "TimeNSType1GalerkinC" << endl;
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  
  N_U = N_BaseFuncts[0];
  
  Orig3 = OrigValues[0]; // u
     
  c4 = coeff[4]; // dot f1
  c5 = coeff[5]; // dot f2
  c6 = coeff[6]; // dot f3
  
  for(i=0;i<N_U;i++)
  {
    test000 = Orig3[i];

    Rhs1[i] += Mult*test000*c4;
    Rhs2[i] += Mult*test000*c5;
    Rhs3[i] += Mult*test000*c6;
  } 
}
void TimeNSType3GalerkinJ3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  Output::print("Nothing has been tested yer: !! :( ");
  ErrThrow("not tested and/or adjusted yet: ");
  
  double **MatrixA11, **MatrixA12, **MatrixA13;
  double **MatrixA21, **MatrixA22, **MatrixA23;
  double **MatrixA31, **MatrixA32, **MatrixA33;
  double val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row;
  double *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row;
  double *Matrix33Row;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U;  // int N_P;
  double c0;
  double u1, u2, u3, u1_x, u1_y, u1_z, u2_x, u2_y, u2_z, u3_x, u3_y, u3_z;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];
  MatrixA31 = LocMatrices[6];
  MatrixA32 = LocMatrices[7];
  MatrixA33 = LocMatrices[8];
  
  N_U = N_BaseFuncts[0];
  //N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_y
  Orig3 = OrigValues[3]; // u
  //Orig4 = OrigValues[4]; // l

  c0 = coeff[0]; // nu
  
  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u2old
  u1_x = param[3]; // u1old
  u2_x = param[4]; // u2old
  u3_x = param[5]; // u2old
  u1_y = param[6]; // u1old
  u2_y = param[7]; // u2old
  u3_y = param[8]; // u2old
  u1_z = param[9]; // u1old
  u2_z = param[10]; // u2old
  u3_z = param[11]; // u2old
  

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix13Row = MatrixA13[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    Matrix23Row = MatrixA23[i];
    Matrix31Row = MatrixA31[i];
    Matrix32Row = MatrixA32[i];
    Matrix33Row = MatrixA33[i];
    
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];
      
      val  = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += u1_x*ansatz000*test000;
      Matrix11Row[j] += Mult * val;

      val  = u1_y*ansatz000*test000;
      Matrix12Row[j] += Mult * val;
      
      val  = u1_z*ansatz000*test000;
      Matrix13Row[j] += Mult * val;
      

      val  = u2_x*ansatz000*test000;
      Matrix21Row[j] += Mult * val;
      
      val  = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += u2_y*ansatz000*test000;
      Matrix22Row[j] += Mult * val;
      
      val  = u2_z*ansatz000*test000;
      Matrix23Row[j] += Mult * val;
      
      
      val  = u3_x*ansatz000*test000;
      Matrix31Row[j] += Mult * val;
      
      val  = u3_y*ansatz000*test000;
      Matrix32Row[j] += Mult * val;
      
      val  = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (u1*ansatz100+u2*ansatz010+u3*ansatz001)*test000;
      val += u3_z*ansatz000*test000;
      Matrix33Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}

// ======================================================================
// Assembling routine for right-hand sides only
// ======================================================================

// ======================================================================
// right-hand side for NSE ONLY
// ======================================================================
void TimeNSRHS3D(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3;
  double test000;
  double *Orig0;
  int i, N_U;
  double c1, c2, c3;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2
  c3 = coeff[3]; // f3

  for(i=0;i<N_U;i++)
  {
    test000 = Orig0[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " "; 
  } // endfor i
}

// ======================================================================
// right-hand side for NSE ONLY 
// ClassicalLES model
// ======================================================================
void TimeNSRHSClassicalLES3D(double Mult, double *coeff,
               double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  Output::print("Nothing has been tested yer: !! :( ");
  ErrThrow("not tested and/or adjusted yet: ");
  
  double *Rhs1, *Rhs2, *Rhs3;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;

  double delta, val1, mu1;
  double D1u1, D2u1, D3u1, D1u2, D2u2, D3u2, D1u3, D2u3, D3u3;
  double a[3][3];
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  D1u1 = param[3]; // D1u1
  D1u2 = param[4]; // D1u2;
  D1u3 = param[5]; // D1u3;
  D2u1 = param[6]; // D2u1
  D2u2 = param[7]; // D2u2;
  D2u3 = param[8]; // D2u2;
  D3u1 = param[9]; // D3u1
  D3u2 = param[10]; // D3u2;
  D3u3 = param[11]; // D3u2;

  // compute Du Du^T
  a[0][0] = D1u1*D1u1 + D2u1*D2u1 + D3u1*D3u1;
  a[0][1] = D1u1*D1u2 + D2u1*D2u2 + D3u1*D3u2;
  a[0][2] = D1u1*D1u3 + D2u1*D2u3 + D3u1*D3u3;
  a[1][0] = a[0][1];
  a[1][1] = D1u2*D1u2 + D2u2*D2u2 + D3u2*D3u2;
  a[1][2] = D1u2*D1u3 + D2u2*D2u3 + D3u2*D3u3;
  a[2][0] = a[0][2];
  a[2][1] = a[1][2];
  a[2][2] = D1u3*D1u3 + D2u3*D2u3 + D3u3*D3u3;
  

  // filter width
  delta =  CharacteristicFilterWidth(hK);
 
  // (delta^2)/(2 gamma) 
  mu1 = 0.5*delta*delta/gamma;
  val1 =  Mult * mu1;

  for(i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    // LES term
    Rhs1[i] += val1*( test100 * a[0][0] + test010 * a[0][1] 
                      + test001 * a[0][2]);
    Rhs2[i] += val1*( test100 * a[1][0] + test010 * a[1][1] 
                      + test001 * a[1][2]);
    Rhs3[i] += val1*( test100 * a[2][0] + test010 * a[2][1] 
                      + test001 * a[2][2]);
  } // endfor i
}

// ======================================================================
// right-hand side for NSE ONLY
// Galdi-Layton model with convolution and auxiliary problem
// ======================================================================
void TimeNSRHSLESModel3D(double Mult, double *coeff,
                         double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs)
{
  Output::print("Nothing has been tested yer: !! :( ");
  ErrThrow("not tested and/or adjusted yet: ");
  
  double *Rhs1, *Rhs2, *Rhs3;
  double test100, test010, test001;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;

  double delta, val1,  mu1;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;
  double gdT11, gdT12, gdT13, gdT22, gdT23, gdT33;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z

  gdT11 = param[0]; // gDeltaT11 or 11 - component of auxiliary problem
  gdT12 = param[1]; // gDeltaT12 or 12 - and 21 - component
  gdT13 = param[2]; // gDeltaT13 or 13 - and 31 - component
  gdT22 = param[3]; // gDeltaT22 or 22 - component of auxiliary problem
  gdT23 = param[4]; // gDeltaT23 or 23 - component of auxiliary problem
  gdT33 = param[5]; // gDeltaT33 or 33 - component of auxiliary problem

  // filter width
  delta =  CharacteristicFilterWidth(hK);
  mu1 = 0.5*delta*delta/gamma;
  val1 =  Mult* mu1;

  for(i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];

    // LES term
    Rhs1[i] += val1*( test100* gdT11 + test010* gdT12 + test001 * gdT13);
    Rhs2[i] += val1*( test100* gdT12 + test010* gdT22 + test001 * gdT23);
    Rhs3[i] += val1*( test100* gdT13 + test010* gdT23 + test001 * gdT33);
  } // endfor i
}

// ======================================================================
// right-hand side for auxiliary problem 
// Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSGL00AuxProblemRHS3D(double Mult, double *coeff,
               double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  Output::print("Nothing has been tested yer: !! :( ");
  ErrThrow("not tested and/or adjusted yet: ");
  
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, *Rhs6, val;
  double mat11, mat12, mat13, mat22, mat23, mat33;
  double test000;
  double *Orig0;
  int i, N_U;
  double D1u1, D2u1, D3u1, D1u2, D2u2, D3u2, D1u3, D2u3, D3u3;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];
  Rhs6 = LocRhs[5];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  D1u1 = param[0]; // D1u1
  D1u2 = param[1]; // D1u2;
  D1u3 = param[2]; // D1u3;
  D2u1 = param[3]; // D2u1
  D2u2 = param[4]; // D2u2;
  D2u3 = param[5]; // D2u2;
  D3u1 = param[6]; // D3u1
  D3u2 = param[7]; // D3u2;
  D3u3 = param[8]; // D3u2;
  
  // compute Du Du^T
  mat11 = D1u1*D1u1 + D2u1*D2u1 + D3u1*D3u1;
  mat12 = D1u1*D1u2 + D2u1*D2u2 + D3u1*D3u2;
  mat13 = D1u1*D1u3 + D2u1*D2u3 + D3u1*D3u3;
  mat22 = D1u2*D1u2 + D2u2*D2u2 + D3u2*D3u2;
  mat23 = D1u2*D1u3 + D2u2*D2u3 + D3u2*D3u3;
  mat33 = D1u3*D1u3 + D2u3*D2u3 + D3u3*D3u3;

  for(i=0;i<N_U;i++)
  {
    test000 = Orig0[i];

    val = Mult*test000;

    Rhs1[i] += val* mat11;
    Rhs2[i] += val* mat12;
    Rhs3[i] += val* mat13;
    Rhs4[i] += val* mat22;
    Rhs5[i] += val* mat23;
    Rhs6[i] += val* mat33;
  } // endfor i
}

// ======================================================================
// right-hand side ONLY for auxiliary problem applied to velocity
// ======================================================================
void TimeNSRHSAuxProblemU(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs)
{
  Output::print("Nothing has been tested yer: !! :( ");
  ErrThrow("not tested and/or adjusted yet: ");
  
  double *Rhs1, *Rhs2, *Rhs3;
  double test000;
  double *Orig0;
  int i, N_U;
  double c1, c2, c3;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  c1 = param[0]; // u1
  c2 = param[1]; // u2
  c3 = param[2]; // u3

  for(i=0;i<N_U;i++)
  {
    test000 = Orig0[i];

    Rhs1[i] += Mult*test000*c1;
    Rhs2[i] += Mult*test000*c2;
    Rhs3[i] += Mult*test000*c3;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " "; 
  } // endfor i
}

// ======================================================================
// right-hand side for additional terms in rhs of small scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_SmallRhs3D(double Mult, double *coeff, 
                           double *param, double hK, 
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  Output::print("Nothing has been tested yer: !! :( ");
  ErrThrow("not tested and/or adjusted yet: ");
  
  double *Rhs1, *Rhs2, *Rhs3;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  double u1, u2, u3,  ho_u1, ho_u2, ho_u3;
  double D1u1, D2u1, D3u1, D1u2, D2u2, D3u2, D1u3, D2u3, D3u3;
  double ho_D1u1, ho_D2u1, ho_D3u1, ho_D1u2, ho_D2u2, ho_D3u2;
  double ho_D1u3, ho_D2u3, ho_D3u3, p, c0;
  int i, N_U;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u_z
  Orig3 = OrigValues[3]; // u

  // large scales
  u1 = param[0];
  u2 = param[1];
  u3 = param[2];
  D1u1 = param[3]; // D1u1
  D1u2 = param[4]; // D1u2;
  D1u3 = param[5]; // D1u3;
  D2u1 = param[6]; // D2u1
  D2u2 = param[7]; // D2u2;
  D2u3 = param[8]; // D2u2;
  D3u1 = param[9]; // D3u1
  D3u2 = param[10]; // D3u2;
  D3u3 = param[11]; // D3u2;
  // small scales
  ho_u1 = param[12];
  ho_u2 = param[13];
  ho_u3 = param[14];
  ho_D1u1 = param[15]; // D1u1
  ho_D1u2 = param[16]; // D1u2;
  ho_D1u3 = param[17]; // D1u3;
  ho_D2u1 = param[18]; // D2u1
  ho_D2u2 = param[19]; // D2u2;
  ho_D2u3 = param[20]; // D2u2;
  ho_D3u1 = param[21]; // D3u1
  ho_D3u2 = param[22]; // D3u2;
  ho_D3u3 = param[23]; // D3u2;
  // pressure
  p = param[24];

  c0 = coeff[0];
  
  for(i=0;i<N_U;i++)
  {
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs1[i] += 2*c0*(D1u1*test100+(D1u2+D2u1)*test010/4+(D1u3+D3u1)*test001/4);
    Rhs1[i] += (u1*D1u1+u2*D2u1+u3*D3u1) * test000;
    Rhs1[i] += (u1*ho_D1u1+u2*ho_D2u1+u3*ho_D3u1) * test000;
    Rhs1[i] += (ho_u1*D1u1+ho_u2*D2u1+ho_u3*D3u1) * test000;
    Rhs1[i] -= p * test100;
    Rhs1[i] *= Mult;
    Rhs2[i] += 2*c0*(D2u2*test010+(D1u2+D2u1)*test100/4+(D2u3+D3u2)*test001/4);
    Rhs2[i] += (u1*D1u2+u2*D2u2+u3*D3u2) * test000;
    Rhs2[i] += (u1*ho_D1u2+u2*ho_D2u2+u3*ho_D3u2) * test000;
    Rhs2[i] += (ho_u1*D1u2+ho_u2*D2u2+ho_u3*D3u2) * test000;
    Rhs2[i] -= p * test010;
    Rhs2[i] *= Mult;
    Rhs3[i] += 2*c0*(D3u3*test001+(D1u3+D3u1)*test100/4+(D2u3+D3u2)*test010/4);
    Rhs3[i] += (u1*D1u3+u2*D2u3+u3*D3u3) * test000;
    Rhs3[i] += (u1*ho_D1u3+u2*ho_D2u3+u3*ho_D3u3) * test000;
    Rhs3[i] += (ho_u1*D1u3+ho_u2*D2u3+ho_u3*D3u3) * test000;
    Rhs3[i] -= p * test001;
    Rhs3[i] *= Mult;
  } // endfor i
}

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// used for : GALERKIN
//            COLETTI (non linear steps)
// ========================================================================
void TimeNSParamsVelo3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old
}

// ========================================================================
// parameters: u1old, u2old,
// used for : COLETTI, Smagorinsky
// ========================================================================
void TimeNSParamsVelo_GradVelo3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  out[12] = in[0]; // x - coordinate for van Driest damping
  out[13] = in[1]; // y - coordinate for van Driest damping
  out[14] = in[2]; // z - coordinate for van Driest damping
}

// ========================================================================
// parameters: gradient(u1), gradient(u2)
// ========================================================================
void TimeNSParamsGradVelo3D(double *in, double *out)
{
  //cout << "GRAD" << endl;
  out[0] = in[6]; // D100(u1old)
  out[1] = in[7]; // D100(u2old)
  out[2] = in[8]; // D100(u3old)
  out[3] = in[9]; // D010(u1old)
  out[4] = in[10]; // D010(u2old)
  out[5] = in[11]; // D010(u3old)
  out[6] = in[12]; // D001(u1old)
  out[7] = in[13]; // D001(u2old)
  out[8] = in[14]; // D001(u3old)
  // cout << in[4] << " E " << in[5] << " " << in[6] << " " << in[7] << endl;
}

// ========================================================================
// parameters: u1old, u2old, u3old
// all partial derivatives
// convolution of u1old, u2old, u3old
// ========================================================================
void TimeNSParamsVelo_GradVelo_ConvVelo3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  out[12] = in[15]; // convolution of u1old
  out[13] = in[16]; // convolution of u2old
  out[14] = in[17]; // convolution of u3old

  out[15] = in[2]; // z - coordinate for van Driest damping
}

// ========================================================================
// used for assembling of term coming from the LES model
// GL00AuxProb and GL00Conv
// parameters used in TimeNSRHSLESModel3D
// ========================================================================
void TimeNSParamsRHSLES3D(double *in, double *out)
{
  // components of convolved tensor
  out[0] = in[3]; // g_d\ast(D1u1*D1u1+D2u1*D2u1+D3u1*D3u1)
  out[1] = in[4]; // g_d\ast(D1u1*D1u2+D2u1*D2u2+D3u1*D3u2)
  out[2] = in[5]; // g_d\ast(D1u1*D1u3+D2u1*D2u3+D3u1*D3u3)
  out[3] = in[6]; // g_d\ast(D1u2*D1u2+D2u2*D2u2+D3u2*D3u2)
  out[4] = in[7]; // g_d\ast(D1u2*D1u3+D2u2*D2u3+D3u2*D3u3)
  out[5] = in[8]; // g_d\ast(D3u1*D3u1+D3u2*D3u2+D3u3*D3u3)

}

// ========================================================================
// used for : classical LES, first iteration step, viscosity type = 4
// ========================================================================
void TimeNSParamsVelo_GradVeloNuT4_3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  // convolution of the solution
  out[12] = in[15]; // g_\delta \ast u1
  out[13] = in[16]; // g_\delta \ast u2
  out[14] = in[17]; // g_\delta \ast u3
}

// ========================================================================
// used for : GL00Convolution, rhs assembling, viscosity type = 4
// ========================================================================
void TimeNSParamsRHSGL00ConvolutionNuT4_3D(double *in, double *out)
{
  // grad u
  out[0] = in[6]; // D1u1
  out[1] = in[7]; // D1u2
  out[2] = in[8]; // D1u3
  out[3] = in[9]; // D2u1
  out[4] = in[10]; // D2u2
  out[5] = in[11]; // D2u3
  out[6] = in[12]; // D3u1
  out[7] = in[13]; // D3u2
  out[8] = in[14]; // D3u3

  // components of convolved tensor or solution of auxiliary problem
  out[9] = in[15]; // g_d\ast(D1u1*D1u1+D2u1*D2u1+D3u1*D3u1)
  out[10] = in[16]; // g_d\ast(D1u1*D1u2+D2u1*D2u2+D3u1*D3u2)
  out[11] = in[17]; // g_d\ast(D1u1*D1u3+D2u1*D2u3+D3u1*D3u3)
  out[12] = in[18]; // g_d\ast(D1u2*D1u2+D2u2*D2u2+D3u2*D3u2)
  out[13] = in[19]; // g_d\ast(D1u2*D1u3+D2u2*D2u3+D3u2*D3u3)
  out[14] = in[20]; // g_d\ast(D3u1*D3u1+D3u2*D3u2+D3u3*D3u3)

  // velocity 
  out[15] = in[3]; // u1
  out[16] = in[4]; // u2
  out[17] = in[5]; // u3

  // convolved velocity
  out[18] = in[21]; //  g_d\ast u1
  out[19] = in[22]; //  g_d\ast u2
  out[20] = in[23]; //  g_d\ast u3
  
}

// ========================================================================
// parameters: 
// used for : GL00AuxProblem, viscosity type = 4
// ========================================================================
void TimeNSParamsGL00AuxProblemNuT4_3D(double *in, double *out)
{

  OutPut("TimeNSParamsGL00AuxProblemNuT43D not implemented " << endl);
  exit(4711);
  // \nabla u
  out[0] = in[4]; // D1u1
  out[1] = in[5]; // D1u2
  out[2] = in[6]; // D2u1
  out[3] = in[7]; // D2u2

  //  solution of auxiliary problem
  out[4] = in[ 8]; // sol11
  out[5] = in[ 9]; // sol12 = sol21
  out[6] = in[10]; // sol22

  // solution
  out[7] = in[2]; // u1
  out[8] = in[3]; // u2

  // convolution of the solution
  out[9] = in[11]; // g_\delta \ast u1
  out[10] = in[12];// g_\delta \ast u2
}

// ========================================================================
// used for VMS, assembling of rhs for small scale equation 
// ========================================================================
void TimeNSParams_VMS_SmallRhs3D(double *in, double *out)
{
  // large scales
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  // small scales
  out[12] = in[15]; // u1old
  out[13] = in[16]; // u2old
  out[14] = in[17]; // u3old

  out[15] = in[18]; // D1u1
  out[16] = in[19]; // D1u2
  out[17] = in[20]; // D1u3
  out[18] = in[21]; // D2u1
  out[19] = in[22]; // D2u2
  out[20] = in[23]; // D2u3
  out[21] = in[24]; // D3u1
  out[22] = in[25]; // D3u2
  out[23] = in[26]; // D3u3
  // large pressure
  out[24] = in[27]; // p
}
// ========================================================================
// parameters: u1old, u2old, G^H
// used for : projection-based VMS
// ========================================================================
void TimeNSParamsVelo_GradVelo_LargeScale3D(double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  out[12] = in[0]; // x - coordinate for van Driest damping
  out[13] = in[1]; // y - coordinate for van Driest damping
  out[14] = in[2]; // z - coordinate for van Driest damping
  
  out[15] = in[15]; // projection space label
}

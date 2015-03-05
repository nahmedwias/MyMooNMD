// ======================================================================
// %W% %G%
//
// common declaration for all 3D convection diffusion problems
// ======================================================================

#include <Database.h>
#include <LinAlg.h>
#include <ConvDiff3D_Routines.h>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>

extern double bound;

/******************************************************************************/
//
// computation of the size of a mesh cell in convection direction
// approximation formula by Tezduyar and Park, CMAME 59, 307 - 325, 1986
//
/******************************************************************************/

double Mesh_size_in_convection_direction(double hK, double b1, double b2, 
					 double b3)
{
    int i;
    double x[8], y[8], z[8], sx, sy, sz, a[64], b[64], den, val, norm_b;

    // already computed for this mesh cell
    if (TDatabase::ParamDB->INTERNAL_HK_CONVECTION >=0)
	return TDatabase::ParamDB->INTERNAL_HK_CONVECTION;

    // tetrahedra
    if (TDatabase::ParamDB->INTERNAL_VERTEX_X[4] == -4711)
    {
	for (i=0;i<4;i++)
	{
	    x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
	    y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
	    z[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Z[i];
	}
        // initialize rhs
	memset(b,0,16*SizeOfDouble);

        // set matrices for computation of the coefficients
        // of the bilinear function
        for (i=0;i<4;i++)
        {
          a[4*i] = 1;
          a[4*i+1] = x[i];
          a[4*i+2] = y[i];
          a[4*i+3] = z[i];
	  b[5*i] = 1;
        }
        // solve system for the coefficients of the bilinear function
	// which is faster ?
        // SolveMultipleSystemsLapack(a,b,4,4,4,4);
        SolveMultipleSystems(a,b,4,4,4,4);
      	
	// compute numerator
	norm_b = sqrt(b1*b1 + b2*b2 + b3*b3);
	// compute denominator 
	den = 0;
	for (i=0;i<4;i++)
	{
	    // value of gradient basis fct. in bary centre
	    // is a constant
	    den += fabs(b1*b[4*i+1]+b2*b[4*i+2]+b3*b[4*i+3]);
	} 
	// return the mesh size in convection direction
	if (den<1e-10)
	{
	    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = hK;
	    return(hK);
	}
	else
	{
	    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = 2*norm_b/den;
	    return(2*norm_b/den);
	}
    }
    else
    { // quadrilateral 
	sx = sy = sz = 0;
	for (i=0;i<8;i++)
	{
	    x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
	    y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
	    z[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Z[i];
	    //OutPut(x[i] <<  " " << y[i] << " ");
	    sx += x[i];
	    sy += y[i];
	    sz += z[i];
	}
	//bary centre
	sx /= 8;
	sy /= 8;
	sz /= 8;
        // initialize rhs
	memset(b,0,64*SizeOfDouble);
	for (i=0;i<8;i++)
	{
	    b[9*i] = 1;
	}
        // set matrices for computation of the coefficients
        // of the bilinear function
	for (i=0;i<8;i++)
	{
	    a[8*i] = 1;
	    a[8*i+1] = x[i];
	    a[8*i+2] = y[i];
	    a[8*i+3] = z[i];
	    a[8*i+4] = x[i]*y[i];
	    a[8*i+5] = x[i]*z[i];
	    a[8*i+6] = y[i]*z[i];
	    a[8*i+7] = x[i]*y[i]*z[i];
	}

        // solve system for the coefficients of the bilinear function
	// which is faster ??
	//SolveMultipleSystemsLapack(a,b,8,8,8,8);
	SolveMultipleSystems(a,b,8,8,8,8);

	// compute numerator
	norm_b = sqrt(b1*b1 + b2*b2 + b3*b3);
	// compute denominator 
	den = 0;
	for (i=0;i<8;i++)
	{
	    // value of gradient basis fct. in bary centre
	    val = b1*(b[8*i+1]+ b[8*i+4]*sy + b[8*i+5]*sz + b[8*i+7]*sy*sz);
	    val += b2*(b[8*i+2]+ b[8*i+4]*sx + b[8*i+6]*sz + b[8*i+7]*sx*sz);
	    val += b3*(b[8*i+3]+ b[8*i+5]*sx + b[8*i+6]*sy + b[8*i+7]*sx*sy);
	    den += fabs(val);
	} 
	// return the mesh size in convection direction
	//OutPut(b1 << " " << b2 << " " << fabs(den) << " " << 2*norm_b/fabs(den) << " " );
	if (den<1e-10)
	{
	    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = hK;
	    return(hK);
	}
	else
	{
	    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = 2*norm_b/den; 	    
	    return(2*norm_b/den);
	}
    }
}

double Mesh_size_in_convection_direction_without_storing(double hK, double b1, double b2, 
					 double b3)
{
    int i;
    double x[8], y[8], z[8], sx, sy, sz, a[64], b[64], den, val, norm_b;

    // tetrahedra
    if (TDatabase::ParamDB->INTERNAL_VERTEX_X[4] == -4711)
    {
	for (i=0;i<4;i++)
	{
	    x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
	    y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
	    z[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Z[i];
	}
        // initialize rhs
	memset(b,0,16*SizeOfDouble);

        // set matrices for computation of the coefficients
        // of the bilinear function
        for (i=0;i<4;i++)
        {
          a[4*i] = 1;
          a[4*i+1] = x[i];
          a[4*i+2] = y[i];
          a[4*i+3] = z[i];
	  b[5*i] = 1;
        }
        // solve system for the coefficients of the bilinear function
	// which is faster ?
        // SolveMultipleSystemsLapack(a,b,4,4,4,4);
        SolveMultipleSystems(a,b,4,4,4,4);
      	
	// compute numerator
	norm_b = sqrt(b1*b1 + b2*b2 + b3*b3);
	// compute denominator 
	den = 0;
	for (i=0;i<4;i++)
	{
	    // value of gradient basis fct. in bary centre
	    // is a constant
	    den += fabs(b1*b[4*i+1]+b2*b[4*i+2]+b3*b[4*i+3]);
	} 
	// return the mesh size in convection direction
	if (den<1e-10)
	{
	    return(hK);
	}
	else
	{
	    return(2*norm_b/den);
	}
    }
    else
    { // quadrilateral 
	sx = sy = sz = 0;
	for (i=0;i<8;i++)
	{
	    x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
	    y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
	    z[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Z[i];
	    //OutPut(x[i] <<  " " << y[i] << " ");
	    sx += x[i];
	    sy += y[i];
	    sz += z[i];
	}
	//bary centre
	sx /= 8;
	sy /= 8;
	sz /= 8;
        // initialize rhs
	memset(b,0,64*SizeOfDouble);
	for (i=0;i<8;i++)
	{
	    b[9*i] = 1;
	}
        // set matrices for computation of the coefficients
        // of the bilinear function
	for (i=0;i<8;i++)
	{
	    a[8*i] = 1;
	    a[8*i+1] = x[i];
	    a[8*i+2] = y[i];
	    a[8*i+3] = z[i];
	    a[8*i+4] = x[i]*y[i];
	    a[8*i+5] = x[i]*z[i];
	    a[8*i+6] = y[i]*z[i];
	    a[8*i+7] = x[i]*y[i]*z[i];
	}

        // solve system for the coefficients of the bilinear function
	// which is faster ??
	//SolveMultipleSystemsLapack(a,b,8,8,8,8);
	SolveMultipleSystems(a,b,8,8,8,8);

	// compute numerator
	norm_b = sqrt(b1*b1 + b2*b2 + b3*b3);
	// compute denominator 
	den = 0;
	for (i=0;i<8;i++)
	{
	    // value of gradient basis fct. in bary centre
	    val = b1*(b[8*i+1]+ b[8*i+4]*sy + b[8*i+5]*sz + b[8*i+7]*sy*sz);
	    val += b2*(b[8*i+2]+ b[8*i+4]*sx + b[8*i+6]*sz + b[8*i+7]*sx*sz);
	    val += b3*(b[8*i+3]+ b[8*i+5]*sx + b[8*i+6]*sy + b[8*i+7]*sx*sy);
	    den += fabs(val);
	} 
	// return the mesh size in convection direction
	//OutPut(b1 << " " << b2 << " " << fabs(den) << " " << 2*norm_b/fabs(den) << " " );
	if (den<1e-10)
	{
	    return(hK);
	}
	else
	{
	    return(2*norm_b/den);
	}
    }
}

double Compute_SDFEM_delta(double hK, double eps, double b1, double b2, double b3, double react, double linfb)
{ 
  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  double alpha,alpha2,tmp1,tmp2,tmp3, delta, norm_b, nu, hk_project;

  if (fabs(eps)<1e-20)
	eps = 1e-20;
 
  switch (TDatabase::ParamDB->SDFEM_TYPE)
  {
    case 0:                     // normal meshes
      if (!TDatabase::ParamDB->SHISHKIN_MESH)
      {
        if(eps < hK*linfb)
          delta = delta0 * hK/linfb;
        else
          delta = delta1 *hK*hK/eps ;
      }
      else                       // delta for SDFEM only in coarse part of Shishkin mesh
      {
        if(hK > bound)
          delta = delta0 * hK/linfb;
        else
          delta = 0;
      }
      break;
    case 1:                    // delta based on 1d Green's formula
      norm_b = sqrt(b1*b1+b2*b2+b3*b3);
      //norm_b = linfb;
      if (norm_b > 0)
      {
        alpha = norm_b*hK/(2*eps);
        delta = hK*(1/tanh(alpha) - 1/alpha)/(2*norm_b);
      }
      else
        delta = 0;
      break;
    case 2:                       // delta based on 1d Green's formula
                                  // for SOLD-papers paper with Petr
      norm_b = sqrt(b1*b1+b2*b2+b3*b3);
      nu = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      hk_project = Mesh_size_in_convection_direction(hK, b1, b2, b3);
      if (norm_b > 0)
      {
        alpha = nu*norm_b*hk_project/(2*eps);
        delta = nu*hk_project*(1/tanh(alpha) - 1/alpha)/(2*norm_b);
      }
      else
        delta = 0;
      break;
      case 3: // for bulk precepetation
      norm_b = sqrt(b1*b1+b2*b2+b3*b3);

      if  (TDatabase::ParamDB->CELL_MEASURE==4)
	  hK = Mesh_size_in_convection_direction(hK, b1, b2, b3);
	  
      alpha = hK*norm_b / (2*eps);
      if(alpha<1e-2)
	  delta  = alpha/(6*norm_b);
      else
	  delta = (1/tanh(alpha)-1/alpha)/(2*norm_b);
     delta *= delta0*hK;
      break;
    case 4:
	// Lube/Rapin (2006)
	// equations with reaction
        // under the assumption that the inverse inequality holds 
	// with nu_inv = 0 (e.g. linear and bilinear elements)
	// polynomial degree = 1
	if  (TDatabase::ParamDB->CELL_MEASURE==4)
	    hK = Mesh_size_in_convection_direction(hK, b1, b2, b3);
	    
        delta = -1;
        
	if (linfb > 0)
	    delta = hK/linfb;

	if (react > 0)
	{
            alpha = 1.0/react;
	    if ((delta < 0) || (alpha<delta))
		delta = alpha;
	}

	if (eps>0)
	{
	    alpha = hK*hK/eps;
	    if ((delta < 0) || (alpha<delta))
		delta = alpha;	    
	}
	
	delta *= delta0;
	if (delta<0)
	    delta = 0;
	break;
     case 5:
	// Lube/Rapin (2006), slightly modified
	// equations with reaction
        // under the assumption that the inverse inequality holds 
	// with nu_inv = 0 (e.g. linear and bilinear elements)
	// polynomial degree = 1
	 if  (TDatabase::ParamDB->CELL_MEASURE==4)
	     hK = Mesh_size_in_convection_direction(hK, b1, b2, b3);
 
        delta = -1;
	norm_b = sqrt(b1*b1+b2*b2+b3*b3);
	if (norm_b > 0)
	    delta = hK/(2*norm_b);

	if (react > 0)
	{
            alpha = 1.0/react;
	    if ((delta < 0) || (alpha<delta))
		delta = alpha;
	}

	if (eps>0)
	{
	    alpha = hK*hK/eps;
	    if ((delta < 0) || (alpha<delta))
		delta = alpha;	    
	}
	
	delta *= delta0;
	if (delta<0)
	    delta = 0;
	break;
      
       case 6:
	// Franca/Valentin (2000)
	// equations with reaction
        // under the assumption that the inverse inequality holds 
	// with nu_inv = 1/3 (e.g. linear and bilinear elements)
	// polynomial degree = 1
        norm_b = sqrt(b1*b1+b2*b2+b3*b3);

      	if  (TDatabase::ParamDB->CELL_MEASURE==4)
	    hK = Mesh_size_in_convection_direction(hK, b1, b2, b3);
	          
	if (react<1e-20)
	    react = 1e20;
	
	alpha  = 6*eps/(hK*hK*react);
	if(alpha <= 1)
	    alpha = 1;
	
        alpha2 = hK*norm_b/(3*eps);
	if(alpha2 <= 1)
	    alpha2 = 1;
   
	delta=1/(react*hK*hK*alpha+6*eps*alpha2);
	delta *= delta0*hK*hK;
        break;
      
      case 7:
	// Codina (2000)
	// equations with reaction
        norm_b = sqrt(b1*b1+b2*b2+b3*b3);
       
	if  (TDatabase::ParamDB->CELL_MEASURE==4)
	    hK = Mesh_size_in_convection_direction(hK, b1, b2, b3);
        
	delta = delta0 * hK * hK/(4*eps+2 * hK * norm_b + hK*hK * react);
       break;
    default :
      OutPut("SDFEM_TYPE "<<TDatabase::ParamDB->SDFEM_TYPE<<
        " not implemented !!!" << endl);
      exit(4711);
  }

  return(delta);
}

// ========================================================================
// compute parameters for SOLD schemes
// stationary cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_z
//        param[4] = ||u^h||_{H^1,K}
//        param[5] = ||R(u^h)||_{L^2,K}
// time-dependent cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_z
//        param[4] = u_old
//        param[5] = u_old_x
//        param[6] = u_old_y
//        param[7] = u_old_z
// ========================================================================

double Compute_SOLD_sigma(double hK, double eps, double b1,
			  double b2, double b3, double c, double f,
			  double linfb, double deltaK, double *param,
			  double residual, int residual_computed, 
			  int time_dependent_problem)
{
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  double u_x,u_y, u_z, norm_u, norm_res, sigma, res, norm_b2, value;
  double b1_orth, b2_orth,  b3_orth, norm_der_u2, linfb_orth, z1, z2, z3, linfz, normz;
  double alpha, beta, gamma, lambda, kappa, omega, rho, norm_b_orth;
  double epsilon= 1e-10, hK_project, y, z;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;

  u_x = param[1];
  u_y = param[2];
  u_z = param[3];
  norm_b2 =  b1*b1+b2*b2+b3*b3;
  norm_der_u2 = u_x*u_x+u_y*u_y+u_z*u_z;
 
  // compute the residual, without Laplacian
  if (!residual_computed)
      // for stationary problems
      res = b1*u_x + b2*u_y + b3*u_z + c*param[0] - f;
  else
      res = residual;

  // compute the parameter for the SOLD scheme
  switch (sold_parameter_type)
  {
    case JSW87:                      // Johnson,Schatz,Wahlbin (1987) (linear)
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2, b3);	
      sigma = sqrt(norm_b2)*hK_project*sqrt(hK_project) - eps;
      if (sigma < 0)
        sigma = 0;
      break;
    case HMM86:                      // Hughes, Mallet, Mizukami (1986)
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      value = b1*u_x + b2*u_y + b3* u_z;   // b \cdot \nabla u^h
                                 // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = value * u_x/ norm_der_u2;
      b2_orth = value * u_y/ norm_der_u2;
      b3_orth = value * u_z/ norm_der_u2;
                                 // ||b_orth||_\infty
      if (fabs(b2_orth)< fabs (b1_orth))
        linfb_orth = fabs(b1_orth);
      else
        linfb_orth = fabs(b2_orth);
      if (linfb_orth < fabs(b3_orth))
	  linfb_orth = fabs(b3_orth);
      // \tau(\b_orth)
      value = Compute_SDFEM_delta(hK, eps, b1_orth, b2_orth, b3_orth, c, linfb_orth);
      // sigma = max ( \tau(b_orth) - \tau(b))
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }
      //OutPut("sigma " << sigma << " orth " << value << " deltaK " << deltaK << endl);
      //OutPut("u_x " << u_x << " u_y " << u_y << endl); 
      value = b1*u_x + b2*u_y + b3*u_z;   // b \cdot nabla u^h
      if (norm_der_u2>0)
                                 // sigma = sigma * residual * (b\cdot u^h)/||\nabla u^h||^2
        sigma *= res*value/norm_der_u2;
      else
        sigma = 0;
      break;
    case TP86_1:
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      alpha = b1*u_x + b2*u_y + b3*u_z;   // b \cdot \nabla u^h
                                 // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = alpha * u_x/ norm_der_u2;
      b2_orth = alpha * u_y/ norm_der_u2;
      b3_orth = alpha * u_z/ norm_der_u2;
      
      rho =  sqrt(b1_orth*b1_orth + b2_orth*b2_orth + b3_orth*b3_orth);
      value = rho/sqrt(norm_b2);
      value = 2*value*(1-value);
      kappa = Mesh_size_in_convection_direction_without_storing(hK,b1_orth,b2_orth,b3_orth);
      lambda = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      sigma = lambda * kappa * value/(2*rho)*res*alpha/norm_der_u2;
      break;
    case TP86_2:
      if (TDatabase::ParamDB->SOLD_U0==0)
      {
	  OutPut("Paramter SOLD_U0 is zero " << endl);
	  exit(4711);
      }
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      alpha = b1*u_x + b2*u_y + b3*u_z;   // b \cdot \nabla u^h
                                 // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = alpha * u_x/ norm_der_u2;
      b2_orth = alpha * u_y/ norm_der_u2;
      b3_orth = alpha * u_z/ norm_der_u2;
      
      rho =  sqrt(b1_orth*b1_orth + b2_orth*b2_orth + b3_orth*b3_orth);
      value = rho/sqrt(norm_b2);
      value = 2*value*(1-value);
      kappa = Mesh_size_in_convection_direction_without_storing(hK,b1_orth,b2_orth,b3_orth);
      lambda = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      sigma = lambda * kappa * kappa* value/(2*rho)*res*alpha/sqrt(norm_der_u2);
      sigma /= TDatabase::ParamDB->SOLD_U0;
      break;
    case GdC88:                      // Galeao, do Carmo (1988)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      if (time_dependent_problem)
      {
	  z1 *= theta1*time_step;
	  z2 *= theta1*time_step;
	  z3 *= theta1*time_step;
      }
      if (fabs(z2)< fabs (z1))
        linfz = fabs(z1);
      else
        linfz = fabs(z2);
      if (linfz < fabs(z3))
	  linfz = fabs(z3);
      
      value = Compute_SDFEM_delta(hK, eps, z1, z2, z3, c, linfz);
      //OutPut(eps << " " << z1 << " " << z2 << " " << c << " " << value << " " << deltaK << endl); 
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }
;
      sigma *= res*res/norm_der_u2;
      break;
    case dCG91:                      // do Carmo, Galeao (1991)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      normz= sqrt(z1*z1+z2*z2+z3*z3);
      if (normz==0)
      {
        sigma = 0;
        break;
      }
      if (time_dependent_problem)
	  normz *= theta1*time_step;
      sigma = deltaK*(sqrt(norm_b2)/normz-1);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= res*res/norm_der_u2;
      break;
    case dCA03:                      // do Carmo, Alvarez (2003)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      if (time_dependent_problem)
      {
	  z1 *= theta1*time_step;
	  z2 *= theta1*time_step;
	  z3 *= theta1*time_step;
      }
      if (fabs(z2)< fabs (z1))
        linfz = fabs(z1);
      else
        linfz = fabs(z2);
      if (linfz < fabs(z3))
	  linfz = fabs(z3);
      value = Compute_SDFEM_delta(hK, eps, z1, z2, z3, c, linfz);
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }
     
      // sigma of [GdC88] is computed, now compute rho
      normz = sqrt(z1*z1+z2*z2+z3*z3);
     
      alpha = normz/sqrt(norm_b2);
      if (alpha >= 1)
      {
	  rho = 1;
      }
     
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2,b3);	
      beta = pow(hK_project,1-alpha*alpha);
      if (beta > 1)
        beta = 1;
     
      gamma = (alpha+beta)/2.0;
      if (gamma > beta)
        gamma = beta;
     
      if (alpha > fabs (res))
        value = alpha;
      else
        value = fabs(res);
      lambda = pow(value,3+alpha/2+alpha*alpha);
      if (0.25+alpha>0.5)
        value = 0.25+alpha;
      else
        value = 0.5;
      lambda /= pow(gamma, value);
      if (lambda >= 1)
      {
	  rho = 1;
      }

      value = (1-lambda)/(1+lambda);
      kappa = pow(fabs(2-lambda),value);
      kappa -= 1;

      value = pow(gamma,2-alpha*alpha);
      omega = alpha*alpha * value/deltaK;

      if ((alpha < 1) && (lambda < 1))
      {
	  rho = pow(omega*sigma,kappa);
      }
      sigma *= rho;
      sigma *= res*res/norm_der_u2;
      break;
    case AS97:                      // Almeida, Silva (1997)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      z3 = res * u_z / norm_der_u2;
      normz= sqrt(z1*z1+z2*z2+z3*z3);
      if (normz==0)
      {
        sigma = 0;
        break;
      }
      if (time_dependent_problem){//OutPut("T");
	  normz *= theta1*time_step;}
      value = b1*u_x + b2*u_y + b3*u_z;

      if (res==0)
        value = 1;
      else
	  value /= res;

      if (value < 1)
        value = 1;
      //else
      //OutPut("val " << value);
 	  
      //OutPut(" res " << res << " tau " << deltaK << " norm b " << sqrt(norm_b2) << " norm z " << normz << " zeta " << value);
      sigma = deltaK*(sqrt(norm_b2)/normz-value);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      //OutPut(" sigma " << sigma << endl);
      sigma *= res*res/norm_der_u2;
      break;
    case C93:                      // Codina (1993)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      value = b1*u_x + b2*u_y + b3*u_z;
      b1_orth = value * u_x/ norm_der_u2;
      b2_orth = value * u_y/ norm_der_u2;
      b3_orth = value * u_z/ norm_der_u2;
      norm_b_orth = sqrt(b1_orth*b1_orth+b2_orth*b2_orth+b3_orth*b3_orth);
      if (norm_b_orth == 0)
      {
        sigma = 0;
        break;
      }
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(hK * norm_b_orth);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * fabs(res)/(2*sqrt(norm_der_u2));
      break;
    case KLR02_1:                      // Knopp, Lube, Rapin (2002)
	//case KLR02_3:                      // same as KLR02_1 with SOLD_S = 0, new version see below
      value = TDatabase::ParamDB->SOLD_S + param[4];
      if (value == 0)
      {
        sigma = 0;
        break;
      }
      alpha = param[5] / value;
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(alpha*hK);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * alpha / 2;
      break;
      case KLR02_3: //pointwise evaluation of residual and norm of u_h
	  if (sqrt(norm_der_u2)<epsilon)
	  {
	      sigma = 0;
	      break;
	  }
	  sigma =  TDatabase::ParamDB->SOLD_CONST * hK *fabs(res)/(2*sqrt(norm_der_u2)) - eps;
	  if (sigma<0)
	      sigma = 0;
	  break;
      case KLR02_4: //pointwise evaluation of residual and norm of u_h with additive constant
	  // THIS IS NOT KLR02_4 FROM John/Knobloch 2007 !!!
	  if (sqrt(norm_der_u2)<epsilon)
	  {
	      sigma = 0;
	      break;
	  }
	  sigma =  TDatabase::ParamDB->SOLD_CONST * hK *fabs(res)/(2*(TDatabase::ParamDB->SOLD_S +
								      sqrt(norm_der_u2))) - eps;
	  if (sigma<0)
	      sigma = 0;
	  break;
     case KLR02_2:                      // similar to KLR02_3 
     case CS99:
	 // norm of convection
      value = sqrt(norm_der_u2);
      if (value == 0)
      {
        sigma = 0;
        break;
      }
      // Q = res/|b|
      alpha = fabs(res) / value;
      // C - (2 eps)/(h Q)
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(alpha*hK);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * alpha / 2;
      break;
      case J90:                      // Johnson (1990)
      alpha = TDatabase::ParamDB->SOLD_CONST;
      kappa = TDatabase::ParamDB->SOLD_POWER;
      sigma = alpha * pow(hK,kappa)*fabs(res)-eps;
      if (sigma < 0)
        sigma = 0;
      break;
      case BE02_1:                     // Burman, Ern 2002
	  alpha = Pi/6;
	  res = res*tanh(res/2.0);
	  value = sqrt(norm_b2)*sqrt(norm_der_u2)+fabs(res);
	  if (value > 0)
	      sigma = deltaK * norm_b2*fabs(res)/value;
	  else
	  {
	      sigma = 0;
	      break;
	  }
	  if (norm_b2>0)
	  {
	      z1 = (1-b1*b1/norm_b2)*u_x - b1*(b2*u_y+ b3*u_z)/norm_b2;
	      z2 = -b2*(b1*u_x + b3*u_z)/norm_b2 + (1-b2*b2/norm_b2)*u_y;
	      z3 = -b3*(b1*u_x + b2*u_y)/norm_b2 + (1-b3*b3/norm_b2)*u_z;
	  }
	  else
	  {
	      sigma = 0;
	      break;
	  }
	  value = fabs(res)+tan(alpha)*sqrt(norm_b2)*sqrt(z1*z1+z2*z2+z3*z3);
	  if (value > 0)
	      sigma *= (sqrt(norm_b2)*sqrt(norm_der_u2) + value)/value;
	  else
	  {
	      sigma = 0;
	  }	  
	  break;
      case BE02_2:                     // modified Burman, Ern 2002
	  value = sqrt(norm_b2)*sqrt(norm_der_u2)+fabs(res);
	  if (value > 0)
	      sigma = deltaK * norm_b2*fabs(res)/value;
	  else
	      sigma = 0;
	  break;
      case BE02_3:                     // Burman, Ern 2002, formula (29)
	  alpha = Pi/6;
	  res = res*tanh(res/2.0);
	  if (norm_b2>0)
	  {
	      z1 = (1-b1*b1/norm_b2)*u_x - b1*(b2*u_y+ b3*u_z)/norm_b2;
	      z2 = -b2*(b1*u_x + b3*u_z)/norm_b2 + (1-b2*b2/norm_b2)*u_y;
	      z3 = -b3*(b1*u_x + b2*u_y)/norm_b2 + (1-b3*b3/norm_b2)*u_z;
	  }
	  else
	  {
	      sigma = 0;
	      break;
	  }
	  value = sqrt(res*res + tan(alpha)*tan(alpha)*norm_b2 * (z1*z1+z2*z2+z3*z3));
	  if (value > 0)
	      sigma = deltaK * norm_b2*fabs(res)/value;
	  else
	      sigma = 0;
	  break;
      case Y_Z_BETA:                     // Tezduyar
	  y = fabs(TDatabase::ParamDB->SOLD_U0);
	  beta = TDatabase::ParamDB->SOLD_POWER;
	  if (y==0)
	  {
	      sigma = 0;
	      break;
	  }
	  if (norm_der_u2==0)
	  {
	      sigma = 0;
	      break;
	  }
	  z = fabs(res)/y;
	  norm_der_u2 /= (y*y);
	  norm_der_u2 = pow(norm_der_u2,beta/2.0-1);
	  //hK_project = Mesh_size_in_convection_direction(hK, b1, b2)/2.0;
	  hK_project = Mesh_size_in_convection_direction_without_storing(hK, u_x, u_y, u_z)/2.0;
	  sigma = TDatabase::ParamDB->SOLD_CONST * z * norm_der_u2 * pow(hK_project,beta);
	  //OutPut(sigma << " ");
	  break;
    case JSW87_1:                      // Johnson,Schatz,Wahlbin (1987) (linear)
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2, b3);	
      sigma = sqrt(norm_b2)*hK_project*sqrt(hK_project) - eps;
      sigma /= theta1*time_step;
      if (sigma < 0)
        sigma = 0;
      break;
      case FEM_TVD:                    // algebraic flux correction
	  sigma = 0;
	  break;
      case GENERAL_SOLD:
	  // TDatabase::ParamDB->SOLD_CONST is the eta from the SOLD2-paper
	  if (norm_der_u2>0)
	      sigma = TDatabase::ParamDB->SOLD_CONST*hK*fabs(res)/(2*sqrt(norm_der_u2));
	  else 
	      sigma = 0.0;
	  break;
      case BH04:                       // Burman, Hansbo 2004, edge stabilization
      case BE05_1:                     // Burman, Hansbo 2004, edge stabilization
      case BE05_2:                     // Burman, Hansbo 2004, edge stabilization
      case LP96:                       // Layton, Polman 1996
      case MH_Kno06:                   // improved Mizukami-Hughes, by Knobloch 2006
      default :
	  OutPut("SOLD type " << sold_parameter_type << " not available" << endl);
	  exit(4711);
  }

  // scaling of sigma accordingly to Knopp, Lube, Rapin (2002)
  if (TDatabase::ParamDB->SOLD_PARAMETER_SCALING)
  {
    value = TDatabase::ParamDB->SOLD_S + param[4];
    if (value == 0)
      sigma = 0;
    else
    {
      alpha = param[5] / value;
      sigma *= alpha*alpha;
    }
  }
  else
    sigma *= TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR;

  return (sigma);
}

void Compute_Crosswind_Plane(double *n, double *v1, double *v2)
{
  double nx, ny, nz, len, nn;

  len = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  if (len==0)
  {
    OutPut("no crosswind plane " << endl);
    exit(4711);
  }

  nx = n[0]/len;
  ny = n[1]/len;
  nz = n[2]/len;

  if ( (fabs(nx)>=0.5) || (fabs(ny)>=0.5))
  {
    nn = sqrt(nx*nx+ny*ny);
    v1[0] = ny/nn;
    v1[1] = -nx/nn;
    v1[2] = 0;
    v2[0] = -v1[1]*nz;
    v2[1] = v1[0]*nz;
    v2[2] = v1[1]*nx-v1[0]*ny;
  }
  else
  {
    nn = sqrt(ny*ny+nz*nz);
    v1[0] = 0;
    v1[1] = -nz/nn;
    v1[2] = ny/nn;
    v2[0] = v1[2]*ny-v1[1]*nz;
    v2[1] = - v1[2]*nx;
    v2[2] = v1[1]*nx;
  }
}


void BilinearAssemble(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0];
  c1 = coeff[1];
  c2 = coeff[2];
  c3 = coeff[3];
  c4 = coeff[4];
  c5 = coeff[5];

  
//   cout << " c0 " << c0 <<" c5 " << c5 << endl;
  
  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs[i] += Mult*test000*c5;
    
//      if(TDatabase::ParamDB->P14<Mult) TDatabase::ParamDB->P14=Mult;    

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001)*test000;
      val += c4*ansatz000*test000;

      val *= Mult;

      MatrixRow[j] += val;
    }                            // endfor j
  }                              // endfor i
}


// SDFEM without 2nd derivatives
void BilinearAssemble_SD(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, c6;
  double delta, bgradv;
 
  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // b_1
  c2 = coeff[2];                 // b_2
  c3 = coeff[3];                 // b_3
  c4 = coeff[4];                 // c
  c5 = coeff[5];                 // f
  c6 = coeff[6];                 // ||b||

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c4, c6);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    bgradv = c1*test100+c2*test010+c3*test001;

    Rhs[i] += Mult*(test000+delta*bgradv)*c5;

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001)*test000;
      val += c4*ansatz000*test000;

      val += delta * (c1*ansatz100+c2*ansatz010+c3*ansatz001
        +c4*ansatz000) * bgradv;

      val *=Mult;

      MatrixRow[j] += val;

    }                            // endfor j
  }                              // endfor i
}


/*
void BilinearAssemble_UPW1(double Mult, double *coeff, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j, N_;
  double c0;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];

  c0 = coeff[0];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = c0*(test10*ansatz10+test01*ansatz01);

      val *= Mult;

      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}
*/

void BilinearAssemble_UPW2(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_;
  double c0, c4, c5;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];

  c0 = coeff[0];
  c4 = coeff[4];
  c5 = coeff[5];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    Rhs[i] += Mult*test000*c5;

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += c4*ansatz000*test000;

      val *= Mult;

      MatrixRow[j] += val;
    }                            // endfor j
  }                              // endfor i
}

// ========================================================================
// SOLD schemes which add isotropic diffusion
// Hughes, Mallet, Mizukami (1986)
// Tezduyar, Park (1986)
// Galeao, do Carmo (1988)
// do Carmo, Galeao (1991)
// do Carmo, Alvarez (2003)
// Almeida, Silva (1997)
// Knopp, Lube, Rapin (2002)
// Johnson (1990)
// Johnson (1992)
// ========================================================================

void BilinearAssemble_SOLD(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, c6;
  double delta, bgradv, sigma, norm_b;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // b_1
  c2 = coeff[2];                 // b_2
  c3 = coeff[3];                 // b_3
  c4 = coeff[4];                 // c
  c5 = coeff[5];                 // f
  c6 = coeff[6];                 // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c4, c6);

  sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, c6, delta, param, 0, 0, 0);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    bgradv = c1*test100+c2*test010+c3*test001;

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = (c0+sigma)*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001)*test000;
      val += c4*ansatz000*test000;

      val += delta * (c1*ansatz100+c2*ansatz010+c3*ansatz001
        +c4*ansatz000) * bgradv;

      val *=Mult;

      MatrixRow[j] += val;

    }                            // endfor j
  }                              // endfor i
}


// ========================================================================
// SOLD schemes which add diffusion only orthogonal to convection
// Johnson,Schatz,Wahlbin (1987) (linear)
// Knopp, Lube, Rapin (2002)
// Codina (1993)
// ========================================================================

void BilinearAssemble_SOLD_Orthogonal(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz000, ansatz100, ansatz010, ansatz001;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5, c6;
  double delta, bgradv, sigma, norm_b;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // b_1
  c2 = coeff[2];                 // b_2
  c3 = coeff[3];                 // b_3
  c4 = coeff[4];                 // c
  c5 = coeff[5];                 // f
  c6 = coeff[6];                 // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c4, c6);

  sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, c6, delta, param, 0, 0, 0);

  norm_b = c1*c1+c2*c2+c3*c3;

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];

    bgradv = c1*test100+c2*test010+c3*test001;

    for(j=0;j<N_;j++)
    {
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      val = c0*(test100*ansatz100+test010*ansatz010+test001*ansatz001);
      val += (c1*ansatz100+c2*ansatz010+c3*ansatz001)*test000;
      val += c4*ansatz000*test000;

      val += delta * (c1*ansatz100+c2*ansatz010+c3*ansatz001
        +c4*ansatz000) * bgradv;

      if (norm_b > 0)
        val += sigma*( (c2*c2+c3*c3)*ansatz100*test100 + (c1*c1+c3*c3)*ansatz010*test010
          + (c1*c1+c2*c2)*ansatz001*test001
          -c1*c3*(ansatz001*test100+ansatz100*test001)
          -c1*c2*(ansatz010*test100+ansatz100*test010)
          -c2*c3*(ansatz001*test010+ansatz010*test001))/norm_b;

      val *=Mult;

      MatrixRow[j] += val;

    }                            // endfor j
  }                              // endfor i
}


// ========================================================================
// parameters:  H1 norm and  norm of residual
// ========================================================================
void DC_CD_Params(double *in, double *out)
{
  out[0] = in[3];                // H1 norm
  out[1] = in[4];                // norm of residual
}


// ========================================================================
// parameters:  partial derivatives
// ========================================================================
void MBE_Params(double *in, double *out)
{
  out[0] = in[3];                // u
  out[1] = in[4];                // u_x
  out[2] = in[5];                // u_y
  out[3] = in[6];                // u_z
}


// ========================================================================
// parameters:  SC_2
// ========================================================================
void SC_2_Params(double *in, double *out)
{
  out[0] = in[3];                // H1 norm
  out[1] = in[4];                // norm of residual
  out[2] = in[5];                // u_x
  out[3] = in[6];                // u_y
  out[4] = in[7];                // u_z
}


// ========================================================================
// parameters:  SOLD
// ========================================================================
void SOLD_Params(double *in, double *out)
{
  out[0] = in[3];                // u
  out[1] = in[4];                // u_x
  out[2] = in[5];                // u_y
  out[3] = in[6];                // u_z
  out[4] = in[7];                // ||u^h||_{H^1,K}
  out[5] = in[8];                // ||R(u^h)||_{L^2,K}
}


// ========================================================================
// write parameters for latex table
// ========================================================================
void ParamsForLatexTable()
{
  /* switch (TDatabase::ParamDB->INTERNAL_MESH_CELL_TYPE)
  {
    case 3: OutPut("$P_");
  break;
    case 4:
  OutPut("$Q_");
  break;
  }
  OutPut(TDatabase::ParamDB->ANSATZ_ORDER << "$&"); */
  if (!TDatabase::ParamDB->SOLD_TYPE)
  {
    OutPut(" &SDFEM " << TDatabase::ParamDB->SDFEM_TYPE <<"& &");
  }
  else
  {
    switch (TDatabase::ParamDB->SOLD_TYPE)
    {
      case 1:
        OutPut("isotr &");
        break;
      case 2:
        OutPut("ortho &");
        break;
      default:
        OutPut("incorrect SOLD_TYPE");
        exit(4711);
    }
    switch (TDatabase::ParamDB->SOLD_PARAMETER_TYPE)
    {
      case 0:
        OutPut("[JSW87] & no para. &");
        break;
      case 1:
        OutPut("[HMM86] & no para. &");
        break;
      case 2:
        break;
      case 3:
        OutPut("[GdC88] & no para. &");
        break;
      case 4:
        OutPut("[dCG91] & no para. &");
        break;
      case 5:
        OutPut("[dCA03] & no para. &");
        break;
      case 6:
        OutPut("[AS97] & no para. &");
        break;
      case 7:
        OutPut("[Cod93] & $C = " << TDatabase::ParamDB->SOLD_CONST << "$&");
        break;
      case 8:
        OutPut("[KLR02] & $S= " << TDatabase::ParamDB->SOLD_S << ", C = " <<
          TDatabase::ParamDB->SOLD_CONST << "$&");
        break;
      case 9:
        OutPut("[Joh92] & $alpha =  " << TDatabase::ParamDB->SOLD_CONST << "$, $kappa = " <<
          TDatabase::ParamDB->SOLD_POWER  << "$&");
        break;
    }
    /*
        if (TDatabase::ParamDB->SOLD_PARAMETER_SCALING)
        {
      OutPut("scale [KLR02] &");
        }
        else
        {
      OutPut("scale " << TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR<< "&");
      }*/
  }
}


void SetSoldParameters(int i)
{
  int max = 16;
  TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR = 1.0;
  TDatabase::ParamDB->SOLD_TYPE = 1;
  if (i>2*max-1)
  {
    i -= 2*max;
    TDatabase::ParamDB->SOLD_TYPE = 2;
  }
  if ((i>max-1)&&(i<2*max))
  {
    i -= max;
    TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR = 0.5;
  }
  /*else
  if ((i>2*max-1))
  {
    i -= 2*max;
    TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR = 1.0;
    TDatabase::ParamDB->SOLD_S = 0;
    TDatabase::ParamDB->SOLD_PARAMETER_SCALING = 1;
    }*/
  switch(i)
  {
    case 0:
      TDatabase::ParamDB->SOLD_TYPE = 0;
      break;
    case 1:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
      break;
    case 2:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 1;
      break;
    case 3:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 3;
      break;
    case 4:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 4;
      break;
    case 5:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 5;
      break;
    case 6:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 6;
      break;
    case 7:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 7;
      TDatabase::ParamDB->SOLD_CONST = 0.25;
      break;
    case 8:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 7;
      TDatabase::ParamDB->SOLD_CONST = 0.5;
      break;
    case 9:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 7;
      TDatabase::ParamDB->SOLD_CONST = 0.75;
      break;
    case 10:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.25;
      TDatabase::ParamDB->SOLD_S = 0.0;
      break;
    case 11:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.5;
      TDatabase::ParamDB->SOLD_S = 0.0;
      break;
    case 12:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.75;
      TDatabase::ParamDB->SOLD_S = 0.0;
      break;
    case 13:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 8;
      TDatabase::ParamDB->SOLD_CONST = 0.5;
      TDatabase::ParamDB->SOLD_S = 1.0;
      break;
    case 14:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 9;
      TDatabase::ParamDB->SOLD_CONST = 1;
      TDatabase::ParamDB->SOLD_POWER = 1.0;
      break;
    case 15:
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 9;
      TDatabase::ParamDB->SOLD_CONST = 1;
      TDatabase::ParamDB->SOLD_POWER = 2.0;
      break;
    default:
      OutPut("Wrong input in SetSoldParameters" << endl);
      exit(4711);
  }
}


void EvaluateResults_old(int *nonlinite, double *scale, double *globmin,
double *globmax, double *bdrmin, double *bdrmax,
double *bdrh11)

{
  int i, rank;
  double supgmax, supgmin, supgbdrmax, supgbdrmin, supgbdrh11;
  double bdrh11min, val;

  supgmax = globmax[0]-1.0;
  supgmin = fabs(globmin[0]);
  supgbdrmax = bdrmax[0]-1.0;
  supgbdrmin = fabs (bdrmin[0]);
  supgbdrh11 = bdrh11[0];
  bdrh11min =  supgbdrh11;

  // set reference values
  for (i=0;i<64;i++)
  {
    // scheme not converged
    if (nonlinite[i] == TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR)
      continue;
    if (bdrh11[i] <  bdrh11min)
      bdrh11min = bdrh11[i];
  }

  for (i=0;i<64;i++)
  {
    // scheme not converged
    if (nonlinite[i] ==TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR)
      continue;
    // recover the parameters of the computation
    SetSoldParameters(i);
    // print the parameters
    ParamsForLatexTable();

    // scale
    OutPut(scale[i] << "&");

    // global minimum
    OutPut(globmin[i] << "&");
    val = fabs(globmin[i]/supgmin);
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank = 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank = 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank = 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank = -1;
          }
          else
          {
            OutPut("$--$ &");
            rank = -5;
          }
        }
      }
    }

    // global maximum
    OutPut(globmax[i] << "&");
    val = fabs((globmax[i]-1)/supgmax);
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank += 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank += 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank += 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank -= 1;
          }
          else
          {
            OutPut("$--$ &");
            rank -= 5;
          }
        }
      }
    }

    // boundary minimum
    OutPut(bdrmin[i] << "&");
    val = fabs(bdrmin[i]/supgbdrmin);
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank += 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank += 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank += 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank += -1;
          }
          else
          {
            OutPut("$--$ &");
            rank += -5;
          }
        }
      }
    }
    // boundary maximum
    OutPut(bdrmax[i] << "&");
    val = fabs((bdrmax[i]-1)/supgbdrmax);
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank += 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank += 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank += 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank += -1;
          }
          else
          {
            OutPut("$--$ &");
            rank += -5;
          }
        }
      }
    }
    /*// boundary H11
    OutPut(bdrh11[i] << "&");
    val = fabs((bdrh11[i]-bdrh11min)/(supgbdrh11-bdrh11min));
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank += 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank += 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank += 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank += -1;
          }
          else
          {
            OutPut("$--$ &");
            rank += -5;
          }
        }
      }
      }*/
    OutPut(rank<< "\\"<<"\\"<<endl);
  }
  return;
}


void EvaluateResults(int *nonlinite, double *scale, double *globmin,
double *globmax, double *bdrmin, double *bdrmax,
double *bdrh11, double *maxprinciple_cells, double *maxprinciple_max,
int ij_end)
{
  int i, rank;
  double supgmax, supgmin, supgbdrmax, supgbdrmin, supgbdrh11;
  double bdrh11min, val, maxprinciple_cell_ref;

  supgmax = globmax[0]-1.0;
  supgmin = fabs(globmin[0]);
  supgbdrmax = bdrmax[0]-1.0;
  supgbdrmin = fabs (bdrmin[0]);
  supgbdrh11 = bdrh11[0];
  bdrh11min =  supgbdrh11;
  maxprinciple_cell_ref = 100 - maxprinciple_cells[0];

  // set reference values
  for (i=0;i<ij_end;i++)
  {
    // scheme not converged
    if (nonlinite[i] == TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR)
      continue;
    if (bdrh11[i] <  bdrh11min)
      bdrh11min = bdrh11[i];
  }

  for (i=0;i<ij_end;i++)
  {
    // scheme not converged
    if (nonlinite[i] ==TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR)
      continue;
    // recover the parameters of the computation
    SetSoldParameters(i);
    // print the parameters
    ParamsForLatexTable();

    // scale
    //OutPut(scale[i] << "&");

    // global minimum
    OutPut(globmin[i] << "&");
    val = fabs(globmin[i]/supgmin);
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank = 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank = 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank = 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank = -1;
          }
          else
          {
            OutPut("$--$ &");
            rank = -5;
          }
        }
      }
    }

    // global maximum
    OutPut(globmax[i] << "&");
    val = fabs((globmax[i]-1)/supgmax);
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank += 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank += 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank += 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank -= 1;
          }
          else
          {
            OutPut("$--$ &");
            rank -= 5;
          }
        }
      }
    }
    // max principle cells
    OutPut(maxprinciple_cells[i] << "&");
    val = fabs(100-maxprinciple_cells[i])/maxprinciple_cell_ref;
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank += 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank += 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank += 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank -= 1;
          }
          else
          {
            OutPut("$--$ &");
            rank -= 5;
          }
        }
      }
    }
    // max principle maximum
    OutPut(maxprinciple_max[i] << "&");
    val = maxprinciple_max[i]/maxprinciple_max[0];
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank += 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank += 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank += 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank -= 1;
          }
          else
          {
            OutPut("$--$ &");
            rank -= 5;
          }
        }
      }
    }

    // boundary minimum
    OutPut(bdrmin[i] << "&");
    val = fabs(bdrmin[i]/supgbdrmin);
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank += 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank += 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank += 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank += -1;
          }
          else
          {
            OutPut("$--$ &");
            rank += -5;
          }
        }
      }
    }
    // boundary maximum
    OutPut(bdrmax[i] << "&");
    val = fabs((bdrmax[i]-1)/supgbdrmax);
    if (val <= 0.1)
    {
      OutPut("$++$ &");
      rank += 2;
    }
    else
    {
      if (val<=0.25)
      {
        OutPut("$+$&");
        rank += 1;
      }
      else
      {
        if (val<=0.5)
        {
          OutPut("$\\circ$ &");
          rank += 0;
        }
        else
        {
          if (val<=0.75)
          {
            OutPut("$-$ &");
            rank += -1;
          }
          else
          {
            OutPut("$--$ &");
            rank += -5;
          }
        }
      }
    }
    /*    // boundary H11
        OutPut(bdrh11[i] << "&");
        val = fabs((bdrh11[i]-bdrh11min)/(supgbdrh11-bdrh11min));
        if (val <= 0.1)
        {
          OutPut("$++$ &");
          rank += 2;
        }
        else
        {
          if (val<=0.25)
          {
            OutPut("$+$&");
            rank += 1;
          }
          else
          {
            if (val<=0.5)
            {
              OutPut("$\circ$ &");
              rank += 0;
            }
            else
            {
              if (val<=0.75)
              {
                OutPut("$-$ &");
                rank += -1;
              }
              else
              {
                OutPut("$--$ &");
                rank += -5;
              }
            }
          }
          }*/
    OutPut(rank<< "\\"<<"\\"<<endl);
  }
  return;
}

// ======================================================================
// @(#)ConvDiff2D.C        1.4 06/27/00
//
// common declaration for all convection diffusion problems
// ======================================================================

#include <Database.h>
#include <FEDatabase2D.h>
#include <FEFunction2D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <MooNMD_Io.h>

#include <stdlib.h>

extern double bound;

/******************************************************************************/
//
// computation of the size of a mesh cell in convection direction
// approximation formula by Tezduyar and Park, CMAME 59, 307 - 325, 1986
//
/******************************************************************************/

double Mesh_size_in_convection_direction(double hK, double b1, double b2)
{
  int i;
  double x[4], y[4], sx, sy, a[16], b[16], den, val, norm_b;

  // already computed for this mesh cell
  if (TDatabase::ParamDB->INTERNAL_HK_CONVECTION >=0)
    return TDatabase::ParamDB->INTERNAL_HK_CONVECTION;

  // triangles
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[3] == -4711)
  {
    for (i=0;i<3;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
    }
    // initialize rhs
    memset(b,0,9*SizeOfDouble);

    // set matrices for computation of the coefficients
    // of the linear function
    for (i=0;i<3;i++)
    {
      a[3*i] = 1;
      a[3*i+1] = x[i];
      a[3*i+2] = y[i];
      b[4*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    SolveMultipleSystems(a,b,3,3,3,3);

    // compute numerator
    norm_b = sqrt(b1*b1 + b2*b2);
    // compute denominator
    den = 0;
    for (i=0;i<3;i++)
    {
      // value of gradient basis fct. in bary centre
      // is a constant
      den += fabs(b1*b[3*i+1]+b2*b[3*i+2]);
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
  {                                               // quadrilateral
    sx = sy = 0;
    for (i=0;i<4;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
      //OutPut(x[i] <<  " " << y[i] << " ");
      sx += x[i];
      sy += y[i];
    }
    //bary centre
    sx /= 4;
    sy /= 4;
    // initialize rhs
    memset(b,0,16*SizeOfDouble);

    // set matrices for computation of the coefficients
    // of the bilinear function
    for (i=0;i<4;i++)
    {
      a[4*i] = 1;
      a[4*i+1] = x[i];
      a[4*i+2] = y[i];
      a[4*i+3] = x[i]*y[i];
      b[5*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    SolveMultipleSystems(a,b,4,4,4,4);

    // compute numerator
    norm_b = sqrt(b1*b1 + b2*b2);
    // compute denominator
    den = 0;
    for (i=0;i<4;i++)
    {
      // value of gradient basis fct. in bary centre
      val = b1*(b[4*i+1] + b[4*i+3] * sy);
      val += b2*(b[4*i+2] + b[4*i+3] * sx);
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


double Mesh_size_in_convection_direction_without_storing(double hK, double b1, double b2)
{
  int i;
  double x[4], y[4], sx, sy, a[16], b[16], den, val, norm_b;

  // triangles
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[3] == -4711)
  {
    for (i=0;i<3;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
    }
    // initialize rhs
    memset(b,0,9*SizeOfDouble);

    // set matrices for computation of the coefficients
    // of the linear function
    for (i=0;i<3;i++)
    {
      a[3*i] = 1;
      a[3*i+1] = x[i];
      a[3*i+2] = y[i];
      b[4*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    SolveMultipleSystems(a,b,3,3,3,3);

    // compute numerator
    norm_b = sqrt(b1*b1 + b2*b2);
    // compute denominator
    den = 0;
    for (i=0;i<3;i++)
    {
      // value of gradient basis fct. in bary centre
      // is a constant
      den += fabs(b1*b[3*i+1]+b2*b[3*i+2]);
    }
    // return the mesh size in convection direction
    if (den<1e-10)
      return(hK);
    else
      return(2*norm_b/den);
  }
  else
  {                                               // quadrilateral
    sx = sy = 0;
    for (i=0;i<4;i++)
    {
      x[i] = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      y[i] = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
      //OutPut(x[i] <<  " " << y[i] << " ");
      sx += x[i];
      sy += y[i];
    }
    //bary centre
    sx /= 4;
    sy /= 4;
    // initialize rhs
    memset(b,0,16*SizeOfDouble);

    // set matrices for computation of the coefficients
    // of the bilinear function
    for (i=0;i<4;i++)
    {
      a[4*i] = 1;
      a[4*i+1] = x[i];
      a[4*i+2] = y[i];
      a[4*i+3] = x[i]*y[i];
      b[5*i] = 1;
    }
    // solve system for the coefficients of the bilinear function
    SolveMultipleSystems(a,b,4,4,4,4);

    // compute numerator
    norm_b = sqrt(b1*b1 + b2*b2);
    // compute denominator
    den = 0;
    for (i=0;i<4;i++)
    {
      // value of gradient basis fct. in bary centre
      val = b1*(b[4*i+1] + b[4*i+3] * sy);
      val += b2*(b[4*i+2] + b[4*i+3] * sx);
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

double Compute_SDFEM_delta(double hK, double eps, double b1, double b2,
double react, double linfb)
{
  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;
  double alpha,alpha2,tmp1,tmp2,tmp3, delta, norm_b, nu, hk_project;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  int i;

  if (fabs(eps)<1e-20)
    eps = 1e-20;

  switch (TDatabase::ParamDB->SDFEM_TYPE)
  {
    case 0:                                       // normal meshes
      if (!TDatabase::ParamDB->SHISHKIN_MESH)
      {
        if(eps < hK*linfb)
	 { delta = delta0 * hK/linfb; }
        else
	 { delta = delta1 *hK*hK/eps ; }
	
// 	cout << "delta " << delta <<endl;
// 	exit(0);
      }
      else                                        // delta for SDFEM only in coarse part of Shishkin mesh
      {
	cout<< " delta for SDFEM only in coarse part of Shishkin mesh " <<endl;
	exit(0);
        if(hK > bound)
          delta = delta0 * hK/linfb;
        else
          delta = 0;
	  }
      break;
    case 1:                                       // delta based on 1d Green's formula
      norm_b = sqrt(b1*b1+b2*b2);
      //norm_b = linfb;
      if (norm_b > 0)
      {
        alpha = norm_b*hK/(2*eps);
        delta = hK*(1/tanh(alpha) - 1/alpha)/(2*norm_b);
      }
      else
        delta = 0;
      break;
    case 2:                                       // delta based on 1d Green's formula
      // for SOLD-papers paper with Petr
      norm_b = sqrt(b1*b1+b2*b2);
      nu = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      hk_project = Mesh_size_in_convection_direction(hK, b1, b2);
      if (norm_b > 0)
      {
        alpha = nu*norm_b*hk_project/(2*eps);
        delta = nu*hk_project*(1/tanh(alpha) - 1/alpha)/(2*norm_b);
      }
      else
      {
        delta = 0;
      }
      break;
    case 3:                                       // for bulk precepetation
      norm_b = sqrt(b1*b1+b2*b2);

      if  (TDatabase::ParamDB->CELL_MEASURE==4)
        hK = Mesh_size_in_convection_direction(hK, b1, b2);

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
        hK = Mesh_size_in_convection_direction(hK, b1, b2);

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
    case 8:
      // Lube/Rapin (2006), slightly modified
      // equations with reaction
      // under the assumption that the inverse inequality holds
      // with nu_inv = 0 (e.g. linear and bilinear elements)
      // polynomial degree = 1
      // !!! case 8: reaction without term 1
      if  (TDatabase::ParamDB->CELL_MEASURE==4)
        hK = Mesh_size_in_convection_direction(hK, b1, b2);

      delta = -1;
      norm_b = sqrt(b1*b1+b2*b2);
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
      norm_b = sqrt(b1*b1+b2*b2);

      if  (TDatabase::ParamDB->CELL_MEASURE==4)
        hK = Mesh_size_in_convection_direction(hK, b1, b2);

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
      norm_b = sqrt(b1*b1+b2*b2);

      if  (TDatabase::ParamDB->CELL_MEASURE==4)
        hK = Mesh_size_in_convection_direction(hK, b1, b2);

      delta = delta0 * hK * hK/(4*eps+2 * hK * norm_b + hK*hK * react);
      break;
    case 9:
	// first estimate in paper with Julia Novo
	// get the unscaled parameters
	b1 /= (theta1*time_step);
	b2 /= (theta1*time_step);
	eps /= (theta1*time_step);
	if  (TDatabase::ParamDB->CELL_MEASURE==4)
	    hK = Mesh_size_in_convection_direction(hK, b1, b2);
	if (eps <= hK)
	{
	    delta = delta0*time_step;
            //delta = delta0*time_step*time_step;
	}
	else
	{
            delta = delta0*time_step;
	    //delta = delta0*time_step*time_step;
	}
	break;
     case 10:
	// second estimate in paper with Julia Novo
	// get the unscaled parameters
        if(eps <= hK)
          delta = delta0 * hK;
        else
          delta = delta0 *hK*hK/eps ;
	delta = delta0 * hK;
        break;
      case 11:
	// for estimate in paper with Julia Novo
	norm_b = sqrt(b1*b1+b2*b2)/time_step;
	delta = delta0 * hK * sqrt(time_step)/norm_b;
	break;
  case 100:
	// take value from piecewise constant field
	// only on finest level available
	i = TDatabase::ParamDB->INTERNAL_LOCAL_DOF;
	delta =  TDatabase::ParamDB->INTERNAL_P1_Array[i];
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
//        param[3] = u_xx
//        param[4] = u_yy
//        param[5] = ||u^h||_{H^1,K}
//        param[6] = ||R(u^h)||_{L^2,K}
// time-dependent cdr equations
// input: param[0] = u
//        param[1] = u_x
//        param[2] = u_y
//        param[3] = u_xx
//        param[4] = u_yy
//        param[5] = u_old
//        param[6] = u_old_x
//        param[7] = u_old_y
//        param[8] = u_old_xx
//        param[9] = u_old_yy
// ========================================================================

double Compute_SOLD_sigma(double hK, double eps, double b1,
double b2, double c, double f,
double linfb, double deltaK, double *param,
double residual, int residual_computed,
int time_dependent_problem)
{
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  double u_x,u_y, norm_u, norm_res, sigma, res, norm_b2, value;
  double b1_orth, b2_orth, norm_der_u2, linfb_orth, z1, z2, linfz, normz;
  double alpha, beta, gamma, lambda, kappa, omega, rho, norm_b_orth;
  double epsilon= 1e-10, hK_project, y, z;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;

  u_x = param[1];
  u_y = param[2];
  norm_b2 =  b1*b1+b2*b2;
  norm_der_u2 = u_x*u_x+u_y*u_y;

  // compute the residual
  if (!residual_computed)
    // for stationary problems
    res = - eps*(param[3]+param[4]) + b1*u_x + b2*u_y + c*param[0] - f;
  else
    res = residual;

  // compute the parameter for the SOLD scheme
  switch (sold_parameter_type)
  {
    case JSW87:                                   // Johnson,Schatz,Wahlbin (1987) (linear)
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2);
      sigma = sqrt(norm_b2)*hK_project*sqrt(hK_project) - eps;
      if (sigma < 0)
        sigma = 0;
      break;
    case HMM86:                                   // Hughes, Mallet, Mizukami (1986)
      if ((norm_b2< epsilon)||(norm_der_u2 < epsilon))
      {
        sigma = 0;
        break;
      }
      value = b1*u_x + b2*u_y;                    // b \cdot \nabla u^h
      // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = value * u_x/ norm_der_u2;
      b2_orth = value * u_y/ norm_der_u2;
      // ||b_orth||_\infty
      if (fabs(b2_orth)< fabs (b1_orth))
        linfb_orth = fabs(b1_orth);
      else
        linfb_orth = fabs(b2_orth);
      // \tau(\b_orth)
      value = Compute_SDFEM_delta(hK, eps, b1_orth, b2_orth, c, linfb_orth);
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
      value = b1*u_x + b2*u_y;                    // b \cdot nabla u^h
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
      alpha = b1*u_x + b2*u_y;                    // b \cdot \nabla u^h
      // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = alpha * u_x/ norm_der_u2;
      b2_orth = alpha * u_y/ norm_der_u2;

      rho =  sqrt(b1_orth*b1_orth + b2_orth*b2_orth);
      value = rho/sqrt(norm_b2);
      value = 2*value*(1-value);
      kappa = Mesh_size_in_convection_direction_without_storing(hK,b1_orth,b2_orth);
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
      alpha = b1*u_x + b2*u_y;                    // b \cdot \nabla u^h
      // (b \cdot \nabla u^h) u_x/||\nabla u^h||^2
      b1_orth = alpha * u_x/ norm_der_u2;
      b2_orth = alpha * u_y/ norm_der_u2;

      rho =  sqrt(b1_orth*b1_orth + b2_orth*b2_orth);
      value = rho/sqrt(norm_b2);
      value = 2*value*(1-value);
      kappa = Mesh_size_in_convection_direction_without_storing(hK,b1_orth,b2_orth);
      lambda = 1.0/TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE;
      sigma = lambda * kappa * kappa* value/(2*rho)*res*alpha/sqrt(norm_der_u2);
      sigma /= TDatabase::ParamDB->SOLD_U0;
      break;
    case GdC88:                                   // Galeao, do Carmo (1988)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      if (time_dependent_problem)
      {
        z1 *= theta1*time_step;
        z2 *= theta1*time_step;
      }
      if (fabs(z2)< fabs (z1))
        linfz = fabs(z1);
      else
        linfz = fabs(z2);
      value = Compute_SDFEM_delta(hK, eps, z1, z2, c, linfz);
      //OutPut(eps << " " << z1 << " " << z2 << " " << c << " " << value << " " << deltaK << endl);
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }
      sigma *= res*res/norm_der_u2;
      break;
    case dCG91:                                   // do Carmo, Galeao (1991)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      normz= sqrt(z1*z1+z2*z2);
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
      //OutPut(deltaK << " " << sqrt(norm_b2) << " " << normz << " " << sigma << endl);
      break;
    case dCA03:                                   // do Carmo, Alvarez (2003)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      if (time_dependent_problem)
      {
        z1 *= theta1*time_step;
        z2 *= theta1*time_step;
      }
      if (fabs(z2)< fabs (z1))
        linfz = fabs(z1);
      else
        linfz = fabs(z2);
      value = Compute_SDFEM_delta(hK, eps, z1, z2, c, linfz);
      if (value > deltaK)
        sigma = value - deltaK;
      else
      {
        sigma = 0;
        break;
      }

      // sigma of [GdC88] is computed, now compute rho
      normz = sqrt(z1*z1+z2*z2);

      alpha = normz/sqrt(norm_b2);
      if (alpha >= 1)
      {
        rho = 1;
      }

      hK_project = Mesh_size_in_convection_direction(hK, b1, b2);
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
    case AS97:                                    // Almeida, Silva (1997)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      z1 = res * u_x / norm_der_u2;
      z2 = res * u_y / norm_der_u2;
      normz= sqrt(z1*z1+z2*z2);
      if (normz==0)
      {
        sigma = 0;
        break;
      }
      if (time_dependent_problem)                 //OutPut("T");
      {
        normz *= theta1*time_step;
      }
      value = b1*u_x + b2*u_y;

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
    case C93:                                     // Codina (1993)
      if (norm_der_u2 == 0)
      {
        sigma = 0;
        break;
      }
      value = b1*u_x + b2*u_y;
      b1_orth = value * u_x/ norm_der_u2;
      b2_orth = value * u_y/ norm_der_u2;
      norm_b_orth = sqrt(b1_orth*b1_orth+b2_orth*b2_orth);
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
    case KLR02_1:                                 // Knopp, Lube, Rapin (2002)
      //case KLR02_3:                      // same as KLR02_1 with SOLD_S = 0, new version see below
      value = TDatabase::ParamDB->SOLD_S + param[5];
      if (value == 0)
      {
        sigma = 0;
        break;
      }
      alpha = param[6] / value;
      sigma = TDatabase::ParamDB->SOLD_CONST - 2*eps/(alpha*hK);
      if (sigma < 0)
      {
        sigma = 0;
        break;
      }
      sigma *= hK * alpha / 2;
      break;
    case KLR02_3:                                 //pointwise evaluation of residual and norm of u_h
      if (sqrt(norm_der_u2)<epsilon)
      {
        sigma = 0;
        break;
      }
      sigma =  TDatabase::ParamDB->SOLD_CONST * hK *fabs(res)/(2*sqrt(norm_der_u2)) - eps;
      if (sigma<0)
        sigma = 0;
      break;
    case KLR02_4:                                 //pointwise evaluation of residual and norm of u_h with additive constant
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
    case KLR02_2:                                 // similar to KLR02_3
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
    case J90:                                     // Johnson (1990)
      alpha = TDatabase::ParamDB->SOLD_CONST;
      kappa = TDatabase::ParamDB->SOLD_POWER;
      sigma = alpha * pow(hK,kappa)*fabs(res)-eps;
      if (sigma < 0)
        sigma = 0;
      break;
    case BE02_1:                                  // Burman, Ern 2002
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
        z1 = (1-b1*b1/norm_b2)*u_x - b1*b2*u_y/norm_b2;
        z2 = -b1*b2*u_x/norm_b2 + (1-b2*b2/norm_b2)*u_y;
      }
      else
      {
        sigma = 0;
        break;
      }
      value = fabs(res)+tan(alpha)*sqrt(norm_b2)*sqrt(z1*z1+z2*z2);
      if (value > 0)
        sigma *= (sqrt(norm_b2)*sqrt(norm_der_u2) + value)/value;
      else
      {
        sigma = 0;
      }
      break;
    case BE02_2:                                  // modified Burman, Ern 2002
      value = sqrt(norm_b2)*sqrt(norm_der_u2)+fabs(res);
      if (value > 0)
        sigma = deltaK * norm_b2*fabs(res)/value;
      else
        sigma = 0;
      break;
    case BE02_3:                                  // Burman, Ern 2002, formula (29)
      alpha = Pi/6;
      res = res*tanh(res/2.0);
      if (norm_b2>0)
      {
        z1 = (1-b1*b1/norm_b2)*u_x - b1*b2*u_y/norm_b2;
        z2 = -b1*b2*u_x/norm_b2 + (1-b2*b2/norm_b2)*u_y;
      }
      else
      {
        sigma = 0;
        break;
      }
      value = sqrt(res*res + tan(alpha)*tan(alpha)*norm_b2 * (z1*z1+z2*z2));
      if (value > 0)
        sigma = deltaK * norm_b2*fabs(res)/value;
      else
        sigma = 0;
      break;
    case Y_Z_BETA:                                // Tezduyar
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
      hK_project = Mesh_size_in_convection_direction_without_storing(hK, u_x, u_y)/2.0;
      sigma = TDatabase::ParamDB->SOLD_CONST * z * norm_der_u2 * pow(hK_project,beta);
      //OutPut(sigma << " ");
      break;
    case JSW87_1:                                 // Johnson,Schatz,Wahlbin (1987) (linear)
      hK_project = Mesh_size_in_convection_direction(hK, b1, b2);
      sigma = sqrt(norm_b2)*hK_project*sqrt(hK_project) - eps;
      sigma /= theta1*time_step;
      if (sigma < 0)
        sigma = 0;
      break;
    case BH04:                                    // Burman, Hansbo 2004, edge stabilization
    case BE05_1:                                  // Burman, Hansbo 2004, edge stabilization
    case BE05_2:                                  // Burman, Hansbo 2004, edge stabilization
    case LP96:                                    // Layton, Polman 1996
    case MH_Kno06:                                // improved Mizukami-Hughes, by Knobloch 2006
    case FEM_TVD:                                 // algebraic flux correction
      sigma = 0;
      break;
    case GENERAL_SOLD:
      // TDatabase::ParamDB->SOLD_CONST is the eta from the SOLD2-paper
      if (norm_der_u2>0)
        sigma = TDatabase::ParamDB->SOLD_CONST*hK*fabs(res)/(2*sqrt(norm_der_u2));
      else
        sigma = 0.0;
      break;
    default :
      OutPut("SOLD type " << sold_parameter_type << " not available" << endl);
      exit(4711);
  }

  // scaling of sigma accordingly to Knopp, Lube, Rapin (2002)
  if (TDatabase::ParamDB->SOLD_PARAMETER_SCALING)
  {
    value = TDatabase::ParamDB->SOLD_S + param[5];
    if (value == 0)
      sigma = 0;
    else
    {
      alpha = param[6] / value;
      sigma *= alpha*alpha;
    }
  }
  else
    sigma *= TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR;

  return (sigma);
}


// Galerkin discretization
void BilinearAssemble(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  // coefficients
  c0 = coeff[0];                                  // eps
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    // test function
    test10 = Orig0[i];                            // xi derivative
    test01 = Orig1[i];                            // eta derivative
    test00 = Orig2[i];                            // function

    // assemble rhs
    // quad_weigth * test_function * f
    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];                        // xi derivative
      ansatz01 = Orig1[j];                        // eta derivative
      ansatz00 = Orig2[j];                        // function

      // assemble viscous term
      // eps (test_x ansatz_x + test_y ansatz_y)
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      // (b_1 ansatz_x + b_2 ansatz_y) test
      val += (c1*ansatz10+c2*ansatz01)*test00;
      // assembel reactive term
      // c  ansatz test
      val += c3*ansatz00*test00;

      // quad weigth
      val *= Mult;

      // update matrix entry
      MatrixRow[j] += val;
    }                                             // endfor j
  }                                               // endfor i
}



// Galerkin discretization
void BilinearAssemble_Axial3D(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, x, r;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  // coefficients
  c0 = coeff[0];                                  // eps
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  x  = coeff[20];
  r  = fabs(x);
  
  if(r<1e-12)
   {
   OutPut("check BilinearAssemble_Axial3D x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as positive points"<<endl);
   }
   
//    cout<< "r " <<  r<< "  c4  " << c4 << endl;
//     exit(0);
   

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    // test function
    test10 = Orig0[i];                            // xi derivative
    test01 = Orig1[i];                            // eta derivative
    test00 = Orig2[i];                            // function

    // assemble rhs
    // quad_weigth * test_function * f
    Rhs[i] += r*Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];                        // xi derivative
      ansatz01 = Orig1[j];                        // eta derivative
      ansatz00 = Orig2[j];                        // function

      // assemble viscous term
      // eps (test_x ansatz_x + test_y ansatz_y)
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      // (b_1 ansatz_x + b_2 ansatz_y) test
      val += (c1*ansatz10+c2*ansatz01)*test00;
      // assembel reactive term
      // c  ansatz test
      val += c3*ansatz00*test00;

      // quad weigth
      val *= Mult;
      
      //axial 3d
      val *= r;
      // update matrix entry
      MatrixRow[j] += val;
    }                                             // endfor j
  }                                               // endfor i
}




void BilinearAssemble_SD(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;
  double delta, bgradv;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

    Rhs[i] += Mult*(test00+delta*bgradv)*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val += delta * (-c0*(ansatz20+ansatz02)
        +c1*ansatz10+c2*ansatz01
        +c3*ansatz00) * bgradv;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
}


void BilinearAssemble_GLS(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;
  double delta, Lu;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    test20 = Orig3[i];
    test02 = Orig4[i];
      
    Lu = (-c0*(test20+test02)+c1*test10+c2*test01 +c3*test00);
    Rhs[i] += Mult*(test00+delta*Lu)*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val += delta * (-c0*(ansatz20+ansatz02)
             +c1*ansatz10+c2*ansatz01
             +c3*ansatz00) * Lu;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
}



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
    }                                             // endfor j
  }                                               // endfor i
}


void BilinearAssemble_UPW2(double Mult, double *coeff, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c3, c4;

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0];
  c3 = coeff[3];
  c4 = coeff[4];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += c3*ansatz00*test00;

      val *= Mult;

      MatrixRow[j] += val;
    }                                             // endfor j
  }                                               // endfor i
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
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;
  double delta, bgradv, sigma, norm_b;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);

  sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, delta, param, 0, 0, 0);

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val = (c0+sigma)*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val += delta * (-c0*(ansatz20+ansatz02)
        +c1*ansatz10+c2*ansatz01
        +c3*ansatz00) * bgradv;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
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
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;
  double delta, bgradv, sigma, norm_b,sigma0;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);
  sigma = Compute_SOLD_sigma(hK, c0, c1, c2, c3, c4, c5, delta, param, 0, 0, 0);
  norm_b = c1*c1+c2*c2;

  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==CS99)
  {
    //sigma0 = sigma-delta*hK*sqrt(norm_b)/2;          // this is kappa_sl in [CS99]
    sigma0 = sigma-delta*norm_b;                  // this is kappa_sl in [CS99]
    if (sigma0<0)
      sigma0 = 0;
    sigma -= sigma0;                              // effective orthogonal diffusion
    c0 += sigma0;                                 // effective isotropic diffusion
  }

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val += delta * (-c0*(ansatz20+ansatz02)
        +c1*ansatz10+c2*ansatz01
        +c3*ansatz00) * bgradv;

      if (norm_b >0)
        val += sigma * (-c2*ansatz10+c1*ansatz01)*(-c2*test10+c1*test01)/norm_b;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
}


// Layton, Polman, SIAM Sci. Comput. 17, 1328 - 1346, 1996
void RhsAssemble_LP96(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;
  double delta, bgradv, u;
  double umin = 0, umax = 1, rho;

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  u = param[0];

  delta = Compute_SDFEM_delta(hK, c0, c1, c2, c3, c5);
  rho = hK* TDatabase::ParamDB->SOLD_CONST;
  for(i=0;i<N_;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    bgradv = c1*test10+c2*test01;
    val = (test00+delta*bgradv)*c4;
    val = val - (MIN(u-umin,0)+MAX(u-umax,0))* test00/rho;
    Rhs[i] += Mult*val;
  }                                               // endfor i
}

/******************************************************************************/
// RhsAssemble_RhsAdjointEnergyEstimate
// assemble rhs for adjoint problem with energy error estimator of Verfuerth
// (2005)
/******************************************************************************/

void RhsAssemble_RhsAdjointEnergyEstimate(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01, test20, test02;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4;
  double test, ansatz, scal, u[5];

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // u_xx
  Orig4 = OrigValues[4]; // u_yy

  c0 = coeff[0]; // eps
  c1 = coeff[1]; // b_1
  c2 = coeff[2]; // b_2
  c3 = coeff[3]; // c
  c4 = coeff[4]; // f
  // current finite element solution of equation
  u[0] = param[0]; // u
  u[1] = param[1]; // u_x
  u[2] = param[2]; // u_y
  u[3] = param[3]; // u_xx
  u[4] = param[4]; // u_yy

  // residual eps Delta u - b* nabla u - c u + f
  ansatz = c0 * (u[3] + u[4]) - c1*u[1] - c2*u[2] - c3*u[0] + c4;

  // compute scaling of strong residual in energy error estimator
  scal = hK*hK/c0;
  if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
  {
    // update weight for energy norm estimator
    if (1.0/TDatabase::ParamDB->INTERNAL_COERCIVITY<scal)
      scal = 1.0/TDatabase::ParamDB->INTERNAL_COERCIVITY; 
  }
  // scal ansatz
  ansatz *= 2 * scal;

  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
    test20 = Orig3[i];
    test02 = Orig4[i];
    // test function
    test = c0 * (test20 + test02) - c1*test10 - c2*test01 - c3*test00;
    // update array for rhs
    val = ansatz * test;
    Rhs[i] += Mult*val;
  }                                               // endfor i
}

/******************************************************************************/
// RhsAssemble_RhsAdjointL2Error
// assemble rhs for adjoint problem with L2 error to prescribed solution
/******************************************************************************/
void RhsAssemble_RhsAdjointL2Error(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs;
  double ansatz, test, u, u_h;
  double *Orig0;
  int i, N_;

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u

  // current finite element solution of equation
  u_h = param[0]; // u

  // prescribed solution
  u = coeff[10];

  // difference in quadrature point
  ansatz = -2 * (u - u_h);

  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    test = Orig0[i];
    Rhs[i] += Mult * ansatz * test;
  }                                               // endfor i
}

/******************************************************************************/
// RhsAssemble_RhsAdjointH1Error
// assemble rhs for adjoint problem with H1-semi norm error to prescribed solution
/******************************************************************************/
void RhsAssemble_RhsAdjointH1Error(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs, val;
  double test10, test01;
  double ansatz_x, ansatz_y;
  double *Orig0, *Orig1;
  int i,j, N_;
  double scal, u[2];

  Rhs = LocRhs[0];
  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y

  // current finite element solution of equation
  u[0] = param[1]; // u_x
  u[1] = param[2]; // u_y

  // error of first partial derivatives
  ansatz_x = coeff[11] - u[0];
  ansatz_y = coeff[12] - u[1]; 
  scal = -2; //savescu
  ansatz_x *= scal;
  ansatz_y *= scal;

  // loop over the basis functions (test functions)
  for(i=0;i<N_;i++)
  {
    // test function
    test10 = Orig0[i];
    test01 = Orig1[i];

    // update array for rhs
    val = ansatz_x * test10 + ansatz_y * test01;
    Rhs[i] += Mult*val;
  }                                               // endfor i
}


/******************************************************************************/
//
// IMPROVED MIZUKAMI-HUGHES METHOD (Knobloch, CMAME 2007)
//
// assembling only Galerkin part
//
/******************************************************************************/

void BilinearAssemble_MH_Kno06(double Mult, double *coeff, double *param,
double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **Matrix, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_;
  double c0, c1, c2, c3, c4, c5;

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];
  Orig3 = OrigValues[3];
  Orig4 = OrigValues[4];

  c0 = coeff[0];                                  // nu
  c1 = coeff[1];                                  // b_1
  c2 = coeff[2];                                  // b_2
  c3 = coeff[3];                                  // c
  c4 = coeff[4];                                  // f
  c5 = coeff[5];                                  // \|b\|_infty

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      ansatz20 = Orig3[j];
      ansatz02 = Orig4[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      val += (c1*ansatz10+c2*ansatz01)*test00;
      val += c3*ansatz00*test00;

      val *=Mult;

      MatrixRow[j] += val;

    }                                             // endfor j
  }                                               // endfor i
}

// ========================================================================
// parameters:  H1 norm and  norm of residual
// ========================================================================
void DC_CD_Params(double *in, double *out)
{
  out[0] = in[2];                                 // H1 norm
  out[1] = in[3];                                 // norm of residual
}


// ========================================================================
// parameters:  partial derivatives
// ========================================================================
void Params_All_Deriv(double *in, double *out)
{
  out[0] = in[2];                                 // u
  out[1] = in[3];                                 // u_x
  out[2] = in[4];                                 // u_y
  out[3] = in[5];                                 // u_xx
  out[4] = in[6];                                 // u_yy
}

// ========================================================================
// parameters:  partial derivatives
// ========================================================================
void Params_Sol(double *in, double *out)
{
  out[0] = in[2];                                 // u
}

// ========================================================================
// parameters:  SC_2
// ========================================================================
void SC_2_Params(double *in, double *out)
{
  out[0] = in[2];                                 // H1 norm
  out[1] = in[3];                                 // norm of residual
  out[2] = in[4];                                 // u_x
  out[3] = in[5];                                 // u_y
}

// ========================================================================
// parameters:  SOLD
// ========================================================================
void SOLD_Params(double *in, double *out)
{
  out[0] = in[2];                                 // u
  out[1] = in[3];                                 // u_x
  out[2] = in[4];                                 // u_y
  out[3] = in[5];                                 // u_xx
  out[4] = in[6];                                 // u_yy
  out[5] = in[7];                                 // ||u^h||_{H^1,K}
  out[6] = in[8];                                 // ||R(u^h)||_{L^2,K}
}

/*************************************************************************/
// estimate of coercivity constant 
// used eg for residual based estimator of Verf"uhrt 2005
/*************************************************************************/
double EstimateCoercivityConstant(TCollection *Coll,
				  CoeffFct2D *Coeffs)
{
  int i, j, N_Cells, N_V;
  double coerc = 4711.0, x, y, xs, ys, *coeff;
  TBaseCell *cell;
  TVertex *vertex;

  coeff = new double[13];
    
  N_Cells = Coll->GetN_Cells();                   // number of mesh cells
  // loop over the cells
  for (i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    // loop over the vertices
    xs = ys = 0;
    for (j=0;j<N_V;j++)
    {
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x, y);
      Coeffs(1, &x, &y, NULL, &coeff);
      // assuming that convection is divergence-free
      if (coeff[3] < coerc)
	  coerc = coeff[3];
      xs += x;
      ys += y;
    }
    xs /= N_V;
    ys /= N_V;
    Coeffs(1, &xs, &ys, NULL, &coeff);
    // assuming that convection is divergence-free
    if (coeff[3] < coerc)
	coerc = coeff[3];
  }

  delete coeff;
  OutPut("coercivity constant (assuming div-free convection): " << coerc << endl);
  return(coerc);
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
    if (!((TDatabase::ParamDB->SOLD_PARAMETER_TYPE>=BH04)
      && (TDatabase::ParamDB->SOLD_PARAMETER_TYPE<=BE05_2)))
      switch (TDatabase::ParamDB->SOLD_TYPE)
      {
        case 1:
          OutPut("iso&");
          break;
        case 2:
          OutPut("ortho&");
          break;
        default:
          OutPut("incorrect SOLD_TYPE");
          exit(4711);
      }
      switch (TDatabase::ParamDB->SOLD_PARAMETER_TYPE)
      {
      case JSW87:
        OutPut("[JSW87] & no para. &");
        break;
      case HMM86:
        OutPut("[HMM86] & no para. &");
        break;
      case TP86_1:
        OutPut("[TP86_1] & no para. &");
        break;
      case TP86_2:
        OutPut("[TP86_2] & no para. &");
        break;
      case GdC88:
        OutPut("[GdC88] & no para. &");
        break;
      case dCG91:
        OutPut("[dCG91] & no para. &");
        break;
      case dCA03:
        OutPut("[dCA03] & no para. &");
        break;
      case AS97:
        OutPut("[AS97] & no para. &");
        break;
      case C93:
        OutPut("[C93]& $C = " << setprecision(2) << TDatabase::ParamDB->SOLD_CONST << "$&");
        break;
      case KLR02_1:
        OutPut("[KLR02_1] & $S= " << setprecision(2) << TDatabase::ParamDB->SOLD_S << ", C = " <<
          TDatabase::ParamDB->SOLD_CONST << "$&");
        break;
      case KLR02_2:
        OutPut("[KLR02_2] & $C = " << setprecision(2) <<
          TDatabase::ParamDB->SOLD_CONST << "$&");
        break;
      case KLR02_3:
        OutPut("[KLR02_3] & $C = "  << setprecision(2) <<
          TDatabase::ParamDB->SOLD_CONST << "$&");
        break;
      case J90:
        OutPut("[J90] & al = " <<  setprecision(1) <<
          TDatabase::ParamDB->SOLD_CONST << " &");
        break;
      case BE02_1:
        OutPut("[BE02_1] & no para. &");
        break;
      case BE02_2:
        OutPut("[BE02_2] & no para. &");
        break;
      case BE02_3:
        OutPut("[BE02_3] & no para. &");
        break;
      case BH04:
        OutPut("[BH04] & C = " <<  setprecision(2) <<
          TDatabase::ParamDB->SOLD_S << " &");
        break;
      case BE05_1:
        OutPut("[BE05_1] & C = " <<  setprecision(2) <<
          TDatabase::ParamDB->SOLD_CONST << " &");
        break;
      case BE05_2:
        OutPut("[BE05_2] & C = " <<  setprecision(2) <<
          TDatabase::ParamDB->SOLD_CONST << " &");
        break;
      case CS99:
        OutPut("[CS99] & $C = " << setprecision(2) <<
          TDatabase::ParamDB->SOLD_CONST << "$&");
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
  // old version
  /*if (i>2*max-1)
  {
    i -= 2*max;
    TDatabase::ParamDB->SOLD_TYPE = 2;
  }
  if ((i>max-1)&&(i<2*max))
  {
    i -= max;
    TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR = 0.5;
    }*/
  // new version
  if ((i>max-1)&&(i<2*max))
  {
    i -= max;
    TDatabase::ParamDB->SOLD_TYPE = 2;
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
    // boundary H11
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
    }
    OutPut(rank<< "\\"<<"\\"<<endl);
  }
  return;
}


void EvaluateResults(int *nonlinite, double *scale, double *globmin,
double *globmax, double *intlayer_width, double *bdrlay_osc,
double *bdrlay_smear, int ij_end)
{
  int i, rank;
  double supgmax, supgmin, supgbdrlay_osc, supgintlayer_width, supgbdrlay_smear;
  double bdrlay_smearmin, val;

  supgmax = globmax[0]-1.0;
  supgmin = fabs(globmin[0]);
  supgbdrlay_osc = bdrlay_osc[0]-1.0;
  supgintlayer_width = fabs(intlayer_width[0]);
  supgbdrlay_smear = bdrlay_smear[0];
  bdrlay_smearmin =  supgbdrlay_smear;

  // set reference values
  for (i=0;i<ij_end;i++)
  {
    // scheme not converged
    if (nonlinite[i] == TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR)
      continue;
    if (bdrlay_smear[i] <  bdrlay_smearmin)
      bdrlay_smearmin = bdrlay_smear[i];
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
    /*val = fabs(globmin[i]/supgmin);
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
    OutPut("$\circ$ &");
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
    }*/

    // global maximum
    OutPut(globmax[i] << "&");
    /*val = fabs((globmax[i]-1)/supgmax);
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
    */
    OutPut(1+globmax[i]-globmin[i] << "&");

    // width of interior layer
    OutPut(intlayer_width[i] << "&");
    /*val = fabs(intlayer_width[i]/supgintlayer_width);
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

    // oscilations in boundary layer
    OutPut(bdrlay_osc[i] << "&");
    /*val = fabs((bdrlay_osc[i]-1)/supgbdrlay_osc);
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

    // smearing in boundary layer
    OutPut(bdrlay_smear[i] << "&");
    /*val = fabs((bdrlay_smear[i]-bdrlay_smearmin)/(supgbdrlay_smear-bdrlay_smearmin));
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
    }
    OutPut(rank<< "\\"<<"\\"<<endl);*/
    OutPut("\\"<<"\\"<<endl);
  }
  return;
}


void EdgeStabilization(TFESpace2D *fespace,
TFEFunction2D *u,
CoeffFct2D *Coeffs,
double *rhs,
int time_dependent,
double *time_step,
TFEFunction2D *old_u)
{
  int i, j, k, ii, N_Cells, *ColInd, *RowPtr, *GlobalNumbers, *BeginIndex;
  int ActiveBound, *DOF, N_Edges, boundedge, locdof;
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  double val[3], val_neigh[3], h, norm_t, x[3], y[3], oldval[3];
  double x0, x1, y0, y1, xs, ys, t1, t2, *coeff, jump, fac0, fac1, fac2;
  double phi0_x, phi0_y, phi1_x, phi1_y, phi2_x, phi2_y, n1, n2, maxjump;
  double sx, sy, tmp, meas, area, rho = 2.0;
  TBaseCell *cell, *neigh;
  TCollection *coll;
  FE2D CurrentElement;
  TJoint *joint;
  TRefDesc *refdesc;
  TVertex *ver0,*ver1;
  const int *TmpEdVer;

  coeff = new double[6];

  // get arrays with the numbering of the dof
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  // get start of dirichlet nodes in dof array
  ActiveBound = fespace->GetActiveBound();
  // get collection and number of cells
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // assign a numbering to the cells
  for(i=0;i<N_Cells;i++)                          // do for all mesh cells
  {                                               // on the finest level
    cell=coll->GetCell(i);
    cell->SetClipBoard(i);
  }                                               // endfor i

  // loop over all cells for computing the edge stabilization
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    h = cell->GetDiameter();
    meas = cell->GetMeasure();
    // pointer to global indices of dof connected with this cell
    DOF = GlobalNumbers + BeginIndex[i];

    // local dofs are arranged as follows
    // local dof 0 on vertex 0 opposite to edge 1
    // local dof 1 on vertex 1 opposite to edge 2
    // local dof 2 on vertex 2 opposite to edge 0

    CurrentElement = fespace->GetFE2D(i, cell);
    if (CurrentElement!=C_P1_2D_T_A)
    {
      if (sold_parameter_type!=BE05_2)
      {
        OutPut("Edge stabilization for element " << CurrentElement <<
          " not implemented !!!"<< endl);
        exit(4711);
      }
      if ((CurrentElement!=C_Q1_2D_Q_A)&&(CurrentElement!=C_Q1_2D_Q_M))
      {
        OutPut("Edge stabilization for element " << CurrentElement <<
          " not implemented !!!"<< endl);
        exit(4711);
      }
    }
    // # of edges
    N_Edges = cell->GetN_Edges();

    sx = sy = 0;
    // compute derivatives for basis functions
    for (j=0;j<N_Edges; j++)
    {
      x[j] = cell->GetVertex(j)->GetX();
      y[j] = cell->GetVertex(j)->GetY();
      sx += x[j];
      sy += y[j];
      //OutPut(x[j] << " " << y[j] << " ");
      u->FindGradientLocal(cell, i, x[j], y[j], val);
      //OutPut("u"<<j << " " << val[0]<<endl);
    }
    sx /= N_Edges;
    sy /= N_Edges;
    //OutPut(endl);
    // compute twice area of triangle
    if (N_Edges==3)
    {
      area = 2*meas;
      phi0_x = (y[1]-y[2])/area;
      phi0_y = (x[2]-x[1])/area;
      phi1_x = (y[2]-y[0])/area;
      phi1_y = (x[0]-x[2])/area;
      phi2_x = (y[0]-y[1])/area;
      phi2_y = (x[1]-x[0])/area;
    }
    else
    {
      // Q1
      OutPut("Implementation not complete"<<endl);
      exit(4711);
      area = meas;
      phi0_x = (y[1]-y[2])/area;
      phi0_y = (x[2]-x[1])/area;
      phi1_x = (y[2]-y[0])/area;
      phi1_y = (x[0]-x[2])/area;
      phi2_x = (y[0]-y[1])/area;
      phi2_y = (x[1]-x[0])/area;
    }

    /* OutPut("0 " << phi0_x << " " << phi0_y << endl);
     OutPut("1 " << phi1_x << " " << phi1_y << endl);
     OutPut("2 " << phi2_x << " " << phi2_y << endl);
    */
    // get refinement descriptor
    refdesc=cell->GetRefDesc();
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);

    // compute gradient of current solution (constant)
    u->FindGradientLocal(cell, i, sx, sy, val);

    // compute maximal normal jump
    if (sold_parameter_type==BH04)
    {
      maxjump = 0;
      for(j=0;j<N_Edges;j++)                      // loop over all edges of cell
      {
        joint=cell->GetJoint(j);
        ver0=cell->GetVertex(TmpEdVer[2*j]);      // get vertices of face j
        ver1=cell->GetVertex(TmpEdVer[2*j+1]);
        x0 = ver0->GetX();                        // coordinates of face j
        y0 = ver0->GetY();
        x1 = ver1->GetX();
        y1 = ver1->GetY();

        // compute tangential
        t1 = x1 - x0;
        t2 = y1 - y0;
        norm_t = sqrt(t1*t1+t2*t2);
        t1 /= norm_t;
        t2 /= norm_t;
        // compute normal
        n1 = -t2;
        n2 = t1;
        // compute solution (including derivative) in midpoint of tangential
        // from point of view of this mesh cell
        xs = (x1+x0)/2;
        ys = (y1+y0)/2;

        // compute solution (including derivative) in midpoint of tangential
        // from point of view of neighbour mesh cell
        // NO ADAPTIVE MESHES ALLOWED
        neigh=joint->GetNeighbour(cell);          // neighbour cell
        if (neigh!=NULL)
        {
          ii =  neigh->GetClipBoard();
          u->FindGradientLocal(neigh, ii, xs, ys, val_neigh);
        }
        else
        {
          // boundary edge
          // continue;
          val_neigh[0] = val_neigh[1] = val_neigh[2] = 0;
        }
        jump = (n1 * val[1] + n2 * val[2]) - (n1 * val_neigh[1] + n2 * val_neigh[2]);
        jump = fabs(jump);
        if (jump > maxjump)
          maxjump = jump;
      }
      //OutPut(" " << maxjump);
    }

    for(j=0;j<N_Edges;j++)                        // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      ver0=cell->GetVertex(TmpEdVer[2*j]);        // get vertices of face j
      ver1=cell->GetVertex(TmpEdVer[2*j+1]);
      x0 = ver0->GetX();                          // coordinates of face j
      y0 = ver0->GetY();
      x1 = ver1->GetX();
      y1 = ver1->GetY();
      //OutPut("ed " << j << " " << x0 << " " << y0 << " ; " << x1 << " " <<y1
      //      << endl);
      // compute tangential
      t1 = x1 - x0;
      t2 = y1 - y0;
      norm_t = sqrt(t1*t1+t2*t2);
      t1 /= norm_t;
      t2 /= norm_t;
      //OutPut(t1 << " " << t2 << " " << t1*t1+t2*t2 << endl);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of this mesh cell
      xs = (x1+x0)/2;
      ys = (y1+y0)/2;
      //u->FindGradientLocal(cell, i, xs, ys, val);
      //OutPut("grad_i " << val[1] << " " << val[2] << endl);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of neighbour mesh cell
      // NO ADAPTIVE MESHES ALLOWED
      neigh=joint->GetNeighbour(cell);            // neighbour cell
      if (neigh!=NULL)
      {
        ii =  neigh->GetClipBoard();
        //OutPut("ii " << ii << endl);
        u->FindGradientLocal(neigh, ii, xs, ys, val_neigh);
        boundedge = 0;
      }
      else
      {
        // boundary edge
        val_neigh[0] = val_neigh[1] = val_neigh[2] = 0;
        boundedge = 1;
      }
      //OutPut("grad_ii " << val_neigh[1] << " " << val_neigh[2] << endl);

      // compute factor which defines the sign
      fac0 = t1 * val[1] + t2 * val[2];
      //OutPut("fac 0 " << fac0 << " ");
      /*if (fac0 > 1e-8)
          fac0 = 1;
      else
      {
          if (fac0<-1e-8)
        fac0 = -1;
          else
        fac0 = 0;
        }*/
      //OutPut("fac0 " << fac0<<endl);
      tmp = tanh(fac0/1.0)/fac0;
      fac0 = tanh(fac0/1.0);
      //OutPut(fac0 << endl);
      // compute nonlinear factor depending on u (Psi_K(u))
      switch (sold_parameter_type)
      {
        case BH04:
          // compute coefficients in (xs,ys)
          Coeffs(1, &xs, &ys, NULL, &coeff);
          fac1 = coeff[0] * TDatabase::ParamDB->SOLD_CONST;
          if (time_dependent)
            fac1 *= time_step[0];
          fac1 += h * TDatabase::ParamDB->SOLD_S;
          fac1 *= h *maxjump/meas;
          break;
        case BE05_1:
          // compute coefficients in (xs,ys)
          Coeffs(1, &xs, &ys, NULL, &coeff);
          // norm of convection
          fac1 = sqrt(coeff[1]*coeff[1]+coeff[2]*coeff[2]);
          if (time_dependent)
            fac1 *= time_step[0];
          //OutPut("conv " << fac1 << " " << h << endl);
          fac1 *= h*h;
          // reaction term
          fac2 = rho * fabs(coeff[3]);
          if (time_dependent)
          {
            fac2 *= time_step[0];
            fac2 += 1.0;
            //OutPut(fac2 << " ");
          }
          fac2 *=h*h*h;
          fac1 += fac2;
          // jump of gradient
          jump = (val[1]-val_neigh[1]) * (val[1]-val_neigh[1]);
          jump += (val[2]-val_neigh[2]) * (val[2]-val_neigh[2]);
          fac1 *= sqrt(jump);
          //OutPut("jump " << sqrt(jump) << endl);
          fac1 *= TDatabase::ParamDB->SOLD_CONST/meas;
          break;

        case BE05_2:
          // this case does not depend from formerly computed values
          // Simpson rule
          // compute coefficients in (x0,y0)
          Coeffs(1, &x0, &y0, NULL, &coeff);
          // pw_linear_rhs
          if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 100)
          {
            if ((x0<0.25)||(x0>0.75)||(y0<0.25)||(y0>0.75))
              coeff[4] = 0;
          }
          // residual
          if (!time_dependent)
          {
            fac1 = coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0]-coeff[4];
          }
          else
          {
            old_u->FindGradientLocal(cell, i, x0, y0, oldval);
            fac1 = val[0] - oldval[0]
              + time_step[0] * (coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0])
              + time_step[1] * (coeff[1]*oldval[1]+ coeff[2]*oldval[2]+coeff[3]*oldval[0])
              - time_step[2] * coeff[5]
              - time_step[3] * coeff[4];
          }
          Coeffs(1, &xs, &ys, NULL, &coeff);
          // pw_linear_rhs
          if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 100)
          {
            if ((sx<0.25)||(sx>0.75)||(sy<0.25)||(sy>0.75))
              coeff[4] = 0;
          }
          // residual
          if (!time_dependent)
          {
            fac1 += 4*(coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0]-coeff[4]);
          }
          else
          {
            old_u->FindGradientLocal(cell, i, xs, ys, oldval);
            fac1 += 4 * (val[0] - oldval[0]
              + time_step[0] * (coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0])
              + time_step[1] * (coeff[1]*oldval[1]+ coeff[2]*oldval[2]+coeff[3]*oldval[0])
              - time_step[2] * coeff[5]
              - time_step[3] * coeff[4]);
          }
          // compute coefficients in (x1,y1)
          Coeffs(1, &x1, &y1, NULL, &coeff);
          // pw_linear_rhs
          if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 100)
          {
            if ((x1<0.25)||(x1>0.75)||(y1<0.25)||(y1>0.75))
              coeff[4] = 0;
          }
          // residual
          if (!time_dependent)
          {
            fac1 += coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0]-coeff[4];
          }
          else
          {
            old_u->FindGradientLocal(cell, i, x1, y1, oldval);
            fac1 += val[0] - oldval[0]
              + time_step[0] * (coeff[1]*val[1]+ coeff[2]*val[2]+coeff[3]*val[0])
              + time_step[1] * (coeff[1]*oldval[1]+ coeff[2]*oldval[2]+coeff[3]*oldval[0])
              - time_step[2] * coeff[5]
              - time_step[3] * coeff[4];
          }

          fac1 = TDatabase::ParamDB->SOLD_CONST * fabs(fac1)/6.0;
          break;
        default :
          OutPut("Edge stabilization " << sold_parameter_type
            << " not implemented !!!" <<endl);
          exit(4711);
      }
      // norm_t is the length of the edge
      //OutPut(tmp*fac1/norm_t <<endl);
      fac1 = fac0*meas*fac1*norm_t;
      //OutPut("a_i " << fac1/norm_t << endl );
      // update the rhs
      switch(j)
      {
        // edge zero, active dof are local 0 and 1
        case 0:
          // local dof 0
          locdof = DOF[0];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi0_x + t2*phi0_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          // local dof 1
          locdof = DOF[1];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi1_x + t2*phi1_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          break;
          // edge one, active dof are local 1 and 2
        case 1:
          // local dof 1
          locdof = DOF[1];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi1_x + t2*phi1_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          // local dof 2
          locdof = DOF[2];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi2_x + t2*phi2_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          break;
          // edge two, active dof are local 0 and 2
        case 2:
          // local dof 0
          locdof = DOF[0];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi0_x + t2*phi0_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          // local dof 2
          locdof = DOF[2];
          // do nothing for Dirichlet boundary
          if (locdof< ActiveBound)
          {
            fac2 = t1*phi2_x + t2*phi2_y;
            fac2 *= fac1;
            rhs[locdof] += fac2;
          }
          break;
      }
    }
  }                                               // loop over cells
  /*
    // loop over all cells for computing the edge stabilization
    for(i=0;i<N_Cells;i++)
    {
      // next cell
      cell = coll->GetCell(i);
      // pointer to global indices of dof connected with this cell
      DOF = GlobalNumbers + BeginIndex[i];

      // local dofs are arranged as follows
      // local dof 0 on vertex 0 opposite to edge 1
  // local dof 1 on vertex 1 opposite to edge 2
  // local dof 2 on vertex 2 opposite to edge 0

  // # of edges
  N_Edges = cell->GetN_Edges();

  // compute derivatives for basis functions
  for (j=0;j<N_Edges; j++)
  {
  x[j] = cell->GetVertex(j)->GetX();
  y[j] = cell->GetVertex(j)->GetY();
  OutPut(x[j] << " " << y[j] << " ori " << rhsori[DOF[j]] << " update " << rhs[DOF[j]] <<
  " diff " <<  rhsori[DOF[j]] - rhs[DOF[j]] << endl);
  }
  }
  */
  delete coeff;

}

// THIS WORKS ONLY FOR P1 !!!
void JumpTermsForAdjointProblemP1(TFESpace2D *fespace,
				TFEFunction2D *u,
				CoeffFct2D *Coeffs,
				BoundCondFunct2D *BoundaryConditions,
				double *rhs)
{
  int i, j, k, ii, N_Cells, *ColInd, *RowPtr, *GlobalNumbers, *BeginIndex;
  int ActiveBound, *DOF, *DOF_n, N_Edges, boundedge, locdof, found;
  int com00, com01, com10, com11, com20, com21, com_other0, com_other1;
  int loc_vert_n, comp;
  double val[3], val_neigh[3], h, norm_t, x[3], y[3], oldval[3];
  double x_n[3], y_n[3], eps = 1e-6;
  double x0, x1, y0, y1, xs, ys, t1, t2, *coeff, jump, fac0, fac1, fac2;
  double phi0_x, phi0_y, phi1_x, phi1_y, phi2_x, phi2_y, n1, n2, maxjump;
  double phi0_n_x, phi0_n_y, phi1_n_x, phi1_n_y, phi2_n_x, phi2_n_y;
  double phi_n_other_x, phi_n_other_y, p0, p1;
  double sx, sy, tmp, meas, area, rho = 2.0, ansatz, test, area_n, meas_n;
  TBaseCell *cell, *neigh;
  TCollection *coll;
  FE2D CurrentElement;
  TJoint *joint;
  TRefDesc *refdesc;
  TVertex *ver0,*ver1;
  BoundCond BdCond;
  TBoundComp *BdComp;
  TBoundEdge *bound_edge;
  TIsoBoundEdge *isobound_edge;
  const int *TmpEdVer;

  coeff = new double[13];

  // get arrays with the numbering of the dof
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  // get start of dirichlet nodes in dof array
  ActiveBound = fespace->GetActiveBound();
  // get collection and number of cells
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // assign a numbering to the cells
  for(i=0;i<N_Cells;i++)                          // do for all mesh cells
  {                                               // on the finest level
    cell=coll->GetCell(i);
    cell->SetClipBoard(i);
  }                                               // endfor i

  // loop over all cells for computing the jump terms
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    h = cell->GetDiameter();
    meas = cell->GetMeasure();
    // pointer to global indices of dof connected with this cell
    DOF = GlobalNumbers + BeginIndex[i];

    // local dofs are arranged as follows
    // local dof 0 on vertex 0 opposite to edge 1
    // local dof 1 on vertex 1 opposite to edge 2
    // local dof 2 on vertex 2 opposite to edge 0

    CurrentElement = fespace->GetFE2D(i, cell);
    if (CurrentElement!=C_P1_2D_T_A)
    {
        OutPut("JumpTermsForAdjointProblem for element " << CurrentElement <<
	       " not implemented !!!"<< endl);
        exit(4711);
    }
    // # of edges
    N_Edges = cell->GetN_Edges();

    sx = sy = 0;
    // compute derivatives for basis functions
    for (j=0;j<N_Edges; j++)
    {
      x[j] = cell->GetVertex(j)->GetX();
      y[j] = cell->GetVertex(j)->GetY();
      sx += x[j];
      sy += y[j];
      //OutPut(x[j] << " " << y[j] << " ");
      u->FindGradientLocal(cell, i, x[j], y[j], val);
      //OutPut("u "<<j << " " << val[1]<<endl);
    }
    sx /= N_Edges;
    sy /= N_Edges;
    //OutPut(endl);
    // compute twice area of triangle
    if (N_Edges==3)
    {
      area = 2 * meas;
      phi0_x = (y[1]-y[2])/area;
      phi0_y = (x[2]-x[1])/area;
      phi1_x = (y[2]-y[0])/area;
      phi1_y = (x[0]-x[2])/area;
      phi2_x = (y[0]-y[1])/area;
      phi2_y = (x[1]-x[0])/area;
    }
    /* OutPut("0 " << phi0_x << " " << phi0_y << endl);
     OutPut("1 " << phi1_x << " " << phi1_y << endl);
     OutPut("2 " << phi2_x << " " << phi2_y << endl);
    */
    // get refinement descriptor
    refdesc=cell->GetRefDesc();
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);

    // compute gradient of current solution (constant)
    u->FindGradientLocal(cell, i, sx, sy, val);

    for(j=0;j<N_Edges;j++)                        // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      ver0=cell->GetVertex(TmpEdVer[2*j]);        // get vertices of face j
      ver1=cell->GetVertex(TmpEdVer[2*j+1]);
      x0 = ver0->GetX();                          // coordinates of face j
      y0 = ver0->GetY();
      x1 = ver1->GetX();
      y1 = ver1->GetY();
      //OutPut(endl << "ed " << j << " " << x0 << " " << y0 << " ; " << x1 << " " <<y1
      //    << endl);
      // compute tangential
      t1 = x1 - x0;
      t2 = y1 - y0;
      norm_t = sqrt(t1*t1+t2*t2);
      t1 /= norm_t;
      t2 /= norm_t;
      // compute normal
      n1 = t2;
      n2 = -t1;
      //OutPut(t1 << " " << t2 << " " << t1*t1+t2*t2 << endl);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of this mesh cell
      xs = (x1+x0)/2;
      ys = (y1+y0)/2;
      //OutPut("locd " << x0 << " " << y0 << " : " 
//	     << x1 << " " << y1 << " : "  << DOF[(j)%3] << 
      //     " " << DOF[(j+1)%3] << endl);
      //u->FindGradientLocal(cell, i, xs, ys, val);
      //OutPut("grad_i " << val[1] << " " << val[2] << endl);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of neighbour mesh cell
      // NO ADAPTIVE MESHES ALLOWED
      neigh=joint->GetNeighbour(cell);            // neighbour cell
      if (neigh!=NULL)
      {
        ii =  neigh->GetClipBoard();
        //OutPut("ii " << ii << endl);
	DOF_n = GlobalNumbers + BeginIndex[ii];
        u->FindGradientLocal(neigh, ii, xs, ys, val_neigh);
        boundedge = 0;
	meas_n = neigh->GetMeasure();
	area_n = 2 * meas_n;
	// vertices of neighbour
	for (k=0;k<N_Edges; k++)
	{
	    x_n[k] = neigh->GetVertex(k)->GetX();
	    y_n[k] = neigh->GetVertex(k)->GetY();
	}
	
	// compute derivatives of basis fcts. in neighbour 
	// cell
	switch(j)
	{
	    // edge 0: vertices 0, 1
	    // local dof 0, 1
	    case 0:
		// find common vertices
		for (k=0;k<N_Edges; k++)
		{
		    found = 0;
		    if ((fabs(x_n[k] - x[0])<eps) && (fabs(y_n[k] - y[0])<eps))
		    {
			com00 = (k+1)%3;
			com01 = (k+2)%3;
			found++;
		    }
		    if ((fabs(x_n[k] - x[1])<eps) && (fabs(y_n[k] - y[1])<eps))
		    {
			com10 = (k+1)%3;
			com11 = (k+2)%3;
			found++;
		    }
		    if (!found)
		    {
			com_other0 = (k+1)%3;
			com_other1 = (k+2)%3;
			loc_vert_n = k;
		    }
		}
		phi0_n_x = (y_n[com00]-y_n[com01])/area_n;
		phi0_n_y = (x_n[com01]-x_n[com00])/area_n;
		phi1_n_x = (y_n[com10]-y_n[com11])/area_n;
		phi1_n_y = (x_n[com11]-x_n[com10])/area_n;
		phi2_n_x = 0;
		phi2_n_y = 0;
		phi_n_other_x = (y_n[com_other0] - y_n[com_other1])/area_n;
		phi_n_other_y = (x_n[com_other1] - x_n[com_other0])/area_n;
		//OutPut(phi0_n_x << " " << phi0_n_y << " done " << endl);
		break;
	    // edge 1: vertices 1, 2
	    // local dof 1, 2
	    case 1:
		// find common vertices
		for (k=0;k<N_Edges; k++)
		{
		    found = 0;
		    if ((fabs(x_n[k] - x[1])<eps) && (fabs(y_n[k] - y[1])<eps))
		    {
			com10 = (k+1)%3;
			com11 = (k+2)%3;
			found++;
			//OutPut(x_n[k] << " " << y_n[k] << endl);
		    }
		    if ((fabs(x_n[k] - x[2])<eps) && (fabs(y_n[k] - y[2])<eps))
		    {
			com20 = (k+1)%3;
			com21 = (k+2)%3;
			found++;
		    }
		    if (!found)
		    {
			com_other0 = (k+1)%3;
			com_other1 = (k+2)%3;
			loc_vert_n = k;
		    }		    
		}
		phi0_n_x = 0;
		phi0_n_y = 0;
		phi1_n_x = (y_n[com10]-y_n[com11])/area_n;
		phi1_n_y = (x_n[com11]-x_n[com10])/area_n;
		phi2_n_x = (y_n[com20]-y_n[com21])/area_n;
		phi2_n_y = (x_n[com21]-x_n[com20])/area_n;
		phi_n_other_x = (y_n[com_other0] - y_n[com_other1])/area_n;
		phi_n_other_y = (x_n[com_other1] - x_n[com_other0])/area_n;
		//OutPut(phi_n_other_x << " " << phi_n_other_y << " other " << endl);
		break;
	    // edge 2: vertices 2, 0
	    // local dof 2, 0
	    case 2:
		// find common vertices
		for (k=0;k<N_Edges; k++)
		{
		    found = 0;
		    if ((fabs(x_n[k] - x[0])<eps) && (fabs(y_n[k] - y[0])<eps))
		    {
			com00 = (k+1)%3;
			com01 = (k+2)%3;
			found++;
			//OutPut(x_n[k] << " " << y_n[k] << endl);
		    }
		    if ((fabs(x_n[k] - x[2])<eps) && (fabs(y_n[k] - y[2])<eps))
		    {
			com20 = (k+1)%3;
			com21 = (k+2)%3;
			found++;
		    }
		    if (!found)
		    {
			com_other0 = (k+1)%3;
			com_other1 = (k+2)%3;
			loc_vert_n = k;
		    }		    
		}
		phi0_n_x = (y_n[com00]-y_n[com01])/area_n;
		phi0_n_y = (x_n[com01]-x_n[com00])/area_n;
		phi1_n_x = 0;
		phi1_n_y = 0;
		phi2_n_x = (y_n[com20]-y_n[com21])/area_n;
		phi2_n_y = (x_n[com21]-x_n[com20])/area_n;
		phi_n_other_x = (y_n[com_other0] - y_n[com_other1])/area_n;
		phi_n_other_y = (x_n[com_other1] - x_n[com_other0])/area_n;
		//OutPut(phi0_n_x << " " << phi0_n_y << " done " << endl);
		break;
	}
      }
      else
      {
        // boundary edge
	if (cell->GetJoint(j)->GetType() == BoundaryEdge)
	{
	    bound_edge = (TBoundEdge *)cell->GetJoint(j);
	    BdComp = bound_edge->GetBoundComp();
            bound_edge->GetParameters(p0, p1);	    
	}
	if (cell->GetJoint(j)->GetType() == IsoBoundEdge)
	{
            isobound_edge = (TIsoBoundEdge *)joint;
            BdComp = isobound_edge->GetBoundComp();
            isobound_edge->GetParameters(p0, p1);
	}
	// get id of the boundary component
	comp=BdComp->GetID();
	BoundaryConditions(comp, (p0+p1)/2.0, BdCond);
	if (BdCond==NEUMANN)
	{
	    boundedge = 2;
	    // continuation of solution
	    val_neigh[0] = val_neigh[1] = val_neigh[2] = 0;
	    // continuation of test functions
	    phi0_n_x = phi0_n_y = phi1_n_x = phi1_n_y = phi2_n_x = phi2_n_y;
	}
	if (BdCond==DIRICHLET)
	{
	    boundedge = 1;
	}
      }
      //OutPut("grad_ii " << val_neigh[1] << " " << val_neigh[2] << endl);
      // do nothing for Dirichlet edges
      if (boundedge==1)
	  continue;
      // compute ansatz factor
      fac0 = val[1] * n1 + val[2] * n2;
      fac1 = val_neigh[1] * n1 + val_neigh[2] * n2;
      //OutPut("jump " << val[1] - val_neigh[1] << " " << val[2] - val_neigh[2] << endl);
      ansatz = fac0 - fac1;
      //OutPut("("<<x0<<"," << y0<<") (" << x1 << "," << y1 << "):  normal jump " << ansatz << endl);
      //OutPut(" a "  << ansatz << " " );
      ansatz *= 2;
      Coeffs(1, &xs, &ys, NULL, &coeff);
      ansatz *= coeff[0]*sqrt(coeff[0]);      
      // length of edge
      ansatz *= norm_t; 
      // weight of error estimator
      fac0 = norm_t/sqrt(coeff[0]);
      if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
      {
	  if (1.0/sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY)<fac0)
	      fac0 = 1.0/sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY); // update weight 
      }
      //OutPut(" fac0 "  << fac0 << " ");
      ansatz *= fac0;
      for (k=0;k<N_Edges; k++)
      {
	  // compute jumps of derivative of test function
 	  switch(k)
	  {
	      case 0:
		  fac0 = phi0_x * n1 + phi0_y * n2;
		  fac1 = phi0_n_x * n1 + phi0_n_y * n2;
		  break;
	      case 1:
		  fac0 = phi1_x * n1 + phi1_y * n2;
		  fac1 = phi1_n_x * n1 + phi1_n_y * n2;
		  break;
	      case 2:
		  fac0 = phi2_x * n1 + phi2_y * n2;
		  fac1 = phi2_n_x * n1 + phi2_n_y * n2;
		  break;
	  }
	  test = fac0 - fac1;
	  //OutPut("("<<x0<<"," << y0<<") (" << x1 << "," << y1 << "):  testjump " << test << endl);
 	  test *= ansatz;
	  //OutPut(k << " t " << test << " " );
	  
	  // update the rhs
	  switch(k)
	  {
	      // dof zero 
	      case 0:
		  // local dof 0
		  locdof = DOF[0];
		  break;
		  // dof one
	      case 1:
		  // local dof 1
		  locdof = DOF[1];
		  break;
		  // dof two
	      case 2:
		  // local dof 2
		  locdof = DOF[2];
		  break;
	  }
	  //OutPut(locdof << " integral ("<<x0<<"," << y0<<") (" << x1 << "," << y1 << "):  integral " << test << endl);
	  rhs[locdof] += test;
	  //OutPut("locd " << locdof << " " <<  rhs[locdof] << endl);
      } // end k (test functions)
      if (boundedge)
	  continue;
      // the test functions opposite the edge
      fac1 = phi_n_other_x * n1 + phi_n_other_y * n2;
      test = - fac1;

      test *= ansatz;
      // compute the dof
      locdof = DOF_n[loc_vert_n];
      //OutPut(locdof << " integral ("<<x0<<"," << y0<<") (" << x1 << "," << y1 << "):  integral_bc " 
      //     << test << "::"<< -phi_n_other_x  << " " << -phi_n_other_y <<endl);
      rhs[locdof] += test;
      //OutPut("locd " << locdof << " " <<  rhs[locdof] << endl);
    } // end j (edges)
  }                                               // loop over cells
  /*
    // loop over all cells for computing the edge stabilization
    for(i=0;i<N_Cells;i++)
    {
      // next cell
      cell = coll->GetCell(i);
      // pointer to global indices of dof connected with this cell
      DOF = GlobalNumbers + BeginIndex[i];

      // local dofs are arranged as follows
      // local dof 0 on vertex 0 opposite to edge 1
  // local dof 1 on vertex 1 opposite to edge 2
  // local dof 2 on vertex 2 opposite to edge 0

  // # of edges
  N_Edges = cell->GetN_Edges();

  // compute derivatives for basis functions
  for (j=0;j<N_Edges; j++)
  {
  x[j] = cell->GetVertex(j)->GetX();
  y[j] = cell->GetVertex(j)->GetY();
  OutPut(x[j] << " " << y[j] << " ori " << rhsori[DOF[j]] << " update " << rhs[DOF[j]] <<
  " diff " <<  rhsori[DOF[j]] - rhs[DOF[j]] << endl);
  }
  }
  */
  delete coeff;

}

// =======================================================================
//
// JumpTermsForAdjointProblem
//
// computes jump terms in rhs of adjoint problem for energy norm 
// error estimator
//
// TFESpace2D *fespace                   -- finite element space
// TFEFunction2D *u                      -- finite element function
// CoeffFct2D *Coeff                     -- coefficients of the equation
// BoundCondFunct2D *BoundaryConditions  -- pointer to function for boundary
//                                          conditions
// double *rhs                           -- array for right hand side
//
// =======================================================================

void JumpTermsForAdjointProblem(TFESpace2D *fespace,
				TFEFunction2D *u,
				CoeffFct2D *Coeff,
				BoundCondFunct2D *BoundaryConditions,
				double *rhs)
 {
  const int MaxN_BaseFunctions2D_Ersatz = 100;

  double hK,w,integrand,edge_par,sigma_par,diffusion;
  int out;
  int i,j,k,l,l1,l2,l3,n,n_neigh,m,r,q,dummy,N_UsedElements,N_LocalUsedElements,ii,jj,ll;
  int N_Cells, N_Points, N_Parameters, N_Points1D, N_Edges, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints, n_rhs = 1, n_fespaces = 1;
  int Used[N_FEs2D];
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TBaseFunct2D *bf;
  FE2D *UsedElements, LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1D;
  BaseFunct2D BaseFunctCell;
  TCollection *Coll;
  TBaseCell *cell;
  TJoint *joint;
  TBoundEdge *boundedge;
  TIsoBoundEdge *isoboundedge;
  int **RhsGlobalNumbers, **RhsBeginIndex;
  int **TestGlobalNumbers, **TestBeginIndex;
  int **AnsatzGlobalNumbers, **AnsatzBeginIndex;
  TFE2D *ele;
  TFEDesc2D *FEDesc_Obj;
  BF2DRefElements bf2Drefelements;
  double *weights, *xi, *eta, *weights1D, *weights_neigh, *xi_neigh, *eta_neigh, *weights1D_neigh;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], X_neigh[MaxN_QuadPoints_2D], Y_neigh[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D], AbsDetjk_neigh[MaxN_QuadPoints_2D],*AbsDetjk1D[4];
  double *Param[MaxN_QuadPoints_2D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux, *aux2, *aux3, *aux4;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  double *Coeffs[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries,*Entries1,*Entries2,*Entries3, *Entries4, *Entries5;
  int *ColInd, *RowPtr;
  int *ColInd1, *RowPtr1,*ColInd2, *RowPtr2, *ColInd3, *RowPtr3;
  int *ColInd4, *RowPtr4, *ColInd5, *RowPtr5;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling, v;
  TBoundComp *BoundComp;
  double t0, t1, t, s, integral, u_val[4], u_val_neigh[4];
  double *u_jump_x, *u_jump_y;
  int comp, dof_ii,dof_jj, found;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  BoundValueFunct2D *BoundaryValue;
  TNodalFunctional2D *nf;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D_Ersatz];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, x1, y0, y1, hE, nx, ny, tx, ty, x, y, val, eps=1e-12;
  double penetration_penalty, friction_parameter;
  double **JointValues, *JointValue, u1_values[3], u2_values[3];
  double delta;
  bool *SecondDer;

  double *Coefficients1D[MaxN_QuadPoints_2D];
  double *Parameters1D[MaxN_QuadPoints_2D];

  double xi1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D], eta1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D];
  //double xietaval_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  //double xideriv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  //double etaderiv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double**** xietaval_ref1D = new double*** [N_BaseFuncts2D];
  double**** xideriv_ref1D = new double*** [N_BaseFuncts2D];
  double**** etaderiv_ref1D = new double*** [N_BaseFuncts2D];
  double *xyval_ref1D[4][MaxN_QuadPoints_1D];
  double *xderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *yderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *X1D[4], *Y1D[4], *X1D_neigh[4], *Y1D_neigh[4];
  RefTrans2D RefTrans;
  int N_DOF;
  double *Values;
  //double value_basefunct_ref1D[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz],value_basefunct_ori[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz];
  //double xderiv_basefunct_ref1D[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz],xderiv_basefunct_ori[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz];
  //double yderiv_basefunct_ref1D[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz],yderiv_basefunct_ori[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz];
  double*** value_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** xderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** yderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double *value_basefunct_ori[6];
  double *xderiv_basefunct_ori[6];
  double *yderiv_basefunct_ori[6];
  double x_pos_ref[6];
  double y_pos_ref[6];
  double x_pos[6];
  double y_pos[6];
  double *value_basefunct_ori_neigh[6];
  double *xderiv_basefunct_ori_neigh[6];
  double *yderiv_basefunct_ori_neigh[6];
  double x_pos_neigh[6];
  double y_pos_neigh[6];
  double dummy2[6];

  int neigh_edge;
  int neigh_N_,N_Neigh;
  double absdet1D_neigh[MaxN_QuadPoints_2D];
  double xi1DNeigh[N_BaseFuncts2D][MaxN_QuadPoints_1D], eta1DNeigh[N_BaseFuncts2D][MaxN_QuadPoints_1D];
  double *X1DNeigh,*Y1DNeigh;
  TBaseCell *neigh;
  FE2D LocalUsedElements_neigh[N_FEs2D], CurrEleNeigh;
  BaseFunct2D BaseFunctNeigh;
  QuadFormula2D QuadFormulaNeigh;
  TQuadFormula2D *qfNeigh;
  QuadFormula1D LineQuadFormulaNeigh;
  TQuadFormula1D *qf1DNeigh;
  int LocN_BF_neigh[N_BaseFuncts2D];
  BaseFunct2D LocBF_neigh[N_BaseFuncts2D];
  int N_Points1DNeigh,N_PointsNeigh;
  double *weights1DNeigh,*zetaNeigh,*weightsNeigh,*xiNeigh,*etaNeigh;
  TFE2D *eleNeigh;
  RefTrans2D RefTransNeigh;
  BF2DRefElements bf2DrefelementsNeigh;
  int *DOF_neigh;
  double xietaval_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double xideriv_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double etaderiv_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double *xyval_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_Neigh1D, *yderiv_Neigh1D, *xyval_Neigh1D;

  double jump_xyval[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_yderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];

  out=1;
  OutPut("JumpTerms"<<endl);

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  // get basis functions from data base
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  // get number of basis functions from data base
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // for the right hand side (n_rhs == 1)
  if(n_rhs)
  {
      // get information on the global degrees of freedom 
      RhsBeginIndex = new int* [n_rhs];
      RhsGlobalNumbers = new int* [n_rhs];
      for(i=0;i<n_rhs;i++)
      {
	  RhsBeginIndex[i] = fespace->GetBeginIndex();
	  RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();
      }                                             // endfor
      // allocate arrays for the local computation of the right hand side
      LocRhs = new double* [n_rhs];
      righthand = new double [n_rhs*MaxN_BaseFunctions2D];
      for(i=0;i<n_rhs;i++)
	  LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs
  
  // no second derivatives necessary
  SecondDer = new bool[1];
  SecondDer[0] = FALSE;

  // for all possible basis functions
  // allocate arrays for values on the actual mesh cells
  for (i=0;i<N_BaseFuncts2D;i++)
  {
    value_basefunct_ref1D[i] = new double* [6];
    xderiv_basefunct_ref1D[i] = new double* [6];
    yderiv_basefunct_ref1D[i] = new double* [6];
    for (j=0;j<6;j++)
    {

      value_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];
      xderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];
      yderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];

      memset( value_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      memset( xderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      memset( yderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );

    }
  }

  // for all possible basis functions
  // allocate arrays for values on the reference mesh cell
  for (i=0;i<N_BaseFuncts2D;i++)
  {
    xietaval_ref1D[i] = new double** [4];
    xideriv_ref1D[i] = new double** [4];
    etaderiv_ref1D[i] = new double** [4];
    for (j=0;j<4;j++)
    {
      xietaval_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      xideriv_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      etaderiv_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      for (n=0;n<MaxN_QuadPoints_1D;n++)
      {
        xietaval_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];
        xideriv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];
        etaderiv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];

        memset( xietaval_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
        memset( xideriv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
        memset( etaderiv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      }
    }
  }

  memset(Used, 0, N_FEs2D*SizeOfInt);

  // for the finite element space
  for(i=0;i<n_fespaces;i++)
  {
    n = fespace->GetN_UsedElements();             /* # used finite elements */
    UsedElements = fespace->GetUsedElements();    /* used finite elements */
    for(j=0;j<n;j++)                              /* for all finite elements */
    {
      CurrentElement = UsedElements[j];
      Used[CurrentElement] = 1;
    }                                             // enfor j
  }                                               // endfor i

  N_UsedElements = 0;                             /* compute number of used elements */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i]) N_UsedElements++;

  UsedElements = new FE2D[N_UsedElements];        /* store used finite elements */
  j=0;                                            /* in array */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i])
  {
    UsedElements[j] = (FE2D)i;
    if (out==2)
	OutPut("element " << UsedElements[j] << endl);
    j++;
  }                                               // endif

  // ########################################################################
  // calculate values of base functions and derivatives on ref element
  // ########################################################################
  if (out==2)
  {
      OutPut("N_UsedElements: " << N_UsedElements << " " << endl); fflush(0);
      OutPut("N_BaseFuncts2D: " << N_BaseFuncts2D << " " << endl);
      OutPut("MaxN_QuadPoints_1D: " << MaxN_QuadPoints_1D << " " << endl);
      OutPut("MaxN_BaseFunctions2D_Ersatz: " << MaxN_BaseFunctions2D_Ersatz << " " << endl);
  }

  for(n=0;n<N_UsedElements;n++)                   // for used finite elements
  {
    CurrentElement = UsedElements[n];
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
    LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1D->GetFormulaData(N_Points1D, weights1D, zeta);
    BaseFunctCell = BaseFuncts[CurrentElement];
                                                  // get base functions
    bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctCell);
    bf2Drefelements = bf->GetRefElement();
    switch(bf2Drefelements)                       // compute coordinates of line quadrature
    {                                             // points in reference cell
      // quadrilateral cell
      case BFUnitSquare :                         // edge 0

        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][0][j] = zeta[j];
          eta1D[BaseFunctCell][0][j] = -1;
          bf->GetDerivatives(D00, zeta[j], -1, xietaval_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D10, zeta[j], -1, xideriv_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D01, zeta[j], -1, etaderiv_ref1D[BaseFunctCell][0][j]);
        }                                         // edge 1
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][1][j] = 1;
          eta1D[BaseFunctCell][1][j] = zeta[j];
          bf->GetDerivatives(D00, 1, zeta[j], xietaval_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D10, 1, zeta[j], xideriv_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D01, 1, zeta[j], etaderiv_ref1D[BaseFunctCell][1][j]);
        }                                         // edge 2
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][2][j] = -zeta[j];
          eta1D[BaseFunctCell][2][j] = 1;
          bf->GetDerivatives(D00, -zeta[j], 1, xietaval_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D10, -zeta[j], 1, xideriv_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D01, -zeta[j], 1, etaderiv_ref1D[BaseFunctCell][2][j]);
        }                                         // edge 3
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][3][j] = -1;
          eta1D[BaseFunctCell][3][j] = -zeta[j];
          bf->GetDerivatives(D00, -1, -zeta[j], xietaval_ref1D[BaseFunctCell][3][j]);
          bf->GetDerivatives(D10, -1, -zeta[j], xideriv_ref1D[BaseFunctCell][3][j]);
          bf->GetDerivatives(D01, -1, -zeta[j], etaderiv_ref1D[BaseFunctCell][3][j]);
        }
        break;

      case BFUnitTriangle :                       // triangular cell

        bf->GetDerivatives(D00, 0, 0, value_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D10, 0, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[0] = 0;
        y_pos_ref[0] = 0;
        bf->GetDerivatives(D00, 1, 0, value_basefunct_ref1D[BaseFunctCell][1]);
        bf->GetDerivatives(D10, 1, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 1, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[1] = 1;
        y_pos_ref[1] = 0;
        bf->GetDerivatives(D00, 0, 1, value_basefunct_ref1D[BaseFunctCell][2]);
        bf->GetDerivatives(D10, 0, 1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0, 1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[2] = 0;
        y_pos_ref[2] = 1;

        bf->GetDerivatives(D00, 0.5, 0, value_basefunct_ref1D[BaseFunctCell][3]);
        bf->GetDerivatives(D10, 0.5, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0.5, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[3] = 0.5;
        y_pos_ref[3] = 0;
        bf->GetDerivatives(D00, 0.5, 0.5, value_basefunct_ref1D[BaseFunctCell][4]);
        bf->GetDerivatives(D10, 0.5, 0.5, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0.5, 0.5, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[4] = 0.5;
        y_pos_ref[4] = 0.5;
        bf->GetDerivatives(D00, 0, 0.5, value_basefunct_ref1D[BaseFunctCell][5]);
        bf->GetDerivatives(D10, 0, 0.5, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0, 0.5, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[5] = 0;
        y_pos_ref[5] = 0.5;

        for (j=0;j<N_Points1D;j++)                // for all quadrature poin
        {
          xi1D[BaseFunctCell][0][j] = (zeta[j]+1)/2;
          eta1D[BaseFunctCell][0][j] = 0;
          bf->GetDerivatives(D00, (zeta[j]+1)/2, 0, xietaval_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D10, (zeta[j]+1)/2, 0, xideriv_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D01, (zeta[j]+1)/2, 0, etaderiv_ref1D[BaseFunctCell][0][j]);
        }                                         // edge 1
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][1][j] = (-zeta[j]+1)/2;
          eta1D[BaseFunctCell][1][j] = (zeta[j]+1)/2;
          bf->GetDerivatives(D00, (-zeta[j]+1)/2, (zeta[j]+1)/2, xietaval_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D10, (-zeta[j]+1)/2, (zeta[j]+1)/2, xideriv_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D01, (-zeta[j]+1)/2, (zeta[j]+1)/2, etaderiv_ref1D[BaseFunctCell][1][j]);
        }                                         // edge 2
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][2][j] = 0;
          eta1D[BaseFunctCell][2][j] = (-zeta[j] +1)/2;
          bf->GetDerivatives(D00, 0, (-zeta[j]+1)/2, xietaval_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D10, 0, (-zeta[j]+1)/2, xideriv_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D01, 0, (-zeta[j]+1)/2, etaderiv_ref1D[BaseFunctCell][2][j]);
        }
        break;
    }
  }                                               // endfor n
  if (out==2)
      OutPut("basefunct" << endl);

  for(l=0;l<6;l++)
  {
    value_basefunct_ori[l] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    value_basefunct_ori_neigh[l] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
  }

  for(m=0;m<4;m++)                                // arrays for coordinates, values and
  {                                               // determinant for 1D quadrature
    X1D[m] = new double[N_Points1D];              // coordinates of edge i
    Y1D[m] = new double[N_Points1D];
                                                  // determinant of affine mapping
    AbsDetjk1D[m] = new double[MaxN_QuadPoints_2D];
    for (j=0;j<N_Points1D;j++)                    // arrays for values in reference cell
    {
      xyval_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      xderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      yderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
    }
  }                                               // endfor m
  // arrays for jump of the gradient
  u_jump_x = new double[N_Points1D]; 
  u_jump_y = new double[N_Points1D]; 

  for (j=0;j<N_Points1D;j++)                      // arrays for values in reference cell
  {
    xyval_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
  }

  // ########################################################################
  // Arrays for Parameters
  // ########################################################################

  if (out==2)
      OutPut("coeff" << endl);
  // 20 <= number of term
  aux2 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Coeffs[j] = aux2 + j*20;

  aux4 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Coefficients1D[j] = aux4 + j*20;

  // ########################################################################
  // prepare loop over cells
  // ########################################################################

  // all spaces use same Coll
  Coll = fespace->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)                          // set clipboard of cells on finest
  {
    cell=Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // ########################################################################
  // loop over all cells
  // ########################################################################
  for(i=0;i<N_Cells;i++)                          // for all cells on the finest level
  {
    cell = Coll->GetCell(i);                      // next cell

    if (out==2)
	OutPut("cell " << i << endl);

  // calculate all needed derivatives of this FE function
  CurrentElement = fespace->GetFE2D(i,cell);  // finite element on cell

  BaseFunctCell = BaseFuncts[CurrentElement]; // basis functions
  N_ = N_BaseFunct[CurrentElement];           // # basis functions
  DOF = RhsGlobalNumbers[0] + RhsBeginIndex[0][i];  // dof of current mesh cell
  if (out==2)
  {
      OutPut("N_ basefct " << N_ << " dof: ");
      for (j=0;j<N_;j++)
	  OutPut(DOF[j] << " ");
      OutPut(endl);
  }
  LocalUsedElements[0] = CurrentElement;
  LocN_BF[0] = N_BaseFunct[CurrentElement];   // local basis functions
  LocBF[0] = BaseFuncts[CurrentElement];
  SecondDer[0] = FALSE;
  RefTrans = TFEDatabase2D::GetOrig(1, LocalUsedElements,
				    Coll, cell, SecondDer,
				    N_Points, xi, eta, weights, X, Y, AbsDetjk);

  // get coefficients of pde
  if(Coeff) Coeff(N_Points, X, Y, Param, Coeffs);
  // prepare 1D quadrature formula
  l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
  if(out==2)
      OutPut("Polynomial degree on cell: " << l << endl);
  LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
  qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  qf1D->GetFormulaData(N_Points1D, weights1D, zeta);

  if(out==2)
  {
      for(j=0;j<N_Points1D; j++)
      {
          OutPut("weights1D["<<j<<"]:" <<  weights1D[j] << endl);
      }
      OutPut(endl);
  }

  // update data base
  TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
      ->MakeRefElementData(LineQuadFormula);
  N_Edges=cell->GetN_Edges();                 // # edges
  if (out==2)
  {
      for(r=0;r<N_Edges;r++)
      {
	  cell->GetVertex(r)->GetCoords(x0, y0);
	  cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
	  OutPut("Local edge r: " << r << " vertex A " << x0 << " " << y0 << " vertex B " << 
		 x1 << " " << y1 << endl);
      }
  }
 

  for(r=0;r<N_Edges;r++)                      // loop over all edges of cell
  {                                           // get original coordinates of edge quad. points
      TFEDatabase2D::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunctCell][r],
				    eta1D[BaseFunctCell][r],
				    X1D[r], Y1D[r], AbsDetjk1D[r]);
      
      for(j=0;j<N_Points1D;j++)                 // get values and derivatives in original cell
      {
          TFEDatabase2D::GetOrigValues(RefTrans, xi1D[BaseFunctCell][r][j],
				       eta1D[BaseFunctCell][r][j],
				       TFEDatabase2D::GetBaseFunct2D(BaseFunctCell),
				       Coll, (TGridCell *)cell,
				       xietaval_ref1D[BaseFunctCell][r][j],
				       xideriv_ref1D[BaseFunctCell][r][j],
				       etaderiv_ref1D[BaseFunctCell][r][j],
				       xyval_ref1D[r][j],
				       xderiv_ref1D[r][j],
				       yderiv_ref1D[r][j]);
      }
      
  }                                           // endfor r
  
  TFEDatabase2D::GetOrigFromRef(RefTrans,6,x_pos_ref,y_pos_ref,x_pos,y_pos,dummy2);
  for(l=0;l<6;l++)
  {
      TFEDatabase2D::GetOrigValues(RefTrans, x_pos_ref[l],
				   y_pos_ref[l],
				   TFEDatabase2D::GetBaseFunct2D(BaseFunctCell),
				   Coll, (TGridCell *)cell,
				   value_basefunct_ref1D[BaseFunctCell][l],
				   xderiv_basefunct_ref1D[BaseFunctCell][l],
				   yderiv_basefunct_ref1D[BaseFunctCell][l],
				   value_basefunct_ori[l],
				   xderiv_basefunct_ori[l],
				   yderiv_basefunct_ori[l]);
      // OutPut("Hallo: x_pos_ref[l]: " << x_pos_ref[l] << "value_basefunct_ref1D[BaseFunctCell][l]: " << value_basefunct_ref1D[BaseFunctCell][l] << endl);
  }

  for(r=0;r<N_Edges;r++)
  {                                           
      // for each edge, get the corresponding neighbour cell.
      neigh=cell->GetJoint(r)->GetNeighbour(cell);
      //#######################################################################//
      // get coefficients on edges
      //only implemented for coeffs that do not depend on the params
      //#######################################################################//
      
      //  if(N_Parameters>0)                // get parameters of equ.
      // Parameters->GetParameters(N_Points1D, Coll, cell, i, xi1D[BaseFunctCell][r], eta1D[BaseFunctCell][r], X1D[r], Y1D[r], Param1D);
      
      if(Coeff) Coeff(N_Points1D, X1D[r], Y1D[r], Param, Coefficients1D);
      //#######################################################################//
      // If there is a neighbour to the edge, do...
      if(neigh)
      {                                         
          // get the number of this neigbbour cell from the clipboard
          q = neigh->GetClipBoard();
	  if (out==2)
	      OutPut("neighbor " << q << endl);
	  // only neighbors with larger number
	  // to ensure that each inner edge is treated only once
          if(1)//i<q)
          {
	      // calculate all needed derivatives of this FE function
	      // finite element on neighbour
	      CurrEleNeigh = fespace->GetFE2D(q,neigh);
	      BaseFunctNeigh = BaseFuncts[CurrEleNeigh];
	      eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh);
	      //BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbour
	      N_Neigh = eleNeigh->GetN_DOF();       // number of basis functions on neighbour
	      // dof of current mesh cell on neighbour cell
	      DOF_neigh = RhsGlobalNumbers[0] + RhsBeginIndex[0][q];
	      
	      if (out==2)
	      {
		  OutPut("neigh N_ basefct " << N_Neigh << " dof: ");
		  for (j=0;j<N_;j++)
		      OutPut(DOF_neigh[j] << " ");
		  OutPut(endl);
	      }
	      LocalUsedElements_neigh[0] = CurrEleNeigh;
	      // local basis functions
	      LocN_BF_neigh[0] = N_BaseFunct[CurrEleNeigh];
	      LocBF_neigh[0] = BaseFuncts[CurrEleNeigh];
	      
	      RefTransNeigh = TFEDatabase2D::GetOrig(1, LocalUsedElements_neigh,
						     Coll, neigh, SecondDer,
						     N_Points, xi_neigh, eta_neigh, 
						     weights_neigh, X_neigh, Y_neigh, AbsDetjk_neigh);
	      
	      // get edge of the neighbour cell which is the edge r
	      neigh_edge=0;
	      while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh)!=cell) neigh_edge ++;
	      if (out==2)
		  OutPut("neigh_edge " << neigh_edge << endl);

	      // arrays for coordinates on neighbour cell
	      for(m=0;m<4;m++)                      
	      {
		  X1D_neigh[m] = new double[N_Points1D];
		  Y1D_neigh[m] = new double[N_Points1D];
	      }
	      
	      // get original coordinates of edge quad. points of neighbour cell
	      TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1D[BaseFunctNeigh][neigh_edge],
					    eta1D[BaseFunctNeigh][neigh_edge],
					    X1D_neigh[neigh_edge], Y1D_neigh[neigh_edge], 
					    AbsDetjk1D[neigh_edge]);
	      
	      // get values and derivatives on original neighbour cell on edge neigh_edge
	      for (j=0;j<N_Points1D;j++)
	      {                                     
		  if(out==2)
		  {  
		      OutPut("X1D[r][j]: " << X1D[r][j] << " Y1D[r][j]: "  <<  Y1D[r][j] <<  
			     " X1D[neigh_edge][j] " << X1D_neigh[neigh_edge][j] << 
			     " Y1D[neigh_edge][j]: " <<  Y1D_neigh[neigh_edge][j] <<  endl);
		  }
	      }
	      
	      if(X1D_neigh[neigh_edge][0] == X1D[r][0] && Y1D_neigh[neigh_edge][0] == Y1D[r][0] )
	      {
		  if(out==2)
		  {
		      OutPut("Quadrature points on neighbour edge in the correct order." << endl);
		  }
		  for (j=0;j<N_Points1D;j++)
		  {
		      TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1D[BaseFunctNeigh][neigh_edge][j],
						   eta1D[BaseFunctNeigh][neigh_edge][j],
						   TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
						   Coll, (TGridCell *)neigh,
						   xietaval_ref1D[BaseFunctNeigh][neigh_edge][j],
						   xideriv_ref1D[BaseFunctNeigh][neigh_edge][j],
						   etaderiv_ref1D[BaseFunctNeigh][neigh_edge][j],
						   xyval_refNeigh1D[j],
						   xderiv_refNeigh1D[j],
						   yderiv_refNeigh1D[j]);
		  }                                   //endfor j
	      }                                     //endif
	      else
	      {
		  if(out==2)
		  {
		      OutPut("Inverse the order of the quadrature points on neighbour edge !" << endl);
		  }
		  for (j=0;j<N_Points1D;j++)
		  {
		      TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1D[BaseFunctNeigh][neigh_edge][j],
						   eta1D[BaseFunctNeigh][neigh_edge][j],
						   TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
						   Coll, (TGridCell *)neigh,
						   xietaval_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
						   xideriv_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
						   etaderiv_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
						   xyval_refNeigh1D[j],
						   xderiv_refNeigh1D[j],
						   yderiv_refNeigh1D[j]);
		  }                                   //endfor j
	      }                                     //endelse
	      
	      TFEDatabase2D::GetOrigFromRef(RefTransNeigh,6,x_pos_ref,y_pos_ref,
					    x_pos_neigh,y_pos_neigh,dummy2);

	      for(l=0;l<6;l++)
	      {
		  TFEDatabase2D::GetOrigValues(RefTrans, x_pos_ref[l],
					       y_pos_ref[l],
					       TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
					       Coll, (TGridCell *)neigh,
					       value_basefunct_ref1D[BaseFunctNeigh][l],
					       xderiv_basefunct_ref1D[BaseFunctNeigh][l],
					       yderiv_basefunct_ref1D[BaseFunctNeigh][l],
					       value_basefunct_ori_neigh[l],
					       xderiv_basefunct_ori_neigh[l],
					       yderiv_basefunct_ori_neigh[l]);
	      }

	      if (out>2)
	      {
		  for(k=0;k<N_; k++)
		  {
		      for(l=0;l<6;l++)
		      {
			  OutPut("basis fct: "<< DOF[k] << 
				 " (x,y)-coordinate: " << x_pos[l] << " " << y_pos[l] << 
				 " value: " <<  value_basefunct_ori[l][k] << endl);
		      }
		  }
	      }
	      
	      for(k=0;k<N_Neigh; k++)
	      {
		  if (out>2)
		  {
		      for(l=0;l<6;l++)
		      {
			  if(out>2)
			  {
			      OutPut("Basisfkt neigh: "<< DOF_neigh[k] <<"  (x,y)-coordinate: " << 
				     x_pos_neigh[l] << " " << y_pos_neigh[l] << " value: " <<  
				     value_basefunct_ori_neigh[l][k] << endl);
			  }
		      }
		  }
	      }                                     //endfor k
	    
	      // compute the jumps of the basis functions
	      // and of their derivatives in the quadrature points on edge r
	      // first for the basis functions of cell i
	      
	      for(k=0;k<N_; k++)
	      {
		  dummy = 0;
		  l=0;
		  // Check if basis function k of cell i is in the FE-Space of neighbour cell q
		  while(l<N_Neigh && dummy == 0)
		  {
		      if(DOF[k] == DOF_neigh[l])dummy=1;
		      l++;
		  }
		  l = l-1;
		  // if basis function k of cell i is in the local FE-Space of neighbour cell q do
		  if(dummy ==1 )
		  {   
                      // assumption: N_Points1D cell =  N_Points1D neighbour !!!!
		      for(j=0;j<N_Points1D;j++)
		      {
			  jump_xyval[j][k] = xyval_ref1D[r][j][k]  -  xyval_refNeigh1D[j][l];
			  jump_xderiv[j][k] = xderiv_ref1D[r][j][k] - xderiv_refNeigh1D[j][l];
			  jump_yderiv[j][k] = yderiv_ref1D[r][j][k] - yderiv_refNeigh1D[j][l];
			  
			  if(out>2)
			  {
			      OutPut("xyval_cell: "<< xyval_ref1D[r][j][k] <<" xyval_Neigh: " << xyval_refNeigh1D[j][l] << " jump= " <<  jump_xyval[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i << " on edge (local): " << r << " in quadrature point: " << j << endl);
			      OutPut("xderiv_cell: "<< xderiv_ref1D[r][j][k] <<  " xderiv_Neigh: " << xderiv_refNeigh1D[j][l] << " jump= " << jump_xderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
			      OutPut("yderiv_cell: "<< yderiv_ref1D[r][j][k] <<  " yderiv_Neigh: " << yderiv_refNeigh1D[j][l] << " jump= " << jump_yderiv[j][k] << " of basefunction: "<< DOF[k]  << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
			      OutPut(endl);
			  }
		      }
		  }
		  //endif

		  // if basis function k of cell i is NOT in the local FE-Space of neighbour cell q 
		  // extend them by zero to cell q
		  if (dummy == 0)
		  {
		      for(j=0;j<N_Points1D;j++)
		      {
			  jump_xyval[j][k]  = xyval_ref1D[r][j][k] ;
			  jump_xderiv[j][k] = xderiv_ref1D[r][j][k];
			  jump_yderiv[j][k] = yderiv_ref1D[r][j][k];
			  
			  if(out>2)
			  {
			      OutPut("No Neighbour: xyval_cell: "<< xyval_ref1D[r][j][k] << " jump= " <<  jump_xyval[j][k] <<  " of basefunction: "<< DOF[k] << " in cell: " << i << " on edge (local): " << r << " in quadrature point: " << j << endl);
			      OutPut("No Neighbour: x-deriv-jump= " <<  jump_xderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
			      OutPut("No Neighbour: y-deriv-jump= " <<  jump_yderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << "\n" << endl);
			      OutPut(endl);
			  }
		      }                                 //endfor j
		  }                                   //endif
	      }                                     //endfor k
	      
	      // now for the basis functions of neighbour cell q
	      // which are not in the local FE-Space of cell i
	      for(l=0;l<N_Neigh; l++)
	      {
		  dummy = 0;
		  k=0;
		  while(k<N_ && dummy == 0 )
		  {
		      if(DOF_neigh[l] == DOF[k]) dummy=1 ;
		      k++;
		  }
		  k=k-1;

		  // if basis function l of neighbour cell q is NOT in the local FE-Space of cell i do
		  // extend them by zero to cell i
		  if( dummy == 0)
		  {		      
		      for(j=0;j<N_Points1D;j++)
		      {
			  jump_xyval[j][l+N_] = -xyval_refNeigh1D[j][l] ;
			  jump_xderiv[j][l+N_]= -xderiv_refNeigh1D[j][l];
			  jump_yderiv[j][l+N_]= -yderiv_refNeigh1D[j][l];
			  if(out>2)
			  {
			      OutPut("Neighbour!!" << "xyval_Neigh: " << xyval_refNeigh1D[j][l] << " jump= " <<  jump_xyval[j][l+N_]<<  " of basefunction: "<< DOF_neigh[l] << " in cell: " << q << " on edge (local): " << r << " in quadrature point: " << j << endl);
			      OutPut("Neighbour!! " << "x-deriv-jump: " << jump_xderiv[j][l+N_] << " of basefunction: "<< DOF_neigh[l] << " in cell: " << q <<  " on edge: " << r << " in quadrature point: " << j << endl);
			      OutPut("Neighbour!! y-deriv-jump= " <<  jump_yderiv[j][l+N_]<< " of basefunction: "<< DOF_neigh[l] << " in cell: " << q <<  " on edge: " << r << " in quadrature point: " << j << "\n" << endl);
			      OutPut(endl);
			  }
		      }                                 //endfor j
		  }                                   //endif
	      }                                     //endfor l
	      
	      // compute jumps of gradient of current solution in the
	      // quadrature points on the edge
	      for(j=0;j<N_Points1D;j++)
	      {
		  // quadrature point s
		  if (out==2)
		      OutPut("quad point " << X1D_neigh[neigh_edge][j] << 
			     " " << Y1D_neigh[neigh_edge][j] << endl);
		  // compute values on cell i
		  u->FindGradientLocal(cell,i,X1D_neigh[neigh_edge][j],Y1D_neigh[neigh_edge][j],
				       u_val);
		  // compute values on cell q
		  u->FindGradientLocal(neigh,q,X1D_neigh[neigh_edge][j],Y1D_neigh[neigh_edge][j],
				       u_val_neigh);
		  // compute jumps (value of cell i - value of q)
		  // same direction as for basis functions
		  u_jump_x[j] = u_val[1] - u_val_neigh[1];
		  u_jump_y[j] = u_val[2] - u_val_neigh[2];
		  
		  if (out==2)
		      OutPut("jump of gradient " << u_jump_x[j] << " " << u_jump_y[j] << endl);
	      }
	      // #################################################################################
	      // Compute the edge integrals with the jumps of the basis functions and their derivatives
	      // #################################################################################
	      
	      // get vertices of the edge
	      // the edge is counter-clockwise wrt cell i
	      cell->GetVertex(r)->GetCoords(x0, y0);
	      cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
	      
	      if(out==2)
		  OutPut("vertex A " << x0 << " " << y0 << " vertex B " << x1 << " " << y1 << endl);
	      // compute length of the edge
	      hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
	      // compute normal vector to this edge (normalized)
	      // this normal is pointing outward of cell i (an inward of cell q)
	      nx = (y1-y0)/hE;
	      ny = (x0-x1)/hE;
	      // tangential normal vector to this edge (normalized)
	      // directed from x0 to x1
	      tx = (x1-x0)/hE;
	      ty = (y1-y0)/hE;
	      
	      // compute weight of the jump term
	      sigma_par = TDatabase::ParamDB->INTERNAL_COERCIVITY;
	      diffusion = Coeffs[0][0];
	      edge_par = hE/sqrt(diffusion);
	      if (sigma_par > 0)
	      {
		  sigma_par = 1/sqrt(sigma_par);
		  if (sigma_par < edge_par)
		      edge_par = sigma_par;
	      }
	      edge_par = 2 * edge_par/sqrt(diffusion);
	      if (out==2){ OutPut("weight of jump term: " << edge_par << " " << Coeffs[0][0] << endl)};
	      
	      ActiveBound = fespace->GetActiveBound();
	      if(out>2)
	      {
		  OutPut("ActiveBound "<< ActiveBound  << endl);
		  for(j=0;j<N_Points1D; j++)
		  {
		      OutPut("Det["<<j<<"]:" <<  AbsDetjk1D[r][j] <<" b1: "<<Coefficients1D[j][1] <<" b2: " << Coefficients1D[j][2] << endl);
		  }
	      }
	  
	      // edge integral: test function from cell
	      if(out==2)
	      {  
		  OutPut("testfct. cell" << endl);
	      }
	      for (ii=0;ii<N_;ii++)                 //ii - test function of cell
	      {
		  // look for 'ii'-th row in all matrices
		  dof_ii = DOF[ii];
		  // Dirichlet node
		  if (dof_ii>=ActiveBound)
		      continue;
		  
		  // initialize edge integral
		  integral=0;
		  // compute edge integral
		  for (j=0;j<N_Points1D;j++)       
		  {
		      if (out==2)
		      {
			  OutPut("(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): normal jump " << u_jump_x[j]*nx + u_jump_y[j]*ny << endl);
		      OutPut("(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): testjump " << jump_xderiv[j][ii]*nx + jump_yderiv[j][ii]*ny << " : " << jump_xderiv[j][ii]*nx << " " <<  jump_yderiv[j][ii]*ny << endl);
		      }
		      integrand = edge_par * diffusion * diffusion 
			  * (u_jump_x[j]*nx + u_jump_y[j]*ny) *
			  (jump_xderiv[j][ii]*nx + jump_yderiv[j][ii]*ny);
		      // hE/2 : determinant of reference transformation
		      w = weights1D[j]*hE/2;
		      integral += w*integrand;        // integral on the edge
		  }
		  if (out==2)
		      OutPut(dof_ii << " integral " << "(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): "
			     << integral << endl);
		  // update rhs
		  rhs[dof_ii] += integral;
	      }                                     // end outer loop over dof (ii)
	      
	      // edge integral: test function from neigh cell q
	      // test functions that belong also to cell i are already treated
	      if(out==2)
	      {  
		  OutPut("testfct. neigh cell" << endl);
	      }
	      for (ii=0;ii<N_;ii++)                 //ii - test function of cell
	      {
		  // look for 'ii'-th row in all matrices
		  dof_ii = DOF_neigh[ii];
		  // Dirichlet node
		  if (dof_ii>=ActiveBound)
		      continue;
		  // check if test function belongs also to mesh cell i
		  dummy = 0;
		  k=0;
		  while(k<N_ && dummy == 0)
		  {
		      if(dof_ii == DOF[k])
			  dummy=1;
		      k++;
		  }
		  if (dummy) 
		      continue;
		  // initialize edge integral
		  integral=0;
		  // compute edge integral
		  for (j=0;j<N_Points1D;j++)       
		  {
		      if (out==2)
		      {
			  OutPut("(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): normal jump " << u_jump_x[j]*nx + u_jump_y[j]*ny << endl);
		      OutPut("(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): testjump " << jump_xderiv[j][N_+ii]*nx + jump_yderiv[j][N_+ii]*ny << " : " << jump_xderiv[j][N_+ii]*nx << " " <<  jump_yderiv[j][N_+ii]*ny << endl);
		      }
		      integrand = edge_par * diffusion * diffusion 
			  *(u_jump_x[j]*nx + u_jump_y[j]*ny) *
			  (jump_xderiv[j][N_+ii]*nx + jump_yderiv[j][N_ + ii]*ny);
		      // hE/2 : determinant of reference transformation
		      w = weights1D[j]*hE/2;
		      integral += w*integrand;        // integral on the edge
		  }
		  if (out==2)
		      OutPut(dof_ii <<" integral b" << "(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): "
			     << integral << " :: " << jump_xderiv[0][N_ + ii]  << " " << jump_yderiv[0][N_ + ii] << endl);
		  // update rhs
		  rhs[dof_ii] += integral;
	      }                                     // end outer loop over dof (ii)
	      for (m=0;m<4;m++)
	      {
		  delete X1D_neigh[m];
		  delete Y1D_neigh[m];
	      }                                     //endfor m
	  } // end if (i<q)
      } // end if(neigh)
  } // end loop over the edges (r)
  } // end loop over the cells (i)

  if (out==2)
      OutPut("free memory " << endl);

  delete SecondDer;

  delete u_jump_x;
  delete u_jump_y;

  delete UsedElements;
  for (i=0;i<4;i++)
  {
    delete X1D[i];
    delete Y1D[i];
    delete AbsDetjk1D[i];
    for (j=0;j<N_Points1D;j++)
    {
      delete xyval_ref1D[i][j];
      delete xderiv_ref1D[i][j];
      delete yderiv_ref1D[i][j];
    }
  }
  for(l=0;l<6;l++)
  {
    delete value_basefunct_ori[l];
    delete xderiv_basefunct_ori[l];
    delete yderiv_basefunct_ori[l];
    delete value_basefunct_ori_neigh[l];
    delete xderiv_basefunct_ori_neigh[l];
    delete yderiv_basefunct_ori_neigh[l];
  }

  for (i=0;i<N_BaseFuncts2D;i++)
  {
     for (j=0;j<6;j++)
    {
      delete value_basefunct_ref1D[i][j];
      delete xderiv_basefunct_ref1D[i][j];
      delete yderiv_basefunct_ref1D[i][j];
    }
     delete value_basefunct_ref1D[i];
     delete xderiv_basefunct_ref1D[i];
     delete yderiv_basefunct_ref1D[i];
  }


  for (i=0;i<N_Points1D;i++)
  {
    delete xyval_refNeigh1D[i];
    delete  xderiv_refNeigh1D[i];
    delete   yderiv_refNeigh1D[i];
  }

  delete aux2;
  delete aux4;

  if(n_rhs)
  {
    delete righthand;
    delete LocRhs;
    delete RhsBeginIndex;
    delete RhsGlobalNumbers;
  }
  for(i=0; i < N_BaseFuncts2D; i++)
  {
    for(j=0; j < 4; j++)
      {
	for(m=0; m < MaxN_QuadPoints_1D; m++)
          {
	    //for(n=0; n < MaxN_BaseFunctions2D_Ersatz; n++)
	    //{
	    delete xietaval_ref1D[i][j][m];
	    delete xideriv_ref1D[i][j][m];
	    delete etaderiv_ref1D[i][j][m];
	    // }
	  }
	delete xietaval_ref1D[i][j];
	delete xideriv_ref1D[i][j];
	delete etaderiv_ref1D[i][j];
      }
    delete xietaval_ref1D[i];
    delete xideriv_ref1D[i];
    delete  etaderiv_ref1D[i];
  }

  delete value_basefunct_ref1D;
  delete xderiv_basefunct_ref1D;
  delete yderiv_basefunct_ref1D;

  delete xietaval_ref1D;
  delete xideriv_ref1D;
  delete etaderiv_ref1D;

  int N_Rows;
  // ####################################################################
  // print jump term
  // ####################################################################
  if(out==2)
  {
    N_Rows = fespace->GetN_DegreesOfFreedom();
    for(i=0;i<N_Rows;i++)
      OutPut(setw(5) << i << setw(20) << rhs[i] << endl);
  }

}                                                 // end of JumpTermsForAdjointProblem

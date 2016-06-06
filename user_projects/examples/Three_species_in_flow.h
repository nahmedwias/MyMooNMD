/*****************************************************************************
 *  @name Three_species_in_flow.h
 *  @brief Example file for a time dependent system of three cdr equations
 *  coupled in the reaction part.
 *
 *  To be used as include example in class Example_CoupledCDR2D.
 *
 *  Use the example with unit square geometry and initial mesh
 *  (e.g. by using the default domain initialisation:
 *   domain.Init("Default_UnitSquare", "UnitSquare"); ).
 *
 *  This example was created by adapting the MooNMD example
 *  "/MooNMD/Examples/TNSE_2D/Bulk_Fallung_Driven_Cavity.h"
 *  to the requirements of the new coupled solver module Coupled_Time_CDR_2D.
 *  The original example was used for the publication
 *  John & Roland 2010: On the impact of the scheme for solving the higher
 *  dimensional equation in coupled population balance systems. Int. J. Numer.
 *  Meth. Engng 82: 1450--1474.
 *
 *  Especially: TODO
 *    - remove particle terms and equations
 *    - remove as much database dependency as is possible
 *    - hard-code the assembling functions trinity
 *    - make flow a stationary, precomputed coefficient function
 *
 *  Note that in addition to the usual functions, which ParMooN examples have to
 *  supply, this also provides a parameter- and assembling functions.
 *  It is due to the MooNMD heritage, that these get passed around the program
 *  as function pointers and thus have to be hardcoded somewhere outside of
 *  class scope. Putting them into the example file seemed the best of the
 *  bad alternatives.
 *
 *  @date Jan 6 2016
 *  @author Clemens Bartsch
 *****************************************************************************/

void ExampleFile()
{
  int range;

  OutPut(" Example: Bulk_Fallung_Driven_Cavity.h ") ;
  OutPut(" inflow (u_infty)" << TDatabase::ParamDB->BULK_u_infty);
  OutPut(" upper lid " << TDatabase::ParamDB->P5);
  OutPut(" left " << (int)TDatabase::ParamDB->P7);
  OutPut(" right " << (int)TDatabase::ParamDB->P8);
  OutPut(" lower " << (int)TDatabase::ParamDB->P9 << endl);

  range = (int)TDatabase::ParamDB->P7;
  if ((range<1)||(range>30))
  {
      OutPut("left boundary out of range !!!"<< endl);
      exit(4711);
  }

  range = (int)TDatabase::ParamDB->P8;
  if ((range<1)||(range>30))
  {
      OutPut("right boundary out of range !!!"<< endl);
      exit(4711);
  }

  range = (int)TDatabase::ParamDB->P9;
  if ((range<1)||(range>30))
  {
      OutPut("lower boundary out of range !!!"<< endl);
      exit(4711);
  }
  // set some parameters
  TDatabase::ParamDB->BULK_D_P_MIN = TDatabase::ParamDB->BULK_D_P_0/TDatabase::ParamDB->BULK_D_P_MAX;
  TDatabase::ParamDB->BULK_c_C_infty = TDatabase::ParamDB->BULK_c_C_infty_sat
      * exp(TDatabase::ParamDB->BULK_C_2/TDatabase::ParamDB->BULK_D_P_0);
  TDatabase::ParamDB->BULK_f_infty = TDatabase::ParamDB->BULK_u_infty/(TDatabase::ParamDB->BULK_C_g
      *TDatabase::ParamDB->BULK_k_g*pow(TDatabase::ParamDB->BULK_D_P_MAX,3)*TDatabase::ParamDB->BULK_l_infty);
  OutPut("BULK d_p_min " << TDatabase::ParamDB->BULK_D_P_MIN  <<
	 " u_infty " << TDatabase::ParamDB->BULK_u_infty <<
	 " c_C_infty " << TDatabase::ParamDB->BULK_c_C_infty <<
	 " f_infty " << TDatabase::ParamDB->BULK_f_infty << endl);
  OutPut("PB_DISC_TYPE = " <<   TDatabase::ParamDB->PB_DISC_TYPE << endl);
  OutPut("PB_TIME_DISC = " <<   TDatabase::ParamDB->PB_TIME_DISC << endl);

  //TDatabase::ParamDB->SAVE_DATA = TRUE;
}
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  int lower = (int)TDatabase::ParamDB->P9;

  if ((i==0)&&((t>lower/32.0)&&(t<(lower+2)/32.0)))
      cond = NEUMANN;
  else
      cond = DIRICHLET;
  // cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0; 
}

void U1BoundValue(int BdComp, double Param, double &value)
{
    int range;
    double y, fac = 1024; // 32^2
    

  switch(BdComp)
  {
     case 0: 
        value = 0;
        break;
     case 1:
	range = (int)TDatabase::ParamDB->P8;
	if ((Param>range/32.0)&&(Param<(range+1)/32.0))
        {
           y = Param;
           value =  6*(y-range/32.0)*(y-(range+1)/32.0)*fac;
	}
	else
	    value = 0;
	break;
	// upper boundary
    case 2: if(Param<0.00001 || Param>0.99999) 
              value = 0;
            else
               value = TDatabase::ParamDB->P5/TDatabase::ParamDB->BULK_u_infty;
               //value = 0;
            break;
    case 3: 
	range = (int)TDatabase::ParamDB->P7;
        y = 1-Param;
	if ((y>range/32.0)&&(y<(range+1)/32.0))
        {
           value = -6*(y-range/32.0)*(y-(range+1)/32.0)*fac;
        }
	else
	    value = 0;
	break;
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
/*  int lower = (int)TDatabase::ParamDB->P9;

  if ((BdComp==0)&&((Param>lower/32.0)&&(Param<(lower+1)/32.0)))
      value = -2;
  else
      value = 0;
*/
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff;
  double l_infty = TDatabase::ParamDB->BULK_l_infty;
  double u_infty = TDatabase::ParamDB->BULK_u_infty;
  double density = TDatabase::ParamDB->BULK_density;
  double dynvisc = TDatabase::ParamDB->BULK_dynamic_viscosity;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = dynvisc/(density*l_infty*u_infty);
    //coeff[0] = 1;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}


// ========================================================================
// definitions for the substance A
// ========================================================================

// initial conditon
void InitialCondition_c_A(double x, double y, double *values)
{
   int range;
   double eps=1e-8;
   double t = TDatabase::TimeDB->CURRENTTIME;
   double t0 = TDatabase::TimeDB->T0;

   if (t<t0)
     values[0] = 0;
   else
   {		
     range = (int)TDatabase::ParamDB->P7;
     if ((fabs(x)<1e-7)&&(y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
       values[0] = 1;
     else
       values[0] = 0;
   }
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_A(int BdComp, double Param, BoundCond &cond)
{
   int range;
   double y,eps=1e-8;

   if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) && 
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
    cond = NEUMANN;
   else
   {      
     if (BdComp==3)
     {
      range = (int)TDatabase::ParamDB->P7;
      y = 1-Param;
      if ((y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
      {
         cond = DIRICHLET;
      }
      else
	  cond = NEUMANN;
   }
   else
     cond = NEUMANN;
  }
}

// value of boundary condition
void BoundValue_c_A(int BdComp, double Param, double &value)
{
   int range;
   double y,eps=1e-8;
   double t = TDatabase::TimeDB->CURRENTTIME;
   double t0 = TDatabase::TimeDB->T0;
  
   if (t<t0)
     value = 0;
   else
   {
   if (BdComp==3)
   {
      range = (int)TDatabase::ParamDB->P7;
      y = 1-Param;
      if ((y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
      {
	  value = 1;
      }
      else
	  value = 0;
   }
   else
     value = 0;
   }
}

void NoCoeffs(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
    return;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double L_infty = TDatabase::ParamDB->BULK_l_infty;
  double U_infty = TDatabase::ParamDB->BULK_u_infty;
  double C_infty = TDatabase::ParamDB->BULK_c_infty;
  double D_A = TDatabase::ParamDB->BULK_D_A;
  double k_r = TDatabase::ParamDB->BULK_k_r;
  double T_infty, eps, c;
  
  T_infty = L_infty/U_infty;
  eps = D_A/(L_infty*U_infty);
  c = k_r*C_infty * L_infty /U_infty;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = param[1]; // u1
    coeff[2] = param[2]; // u2
    coeff[3] = c * param[0];
    coeff[4] = 0; 
    //OutPut(param[0] << " " << param[1] << " " << param[2] << endl); 
  }
}

// ========================================================================
// definitions for the substance B
// ========================================================================

// initial conditon
void InitialCondition_c_B(double x, double y, double *values)
{
   int range;
   double eps = 1e-8;
   double t = TDatabase::TimeDB->CURRENTTIME;
   double t0 = TDatabase::TimeDB->T0;

   if (t<t0)
     values[0] = 0;
    else
    {
    range = (int)TDatabase::ParamDB->P8;
    if ((fabs(1-x) <1e-7)&&(y>=range/32.0-eps)&&(y<=(range+1)/32.0+eps))
	values[0] = 1;
    else
	values[0] = 0;
}	
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_B(int BdComp, double Param, BoundCond &cond)
{
   int range;
   double eps = 1e-8;

  if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) && 
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
    cond = NEUMANN;
   else
   {
     if (BdComp==1)
     {
       range = (int)TDatabase::ParamDB->P8;
       if ((Param>=range/32.0-eps)&&(Param<=(range+1)/32.0+eps))
       {
         cond = DIRICHLET;
       }
       else
         cond = NEUMANN;
     }
     else
       cond = NEUMANN;
   }
}

// value of boundary condition
void BoundValue_c_B(int BdComp, double Param, double &value)
{
   int range;
   double t = TDatabase::TimeDB->CURRENTTIME;
   double t0 = TDatabase::TimeDB->T0;
  
   if (t<t0)
     value = 0;
   else
   {   
   if (BdComp==1)
   {
      range = (int)TDatabase::ParamDB->P8;
      if ((Param>=range/32.0)&&(Param<=(range+1)/32.0))
      {
	  value = 1;
      }
      else
	  value = 0;
   }
   else
     value = 0;
   }
}

// ========================================================================
// definitions for the substance C
// ========================================================================

// initial conditon
void InitialCondition_c_C(double x, double y, double *values)
{
  values[0] = 0;	
}

// kind of boundary condition (for FE space needed)
void BoundCondition_c_C(int BdComp, double Param, BoundCond &cond)
{
    int range;
    double y;

   cond = NEUMANN;

/*
   if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) && 
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
     cond = NEUMANN;
   else
   {
    cond = NEUMANN;

    switch(BdComp)
    {
      case 1:
        range = (int)TDatabase::ParamDB->P8;
        if ((Param>range/32.0)&&(Param<(range+1)/32.0))
        {
          cond = DIRICHLET;
        }
      break;
      case 3:
        range = (int)TDatabase::ParamDB->P7;
        y = 1-Param;
        if ((y>range/32.0)&&(y<(range+1)/32.0))
        {
          cond = DIRICHLET;
        }
      break;
    }
  }
*/
}

// value of boundary condition
void BoundValue_c_C(int BdComp, double Param, double &value)
{
   value = 0;
}

// param[2] : c_A
// param[3] : c_B
// param[4] : c_C (old)
// param[5] : r_g
void BilinearCoeffs_Cc(int n_points, double *X, double *Y,
		    double **parameters, double **coeffs)
{
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double L_infty = TDatabase::ParamDB->BULK_l_infty;
  double U_infty = TDatabase::ParamDB->BULK_u_infty;
  double C_infty = TDatabase::ParamDB->BULK_c_infty;
  double c_C_infty_sat =  TDatabase::ParamDB->BULK_c_C_infty_sat;
  double C_g = TDatabase::ParamDB->BULK_C_g;
  double C_nuc = 15.33;
  double C_2 = TDatabase::ParamDB->BULK_C_2;
  double D_A = TDatabase::ParamDB->BULK_D_A;
  double d_p_0 = TDatabase::ParamDB->BULK_D_P_0;
  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
  double k_g = TDatabase::ParamDB->BULK_k_g;
  double k_r = TDatabase::ParamDB->BULK_k_r;
  double k_nuc =  TDatabase::ParamDB->BULK_k_nuc;
  double eps,  B_C_c, T_infty, lambda_chem, lambda_nuc;
  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
  double c_C_infty = TDatabase::ParamDB->BULK_c_C_infty;
  double r_g;

  // compute derived quantities of the model
  T_infty = L_infty/U_infty;

  // compute coefficients of the equation
  eps = D_A/(L_infty*U_infty);
  lambda_chem = k_r*C_infty*C_infty*L_infty /(U_infty*c_C_infty);
  lambda_nuc = C_nuc*k_nuc*pow(d_p_0,3)*L_infty*pow(c_C_infty,4)/U_infty;
  //OutPut("lambda_chem " << lambda_chem << "  lambda_nuc " << lambda_nuc << " " << endl);
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    if (param[4] < 1)
      B_C_c = 0;
    else
      B_C_c = pow(param[4] - 1,5);

    coeff[0] = eps;
    coeff[1] = param[0]; // u1
    coeff[2] = param[1]; // u2
    coeff[3] = 0;
    // following script from 05/08/23 for r_g
    r_g = param[4] - c_C_infty_sat/c_C_infty;
    // with d_p -- dependent definition of G(c_C,d_p)
    if (TDatabase::ParamDB->BULK_GROWTH_RATE == 2)
       r_g = 1;
    // with d_p - dependent definition of G(c_C,d_p)
     coeff[4] = lambda_chem*param[2]*param[3] - lambda_nuc*B_C_c 
	- r_g*param[5];//r_chem - r_nuc - r_g    
  }
}

// CB 2016/02/03: We're not using MOM at the moment, so all the functions
// necessary for its assembling are commented out here.

// /******************************************************************************/
// // MOM
// /******************************************************************************/
//
//void BoundCondition_mom(int BdComp, double Param, BoundCond &cond)
//{
//    int range;
//    double y;
//
//    cond = NEUMANN;
//
//    switch(BdComp)
//    {
//	case 1:
//	    range = (int)TDatabase::ParamDB->P8;
//	    if ((Param>range/32.0)&&(Param<(range+1)/32.0))
//	    {
//		cond = DIRICHLET;
//	    }
//	    break;
//	case 3:
//	    range = (int)TDatabase::ParamDB->P7;
//	    y = 1-Param;
//	    if ((y>range/32.0)&&(y<(range+1)/32.0))
//	    {
//		cond = DIRICHLET;
//	    }
//	    break;
//    }
//}
//
//// value of boundary condition
//void BoundValue_mom(int BdComp, double Param, double &value)
//{
//   value = 0;
//}
//
//// initial conditon
//void InitialCondition_mom(double x, double y, double *values)
//{
//  values[0] = 0;
//}
//void BilinearCoeffs_mom(int n_points, double *X, double *Y,
//		    double **parameters, double **coeffs)
//{
//  int i, k = TDatabase::ParamDB->INTERNAL_MOMENT;
//  double *coeff, *param;
//  double x, y;
//  double t = TDatabase::TimeDB->CURRENTTIME;
//  double l_infty = TDatabase::ParamDB->BULK_l_infty;
//  double u_infty = TDatabase::ParamDB->BULK_u_infty;
//  double c_C_infty = TDatabase::ParamDB->BULK_c_infty;
//  double c_C_infty_sat = TDatabase::ParamDB->BULK_c_C_infty_sat;
//  double f_infty = TDatabase::ParamDB->BULK_f_infty;
//  double d_p_min = TDatabase::ParamDB->BULK_D_P_MIN;
//  double d_p_max = TDatabase::ParamDB->BULK_D_P_MAX;
//  double k_g = TDatabase::ParamDB->BULK_k_g;
//  double k_nuc = TDatabase::ParamDB->BULK_k_nuc;
//  double eps, factor_G, factor_G0, G_c_C, c_d_p, B_c_C, f_d_p_min;
//
//  // this is just for testing
//  eps = 1e-10;
//  // computed model constants
//  factor_G = k_g*c_C_infty*l_infty/(u_infty*d_p_max);
//  factor_G0 = k_g*c_C_infty*f_infty;
//
//  for(i=0;i<n_points;i++)
//  {
//    coeff = coeffs[i];
//    param = parameters[i];
//
//    if (TDatabase::ParamDB->BULK_GROWTH_RATE==2)
//    {
//	OutPut("MOM not for BULK_GROWTH_RATE==2 implemented !!!" << endl);
//	exit(4711);
//    }
//    else
//    {
//       G_c_C = factor_G0*(param[2] - c_C_infty_sat/c_C_infty);
//       c_d_p = factor_G*(param[2] - c_C_infty_sat/c_C_infty);
//    }
//
//    if (G_c_C > 1e-10)
//    {
//        // compute rate of nucleation
//        B_c_C = k_nuc*pow(c_C_infty*(param[2] - 1),5);
//        // truncate negative values
//        if (B_c_C < 0)
//	    B_c_C = 0;
//	f_d_p_min = B_c_C / G_c_C;
//    }
//    else
//	f_d_p_min = 0;
//
//    coeff[0] = eps;
//    coeff[1] = param[0]; // u1
//    coeff[2] = param[1]; // u2
//    coeff[3] = 0;
//    if (k>0)
//	coeff[4] = c_d_p * f_d_p_min * pow(d_p_min,k) +  c_d_p * k * param[3];
//    else
//	coeff[4] = c_d_p * f_d_p_min;
//  }
//}
//
////#endif

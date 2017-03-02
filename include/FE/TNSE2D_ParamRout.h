// ======================================================================
// @(#)TNSE2D_ParamRout.h        1.2 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#ifndef __TNSE2D_PARAMROUT__
#define __TNSE2D_PARAMROUT__

#include <Constants.h>
#include <Enumerations.h>

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void TimeNSParams2(double *in, double *out);

// ========================================================================
// parameters: u1old, u2old, u1_previous time, u2_previous time
// ========================================================================
void TimeNSParamsVelo_GradVelo(double *in, double *out);


// ========================================================================
// coletti, without g_\delta \ast u, in ALE 
// ========================================================================
void TimeNSParamsVelo_GradVelo_ALE(double *in, double *out);


// ========================================================================
// velocity, gradient and convolution of velocity
// ========================================================================

void TimeNSParamsVelo_GradVelo_ConvVelo(double *in, double *out);

// ========================================================================
// Coletti,
// turbulent viscosity \|u - g_\delta \ast u\|_2
// ========================================================================
void TimeNSParamsVelo_GradVeloNuT4(double *in, double *out);
//======================================================================
// Galdi/Layton with convolution, without g_\delta \ast u 
// ========================================================================
void TimeNSParamsGL00Convolution(double *in, double *out);
// ========================================================================
// Galdi/Layton with convolution, rhs assembling, without g_\delta \ast u
// ========================================================================
void TimeNSParamsRHSLES(double *in, double *out);
// ========================================================================
// Galdi/Layton with convolution, rhs assembling, 
// turbulent viscosity \|u - g_\delta \ast u\|_2
// ========================================================================
void TimeNSParamsRHSGL00ConvolutionNuT4(double *in, double *out);
// ========================================================================
// parameters for Galdi/Layton model with auxiliary problem
// ========================================================================

// ========================================================================
// parameters for Galdi/Layton model with auxiliary problem
// turbulent viscosity \|u - g_\delta \ast u\|_2
// ========================================================================

void TimeNSParamsGL00AuxProblemNuT4(double *in, double *out);

void TimeNSParamsGL00AuxProblemPaper2(double *in, double *out);

// ========================================================================
// parameters: gradient(u1), gradient(u2)
// ========================================================================
void TimeNSParamsGrad(double *in, double *out);

// ========================================================================
// parameters for VMS
// ========================================================================
void TimeNSParams_VMS_SmallRhs2D(double *in, double *out);

void TimeNSParams_VMS_LargeRhs2D(double *in, double *out);

// ========================================================================
// parameters: low order : u1old, u2old, higher order  : u1old, u2old
// ========================================================================
void TimeNSParams_NLGalerkin_VMS_1_2D(double *in, double *out);
// ========================================================================
// right-hand side ONLY, defect correction type 1, u2
// ========================================================================

void TimeNSParamsVelo_GradVeloOld2(double *in, double *out);
// ========================================================================
// parameters: x, y, u1old, u2old
// ========================================================================

void TimeNSParamsVeloPos(double *in, double *out);
// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void MovingTNSParams(double *in, double *out);

// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
//  for axial symmetric case
// ======================================================================
void MovingTNSParams_Axial3D(double *in, double *out);

// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
//  for axial symmetric case
// ======================================================================
void MovingTNSParams_Axial3D_HeatLine(double *in, double *out);

// ======================================================================
// parameters for heatline   
// ======================================================================
void  ParamsFct_HeatLine(double *in, double *out);

// ========================================================================
// parameters: u1old, u2old, u1 previous, u2 previous
// used for : SUPG
// ========================================================================
void TimeNSParams4(double *in, double *out);

#endif






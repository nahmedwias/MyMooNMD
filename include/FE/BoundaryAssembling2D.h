/** ************************************************************************ 
*
* @brief     Common declaration for boundary integrals
*     
* 
* @author    Alfonso Caiazzo & Laura Blank
* @date      18.05.16
 ************************************************************************  */

#ifndef __BoundaryAssembling2D__
#define __BoundaryAssembling2D__

#include <Enumerations.h>
#include <FESpace2D.h>
#include <FEFunction2D.h>

/**
   @brief compute (pressure(given),v.n)_{[compBC]}
   @param[in] compBC: the boundary component to integrate on
   @param[in] mult: given multiplicative factor
   @param[in] pFunct: the boundary pressure (as a finite element function)
 */
void  BoundaryAssemble_on_rhs_g_v_n(double **rhs, 
			const TFESpace2D *U_Space, 
			TFEFunction2D *pFunct,
			int compBC,
			double mult
			 );
  
// ======================================================================
// Type 1, Standard Galerkin for Brinkman in [p div u] formulation
// ======================================================================
//void BrinkmanType1Galerkin(double Mult, double *coeff,
//                   double *param, double hK,
//                   double **OrigValues, int *N_BaseFuncts,
//                   double ***LocMatrices, double **LocRhs);
// ======================================================================

#endif

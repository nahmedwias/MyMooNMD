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

class BoundaryAssembling2D
{
  public:
    void BoundaryAssemble_on_rhs_g_v_n(double **rhs,
                                       const TFESpace2D *U_Space,
                                       TFEFunction2D *pFunct,
                                       int compBC,
                                       double mult
                                       ) ;


    //    void GetQuadFormulaData(int FEId, int &nQuadPoints, QuadFormula1D &LineQuadFormula,double &quadPoints, double &quadWeights)
//    {
//        // get a quadrature formula good enough for the velocity FE space
//        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//        // get the type of required quadrature (include/FE/Enumerations.h)
//        QuadFormula1D LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
//        // initialize points and weights of quadrature
//        TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
//        
//        ///@todo rewrite the GetFormulaData using (e.g.) vector<> class
//        qf1->GetFormulaData(int &nQuadPoints, double &quadWeights, double &quadPoints);
//        TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
//    }
    
    
 ///   void BundaryAssemble(std::vector<TJoint*> joints, rhs, );
  protected:
    void GetQuadFormulaData(int degree,
			    std::vector<double> &P,
			    std::vector<double> &W);

    std::vector<TJoint*> get_joints_of_component(int boundary_component);

    ///@brief type of quadrature used for line-integration
    QuadFormula1D LineQuadFormula;


};
/**
   @brief compute (pressure(given),v.n)_{[compBC]}
   @param[in] compBC: the boundary component to integrate on
   @param[in] mult: given multiplicative factor
   @param[in] pFunct: the boundary pressure (as a finite element function)
 */

// ======================================================================
// Type 1, Standard Galerkin for Brinkman in [p div u] formulation
// ======================================================================
//void BrinkmanType1Galerkin(double Mult, double *coeff,
//                   double *param, double hK,
//                   double **OrigValues, int *N_BaseFuncts,
//                   double ***LocMatrices, double **LocRhs);
// ======================================================================

#endif

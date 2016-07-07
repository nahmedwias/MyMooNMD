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
#include <BoundEdge.h>
#include <Collection.h>
#include <Database.h>
#include <BlockFEMatrix.h>

class BoundaryAssembling2D
{
public:
    void rhs_g_v_n(double **rhs,
		   const TFESpace2D *U_Space,
		   TFEFunction2D *pFunct,
		   int compBC,
		   double mult
		   ) ;
    void rhs_g_v_n(double **rhs,
		   const TFESpace2D *U_Space,
		   TFEFunction2D *given_data,
		   std::vector<TBoundEdge*> &edge,
		   double mult
		   );
    
    void rhs_g_v(double **rhs,
                   const TFESpace2D *U_Space,
                   TFEFunction2D *pFunct,
                   int compBC,
                   double mult
                   ) ;
    void rhs_g_v(double **rhs,
                   const TFESpace2D *U_Space,
                   TFEFunction2D *given_boundary_data,
                   std::vector<TBoundEdge*> &edge,
                   double mult
                   );
    
    void matrix_v_n_v_n(BlockFEMatrix &M,
			const TFESpace2D *U_Space,
			int boundary_component_id,
			double mult
			);
    
    void matrix_v_n_v_n(BlockFEMatrix &M,
			const TFESpace2D *U_Space,
			std::vector<TBoundEdge*> &edge,
			double mult);
    
    void matrix_gradv_n_v(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          int boundary_component_id,
                          double mult
                          );
    
    void matrix_gradv_n_v(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          std::vector<TBoundEdge*> &edge,
                          double mult);
    
    void matrix_u_v(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          int boundary_component_id,
                          double mult
                          );
    
    void matrix_u_v(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          std::vector<TBoundEdge*> &edge,
                          double mult);
    
    
    void matrix_p_v_n(BlockFEMatrix &M,
                    const TFESpace2D *U_Space,
                    const TFESpace2D *P_Space,
                    int boundary_component_id,
                    double mult
                    );
    
    void matrix_p_v_n(BlockFEMatrix &M,
                    const TFESpace2D *U_Space,
                    const TFESpace2D *P_Space,
                    std::vector<TBoundEdge*> &edge,
                    double mult);
    
    //todo for symmetry (nitsche): ((u-g)n,q) and -((u-g),dnvn)

    
protected:
    void get_quadrature_formula_data(int degree,
                                     std::vector<double> &P,
                                     std::vector<double> &W);
    
    void get_original_values(FE2D FEId, int joint_id, TBaseCell *cell,
			     std::vector<double> quadPoints, int BaseVectDim,
			     std::vector< std::vector<double> > &v00,
			     std::vector< std::vector<double> > &v10,
                             std::vector< std::vector<double> > &v01);
    
    
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

#endif

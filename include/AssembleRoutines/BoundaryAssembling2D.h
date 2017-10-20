/** ************************************************************************
 *
 * @brief     Assemble boundary integrals to the matrix and add their values to the rhs
 *
 *
 * @author    Alfonso Caiazzo & Laura Blank
 * @date      18.05.16
 ************************************************************************  */

#ifndef __BoundaryAssembling2D__
#define __BoundaryAssembling2D__

#include <Enumerations.h>
#include <FESpace2D.h>
#include <BoundEdge.h>
#include <Database.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>

class BoundaryAssembling2D
{
public:
    
    /** @brief assemble integral (given 1D function, v \cdot n)_{L^2([boundary_component_id])}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     @param[in] given_boundary_data: the boundary pressure (as a finite element function)
     */
    void rhs_g_v_n(BlockVector &rhs,
                   const TFESpace2D *U_Space,
                   BoundValueFunct2D *given_boundary_data,
                   int boundary_component_id,
                   double mult
                   ) ;
    void rhs_g_v_n(BlockVector &rhs,
                   const TFESpace2D *U_Space,
                   BoundValueFunct2D *given_boundary_data,
                   std::vector<TBoundEdge*> &edge,
                   double mult
                   );
    
    /** @brief assemble integral (u \cdot n, v \cdot n)_{L^2([boundary_component_id])}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void matrix_u_n_v_n(BlockFEMatrix &M,
                        const TFESpace2D *U_Space,
                        int boundary_component_id,
                        double mult
                        );
    
    void matrix_u_n_v_n(BlockFEMatrix &M,
                        const TFESpace2D *U_Space,
                        std::vector<TBoundEdge*> &edge,
                        double mult);
    
    /** @brief assemble integral (u, v)_{L^2([boundary_component_id])}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     @param[in] rescale_by_h: true: divide by length of the edges
     false: do not divide by length of the edges
     */
    void matrix_u_v(BlockFEMatrix &M,
                    const TFESpace2D *U_Space,
                    int boundary_component_id,
                    double mult,
                    bool rescale_by_h
                    );
    
    void matrix_u_v(BlockFEMatrix &M,
                    const TFESpace2D *U_Space,
                    std::vector<TBoundEdge*> &edge,
                    double mult,
                    bool rescale_by_h);
    
    
    /** @brief assemble integral (given 2D function, v)_{L^2([boundary_component_id])}
     @param[in] [given_boundary_data1,given_boundary_data2]: e.g., the 2D Dirichlet boundary velocity (as a finite element function)
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     @param[in] rescale_by_h: true: divide by length of the edges
     false: do not divide by length of the edges
     */
    
    void rhs_uD_v(BlockVector &rhs,
                 const TFESpace2D *U_Space,
                 BoundValueFunct2D *given_boundary_data1,
                 BoundValueFunct2D *given_boundary_data2,
                 int boundary_component_id,
                 double mult,
                 bool rescale_by_h
                 ) ;
    void rhs_uD_v(BlockVector &rhs,
                 const TFESpace2D *U_Space,
                 BoundValueFunct2D *given_boundary_data1,
                 BoundValueFunct2D *given_boundary_data2,
                 std::vector<TBoundEdge*> &edge,
                 double mult,
                 bool rescale_by_h
                 );

    
    /** @brief assemble integral (\nabla u \cdot n, v)_{L^2([boundary_component_id])}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    
    void matrix_gradu_n_v(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          int boundary_component_id,
                          double mult
                          );
    
    void matrix_gradu_n_v(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          std::vector<TBoundEdge*> &edge,
                          double mult,
                          int boundary_component_id);
    
    
    /** @brief assemble integral (\nabla v \cdot n, u)_{L^2([boundary_component_id])}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void matrix_gradv_n_u(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          int boundary_component_id,
                          double mult
                          );
    
    void matrix_gradv_n_u(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          std::vector<TBoundEdge*> &edge,
                          double mult);

    
    /** @brief assemble integral (\nabla v \cdot n, given 2D function)_{L^2([boundary_component_id])}
     @param[in] [given_boundary_data1, given_boundary_data2]: e.g., the 2D Dirichlet boundary velocity (as a finite element function)
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void rhs_gradv_n_uD(BlockVector &rhs,
                        const TFESpace2D *U_Space,
                        BoundValueFunct2D *given_boundary_data1,
                        BoundValueFunct2D *given_boundary_data2,
                        int boundary_component_id,
                        double mult
                        ) ;
    void rhs_gradv_n_uD(BlockVector &rhs,
                        const TFESpace2D *U_Space,
                        BoundValueFunct2D *given_boundary_data1,
                        BoundValueFunct2D *given_boundary_data2,
                        std::vector<TBoundEdge*> &edge,
                        double mult
                        );
    
    /** @brief assemble integral (p, v \cdot n)_{L^2([boundary_component_id])} into the matrix
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
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

    
    /** @brief assemble integral (q, u \cdot n)_{L^2([boundary_component_id])} into the matrix
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void matrix_q_u_n(BlockFEMatrix &M,
                      const TFESpace2D *U_Space,
                      const TFESpace2D *P_Space,
                      int boundary_component_id,
                      double mult
                      );
    
    void matrix_q_u_n(BlockFEMatrix &M,
                      const TFESpace2D *U_Space,
                      const TFESpace2D *P_Space,
                      std::vector<TBoundEdge*> &edge,
                      double mult);
    
    
    /** @brief assemble integral (q, given 2D function \cdot n )_{L^2([boundary_component_id])}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     @param[in] given_boundary_data: the boundary velocity (as a finite element function)
     */
    void rhs_q_uD_n(BlockVector &rhs,
                   const TFESpace2D *U_Space,
                   const TFESpace2D *P_Space,
                   BoundValueFunct2D *given_boundary_data1,
                   BoundValueFunct2D *given_boundary_data2,
                   int boundary_component_id,
                   double mult
                   ) ;
    void rhs_q_uD_n(BlockVector &rhs,
                   const TFESpace2D *U_Space,
                   const TFESpace2D *P_Space,
                   BoundValueFunct2D *given_boundary_data1,
                   BoundValueFunct2D *given_boundary_data2,
                   std::vector<TBoundEdge*> &edge,
                   double mult
                   );

    

    
    
protected:
    
    /** @brief type of quadrature formula used for line-integration
     */
    QuadFormula1D LineQuadFormula;
    
    /** @brief Get the quadrature points and weights according to LineQuadFormula = Quadrature formula for specified degree. LineQuadFormula has to be set before the call of this function
     */
    void get_quadrature_formula_data(std::vector<double> &P,
                                     std::vector<double> &W);
    
    /** @brief access to the coordinates and first order partial derivatives of the solution (e.g. $\boldsymbol{u}$ or $p$) on joint with joint_id in reference element and transformation to the actual element with output of the coordinates [u00] and partial first order derivatives [u10],[u01] of the actual solution (e.g. u or p)
     */
    void get_original_values(FE2D FEId, int joint_id, TBaseCell *cell,
                             std::vector<double> quadPoints, int BaseVectDim,
                             std::vector< std::vector<double> > &u00,
                             std::vector< std::vector<double> > &u10,
                             std::vector< std::vector<double> > &u01);
};

#endif


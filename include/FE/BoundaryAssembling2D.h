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
    
    /** @brief integral (pressue(given 1D function), v)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor
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
    
    /** @brief integral (given 2D function, v)_{[boundary_component_id]}
     @param[in] [given_boundary_data1,given_boundary_data2]: the e.g. 2D Dirichlet boundary velocity (as a finite element function)
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor
     @param[in] rescale_by_h: true: divide by length of the edges
     false: do not divide by length of the edges
     */
    void rhs_g_v(BlockVector &rhs,
                 const TFESpace2D *U_Space,
                 BoundValueFunct2D *given_boundary_data1,
                 BoundValueFunct2D *given_boundary_data2,
                 int boundary_component_id,
                 double mult,
                 bool rescale_by_h
                 ) ;
    void rhs_g_v(BlockVector &rhs,
                 const TFESpace2D *U_Space,
                 BoundValueFunct2D *given_boundary_data1,
                 BoundValueFunct2D *given_boundary_data2,
                 std::vector<TBoundEdge*> &edge,
                 double mult,
                 bool rescale_by_h
                 );
    
    /** @brief integral (u \cdot n, v \cdot n)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor
     */
    void matrix_v_n_v_n(BlockFEMatrix &M,
                        const TFESpace2D *U_Space,
                        int boundary_component_id,
                        double mult
                        );
    
    void matrix_v_n_v_n(BlockFEMatrix &M,
                        const TFESpace2D *U_Space,
                        std::vector<TBoundEdge*> &edge,
                        double mult);
    
    /** @brief integral (\nabla u \cdot n, v)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor
     */
    void matrix_gradv_n_v(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          int boundary_component_id,
                          double mult
                          );
    
    void matrix_gradv_n_v(BlockFEMatrix &M,
                          const TFESpace2D *U_Space,
                          std::vector<TBoundEdge*> &edge,
                          double mult);

    /** @brief integral ((\nabla u \cdot n).t, v.t)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor
     */
    void matrix_gradv_n_v_tangential(BlockFEMatrix &M,
				     const TFESpace2D *U_Space,
				     int boundary_component_id,
				     double mult
				     );
    
    void matrix_gradv_n_v_tangential(BlockFEMatrix &M,
				     const TFESpace2D *U_Space,
				     std::vector<TBoundEdge*> &edge,
				     double mult);
    
    /** @brief integral (u, v)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor
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

    /** @brief integral |u.n|_(u, v)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor
     @param[in] rescale_by_h: true: divide by length of the edges
     false: do not divide by length of the edges
     */
    void matrix_u_v_backflow_stab(BlockFEMatrix &M,
				  const TFESpace2D *U_Space,
				  std::vector< TFEFunction2D* >& u_conv,
				  int boundary_component_id,
				  double mult);
    
    void matrix_u_v_backflow_stab(BlockFEMatrix &M,
				  const TFESpace2D *U_Space,
				  std::vector< TFEFunction2D* >& u_conv,
				  std::vector<TBoundEdge*> &edge,
				  double mult);

     /** @brief integral |u.n|_(dt(u), dt(v))_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor
     @param[in] rescale_by_h: true: divide by length of the edges
     false: do not divide by length of the edges
     */
    void matrix_dtu_dtv_backflow_stab(BlockFEMatrix &M,
				  const TFESpace2D *U_Space,
				  std::vector< TFEFunction2D* >& u_conv,
				  int boundary_component_id,
				  double mult);
    
    void matrix_dtu_dtv_backflow_stab(BlockFEMatrix &M,
				  const TFESpace2D *U_Space,
				  std::vector< TFEFunction2D* >& u_conv,
				  std::vector<TBoundEdge*> &edge,
				  double mult);
    
    /** @brief assemble integral (p, v \cdot n)_{[boundary_component_id]} into the matrix
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor
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
    
    //todo: for symmetry (nitsche): ((u-g)n,q) and -((u-g),dnvn)
    
protected:
    
    /** @brief type of quadrature used for line-integration
     */
    QuadFormula1D LineQuadFormula;
    
    /** @brief Get the quadrature points and weights according to LineQuadFormula = Quadrature formula for specified degree), LineQuadFormula has to be set before the call of this function
     */
    void get_quadrature_formula_data(std::vector<double> &P,
                                     std::vector<double> &W);
    
    /** @brief access to the coordinates and first order partial derivatives of the solution (e.g. u or p) on joint with joint_id in reference element and transformation to the actual element with output of the coordinates [u00] and partial first order derivatives [u10],[u01] of the actual solution (e.g. u or p)
     */
    void get_original_values(FE2D FEId, int joint_id, TBaseCell *cell,
                             std::vector<double> quadPoints, int BaseVectDim,
                             std::vector< std::vector<double> > &u00,
                             std::vector< std::vector<double> > &u10,
                             std::vector< std::vector<double> > &u01);
};

#endif


/** ************************************************************************
 * @brief     Assemble boundary integrals to the matrix and add their values to the rhs
 * @author    Alfonso Caiazzo & Laura Blank
 * @date      28.10.16
 ************************************************************************  */

#ifndef __BoundaryAssembling3D__
#define __BoundaryAssembling3D__

#include <Enumerations.h>
#include <FESpace3D.h>
#include <Database.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Point.h>
#include <BoundFace.h>

class BoundaryAssembling3D
{
public:
    
    /** 
      @brief assemble (g, v \cdot n)_{L^2([boundary_component_id])}
      @param[in] boundary_component_id: the boundary component to integrate on
      @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
      @param[in] given_boundary_data: the boundary pressure (as a finite element function)
      @param[in] boundaryCells: the list of all boundary cells 
     */
  void rhs_g_v_n(BlockVector &rhs,
		 const TFESpace3D *U_Space,
		 BoundValueFunct3D *given_boundary_data,
		 std::vector<TBaseCell*> &boundaryCells,
		 int bd_component_id, double mult
		 );


  void rhs_g_v_n(BlockVector &rhs,
		 const TFESpace3D *U_Space,
		 BoundValueFunct3D *given_boundary_data,
		 std::vector<TBoundFace*> boundaryFaceList,
		 int componentID,
		 double mult);



  
    /** @brief integral (u, v)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     @param[in] rescale_by_h: true: divide by length of the edges
     false: do not divide by length of the edges
     */
    void matrix_u_v(BlockFEMatrix &M,
                    const TFESpace3D *U_Space,
                    std::vector<TBoundFace*> boundaryFaceList,
                    int componentID,
                    double mult,
                    bool rescale_by_h);

    
    /** @brief integral (given 3D function, v)_{L^2([boundary_component_id])}
     @param[in] [given_boundary_data1,given_boundary_data2]: the e.g. 3D Dirichlet boundary velocity (as a finite element function)
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     @param[in] rescale_by_h: true: divide by length of the edges
     false: do not divide by length of the edges
     */
    void rhs_uD_v(BlockVector &rhs,
                 const TFESpace3D *U_Space,
                 BoundValueFunct3D *given_boundary_data1,
                 BoundValueFunct3D *given_boundary_data2,
                 BoundValueFunct3D *given_boundary_data3,
                 std::vector<TBoundFace*> boundaryFaceList,
                 int componentID,
                 double mult,
                 bool rescale_by_h);

   
    /** @brief integral (\nabla v \cdot n, u)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void matrix_gradu_n_v(BlockFEMatrix &M,
                          const TFESpace3D *U_Space,
                          std::vector<TBoundFace*> boundaryFaceList,
                          int componentID,
                          double mult);
    
    
    /** @brief integral (\nabla v \cdot n, u)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void matrix_gradv_n_u(BlockFEMatrix &M,
                          const TFESpace3D *U_Space,
                          std::vector<TBoundFace*> boundaryFaceList,
                          int componentID,
                          double mult);
    
    
    /** @brief integral (\nabla v \cdot n, given 3D function)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void rhs_gradv_n_uD(BlockVector &rhs,
                        const TFESpace3D *U_Space,
                        BoundValueFunct3D *given_boundary_data1,
                        BoundValueFunct3D *given_boundary_data2,
                        BoundValueFunct3D *given_boundary_data3,
                        std::vector<TBoundFace*> boundaryFaceList,
                        int componentID,
                        double mult);


    /** @brief integral (\nabla v \cdot n, given 3D function)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void matrix_p_v_n(BlockFEMatrix &M,
                      const TFESpace3D *U_Space,
                      const TFESpace3D *P_Space,
                      std::vector<TBoundFace*> boundaryFaceList,
                      int componentID,
                      double mult);
    
    
    /** @brief integral (\nabla v \cdot n, given 3D function)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void matrix_q_u_n(BlockFEMatrix &M,
                      const TFESpace3D *U_Space,
                      const TFESpace3D *P_Space,
                      std::vector<TBoundFace*> boundaryFaceList,
                      int componentID,
                      double mult);
    
    
    /** @brief integral (\nabla v \cdot n, given 3D function)_{[boundary_component_id]}
     @param[in] boundary_component_id: the boundary component to integrate on
     @param[in] mult: given multiplicative factor (e.g., viscosity or a penalty)
     */
    void rhs_q_uD_n( BlockVector &rhs,
                    const TFESpace3D *V_Space,
                    const TFESpace3D *P_Space,
                    BoundValueFunct3D *given_boundary_data1,
                    BoundValueFunct3D *given_boundary_data2,
                    BoundValueFunct3D *given_boundary_data3,
                    std::vector<TBoundFace*> boundaryFaceList,
                    int componentID,
                    double mult);

    
    /** @brief Get the quadrature points and weights (according to FaceQuadFormula = Quadrature formula for specified degree which is computed inside this member function) on a joint
     @param[in] *fespace: Finite element space
     @param[in] *cell: Cell
     @param[in] m: Joint ID
     @param[out]  qWeights: Quadrature weights
     @param[out] qPointsT: Quadrature points (Gauss points)
     @param[out] values: Values of the FE basis functions on the Gauss points
     */
    void getQuadratureData(const TFESpace3D *fespace,
                           TBaseCell *cell,
                           int m,
                           std::vector<double>& qWeights,
                           std::vector<double>& qPointsT,
                           std::vector<double>& qPointsS,
                           std::vector< std::vector<double> >& values);
    
    
    /** @brief Get the quadrature points and weights (according to FaceQuadFormula = Quadrature formula for specified degree which is computed inside this member function) and corresponding values of the basis functions on a joint
     @param[in] *fespace: Finite element space
     @param[in] *cell: Cell
     @param[in] m: Joint ID
     @param[out]  qWeights: Quadrature weights
     @param[out] qPointsT: Quadrature points (Gauss points)
     @param[out] values: Values of the FE basis functions on the Gauss points
          @param[out] basisFunctionsValues_derivative_x: Values of the partial derivatives of the FE basis functions in x-direction on the Gauss points
          @param[out] basisFunctionsValues_derivative_y: Values of the partial derivatives of the FE basis functions in y-direction on the Gauss points
          @param[out] basisFunctionsValues_derivative_z: Values of the partial derivatives of the FE basis functions in z-direction on the Gauss points
     */
    void getQuadratureDataIncludingFirstDerivatives(const TFESpace3D *fespace,TBaseCell *cell,
                                                    int m,
                                                    std::vector<double>& qWeights,
                                                    std::vector<double>& qPointsT,
                                                    std::vector<double>& qPointsS,
                                                    std::vector< std::vector<double> >& basisFunctionsValues,
                                                    std::vector< std::vector<double> >& basisFunctionsValues_derivative_x,
                                                    std::vector< std::vector<double> >& basisFunctionsValues_derivative_y,
                                                    std::vector< std::vector<double> >& basisFunctionsValues_derivative_z);
    
    
    /** @brief For a joint in a cell with joint_id=m compute the normal and the Determinant (later used for transformation to reference element and backwards)
     @param[in] *cell: Cell
     @param[in] m: Joint ID
     @param[out] normal: Normal vector on joint with joint ID m
     @param[out] transformationDeterminant: Determinant of the joint with joint ID m
     */
    void computeNormalAndTransformationData(TBaseCell *cell,
                                            int m,
                                            std::vector<double>& normal,
                                            double& transformationDeterminant);
    
    
//    /** @brief For a joint in a cell with joint_id=m compute the number of vertices
//     @param[in] *cell: cell
//     @param[in] m: joint ID
//     */
//    void compute_h(TBaseCell *cell,
//                   int m,
//                   double &h);
    
    
    /**
      @brief Assemble the terms needed for Nitsche BC on a given boundary
    */
    void nitsche_bc(BlockFEMatrix &s_matrix,BlockVector &s_rhs,
        const TFESpace3D * v_space, const TFESpace3D *p_space,
        BoundValueFunct3D * U1, BoundValueFunct3D *U2, BoundValueFunct3D *U3,
        std::vector<TBoundFace*> boundaryFaceList,
        int bd_comp, double gamma, double mu,
        int sym_u, int sym_p);



    /** @brief For a joint in a cell with joint_id=m compute the number of vertices
     @param[in] *cell: cell
     @param[in] m: joint ID
     */
    int getNumberOfFaceVertices(TBaseCell *cell, int m);
    
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
    void get_original_values(FE3D FEId, int joint_id, TBaseCell *cell,
                             std::vector<double> quadPoints, int BaseVectDim,
                             std::vector< std::vector<double> > &u000,
                             std::vector< std::vector<double> > &u100,
                             std::vector< std::vector<double> > &u010,
                             std::vector< std::vector<double> > &u001);
};

#endif


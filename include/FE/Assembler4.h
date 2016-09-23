/** ***************************************************************************
 *
 * @name   Assembler
 * @brief  Class with the purpose of performing all operations
 *         needed for assembling a FE linear system
 *
 *
 * @author     Alfonso Caiazzo and Laura Blank
 * @date       16.08.2016
 *
 ******************************************************************************/

#pragma once

#include <Constants.h>
#include <FEDatabase2D.h>
#include <LocalAssembling2D.h>
#include <BlockVector.h>
#include <BlockFEMatrix.h>
#include <Example2D.h>

#ifdef __3D__
#include <Aux2D3D.h>
#endif

/**
 * This class should eventually replace the function Assemble2D(...)
 * Currently it works only for NStype = 14, and it should be extended to
 * other types of BlockMatrix. 
 * 
 */

class Assembler4{
    
protected:
    
public:
    
    /** @brief constructor (taken from LocalAssembling)
     *
     * In the current version, the Assembler object is created using
     * a constructor analogous to LocalAssembling.
     * This is just a temporary solution and should probably done differently
     * Note: the LocalAssembling2D is only needed within a function of this class
     */
    Assembler4(LocalAssembling2D_type type, TFEFunction2D **fefunctions2d,
               CoeffFct2D *coeffs);
    
    ///@brief The variational form to be assembled
    LocalAssembling2D la;
    ///@brief The collection on which we want to assemble
    TCollection *Coll;
    
    ///@brief Vector for handling hanging nodes
    std::vector< std::vector<double> > hangingEntries, hangingRhs;
    
    ///@brief blocks to be assembled
    std::vector< TSquareMatrix2D* > square_matrices;
    std::vector< TMatrix2D* > rectangular_matrices;
    std::vector< double* > rhs_blocks;
    
    /// @brief specify sizes of the linear system
    int n_square_matrices;
    int n_rectangular_matrices;
    int n_rhs_blocks;
    int n_all_matrices;
    
    int maximum_number_base_function;
    
    /**
     @brief Set the variables of the class (hanging nodes, allocate pointers)
     
     This is a temporary solution.
     */
    void init(BlockFEMatrix &M,BlockVector &b_rhs,
              std::vector<const TFESpace2D*>& fespaces,
              std::vector<const TFESpace2D*>& ferhs);
    
    /** a function from a finite element space */
    void Assemble2D(BlockFEMatrix &M,BlockVector &b_rhs,
                    std::vector <const TFESpace2D*>& fespaces,
                    std::vector <const TFESpace2D*>& ferhs,
                    const Example2D& example,
                    int AssemblePhaseID=-1
                    );

    /** 
	@brief assemble the selected local form on the cell i
	@param[in] fespaces: the array of fe spaces of the matrix blocks
	@param[in] i: the cell number
	@param[inout] LocMatrices: pointer to local matrices
	@param[inout] LocRhs: pointer to local rhs
    **/
    void assemble_local_system(std::vector <const TFESpace2D*>& fespaces,
                               int i, double ***LocMatrices,double **LocRhs);

    /** 
	@brief assemble the boundary conditions
	@param[in] fespace: the fe space of the considered rhs block
	@param[in] example: the example class where the BC (values) are specified
	@param[in] cell, i: the cell
	@param[in] j: the considered block of the rhs
	@param[inout] RHS: pointer to rhs_blocks[j]
	@param[in] DOF: pointer to global DOF: DOF = ferhs[j]->GetGlobalDOF(i)
    */
    void impose_boundary_conditions(const TFESpace2D *fespace,
				    const Example2D& example,
				    TBaseCell* cell,
				    int i, int j,
				    double *RHS,
				    int *DOF);
 
    /** 
	@brief add the local contribution to the global matrix
	@param[in] i: the cell
	@param[in] LocMatrices: pointer to local square matrices
	@param[in] Matrices: pointer to local rectangular matrices
    */
    void add_local_to_global_matrix(int i,
                                double ***LocMatrices,
                                double **Matrices);
    
    /** 
	@brief add the local contribution to the global rhs
	@param[in] i: the cell
	@param[in] ferhs: fe spaces of rhs blocks
	@param[in] LocRhs: pointer to rhs
	@param[in] example: used for boundary conditions (boundary values)
    */
    void add_local_to_global_rhs(int i,
                             std::vector<const TFESpace2D*>& ferhs,
                             double **LocRhs,
                             const Example2D& example);

    /**
       @brief handle hanging nodes in matrix and rhs
       @param[in] ferhs spaces of the rhs blocks

       Note: the matrix fe spaces are taken from the matrices 
     */
    void handle_hanging_nodes(std::vector<const TFESpace2D*>& ferhs);
 
    
    ~Assembler4();
};

/** ***************************************************************************
 *
 * @name   2D_Assembler
 * @brief  Assemble
 *
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
 @brief Todo...
 */

class Assembler4{
    
protected:
    
public:
    
    /** @brief constructor (taken from LocalAssembling)
     *
     * In the current version, the Assembler object is created using
     * a constructor analogous to LocalAssembling.
     * This is just a temporary solution and should probably done differently
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
    
    // new
    int n_square_matrices;
    int n_rectangular_matrices;
    int n_rhs_blocks;
    int N_AllMatrices;
    
    /** a function from a finite element space */
    void Assemble2D(BlockFEMatrix &M,BlockVector &b_rhs,
                    std::vector <const TFESpace2D*>& fespaces,
                    std::vector <const TFESpace2D*>& ferhs,
                    const Example2D& example,
                    int AssemblePhaseID=-1
                    );
    
    /**
     @brief Set the variables of the class (hanging nodes, allocate pointers)
     
     This is a temporary solution.
     */
    void init(BlockFEMatrix &M,BlockVector &b_rhs,
              std::vector<const TFESpace2D*>& fespaces,
              std::vector<const TFESpace2D*>& ferhs);
    
    
    void impose_boundary_conditions(const TFESpace2D *fespace,
                                  const Example2D& example,
                                  TBaseCell* cell,
                                  int i, int j,
                                  double *RHS,
                                  int *DOF);
    
    
    void handle_hanging_nodes(std::vector<const TFESpace2D*>& ferhs);

    // new
    void add_local_to_global_matrix(int i,
                                double ***LocMatrices,
                                double **Matrices);
    // new
    void add_local_to_global_rhs(int i,
                             std::vector<const TFESpace2D*>& ferhs,
                             double *righthand,
                             const Example2D& example);
    
    void assemble_local_system(std::vector <const TFESpace2D*>& fespaces,
                               int i, double ***LocMatrices,double **LocRhs);
    
    
    ~Assembler4();
    
    
    
    
    
};

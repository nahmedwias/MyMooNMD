/** ***************************************************************************
 *
 * @name       CD2D
 * @brief      store everything needed to solve a convection-diffusion-reaction
 *             (cdr) problem
 *
 *             Store matrix, right hand side, FE spaces, FE functions and 
 *             the solution vector of a convection-diffusion problem. This 
 *             wraps up everything which is necessary to solve a convection 
 *             diffusion problem in 2D.
 *
 * @author     Ulrich Wilbrandt
 * @date       06.09.13
 *
 ******************************************************************************/

#ifndef __CD2D_H__
#define __CD2D_H__

#include <FEFunction2D.h>
#include <BlockMatrixCD2D.h>
#include <Example_CD2D.h>
#include <vector>

class CD2D
{
  protected:
    /** @brief the system matrix (here one block)
     * 
     * More entries in this vector only for multigrid.
     */
    std::vector<BlockMatrixCD2D*> matrix;
    
    /** @brief the right hand side vector 
     * 
     * More entries in this vector only for multigrid.
     */
    std::vector<double*> rhs;
    
    /** @brief Finite Element function 
     * 
     * The finite element function knows its space and finite element vector so
     * that those two are not explicitly stored in this class.
     * 
     * More entries in this vector only for multigrid. 
     */
    std::vector<TFEFunction2D*> function;
    
    /** @brief Definition of the used example */
    const Example_CD2D* example;
    
    /** @brief a multigrid object which is set to nullptr in case it is not 
     *         needed
     */
    TMultiGrid2D * multigrid;
    
    /** @brief set parameters in database
     * 
     * This functions checks if the parameters in the database are meaningful 
     * and resets them otherwise. The hope is that after calling this function
     * this class is fully functional. 
     * 
     * If some parameters are set to unsupported values, an error occurs and 
     * throws an exception.
     */
    void set_parameters();
    
  public:
    /** @brief constructor 
     * 
     * The domain must have been refined a couple of times already. On the finest
     * level the finite element spaces and functions as well as matrices, 
     * solution and right hand side vectors are initialized. 
     */
    CD2D(TDomain *domain, const Example_CD2D* _example = NULL);
    
    /** @brief standard destructor */
    ~CD2D();
    
    /** @brief assemble matrix, 
     * 
     * depending on 'TDatabase::ParamDB->DISCTYPE' different (local) assembling 
     * routines are used. Also in case of multigrid the matrices on all grids are
     * assembled.
     */
    void assemble();
    
    /** @brief solve the system */
    void solve();
    
    /** 
     * @brief measure errors and write pictures 
     * 
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk files are created.
     * 
     * @param i suffix for output file name, -1 means no suffix
     */
    void output(int i = -1);
    
    // getters and setters
    BlockMatrixCD2D* getMatrix() const
    { return matrix[0]; }
    double* getRhs() const
    { return rhs[0]; }
    TFEFunction2D *get_function() const
    { return function[0]; }
    TFESpace2D * getSpace() const
    { return function[0]->GetFESpace2D(); }
    double * getSolution() const
    { return function[0]->GetValues(); }
    unsigned int getSize() const
    { return function[0]->GetLength(); }
    const Example_CD2D* getExample() const
    { return example; }
};

#endif // __CD2D_H__

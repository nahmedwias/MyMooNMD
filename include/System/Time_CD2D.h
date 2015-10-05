/** ************************************************************************
 * 
 * @name         Time_CD2D
 * @brief        store everything needed to solve a time dependent convection 
 *               diffusion reaction problem
 *               Stores matrix, right hand side, FE spaces, FE functions 
 *               and the solution vector
 *               
 * @author       Naveed Ahmed, Ulrich Wilbrandt, Clemens Bartsch
 * @History      16.09.2015
**************************************************************************/

#ifndef __Time_CD2D__
#define __Time_CD2D__

#include <FEFunction2D.h>
#include <BlockMatrixCD2D.h>
#include <BlockVector.h>
#include <Example_CD2D.h>

#include <vector>
#include <deque>

class Time_CD2D
{
  protected:
    /** @brief store a complete system on a paticular grid.
     * 
     * This combines a matrix, rhs, solution, spaces and functions 
     * needed to describe a Time CDR problem in 2D
     */
    struct System_per_grid
    {
      /** @brief Finite element space */
      TFESpace2D fe_space;
      /** @brief Stiffness Matrix */
      BlockMatrixCD2D Stiff_matrix;
      /** @brief Mass matrix */
      BlockMatrixCD2D Mass_Matrix;
      /** @brief right hand side vector */
      BlockVector rhs;
      /** @brief solution vector */
      BlockVector solution;
      /** @brief Finite element function */
      TFEFunction2D fe_function;
      
      /** @brief constructor*/
      System_per_grid(const Example_CD2D& example, TCollection& coll);
    };
    
    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    const Example_CD2D& example;
    
    /** @brief a multigrid object which is set to nullptr in case it is not 
     *         needed
     */
    std::shared_ptr<TMultiGrid2D> multigrid;
    
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
    
    /** @brief old right hand side vectior 
     * this will be used to save the right hand side from the 
     * previous time step that will be used for different 
     * time stepping schemes
     */
    BlockVector old_rhs;
    
    /** @brief store the errors to compute accumulated error norms */
    std::vector<double> errors;
    
  public:
    /** @brief constructor
     * This constructor calls the other constructor creating an Example_CD2D
     * object. 
     */
    Time_CD2D(const TDomain& domain, int reference_id = -4711);
    
    /** @brief constructor 
     * 
     * The domain must have been refined a couple of times already. On the 
     * finest level the finite element spaces and functions as well as 
     * matrices, solution and right hand side vectors are initialized. 
     * 
     * The reference_id can be used if only the cells with the give reference_id
     * should be used. The default implies all cells.
     */
    Time_CD2D(const TDomain& domain, const Example_CD2D& ex, 
              int reference_id = -4711);
    
    /** @brief Assemble all the matrices before the time iterations
     * 
     * This includes the assembling of: Stiff_matrix, Mass_Matrix, 
     * (additional matrixK in case of SUPG stabilization), rhs
     */
    void assemble_initial_time();
    
    /** @brief preparing the right hand side, 
     * this function will prepare the right during the time iterations
     * that will be used for solving the system
     * NOTE: check whether needed or delete
     */
    // void prepare_rhs();
    /** @brief assemble the matrices
     * this function will assemble the stiffness matrix and rhs
     * In addition the system matrix and the rhs which passes to the solver 
     * are also prepared within the function
     */
    void assemble();
    
    /** @brief solve the system
     */
    void solve();
    
    /** @brief measure errors and write solution
     * 
     */
    void output(int m, int& imgage);
};

#endif
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
#include <BlockVector.h>
#include <Example_CD2D.h>
#include <vector>
#include <deque>

class CD2D
{
  protected:
    
    /** @brief store a complete system on a particular grid
     * 
     * This combines a matrix, rhs, solution, space and function needed to 
     * describe one convection-diffusion-reaction problem in 2D.
     */
    struct System_per_grid
    {
      /** @brief Finite Element space */
      TFESpace2D fe_space;
      /** @brief the system matrix */
      BlockMatrixCD2D matrix;
      /** @brief the right hand side vector */
      BlockVector rhs;
      /** @brief solution vector with one component. */
      BlockVector solution;
      /** @brief Finite Element function */
      TFEFunction2D fe_function;
      
      /** @brief constructor */
      System_per_grid( const Example_CD2D& example, TCollection& coll );
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
    
  public:
    
    /** @brief constructor 
     * 
     * This constructor calls the other constructor creating an Example_CD2D
     * object for you. See there for more documentation.
     */
    CD2D(const TDomain& domain, int reference_id = -4711);
    
    /** @brief constructor 
     * 
     * The domain must have been refined a couple of times already. On the 
     * finest level the finite element spaces and functions as well as 
     * matrices, solution and right hand side vectors are initialized. 
     * 
     * The reference_id can be used if only the cells with the give reference_id
     * should be used. The default implies all cells.
     */
    CD2D(const TDomain& domain, const Example_CD2D&, int reference_id = -4711);
    
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
    const BlockMatrixCD2D & get_matrix() const
    { return this->systems.front().matrix; }
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    const TFEFunction2D & get_function() const
    { return this->systems.front().fe_function; }
    const TFESpace2D & get_space() const
    { return this->systems.front().fe_space; }
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    const Example_CD2D& get_example() const
    { return example; }
};

#endif // __CD2D_H__

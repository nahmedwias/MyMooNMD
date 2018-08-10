/** ************************************************************************
 * 
 * @name         Time_CD2D_POD
 * @brief        store everything needed to solve a time dependent convection 
 *               diffusion reaction problem
 *               Stores matrix, right hand side, FE spaces, FE functions 
 *               and the solution vector
 *               
 * @author       Naveed Ahmed, Ulrich Wilbrandt, Clemens Bartsch
 * @History      16.09.2015
**************************************************************************/

#ifndef __Time_CD2D_POD
#define __Time_CD2D_POD

#include <FEFunction2D.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Example_TimeCD2D.h>
#include <Multigrid.h>
#include <Domain.h>
#include <PostProcessing2D.h>
#include <Solver.h>
#include <POD.h>

#include <vector>
#include <deque>

class LocalAssembling2D; //forward declaration

class Time_CD2D_POD
{
  protected:
    /** @brief store a complete system on a particular grid.
     * 
     * This combines a matrix, rhs, solution, spaces and functions 
     * needed to describe a Time CDR problem in 2D
     */
    struct System_per_grid
    {
      /** @brief Finite element space */
      TFESpace2D fe_space;
      /** @brief Mass matrix */
      BlockFEMatrix gramian_matrix;
      /** @brief solution vector */
      BlockVector pod_mode;
      /** @brief Finite element function */
      TFEFunction2D fe_function;

      /** @brief constructor*/
      System_per_grid(const Example_TimeCD2D& example, TCollection& coll);

      //! Delete copy constructor.
      System_per_grid(const System_per_grid&) = delete;

      //! Delete move constructor.
      System_per_grid(System_per_grid&&) = delete;

      //! Delete copy assignment operator.
      System_per_grid& operator=(const System_per_grid&) = delete;

      //! Delete move assignment operator.
      System_per_grid& operator=(System_per_grid&&) = delete;

      //! Default destructor. Most likely causes memory leaks.
      ~System_per_grid() = default;

    };
    
    /** @brief a local parameter database which controls this class
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * default ParMooN parameters). Note that this usually does not include
     * other parameters such as solver parameters. Those are only in the
     * Solver object.
     */
    ParameterDatabase db;

    /**  @brief a POD object which will compute POD basis and write it to file
     */
    POD pod;

    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger than that
     * only in case of multigrid (when it holds as many systems as there are
     * multigrid levels).
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    const Example_TimeCD2D example;
    
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

    /// @brief class for handling (time dependent) output 
    PostProcessing2D timeDependentOutput;
    
  public:
    /** @brief constructor
     * This constructor calls the other constructor creating an Example_CD2D
     * object. 
     */
    Time_CD2D_POD(const TDomain& domain, const ParameterDatabase& param_db,
                int reference_id = -4711);
    
    /** @brief constructor 
     * 
     * The domain must have been refined a couple of times already. On the 
     * finest level the finite element spaces and functions as well as 
     * matrices, solution and right hand side vectors are initialized. 
     * 
     * The reference_id can be used if only the cells with the give reference_id
     * should be used. The default implies all cells.
     */
    Time_CD2D_POD(const TDomain& domain, const ParameterDatabase& param_db,
                const Example_TimeCD2D& ex, int reference_id = -4711);
    
    /** @brief Assemble all the matrices before the time iterations
     * 
     * This includes the assembling of: Stiff_matrix, Mass_Matrix, 
     * (additional matrixK in case of SUPG stabilization), rhs
     */
    void assemble_gramian();
    
    
    /// @brief measure errors and write solution
    void output();
     // getters and setters
    const Example_TimeCD2D& get_example() const
    { return example; }
    const TFEFunction2D & get_function() const
    { return this->systems.front().fe_function; }
    TFEFunction2D & get_function()
    { return this->systems.front().fe_function; }
    const BlockFEMatrix & get_gramian() const
    { return this->systems.front().gramian_matrix; }
    const BlockVector & get_pod_mode() const
    { return this->systems.front().pod_mode; }
    const TFESpace2D & get_space() const
    { return this->systems.front().fe_space; }
    const ParameterDatabase & get_db() const
    { return db; }

  private:

    /**
     * This wraps up the tedious call to Assemble2D and is only used
     * to avoid code duping. Is really not written very sophisticated,
     * use it with care.
     * @param block_mat should be one system's stiffness or mass matrix.
     * @param la_gramian A fittingly constructed LocalAssemble2D object which
     * is responsible for the assembling of the gramian matrix for the
     * computation of the POD basis. It has to set to represent mass matrix
     * (L2 inner product) or stiffness matrix (H1 inner product) without the
     * viscosity parameter.
     */
    void call_assembling_routine(Time_CD2D_POD::System_per_grid& system,
                                 LocalAssembling2D& la_gramian);
};

#endif

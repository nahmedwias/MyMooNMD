/** ************************************************************************ 
*
* @class     Darcy2D
* @brief     stores the information of a 2D Darcy system matrix 
* @author    Ulrich Wilbrandt
* @date      15.03.15
 ************************************************************************  */


#ifndef __SYSTEMDARCY2D__
#define __SYSTEMDARCY2D__

#include <ColoredBlockFEMatrix.h>
#include <Example_Darcy2D.h>
#include <FEFunction2D.h>
#include <BlockVector.h>

#include <deque>
#include <array>

/**class for 2D scalar system matrix */
class Darcy2D
{
  protected:
    
    /** @brief store a complete system on a particular grid
     * 
     * This combines a matrix, rhs, solution, spaces and functions needed to 
     * describe one Darcy problem in 2D.
     */
    struct System_per_grid
    {
      /** @brief Finite Element space for the velocity */
      TFESpace2D velocity_space;
      /** @brief Finite Element space for the pressure */
      TFESpace2D pressure_space;
      /** @brief the system matrix (here one block) 
       *  [ A  BT ]
       *  [ B  C  ]
       */
      ColoredBlockFEMatrix matrix;
      /** @brief the right hand side vector */
      BlockVector rhs;
      /** @brief solution vector with two components. */
      BlockVector solution;
      /** @brief Finite Element function for velocity */
      TFEFunction2D u;
      /** @brief Finite Element function for pressure */
      TFEFunction2D p;
      
      /** @brief constructor */
      System_per_grid(const Example_Darcy2D& example, TCollection& coll);
    };
    
    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    const Example_Darcy2D& example;
    
    /** @brief a multigrid object which is set to nullptr in case it is not 
     *         needed
     * 
     * @note multigrid for Darcy type problems is not yet implemented
     */
    std::shared_ptr<TMultiGrid2D> multigrid;
    
    /** @brief store the errors to access them from outside this class
     * 
     * This array is filled during a call to Darcy2D::output if 
     * TDatabase::ParamDB->MEASURE_ERRORS is set to true. The exact solution is
     * taken from Darcy2D::example. If that example does not provide an exact 
     * solution, typically it is set to be zero, so that this array contains
     * the norms of the solution instead of the error.
     * 
     * The errors are stored in the following order: 
     * 
     *  - L2 error of velocity
     *  - L2 error of divergence of velocity
     *  - H1-semi error of velocity
     *  - L2 error of pressure
     *  - H1-semi error of pressure
     */
    std::array<double, 5> errors;
    
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
     * This constructor calls the other constructor creating an Example_Darcy2D
     * object for you. See there for more documentation.
     */
    Darcy2D(const TDomain& domain, int reference_id = -4711);
    
    /** @brief constructor 
     * 
     * The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and 
     * functions as well as matrices, solution and right hand side vectors are 
     * initialized. 
     * 
     * The reference_id can be used if only the cells with the give reference_id
     * should be used. The default implies all cells.
     */
    Darcy2D(const TDomain& domain, const Example_Darcy2D& ex, 
            int reference_id = -4711);
    
    /** @brief standard destructor */
    ~Darcy2D();
    
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
     * The current errors will be printed out if 
     * TDatabase::ParamDB->MEASURE_ERRORS is true. The errors will be stored
     * and are then accessible via the get*Error methods. If desired, further 
     * output, e.g., vtk files are created. 
     * 
     * @param i suffix for output file name, -1 means no suffix
     */
    void output(int i = -1);
    
    /// @name return computed errors
    ///
    /// You have to call Darcy2D::output for any of these to return a 
    /// meaningful value.
    //@{
    /// @brief return the computed L2 error of the velocity
    double getL2VelocityError() const;
    /// @brief return the computed L2 error of the divergence of the velocity
    double getL2DivergenceError() const;
    /// @brief return the computed H1-semi error of the velocity
    double getH1SemiVelocityError() const;
    /// @brief return the computed L2 error of the pressure
    double getL2PressureError() const;
    /// @brief return the computed L2 error of the pressure
    double getH1SemiPressureError() const;
    //@}
    
    // getters and setters
    const ColoredBlockFEMatrix & get_matrix() const
    { return this->systems.front().matrix; }
    ColoredBlockFEMatrix & get_matrix()
    { return this->systems.front().matrix; }
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    BlockVector & get_rhs()
    { return this->systems.front().rhs; }
    const TFEFunction2D & get_velocity() const
    { return this->systems.front().u; }
    const TFEFunction2D & get_pressure() const
    { return this->systems.front().p; }
    const TFESpace2D & get_velocity_space() const
    { return this->systems.front().velocity_space; }
    const TFESpace2D & get_pressure_space() const
    { return this->systems.front().pressure_space; }
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    BlockVector & get_solution()
    { return this->systems.front().solution; }
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    const Example_Darcy2D& get_example() const
    { return example; }
};

#endif // __SYSTEMMATDARCY2D__

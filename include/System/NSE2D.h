/** ***************************************************************************
 *
 * @name   NSE2D
 * @brief  store everything needed to solve a Navier-Stokes problem
 *
 *         Store matrix, right hand side, FE spaces, FE functions and the 
 *         solution vector of a Stokes problem. This wraps up everything which 
 *         is necessary to solve a Stokes problem in 2D.
 *
 * @author     Ulrich Wilbrandt
 * @date       06.09.13
 *
 ******************************************************************************/

#ifndef NSE2D_H_
#define NSE2D_H_

#include <FEVectFunct2D.h>

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Residuals.h>

#include <Example_NSE2D.h>

#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <NSE_MGLevel14.h>
#include <MultiGrid2D.h>

#include <MainUtilities.h> // FixedSizeQueue

#include <utility>
#include <array>


class NSE2D
{
  enum class Matrix{Type14, Type1, Type2 , Type3, Type4};

  protected:
    
    /** @brief store a complete system on a particular grid
     * 
     * This combines a matrix, rhs, solution, spaces and functions needed to 
     * describe one stationary Navier-Stokes problem in 2D.
     */
    struct System_per_grid
    {
      /** @brief Finite Element space for the velocity */
      TFESpace2D velocity_space;
      /** @brief Finite Element space for the pressure */
      TFESpace2D pressure_space;

      /** The system matrix. */
      BlockFEMatrix matrix;

      /** @brief the right hand side vector */
      BlockVector rhs;
      /** @brief solution vector with two components. */
      BlockVector solution;
      /** @brief Finite Element function for velocity */
      TFEVectFunct2D u;
      /** @brief Finite Element function for pressure */
      TFEFunction2D p;
      
      /** @brief constructor */
      System_per_grid(const Example_NSE2D& example, TCollection& coll, 
                      std:: pair <int,int> velocity_pressure_orders,
                      NSE2D::Matrix type);

      /**
       * Special member functions mostly deleted,
       * for struct takes ownership of the bad
       * classes TFEFunction2D, TFEVectFunct2D and TFESpace2D.
       */
      //! Delete copy constructor.
      System_per_grid(const System_per_grid&) = delete;

      //! Delete move constructor.
      System_per_grid(System_per_grid&&) = delete;

      //! Delete copy assignment operator.
      System_per_grid& operator=(const System_per_grid&) = delete;

      //! Delete move assignment operator.
      System_per_grid& operator=(System_per_grid&&) = delete;

      //! Default destructor.
      ~System_per_grid() = default;
    };
    
    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    const Example_NSE2D example;
    
    /** @brief a multigrid object which is set to nullptr in case it is not 
     *         needed
     */
    std::shared_ptr<TNSE_MultiGrid> multigrid;
    /// This sorry thing is needed for multigrid with NSTypes 1 or 3, where
    /// transposed blocks are not stored explicitely...sad but true.
    std::vector<std::shared_ptr<TStructure>> transposed_B_structures_;
    
    //! @brief An array to store the current defect.
    BlockVector defect;

    ///@brief The norms of residuals from up to 10 previous iterations
    FixedSizeQueue<10, Residuals> oldResiduals;

    /** @brief store the initial residual so that the nonlinear iteration can 
     *         be stopped as soon as a desired reduction is achieved
     */
    double initial_residual;
    
    /** @brief Errors to be accesed from outside the class
     * The array is filled during the function call NSE2D::output()
     * Currently, the errors store the L2 and H1 errors of the velocity
     * and pressure
     */
    std::array<double, int(4)> errors;
    
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
    
    /** @brief set the velocity and pressure orders
     *
     * This function sets the corresponding velocity and
     * pressure orders. The pressure order is set if it is
     * not specified by the readin file. Default is -4711
     */
    void get_velocity_pressure_orders(std::pair <int,int>
                   &velocity_pressure_orders);
    
  public:
    
    /** @brief constructor 
     * 
     * This constructor calls the other constructor creating an Example_NSE2D
     * object for you. See there for more documentation.
     */
    NSE2D(const TDomain& domain, int reference_id = -4711);
    
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
    NSE2D(const TDomain & domain, const Example_NSE2D & _example,
          unsigned int reference_id = -4711);
    
    /** @brief standard destructor */
    ~NSE2D();
    
    /** @brief assemble matrix, 
     * 
     * This assembles everything which is not related to the nonlinear term.
     * I.e. it assembles a Stokes matrix.
     */
    void assemble();
    
    /** @brief assemble nonlinear term
     * 
     * The matrix blocks to which the nonlinear term contributes are reset to 
     * zero and then completely reassembled, including the linear and nonlinear
     * terms. If this->assemble() has been called before, the matrix is now set
     * up correctly. 
     */
    void assemble_nonlinear_term();
    
    /** @brief solve the system */
    void solve();
    
    /** 
     * @brief Compute the defect Ax-b and store its norm in NSE2D::oldResiduals
     * 
     * where A is the current matrix, x is the current solution and b is the 
     * right hand side. Call this function after assembling the nonlinear
     * matrix with the current solution. 
     */
    void computeNormsOfResiduals();
    
    /** @brief check if one of the stopping criteria is fulfilled
     * 
     * either converged, maximun number of iterations reached, or slow 
     * convergence
     * 
     * @param iteration_counter current iterate
     */
    bool stopIt(unsigned int iteration_counter);
    
    /** 
     * @brief measure errors and write pictures 
     * 
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk files are created.
     * 
     * @param i suffix for output file name, -1 means no suffix
     */
    void output(int i = -1);
    
    /**
   * @brief initialize multigrid levels for different NSTYPE's
   */
    TNSE_MGLevel* mg_levels(int i, System_per_grid& s);

    /** @brief multigrid solver */
    void mg_solver();
    
    // getters and setters
    const BlockFEMatrix & get_matrix() const
    { return this->systems.front().matrix; }
    BlockFEMatrix & get_matrix()
    { return this->systems.front().matrix; }
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    BlockVector & get_rhs()
    { return this->systems.front().rhs; }
    const TFEVectFunct2D & get_velocity() const
    { return this->systems.front().u; }
    TFEVectFunct2D & get_velocity()
    { return this->systems.front().u; }
    // try not to use this as it is not const
    TFEFunction2D *get_velocity_component(int i)
    { return (i==0) ? this->systems.front().u.GetComponent(0)
                    : this->systems.front().u.GetComponent(1); }
    const TFEFunction2D & get_pressure() const
    { return this->systems.front().p; }
    TFEFunction2D & get_pressure()
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
    const Example_NSE2D & get_example() const
    { return example; }
    /// @brief get the current residuals  (updated in NSE2D::normOfResidual)
    const Residuals& getResiduals() const;
    /// @brief get the current impuls residual (updated in NSE2D::normOfResidual)
    double getImpulsResidual() const;
    /// @brief get the current mass residual (updated in NSE2D::normOfResidual)
    double getMassResidual() const;
    /// @brief get the current residual (updated in NSE2D::normOfResidual)
    double getFullResidual() const;
    /// @brief return the computed errors
    std::array<double, int(4)> get_errors() const;
};



#endif /* NSE2D_H_ */

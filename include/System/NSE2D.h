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
#include <Example_NSE2D.h>
#include <BlockMatrixNSE2D.h>
#include <MultiGrid2D.h>
#include <MainUtilities.h> // FixedSizeQueue


#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <NSE_MGLevel14.h>
#include <utility>

class NSE2D
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
      /** @brief the system matrix (depends strongly on 
       *         TDatabase::ParamDB->NSTYPE)
       *  [ A11  A12  B1T ]
       *  [ A21  A22  B2T ]
       *  [ B1   B2   C   ]
       */
      BlockMatrixNSE2D matrix;
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
                      std:: pair <int,int> velocity_pressure_orders);
    };
    
    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    const Example_NSE2D & example;
    
    /** @brief a multigrid object which is set to nullptr in case it is not 
     *         needed
     */
    std::shared_ptr<TNSE_MultiGrid> multigrid;
    
    /** @brief an array to store defect, so that we don't have to reallocate
     *         so often
     */
    BlockVector defect;

    /** @brief stores the norms of the residuals of previous iterations.
     * The default length is 10
     */
    std::vector<double> norms_of_residuals;
    
  public:
    /**
     * @brief a simple struct storing one set of residuals
     * 
     * The full residual is the \f$\ell^2\f$-norm of the vector \f$Ax-b\f$ 
     * where \f$A\f$ is the current matrix and \f$b\f$ the right hand side. It
     * is composed of two parts, the impuls and the residual.
     * 
     * If not default constructed it holds 
     *     fullResidual*fullResidual = impulsResidual*impulsResidual
     *                                 +massResidual*massResidual
     */
    struct Residuals
    {
      /// @brief the impuls residual
      double impulsResidual;
      /// @brief the mass residual
      double massResidual;
      /// @brief the fulf residual
      double fullResidual;
      ///@brief standard constructor, initialize with large numbers
      Residuals();
      /// @brief constructor given the \e square of the impuls and mass 
      /// residuals
      Residuals(double imR, double maR);
      /// @brief write out the three numbers to a stream.
      friend std::ostream& operator<<(std::ostream& s, const Residuals& n);
    };
  protected:
    
    /**
     * @brief store the norms of residuals from previous iterations 
     */
    FixedSizeQueue<10, Residuals> oldResiduals;

    /** @brief store the initial residual so that the nonlinear iteration can 
     *         be stopped as soon as a desired reduction is achieved
     */
    double initial_residual;
    
    /** @brief set the velocity and pressure orders
     * 
     * This function sets the corresponding velocity and 
     * pressure orders. The pressure order is set if it is
     * not specified by the readin file. Default is -4711
     */
    void get_velocity_pressure_orders(std::pair <int,int> &velocity_pressure_orders);
    
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
    void normOfResidual();
    
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
    /**
   * @brief multigrid solver
   */
    void mg_solver();
    
    // getters and setters
    const BlockMatrixNSE2D & get_matrix() const
    { return this->systems.front().matrix; }
    BlockMatrixNSE2D & get_matrix()
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
};



#endif /* NSE2D_H_ */

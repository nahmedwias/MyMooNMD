/** ************************************************************************
 * 
 * @name         Time_NSE2D
 * @brief        store everything needed to solve a time dependent convection 
 *               diffusion reaction problem
 *               Stores matrix, right hand side, FE spaces, FE functions 
 *               and the solution vector
 *               
 * @author       Naveed Ahmed
 * @History      16.09.2015
**************************************************************************/

#ifndef __TIME_NSE2D__
#define __TIME_NSE2D__

#include <FEVectFunct2D.h>

#include <BlockFEMatrix.h>
#include <BlockVector.h>

#include <FESpace2D.h>
#include <Example_NSE2D.h>

#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <NSE_MGLevel14.h>

#include <MainUtilities.h>

#include <vector>
#include <deque>
#include <utility>


class Time_NSE2D
{
  enum class Matrix{Type14, Type1, Type2, Type3, Type4};
  protected:
    /** @brief store a complete system on a paticular grid.
     * 
     * This combines a matrix, rhs, solution, spaces and functions 
     * needed to describe a Time Navier-Stokes problem in 2D
     */
    struct System_per_grid
    {
      /** @brief Finite element space for velocity*/
      TFESpace2D velocity_space;
      /** @brief Finite element space for pressure*/
      TFESpace2D pressure_space;
      
      /** @brief system matrix
       *  [ A11  A12  B1 ]
       *  [ A21  A22  B2 ]
       *  [ B3   B4   C  ]
      */
      BlockFEMatrix matrix;
      /** @brief Mass matrix
       *  [ M11  0    0 ]
       *  [ 0    M22  0 ]
       *  [ 0    0    0 ]
       */
      BlockFEMatrix Mass_Matrix;
      /** @brief right hand side vector*/
      BlockVector rhs;
      /** @brief solution vector*/
      BlockVector solution;
      /** @brief Finite element function for velocity*/
      TFEVectFunct2D u;
      /** @brief Finite element function for pressure*/
      TFEFunction2D p;
      
      /** @brief constructor*/
      System_per_grid(const Example_NSE2D& example, TCollection& coll, 
                      std::pair<int,int> order, Time_NSE2D::Matrix type);
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

    protected:
      /**
       * @brief store the square root of the residual from previous iteration
       */
      double oldResidual;
      /** @brief store the initial residual so that the nonlinear iteration can 
       *         be stopped as soon as a desired reduction is achieved
       */
      double initial_residual;
      
      /** @brief right hand side vector from previous time step (on finest mesh)*/
      BlockVector old_rhs;
      
      BlockVector old_solution;
      
      /** @brief store errors  */
      std::vector<double> errors;
    
      /** old time step length used to scale the pressure blocks*/
      double oldtau;
      
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
      
      /** @brief get velocity and pressure space*/
      void get_velocity_pressure_orders(std::pair <int,int> &velocity_pressure_orders);
  public:
    /** @brief constructor
     * This constructor calls the other constructor creating an Example_CD2D
     * object. 
     */
    Time_NSE2D(const TDomain& domain, int reference_id = -4711);
    /** @brief constructor 
     * 
     * The domain must have been refined a couple of times already. On the 
     * finest level the finite element spaces and functions as well as 
     * matrices, solution and right hand side vectors are initialized. 
     * 
     * The reference_id can be used if only the cells with the give reference_id
     * should be used. The default implies all cells.
     */
    Time_NSE2D(const TDomain& domain, const Example_NSE2D& ex, 
	       int reference_id = -4711);
    
    /** @brief Assemble all the matrices and rhs before the time iterations
     * 
     * This includes the assembling of: Stiff_matrix, Mass_Matrix, 
     * (additional matrixK in case of SUPG stabilization), rhs
     */
    void assemble_initial_time();
    
    /** @brief 
     * 1. assembling the right hand side
     * 2. scaling of the B-blocks due to time stepping
     * this function will prepare the right hand side during the time 
     * discretization but should be outside the nonlinear loop
     */
    void assemble_rhs();
    
    /** @brief Assemble the right hand side only
     *  @param la: LocalAssembling2D object
     *  @param rhs: rhs vector to be assembled
    */
    void AssembleRhs();
    
    /** @brief Assemble the system matrix
     * This function will prepare the system which will be 
     * used for solvers
     */
    void assemble_system();
    
    /** @brief assemble nonlinear term
     * 
     * The matrix blocks to which the nonlinear term contributes are reset to 
     * zero and then completely reassembled, including the linear and nonlinear
     * terms. If this->assemble() has been called before, the matrix is now set
     * up correctly. 
     */
    void assemble_nonlinear_term();
    
    /** @brief check if one of the stopping criteria is fulfilled
     * 
     * either converged, maximun number of iterations reached, or slow 
     * convergence
     * 
     * @param it_counter current iterate
     */
    bool stopIte(unsigned int it_counter);
    
    /** @brief solve the system */
    void solve();
    
    /** descale matrices
     * This function will descale all A-blocks which were scaled
     * during the function call Time_NSE2D::assemble_system():
     */
    void deScaleMatrices();
    /**
     * @brief initialize multigrid levels for different NSTYPE's
     */
    TNSE_MGLevel* mg_levels(int i, System_per_grid& s);
    /** @brief multigrid solver */
    void mg_solver();
    
    /** @brief 
     * compute errors and write solution
     */
    void output(int m, int &image);
    
    // getters and setters
    /*const BlockMatrixNSE2D & get_matrix() const
    { return this->systems.front().matrix; }*/
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    const TFEVectFunct2D & get_velocity() const
    { return this->systems.front().u; }
    // try not to use this as it is not const
    TFEFunction2D *get_velocity_component(int i)
    { return (i==0) ? this->systems.front().u.GetComponent(0)
                    : this->systems.front().u.GetComponent(1); }
    const TFEFunction2D & get_pressure() const
    { return this->systems.front().p; }
    const TFESpace2D & get_velocity_space() const
    { return this->systems.front().velocity_space; }
    const TFESpace2D & get_pressure_space() const
    { return this->systems.front().pressure_space; }
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    const Example_NSE2D & get_example() const
    { return example; }
};

#endif // __TIME_NSE2D__

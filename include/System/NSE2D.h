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
#include <SystemMatNSE2D.h>
#include <Example_CD2D.h>
#include <vector>
#include <MultiGrid2D.h>


#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <NSE_MGLevel14.h>

class NSE2D
{
  protected:
    /** @brief the system matrix
     * 
     * This matrix consists of several blocks. What exactly these blocks mean
     * and how many there are depends on TDatabase::ParamDB->NSTYPE.
     * 
     * More entries in this vector only for multigrid.
     */
    std::vector<TSystemMatNSE2D*> matrix;
    
    /** @brief the right hand side vector 
     * 
     * More entries in this vector only for multigrid.
     */
    std::vector<double*> rhs;
    
    /** Finite Element functions for velocity (vector) and pressure and the 
     * velocity components (scalars)
     * 
     * The finite element functions know their spaces and finite element vectors
     * so that those two are not explicitly stored in this class.
     * 
     * More entries in these vectors only for multigrid. 
     */
    std::vector<TFEVectFunct2D *> u;
    std::vector<TFEFunction2D *> p;
    std::vector<TFEFunction2D *> u1;
    std::vector<TFEFunction2D *> u2;
    
    /** @brief Definition of the used example */
    const Example_NSE2D* example;
    
    /** @brief a multigrid object which is set to nullptr in case it is not 
     *         needed
     */
    TNSE_MultiGrid * multigrid;
    
    /** @brief an array to store defect, so that we don't have to reallocate
     *         so often
     */
    std::vector<double> defect;

    /** @brief stores the norms of the residuals of previous iterations.
     * The default length is 10
     */
    std::vector<double> norms_of_residuals;

    /** @brief store the initial residual so that the nonlinear iteration can 
     *         be stopped as soon as a desired reduction is achieved
     */
    double initial_residual;
    
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
    NSE2D(TDomain *domain, const Example_NSE2D* _example = NULL);
    
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
     * @brief Compute the residual Ax-b 
     * 
     * where A is the current matrix, x is the current solution and b is the 
     * right hand side. Call this function after assembling the nonlinear
     * matrix with the current solution. 
     */
    double normOfResidual();
    
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
    TNSE_MGLevel* mg_levels(int i, int index);
    /**
   * @brief multigrid solver
   */
    void mg_solver();
    
    // getters and setters
    TSystemMatNSE2D* getMatrix() const
    { return matrix[0]; }
    double* getRhs() const
    { return rhs[0]; }
    TFEVectFunct2D *get_velocity() const
    { return u[0]; }
    TFEFunction2D *get_velocity_component(int i) const
    { return (i==0) ? u1[0] : u2[0]; }
    TFEFunction2D *get_pressure() const
    { return p[0]; }
    TFESpace2D * get_velocity_space() const
    { return u[0]->GetFESpace2D(); }
    TFESpace2D * get_pressure_space() const
    { return p[0]->GetFESpace2D(); }
    double * get_solution() const
    { return u[0]->GetValues(); }
    unsigned int get_size() const
    { return 2*u1[0]->GetLength() + p[0]->GetLength(); }
    const Example_NSE2D* get_example() const
    { return example; }
};



#endif /* NSE2D_H_ */

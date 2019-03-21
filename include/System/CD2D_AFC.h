/** ***************************************************************************
 *
 * @class       CD2D_AFC
 * @brief      store everything needed to solve convection-diffusion-reaction
 *             (cdr) problem using Algebraic Flux Correction (AFC). This is a
 *             derived class of CD2D.
 *
 * This wraps up everything which is necessary to solve a convection diffusion 
 * problem using AFC in 2D.
 *
 * @author     Abhinav Jha
 * @date       06.09.13
 *
 ******************************************************************************/

#ifndef __CD2D_AFC_H__
#define __CD2D_AFC_H__

#include <CD2D.h>
#include <anderson.h>

class Multigrid;

class CD2D_AFC : public CD2D
{
   protected:
    
    /* 
     * NOTE: There are two AFC members defined in CD2D class that we call here.
     * The reason for that is because in do_algebraic_flux_correction(), we need
     * those members and the loop over s requires them to be the member of system
     * hence they need to be defined there.
     */
    
    
    /** @brief set parameters in database
     * 
     * Here we set the AFC parameters.  And an error is thrown if the parameters
     * don't match
     */
    void set_AFC_parameters();
    
    
    /** @brief an interpolant of the exact solution onto the fe space
     *
     * This is required if we want to compute d_h(u_h; , )
     */
    std::vector<double> exact_interpolant;
    
  public:
    
    /** @brief constructor 
     * 
     * This constructor first calls the base class constructor and then assigns the values
     *
     * @param[in] domain The readily treated (refined/partitioned...) domain 
     *                   object. Must not go out of scope before CD2D does!
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] reference_id The cell reference id, of which cells to create
     *                         the TCollection.
     */
    CD2D_AFC(const TDomain& domain, const ParameterDatabase& param_db,
         int reference_id = -4711);
    
    
    /** @brief assemble matrix, 
     * 
     * depending on 'this->db["space_discretization_type]' different (local)
     * assembling routines are used. Also in case of multigrid the matrices
     * on all grids are assembled.
     * 
     * @param[in] iteration This is used mainly to assemble the matrix for Fixed_point_RHS
     *                      only once.
     */
    void assemble(const int iteration);
    
    /** @brief solve the system 
     * 
     * @param[in] iteration current iteration
     */
    /// @brief returns whether or not the method is converged, unimportant for linear discretizations
    bool solve(const int iteration);
    
    /** 
     * @brief measure errors and write pictures 
     * 
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk or case files are created.
     * 
     * @param i suffix for vtk output file name, -1 means no suffix
     */
    void output(int i = -1);
    
    
 // Special member functions. Disable copy/move, set destructor to default.
    // Will be changed only when the underlying classes follow rule of 0/5.

    //! Delete copy constructor.
    CD2D_AFC(const CD2D_AFC&) = delete;

    //! Delete move constructor.
    CD2D_AFC(CD2D_AFC&&) = delete;

    //! Delete copy assignment operator.
    CD2D_AFC& operator=(const CD2D_AFC&) = delete;

    //! Delete move assignment operator.
    CD2D_AFC& operator=(CD2D_AFC&&) = delete;
    
    //Dynamic Damping from [JK08]
    void dynamic_damping(const int iteration);
        
    //Anderson acceleration from [WN11]
    void anderson_acceleration_damping(int N_Unknowns, int iteration, 
     std::list<std::vector<double>> & solAnderson,
     std::list<std::vector<double>> & deltaAnderson);
    
    //! Destructor. Still leaks memory (esp. because of the multigrid objeect).
    ~CD2D_AFC()= default;

  private:
    /**
     * Apply an algebraic flux correction scheme to the assembled matrix.
     * Should be called within the assemble routine, after the actual assembling
     * has been performed with the INTERNAL_FULL_MATRIX_STRUCTURE switch on.
     *
     * Which afc algorithm is performed is determined by switching over
     * ALGEBRAIC_FLUX_CORRECTION.
     */
    
    /* residual: norm of current residual r_k+1
     * residual_old: norm of previous residue r_k */
     double residual, residual_old; 
   
    /* t: time taken to complete one iteration */
     double t;
    /*
     * rhs_flag: Use to denote when one it takes one iteration of Newton as well as one iteration of Fixed_point_RHS more than 2 seconds
     * newton_flag: Use to denote when iteration moves from Newton to Fixed_point_RHS and hence increase the tolerance
     * up_param: The tolerance taken to change fron one schem to another
     * 
     */
    void do_algebraic_flux_correction(const int iteration, const int check_fpr);
    /*
     * time_total: time taken to compute solve the problem
     */    
    double time_total;

    
    /*
     * rejected_steps: rejected iteration steps
     */
    int rejected_steps;

    int newton_iterate, newton_flag;

     
    BlockVector alphas_x_i;
    BlockVector old_solution;
    //Copy the RHS after the first iteration is done so that it can be used for Fixed_point_RHS
    BlockVector rhs_copy;
    //checks whether we have Fixed_point_RHS and iteration>1
    int is_not_afc_fixed_point_rhs;
    //copy for the AFC steady state function
    BlockFEMatrix matrix_copy;

    
    //storage of anderson previous updates
    std::list<std::vector<double>>  solAnderson;
    std::list<std::vector<double>>  deltaAnderson;
   

    
    //number of steps in anderson acceleration
    int kAnderson;
};

#endif // __CD2D_H__

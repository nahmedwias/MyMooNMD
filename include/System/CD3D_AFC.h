/** ***************************************************************************
 *
 * @name       CD3D_AFC
 * @brief      store everything needed to solve a convection-diffusion-reaction
 *             (cdr) problem usign Algebraic Flux Correction (AFC)
 *
 *             Store matrix, right hand side, FE spaces, FE functions and 
 *             the solution vector of a convection-diffusion problem. This 
 *             wraps up everything which is necessary to solve a convection 
 *             diffusion problem in 3D.
 *
 *             Here are some tasks for class CD3D
 *             \todo enable different discretizations/stabilizations
 *             \todo enable direct solver in MPI
 *             \todo test multigrid with all smoothers (so far: only SSOR)
 *
 * @author     Ulrich Wilbrandt, Clemens Bartsch
 * @date       09.06.15, 2015/10/20
 *
 ******************************************************************************/

#ifndef __CD3D_AFC_H__
#define __CD3D_AFC_H__

#include <CD3D.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <ParameterDatabase.h>
#include <Solver.h>

#include <PostProcessing3D.h>

#include <FEFunction3D.h>
#include <Example_CD3D.h>
#include <anderson.h>
#include <PostProcessing3D.h>

#include <vector>
#include <deque>
#include <array>
#include <list>

class LocalAssembling3D; //forward declaration
class Example_CD3D;


class CD3D_AFC : public CD3D
{
  protected:
    /** @brief store a complete system on a particular grid
     * 
     * This combines a matrix, rhs, solution, space and function needed to
     * describe one convection-diffusion-reaction problem in 3D.
     * In MPI case the parallel infrastructure is stored, too.
     */
    
     /** @brief an interpolant of the exact solution onto the fe space */
    std::vector<double> exact_interpolant;

  public:

    /** @brief The standard constructor, can be used for multigrid and non-multigrid.
     *
     * It is user responsibility to call this constructor with correct TCollections.
     * Here is a "safe" way to do it.
     *
     * When not using multigrid, this must be of length 1,
     * containing a reference to a collection gained e.g. by calling
     * domain.GetCollection(It_Finest, 0) in the main program.
     *
     * In multigrid case, the collections vector must contain TCollections gained
     * by calling domain.GetCollection(It_Finest, 0) in the main program repeatedly,
     * directly after the corresponding refinement step (and after Domain_Crop in MPI
     * case).
     *
     *
     * @param[in] collections A hierarchy of cell collections used for multigrid solve,
     * or just one collection in non-multigrid case. Ordered by fineness of the grid -
     * finest collection first!
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] example The example which is to be calculated.
     * @param[in] maxSubDomainPerDof Only in MPI case! the maximal number of
     * processes which share a single dof. The value is calculated in the
     * main program and handed down to the FESpaces. Getting rid of this
     * construction is a TODO .
     */
#ifdef _MPI
    CD3D_AFC(std::list<TCollection* > collections,
         const ParameterDatabase & param_db, const Example_CD3D& example,
         int maxSubDomainPerDof);
#else
    CD3D_AFC(std::list<TCollection* > collections, 
         const ParameterDatabase & param_db, const Example_CD3D& _example);
#endif
    
    /** @brief Assemble the system matrix resp. matrices in multigrid case.
     * 
     * Depending on 'this->db["space_discretization_type]' different (local)
     * assembling routines are supposed to be used - so far only standard
     * "galerkin" ("1") is supported.
     */
    void assemble(const int iteration);
    
    /** @brief Solve the system.
     *
     * Have a look at the code to see which solvers are already enabled.
     * More are to come.
     */
    bool solve(const int iteration);
    
    /** 
     * @brief Measure errors and write nice pictures.
     * 
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk files are created.
     *
     * @param i Integer suffix for output file name. -1 means no suffix.
     */
    void output(int i = -1);
    


    // getters and setters

    
    // Declaration of special member functions - delete all but destructor.
    // TODO CB Apply rule of zero/five as soon as underlying classes do.

    /** Delete copy constructor. This class may not be copied,
     * until its member classes fulfil the rule of 0/5. */
    CD3D_AFC(const CD3D_AFC&) = delete;

    /** Delete move constructor. This class may not be moved yet. */
    CD3D_AFC(CD3D_AFC&&) = delete;

    /** Delete copy assignment. This class may not be copied yet. */
    CD3D_AFC& operator=(const CD3D_AFC&) = delete;

    /** Delete move assignment. This class may not be moved yet. */
    CD3D_AFC& operator=(CD3D_AFC&&) = delete;
    
    //Dynamic Damping from [JK08]
    void dynamic_damping(const int iteration);
    
    void anderson_acceleration_damping(int N_Unknowns, int iteration, 
     std::list<std::vector<double>> & solAnderson,
     std::list<std::vector<double>> & deltaAnderson);

    //! Default destructor. Does most likely cause memory leaks.
    ~CD3D_AFC() = default;

  private:
    /**
     * Scramble together the parameters which Assemble3D needs and call it.
     * Is only put here to keep the code of assemble() slender.
     * assemble() should take care of the right choice of the LocalAssembling
     * object and whether it fits the system matrix' block structure
     * (which should always be the case in CD3D...).
     * @param s The sytem where rhs and stiffness matrix are to be assembled.
     * @param la_stiff The local assembling object of choice.
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
     * time_rhs: Time taken to compute one Fixed_point_RHS iteration
     */    
    double time_total;
    
    /*
     * rejected_steps: rejected iteration steps
     */
    int rejected_steps;
        
    /*
     * rhs_flag, newton_flag: Used for finding if time taken by one step of Newton as well as one step of Fixed_point_RHS is same
     * newton_iterate: Count the number of Newton iteration
     *rhs_iterate: Count the number of Fixed_point_RHS iteration     
     */
    int newton_flag, newton_iterate;

     
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
        
    void call_assembling_routine(SystemPerGrid& s, LocalAssembling3D& la);

};

#endif // __CD3D_H__

/** ************************************************************************
 *
 * @name         Coupled_Time_CDR_2D
 * @brief        Store everything needed to solve a time dependent
 *               system of convection diffusion reaction problems
 *               which are coupled in the reaction term.
 *               So far there is only one way to decouple the system,
 *               which is by inserting the last known values and put them
 *               onto the right hand side (Richardson iteration,
 *               "linear_decoupled").
 *
 * @author       Clemens Bartsch
 * @date         6.1.2016
**************************************************************************/

#ifndef USER_PROJECTS_INC_COUPLED_TIME_CDR_2D_H_
#define USER_PROJECTS_INC_COUPLED_TIME_CDR_2D_H_


#include <CoupledCDR_2D.h>
#include <Example_TimeCoupledCDR2D.h>
#include <ParameterDatabase.h>

// forward declarations
class Time_CD2D;
class ReactionCoupling;

/*!
 * @brief Control structure for a system of coupled time-dependent
 * CDR equations, which manages the coupled and uncoupled parts.
 *
 * For every time-step the system is decoupled with a simple
 * Richardson-style fixed point iteration so far. Uncoupled CD parts
 * and the extra terms on the right hand side due to this de-coupling
 * are managed seperately.
 */
class Coupled_Time_CDR_2D {

  public:

    /**
     * Standard constructor which gets a domain and an example object.
     * @param domain The domain to compute on.
     * @param example The example to be executed.
     */
    Coupled_Time_CDR_2D(const TDomain& domain,  const ParameterDatabase& db,
                        const Example_TimeCoupledCDR2D& example);

    /**
     * Assemble all matrices and rhs at the 0th time step.
     */
    void assemble_initial_time();


    /**
     * Assembles the decoupled equations for the current time step.
     * Note that the assembling and adding of the coupled parts is done in the
     * loop within the solving routine.
     */
    void assemble_uncoupled_part();



    /**
     * Solve the system. Contains a loop in which the coupled part is assembled,
     * put to the right hand side and the resulting equations is solved.
     * The coupled part may depend on further FE functions (given in further_functions),
     * and the velocity field (which is not used currently).
     */
    void couple_and_solve(
        const TFEVectFunct2D* velocity_field,
        std::vector<const TFEFunction2D*> further_functions = {});

    /**
     * Assembling functions for the case of coupling with a velocity field
     * (and other functions, e.g. a particle size distribution from a PBE).
     */
    void assemble_initial_time(const TFEVectFunct2D* velocity_field,
                               std::vector<const TFEFunction2D*> further_functions = {});

    /// We assume that the "uncoupled" (CD) parts do only depend on the velocity
    /// field, but neither on each other nor on other FE functions .
    void assemble_uncoupled_part(const TFEVectFunct2D* velocity_field);

    /**
     * Produce some output.
     */
    void output();

    //Deletion of special member functions.
    //! Delete copy constructor.
    Coupled_Time_CDR_2D(const Coupled_Time_CDR_2D&) = delete;

    //! Delete move constructor.
    Coupled_Time_CDR_2D(Coupled_Time_CDR_2D&&) = delete;

    //! Delete copy assignment operator.
    Coupled_Time_CDR_2D& operator=(const Coupled_Time_CDR_2D&) = delete;

    //! Delete move assignment operator
    Coupled_Time_CDR_2D& operator=(Coupled_Time_CDR_2D&&) = delete;

    //! Default destructor. Most likely causes memory leaks.
    ~Coupled_Time_CDR_2D() = default;


  private:

    /*! @brief The number of CDR equations forming the system.*/
    size_t nEquations_;

    /*! @brief The involved CD problems - without coupling part! */
    std::vector<std::shared_ptr<Time_CD2D>> cdProblems_;

    /*! @brief The reaction parts which belong to the equation and involve the coupling.*/
    std::vector<std::shared_ptr<ReactionCoupling>> coupledParts_;

    /*! @brief The used example. */
    Example_TimeCoupledCDR2D example_;

    /** A database object holding parameters for controlling the object. */
    ParameterDatabase db_;

    /// Decide whether stopping criterion is hit and solve loop can be broken.
    /// If so, also gives some nice print out.
    /// @param[in] step The current iteration number in the solving loop
    /// @param[out] residuals The current vector of residuals per single equation.
    /// @return Whether to break the solve loop as one of the breaking criteria is hit.
    bool break_iteration(size_t step, std::vector<double> residuals);

};



#endif /* USER_PROJECTS_INC_COUPLED_TIME_CDR_2D_H_ */

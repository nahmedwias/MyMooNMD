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

#include <Example_CoupledCDR2D.h>
#include <CoupledCDR_2D.h>

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
    Coupled_Time_CDR_2D(const TDomain& domain,
                        const Example_CoupledCDR2D& example);

    /**
     * Assemble all matrices and rhs at the 0th time step.
     */
    void assemble_initial_time();

    /**
     * Produce some output.
     * @param[in] image The number of the image, will be used for the naming
     * of the pictures.
     */
    void output(int& image);

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


  protected:

    /*! @brief The number of CDR equations forming the system.*/
    size_t nEquations_;

    /*! @brief The involved CD problems - without coupling part! */
    std::vector<std::shared_ptr<Time_CD2D>> cdProblems_;

    /*! @brief The reaction parts which belong to the equation and involve the coupling.*/
    std::vector<std::shared_ptr<ReactionCoupling>> coupledParts_;

    /*! @brief The used example. */
    Example_CoupledCDR2D example_;

};



#endif /* USER_PROJECTS_INC_COUPLED_TIME_CDR_2D_H_ */

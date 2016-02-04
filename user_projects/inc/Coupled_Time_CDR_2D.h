/** ************************************************************************
 *
 * @name         Coupled_Time_CDR_2D
 * @brief        Store everything needed to solve a time dependent
 *               system of convection diffusion reaction problems
 *               which are coupled in the reaction term.
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

    /*! The strategy for solving the coupled system.
     * (The type is implemented in CoupledCDR_2D) */
    CoupledCDR_2D::SolvingStrategy const strategy_;




};



#endif /* USER_PROJECTS_INC_COUPLED_TIME_CDR_2D_H_ */

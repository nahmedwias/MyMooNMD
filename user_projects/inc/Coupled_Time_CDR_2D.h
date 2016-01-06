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

  public:
    //Declaration of special member functions - rule of zero.
    /**
     * Note that members cdProblems_ and coupledParts_ are only shared pointers
     * so far. So a copied/moved object of class Coupled_Time_CDR_2D will just
     * share ownership of the CD2D and ReactionCoupling objects which
     * the entries of cdProblems_ and coupledParts_ point to.
     */
    // Default copy constructor. Shallow copy!
    Coupled_Time_CDR_2D(const Coupled_Time_CDR_2D&) = default;

    //! Default move constructor.
    Coupled_Time_CDR_2D(Coupled_Time_CDR_2D&&) = default;

    //! Default copy assignment operator. Shallow copy!
    Coupled_Time_CDR_2D& operator=(const Coupled_Time_CDR_2D&) = default;

    //! Default move assignment operator
    Coupled_Time_CDR_2D& operator=(Coupled_Time_CDR_2D&&) = default;

    //! Default destructor.
    ~Coupled_Time_CDR_2D() = default;


};



#endif /* USER_PROJECTS_INC_COUPLED_TIME_CDR_2D_H_ */

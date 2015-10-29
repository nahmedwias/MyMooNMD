/*****************************************************************************
 *  @name CDR_2D_System.h
 *	@brief The class holds a system of strongly coupled stationary
 *	Convection-Diffusion-Reaction equations in 2D and provides
 *	methods to assemble and solve the system.
 *
 *	So far only one strategy to solve the system is implemented.
 *
 *  @date May 6, 2015
 *  @author Clemens Bartsch
 *****************************************************************************/

#ifndef CDR_2D_SYSTEM_H_
#define CDR_2D_SYSTEM_H_

#include <BlockMatrixCD2D.h>
#include <Example_CoupledCDR2D.h>
#include <vector>
#include <memory>

//Forward declarations.
class ReactionCoupling;
class CD2D;
class TDomain;




/*!
 * @brief Control structure for a system of coupled CDR equations, which manages the coupled and uncoupled parts.
 *
 * The idea is to provide a control structure, while the system itself is divided into
 * pure convection-diffusion equations and a structure which takes care of the coupled reaction part.
 * The class CDR_2D_System knows which strategy is used to calculate a solution of the system
 * and how to proceed in every case.
 */
class CDR_2D_System {

public:
	/*!
	 * Three different solution strategies, determines how the systems are assembled and solved.
	 * newton_monolithic - solve the entire linearized system at once by Newton iteration
	 * newton_decoupled  - decouple the systems and solve each seperately by Newton iteration
	 * linear_decoupled  - decouple the systems, put reaction term to rhs (with last known solution)
	 * 					   and solve the linear systems seperately
	 * none				 - Because an enum class requires a default member - constructor throws exception,
	 * 					   when this is chosen.
	 */
	enum class SolvingStrategy{
		newton_monolithic, newton_decoupled, linear_decoupled, none
	};

public:
	/*! @brief Constructor; can be used for multigrid in non-multigrid case.
	 *  @param[in] domain The geometry domain - must be refined.
	 *  @param[in] exam The used example.
	 *  @param[in] strat The desired solving strategy. Defaults to 'none'.
	 */
	CDR_2D_System(const TDomain& domain, const Example_CoupledCDR2D& exam,
			SolvingStrategy strat = SolvingStrategy::none);

	/*! @brief assemble matrix,
	 *
	 * Chooses appropriate assembling subroutine for the chosen solving strategy.
	 */
	void assembleCDPart();

	/*! @brief Solve the system.
	 * Chooses appropriate solving subroutine for the chosen solving strategy and
	 * control parameters from input file.
	 */
	void solve();

	/*! @brief Measure errors and write pictures. */
	void output();


	//Declaration of special member functions - rule of zero.
	/**
	 * Note that members cdProblems_ and coupledParts_ are only shared pointers
	 * so far. So a copied/moved object of class CDR_2D_System will just
	 * share ownership of the CD2D and ReactionCoupling objects which
	 * the entries of cdProblems_ and coupledParts_ point to.
	 */
    // Default copy constructor. Shallow copy!
	CDR_2D_System(const CDR_2D_System&) = default;

    //! Default move constructor.
	CDR_2D_System(CDR_2D_System&&) = default;

    //! Default copy assignment operator. Shallow copy!
	CDR_2D_System& operator=(const CDR_2D_System&) = default;

    //! Default move assignment operator
	CDR_2D_System& operator=(CDR_2D_System&&) = default;

    //! Default destructor.
    ~CDR_2D_System() = default;

protected:

	/*! @brief The number of CDR equations involved. Put it here for the moment.*/
	size_t nEquations_;

	/*! @brief The involved CD problems - make sure these are without coupling part! */
	std::vector<std::shared_ptr<CD2D>> cdProblems_;

	/*! @brief The reaction parts which belong to the equation and involve the coupling.*/
	std::vector<std::shared_ptr<ReactionCoupling>> coupledParts_;

	/*! @brief The used example. Which one is used is determined by the input file.*/
	Example_CoupledCDR2D example_;

	/*! @brief The strategy for solving the coupled system. */
	SolvingStrategy const strategy_;

private:
	/*! @brief Assemble the coupled terms.
	 *
	 * Since these have to be reassembled multiple times, this method is not part of assemble().
	 */
	void assembleCoupledTerms();

};

#endif /*  CDR_2D_SYSTEM_H_ */

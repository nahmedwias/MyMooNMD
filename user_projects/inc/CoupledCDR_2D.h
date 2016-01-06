/*****************************************************************************
 *  @name CoupledCDR_2D.h
 *	@brief The class holds a system of strongly coupled stationary
 *	Convection-Diffusion-Reaction equations in 2D and provides
 *	methods to assemble and solve the system.
 *
 *	So far only one strategy to solve the system is implemented.
 *
 *  @date May 6, 2015
 *  @author Clemens Bartsch
 *****************************************************************************/

#ifndef CoupledCDR_2D_H_
#define CoupledCDR_2D_H_

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
 * The class CoupledCDR_2D knows which strategy is used to calculate a solution of the system
 * and how to proceed in every case.
 */
class CoupledCDR_2D {

public:
	/*!
	 * Three different solution strategies, determines how
	 * the systems are assembled and solved.
	 *
	 * This is still dreams of the future - currently we
	 * are only aiming at the "linear_decoupled" strategy.
	 *
	 * newton_monolithic - solve the entire linearized system at once by Newton iteration
	 * newton_decoupled  - decouple the systems and solve each seperately by Newton iteration
	 * linear_decoupled  - decouple the systems, put reaction term to rhs (with last known solution)
	 * 					           and solve the linear systems seperately
	 * none				       - default member - constructor throws exception, when this is chosen.
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
	CoupledCDR_2D(const TDomain& domain, const Example_CoupledCDR2D& exam,
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
	 * so far. So a copied/moved object of class CoupledCDR_2D will just
	 * share ownership of the CD2D and ReactionCoupling objects which
	 * the entries of cdProblems_ and coupledParts_ point to.
	 */
    // Default copy constructor. Shallow copy!
	CoupledCDR_2D(const CoupledCDR_2D&) = default;

    //! Default move constructor.
	CoupledCDR_2D(CoupledCDR_2D&&) = default;

    //! Default copy assignment operator. Shallow copy!
	CoupledCDR_2D& operator=(const CoupledCDR_2D&) = default;

    //! Default move assignment operator
	CoupledCDR_2D& operator=(CoupledCDR_2D&&) = default;

    //! Default destructor.
    ~CoupledCDR_2D() = default;

protected:

	/*! @brief The number of CDR equations involved. Put it here for the moment.*/
	size_t nEquations_;

	/*! @brief The involved CD problems - make sure these are without coupling part! */
	std::vector<std::shared_ptr<CD2D>> cdProblems_;

	/*! @brief The reaction parts which belong to the equation and involve the coupling.*/
	std::vector<std::shared_ptr<ReactionCoupling>> coupledParts_;

	/*! @brief The used example.*/
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

#endif /*  CoupledCDR_2D_H_ */

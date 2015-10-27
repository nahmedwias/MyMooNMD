#ifndef _EXAMPLE_COUPLEDCDR2D_
#define _EXAMPLE_COUPLEDCDR2D_

#include <memory>
#include <Example_CD2D.h>
//for the typedef CouplingFunction.
#include <Constants.h>


class Example_CoupledCDR2D : public Example2D
{


public:

	//! @brief Constructor. Refers to the input file to chose the correct example.
	Example_CoupledCDR2D();

	//! Destructor. Does nothing special.
	~Example_CoupledCDR2D(){};

	//! Return function pointer to coefficient function of equation nr. equationIndex.
	CoeffFct2D* getCoeffFct(size_t equationIndex) const;

	//! Return function pointer to assembling function of equation nr. equationIndex.
	AssembleFctParam2D* getAssemblingFct(size_t equationIndex) const;

	//! Return function pointer to parameter function of equation nr. equationIndex.
	ParamFct* getParamFct(size_t equationIndex) const;

	/*! Get number of equations
	 * @return The value of nEquations.*/
	size_t getNEquations() const{
		return nEquations_;
	}

	/**
	 * @brief Get one of the underlying CD(R) examples without the coupling term.
	 * TODO CB 27/05/2015 Adapt this to fit all future examples.
	 * @param n 0 or 1 gives c1 or c2 respectively
	 * @return example
	 */
	Example_CD2D* getDecoupledExample(size_t n) const;


private:
	//! Number of contained equations.
	size_t nEquations_;

	//!@brief List of the coefficient function pointers used in the assembling of the uncoupled (convection-diffusion) parts.
	// NOTE: This replaces usage of problem_coefficients from Base class. Imho best of the bad alternatives
	std::vector<CoeffFct2D*> bilinCoeffs_;
	//!@brief List of the coefficient function pointers used in the assembling of the rhs part (linear-decoupled strategy).
	std::vector<AssembleFctParam2D*> rhsAssemblingFunctions_;
	//!@brief List of the parameter function pointers used in the assembling of the rhs part (linear-decoupled strategy).
	std::vector<ParamFct*> parameterFunctions_;


};


#endif

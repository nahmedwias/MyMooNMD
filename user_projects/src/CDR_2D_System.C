/*! ****************************************************************************
 *  @name CDR_2D_System.C
 *	@brief Implements functions of class CDR_2D_System declared in CDR_2D_System.h.
 *
 *  @date May 8, 2015; transferred to ParMooN May 28, 2015
 *  @author Clemens Bartsch
 *****************************************************************************/

#include <AuxParam2D.h>
//#include <MainUtilities.h>
#include <LinAlg.h>

//#include <ItMethod.h>
//#include <MultiGridIte.h>
//#include <FixedPointIte.h>
//#include <FgmresIte.h>
//#include <LocalProjection.h>
#include <Database.h>

#include <CD2D.h>

#include <CDR_2D_System.h>
#include <Example_CoupledCDR2D.h>
#include <CoupledReaction.h>

/*! @brief Standard constructor.*/
CDR_2D_System::CDR_2D_System(const TDomain& domain, const Example_CoupledCDR2D& exam,
		SolvingStrategy strat) :
		example_(exam), strategy_(strat)
{
	/********** Preparation ********************/
	// if strategy is "none" throw exception
	if(strategy_==SolvingStrategy::none){
	    std::stringstream errMsg;
	    errMsg << "No solving strategy chosen for CDR_2D_System.";
	    throw std::runtime_error(errMsg.str().c_str());
	}

	//Everything seems to be all right - start constructing.
	nEquations_ = example_.getNEquations();

	//Code for linearized decoupled solving strategy.
	if(strategy_== SolvingStrategy::linear_decoupled){
		/********** The list of CD problems.  ********************/
		cdProblems_.reserve(nEquations_); //reserve space
		//loop over all equations
		for(size_t index = 0; index<nEquations_;++index){
			//Construct a CD2D problem from the current CD part...
			std::shared_ptr<CD2D> currProblemCD(new CD2D(domain,example_.getDecoupledExample(index)));
			//...and attach it to the list.
			cdProblems_.push_back(currProblemCD);

		}
		/********** The list of coupled parts. ********************/
		//loop over all equations and construct the coupled parts objects
		for(size_t index = 0; index<nEquations_;++index){
			//Get a pointer to the assembling function.
			AssembleFctParam2D* currAssemblingFunction = example_.getAssemblingFct(index);
			//Get a pointer to the param (in-out) function.
			ParamFct* currParamFunction = example_.getParamFct(index);
			// Get a pointer to the rhs space (all spaces the same, ansatz=test (=rhs) )
			const TFESpace2D& currSpace = cdProblems_[index]->get_space();

			//Construct a CoupledReaction from the current coupling function and attach it to the list.
			std::shared_ptr<CoupledReaction> currCoupledPart(
					new CoupledReaction(strategy_, currAssemblingFunction,
							currParamFunction, nEquations_, currSpace));
			coupledParts_.push_back(currCoupledPart);
		}
	} else {
		ErrMsg("Unknown or unimplemented solving strategy chosen.");
	}
}

/*!
 * Assembles the CD part of the system, excluding the coupled part.
 */
void CDR_2D_System::assembleCDPart(){

	/**************** Setting up the CD problems. ***************/
	for(size_t index = 0; index<nEquations_;++index){
		//assemble matrix and rhs of CD part
		cdProblems_[index]->assemble();
	}
}

void CDR_2D_System::solve(){
	switch(strategy_) {

	case SolvingStrategy::linear_decoupled:{

		// Put up an array of pointers to the solutions of previous iteration
		TFEFunction2D** previousSolutions = new TFEFunction2D*[nEquations_];

		// Store the original right hand sides of the CDR Equations without coupling.
		std::vector<BlockVector> originalRightHandSides;
		for (size_t equation = 0; equation<nEquations_;++equation){
			BlockVector vector(cdProblems_[equation]->get_rhs());
			originalRightHandSides.push_back(vector);
		}


		// while(Abbruchbedingung nicht erfuellt) - beginne einfach mit einer festen Anzahl Iterationen.
		for (size_t steps = 0; steps<10;steps++){
			//Fill pointers to available solutions into previousSolutions array
			for(size_t i =0;i<nEquations_;++i){
				previousSolutions[i]=&cdProblems_[i]->get_function();
			}
			//loop over the equations
			for (size_t equation = 0; equation<nEquations_;++equation){

				// assemble the coupling term
				coupledParts_[equation]->assembleLinearDecoupled(previousSolutions);


				// add coupled rhs to rhs of the uncoupled equation
				cdProblems_[equation]->get_rhs().add_scaled(coupledParts_[equation]->getRightHandSide(),1);

				// solve equation with the new right hand side
				cdProblems_[equation]->solve();

				// Set back to original rhs.
				cdProblems_[equation]->get_rhs().copy(originalRightHandSides.at(equation).get_entries());

			}//end loop over equations

		}//endwhile bzw. endfor
		delete[] previousSolutions; previousSolutions = nullptr; //just delete the pointers array
		break;
	}
	case SolvingStrategy::newton_decoupled:
	case SolvingStrategy::newton_monolithic:
	case SolvingStrategy::none:
	default:
		OutPut("No or unsupported solving strategy chosen.");
		exit(0);
	}
}

void CDR_2D_System::output(){
	//Let the work be done by the CDR-Objects for now.
	for (size_t equation = 0; equation<nEquations_;++equation){
		cdProblems_[equation]->output(equation);
	}
}

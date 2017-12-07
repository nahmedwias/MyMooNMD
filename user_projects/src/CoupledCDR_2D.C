/*! ****************************************************************************
 *  @name CoupledCDR_2D.C
 *	@brief Implements functions of class CoupledCDR_2D declared in CoupledCDR_2D.h.
 *
 *  @date May 8, 2015; transferred to ParMooN May 28, 2015
 *  @author Clemens Bartsch
 *****************************************************************************/

#include <AuxParam2D.h>
#include <LinAlg.h>

#include <Database.h>

#include <CD2D.h>

#include <CoupledCDR_2D.h>
#include <Example_CoupledCDR2D.h>
#include <ParameterDatabase.h>
#include <ReactionCoupling.h>

/*! @brief Standard constructor.*/
CoupledCDR_2D::CoupledCDR_2D(
    const TDomain& domain, const ParameterDatabase& db,
    const Example_CoupledCDR2D& exam, SolvingStrategy strat) :
		example_(exam), strategy_(strat)
{
	/********** Preparation ********************/
	// if strategy is "none" throw exception
	if(strategy_==SolvingStrategy::none){
	    std::stringstream errMsg;
	    errMsg << "No solving strategy chosen for CoupledCDR_2D.";
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
		  ParameterDatabase cd2d_db = ParameterDatabase::parmoon_default_database();
		  cd2d_db.merge(db, false);

		  cd2d_db["output_basename"] = cd2d_db["output_basename"].get<std::string>() + std::string("_species_" + index);


			std::shared_ptr<CD2D> currProblemCD(
			    new CD2D(domain, cd2d_db, example_.getDecoupledExample(index)));

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

			//Construct a ReactionCoupling from the current coupling function and attach it to the list.
			std::shared_ptr<ReactionCoupling> currCoupledPart(
					new ReactionCoupling(strategy_, currAssemblingFunction,
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
void CoupledCDR_2D::assembleCDPart(){

	/**************** Setting up the CD problems. ***************/
	for(size_t index = 0; index<nEquations_;++index){
		//assemble matrix and rhs of CD part
		cdProblems_[index]->assemble();
	}
}

void CoupledCDR_2D::solve(){
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

void CoupledCDR_2D::output(){
	//Let the work be done by the CDR-Objects for now.
	for (size_t equation = 0; equation<nEquations_;++equation){
		cdProblems_[equation]->output(equation);
	}
}

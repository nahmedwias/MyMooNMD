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

/*! @brief Constructor to use normally.*/
CDR_2D_System::CDR_2D_System(TDomain* domain, Example_CoupledCDR2D* exam,
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
	nEquations_ = example_->getNEquations();

	/********** The list of CD problems.  ********************/
	//Construct the list of underlying CD problems.
	cdProblems_.reserve(nEquations_); //reserve space
	//loop over all equations
	for(size_t index = 0; index<nEquations_;++index){
		// Get the CD part from current equation of the example
		// FIXME This is hard coded for constant_function example. Does this work with other examples?
		Example_CD2D* currExampleCD = example_->getDecoupledExample(index);

		//Construct a CDR_2D problem from the current CD part and attach it to the list.
		CD2D* currProblemCD = new CD2D(domain,currExampleCD);
		cdProblems_.push_back(currProblemCD);
	}

	/********** The list of coupled parts. ********************/
	//loop over all equations and construct the coupled parts objects
	for(size_t index = 0; index<nEquations_;++index){
		//Get a pointer to the assembling function.
		AssembleFctParam2D* currAssemblingFunction = example_->getAssemblingFct(index);
		//Get a pointer to the param (in-out) function.
		ParamFct* currParamFunction = example_->getParamFct(index);
		// Get a pointer to the rhs space (all spaces the same, ansatz=test (=rhs) )
		TFESpace2D* currSpace = cdProblems_[index]->getSpace();

		//Construct a CoupledReaction from the current coupling function and attach it to the list.
		CoupledReaction* currCoupledPart = new CoupledReaction(currAssemblingFunction, currParamFunction ,nEquations_,currSpace);
		coupledParts_.push_back(currCoupledPart);
	}
}

/*! @brief Standard destructor. */
CDR_2D_System::~CDR_2D_System(){
	//delete matrix;
	//delete rhs;
	//delete solution;
	for (auto cd : cdProblems_) delete cd;
	for (auto cp : coupledParts_) delete cp;
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
		std::vector<double*> originalRightHandSides(nEquations_);
		for (size_t equation = 0; equation<nEquations_;++equation){
			//Copy the entries from the rhs of the cdrProblems to the rhsList
			int size = cdProblems_[equation]->getSize();
			originalRightHandSides[equation] = new double[size];
			Dcopy(size, cdProblems_[equation]->getRhs(), originalRightHandSides[equation]);
		}


		// while(Abbruchbedingung nicht erfuellt) - beginne einfach mit einer festen Anzahl Iterationen.
		for (size_t steps = 0; steps<10;steps++){
			//Fill pointers to available solutions into previousSolutions array
			for(size_t i =0;i<nEquations_;++i){
				previousSolutions[i]=cdProblems_[i]->get_function();
			}
			//loop over the equations
			for (size_t equation = 0; equation<nEquations_;++equation){
				//store size of the current problem (number of dofs)
				int size = cdProblems_[equation]->getSize();
				// assembliere den Kopplungsterme
				coupledParts_[equation]->assembleLinearDecoupled(previousSolutions);
				// add coupled rhs to rhs of the uncoupled equation (using BLAS 1 method Daxpy)
				Daxpy(size,1,coupledParts_[equation]->getRightHandSide(),cdProblems_[equation]->getRhs());
				// Loese die Gleichung mit der neuen rechten Seite
				cdProblems_[equation]->solve();

				// Set back to original rhs. TODO This involves copying and takes long supposedly.
				Dcopy(cdProblems_[equation]->getSize(), originalRightHandSides[equation],  cdProblems_[equation]->getRhs() );


			}//end loop over equations

		}//endwhile bzw. endfor
		delete[] previousSolutions; previousSolutions = nullptr; //just delete the pointers array
		for(auto rhs : originalRightHandSides) delete[] rhs;
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

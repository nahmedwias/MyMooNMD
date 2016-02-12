/**
 *  @brief Implement class Coupled_Time_CDR_2D declared in Coupled_Time_CDR_2D.h
 */

#include <Coupled_Time_CDR_2D.h>
#include <CoupledCDR_2D.h> //for SolvingStrategy type
#include <ReactionCoupling.h>
#include <Time_CD2D.h>


Coupled_Time_CDR_2D::Coupled_Time_CDR_2D(
    const TDomain& domain, const Example_CoupledCDR2D& example)
: example_(example)
{
  nEquations_ = example_.getNEquations();

  /* ********* The list of Time_CD2D problems.  ********************/
  cdProblems_.reserve(nEquations_); //reserve space

  // for each equation build a decoupled example and store it
  for(size_t index = 0; index<nEquations_;++index)
  {
    auto tcd_obj = std::make_shared<Time_CD2D>(domain,example_.getDecoupledExample(index));
    cdProblems_.push_back(tcd_obj);

  }
  /* ********* The list of coupled parts. ********************/
  //loop over all equations and construct the coupled parts objects
  for(size_t index = 0; index<nEquations_;++index){
    //Get a pointer to the assembling function.
    AssembleFctParam2D* currAssemblingFunction = example_.getAssemblingFct(index);
    //Get a pointer to the param (in-out) function.
    ParamFct* currParamFunction = example_.getParamFct(index);
    // Get a pointer to the rhs space (all spaces the same, ansatz=test (=rhs) )
    const TFESpace2D& currSpace = cdProblems_[index]->get_space();

    //strategy is a dummy here, as we use only one
    auto strategy = CoupledCDR_2D::SolvingStrategy::linear_decoupled;

    //Construct a ReactionCoupling object and attach it to the list.
    auto coupled_part =
        std::make_shared<ReactionCoupling>(
            strategy, currAssemblingFunction,
            currParamFunction, nEquations_, currSpace);

    coupledParts_.push_back(coupled_part);
  }

}

void Coupled_Time_CDR_2D::assemble_initial_time()
{
  for (auto tcd : cdProblems_)
  {//ask all decoupled problems politely to assemble their state at initial time
    tcd->assemble_initial_time();
  }
}

void Coupled_Time_CDR_2D::assemble_uncoupled_part()
{
  for (auto tcd : cdProblems_)
  {
    //ask all decoupled problems politely to assemble their state at current time
    tcd->assemble();
  }
}

void Coupled_Time_CDR_2D::couple_and_solve()
{
  Output::print("Call to Couple and solve.");
  // ///////////////// COPY & PASTE FROM STATIONARY ///////////////////////////
    // Put up an array of pointers to the solutions of previous iteration
    TFEFunction2D** previousSolutions = new TFEFunction2D*[nEquations_];

    // Store the original right hand sides of the CDR Equations without coupling.
    std::vector<BlockVector> originalRightHandSides;
    for (size_t equation = 0; equation<nEquations_;++equation){
      BlockVector vector(cdProblems_[equation]->get_rhs());
      originalRightHandSides.push_back(vector);
    }

    // while(Abbruchbedingung nicht erfuellt) - beginne einfach mit einer festen Anzahl Iterationen.
    for (size_t steps = 0; steps<5;steps++){
      //Fill pointers to available solutions into previousSolutions array
      for(size_t i =0;i<nEquations_;++i){
        previousSolutions[i]=&cdProblems_[i]->get_function();
      }
      //loop over the equations
      for (size_t equation = 0; equation<nEquations_;++equation){
        Output::print("Step ", steps, " Equation ", equation);

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

    //descale the stiffness matrices of the problems, which also updates old_Au
    for (auto cd : cdProblems_){
      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      double theta_1 = TDatabase::TimeDB->THETA1;
      cd->descale_stiffness(tau, theta_1);
    }

    delete[] previousSolutions; //just delete the pointers array

    // ///////////////// END COPY & PASTE FROM STATIONARY ///////////////////////////
}

void Coupled_Time_CDR_2D::output(int& image){
  //Let the work be done by the Time_CD objects for now.
  for (int equation = 0; equation<nEquations_;++equation){
//  TODO give meaningful names to the different pictures
    cdProblems_[equation]->output(equation, image);
  }
}

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

void Coupled_Time_CDR_2D::output(int& image){
  //Let the work be done by the Time_CD objects for now.
  for (size_t equation = 0; equation<nEquations_;++equation){
    cdProblems_[equation]->output(equation, image);
  }
}

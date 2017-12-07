/**
 *  @brief Implement class Coupled_Time_CDR_2D declared in Coupled_Time_CDR_2D.h
 */

#include <Coupled_Time_CDR_2D.h>
#include <CoupledCDR_2D.h> //for SolvingStrategy type
#include <FEFunctionInterpolator.h>
#include <LinAlg.h> //in the future: use an actual BLAS instead!
#include <ParameterDatabase.h>
#include <ReactionCoupling.h>
#include <Time_CD2D.h>

#include <algorithm>


ParameterDatabase get_default_Coupled_Time_CDR_2D_parameters()
{
  ParameterDatabase db("TCDRE System Database");

  db.add("tcdre_system_solve_maxit", (size_t) 5,
         "Maximum n of iterations in solve loop.", (size_t) 0, (size_t) 20);

  db.add("tcdre_system_solve_epsilon", 1e-10,
         "Target residual to break solver loop. "
         "Maximum over all equations.", 0.0, 1.0);

  return db;
}


Coupled_Time_CDR_2D::Coupled_Time_CDR_2D(
    const TDomain& domain, const ParameterDatabase& input_db,
    const Example_TimeCoupledCDR2D& example)
: example_(example), db_(get_default_Coupled_Time_CDR_2D_parameters())
{
  db_.merge(input_db, false);
  nEquations_ = example_.getNEquations();

  /* ********* The list of Time_CD2D problems.  ********************/
  cdProblems_.reserve(nEquations_); //reserve space

  // for each equation build a decoupled example and store it
  for(size_t index = 0; index<nEquations_;++index)
  {
    ParameterDatabase tcd2d_db = ParameterDatabase::parmoon_default_database();
    tcd2d_db.merge(input_db, true);
    //TODO this is awful about the new Database!
    tcd2d_db["output_basename"].set_range<std::string>(
        { tcd2d_db["output_basename"].get<std::string>()
      + std::string("_species_" + std::to_string(index) + std::string(".")),
        tcd2d_db["output_basename"].get<std::string>()});
    tcd2d_db["output_basename"] = tcd2d_db["output_basename"].get<std::string>()
        + std::string("_species_" + std::to_string(index) + std::string("."));

    auto tcd_obj = std::make_shared<Time_CD2D>(domain, tcd2d_db, example_.getDecoupledExample(index));
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
  // Put up an array of pointers to the solutions of previous iteration
  TFEFunction2D** previousSolutions = new TFEFunction2D*[nEquations_];

  // Store the original right hand sides of the CDR Equations without coupling.
  std::vector<BlockVector> originalRightHandSides;
  for (size_t equation = 0; equation<nEquations_;++equation){
    BlockVector vector(cdProblems_[equation]->get_rhs());
    originalRightHandSides.push_back(vector);
  }

  // A vector keeping track of the residuals before solving
  std::vector<double> resids( nEquations_ , 10e10 );

  for (size_t step = 1 ; true ; ++step){

    //Fill pointers to available solutions into previousSolutions array
    for(size_t i =0;i<nEquations_;++i){
      previousSolutions[i]=&cdProblems_[i]->get_function();
    }

    //loop over the equations
    for (size_t equation = 0; equation<nEquations_;++equation){

      // assemble the coupling term
      coupledParts_[equation]
                    ->assembleLinearDecoupled(previousSolutions);

      // add coupled rhs to rhs of the uncoupled equation
      cdProblems_[equation]->get_rhs().add_scaled(coupledParts_[equation]->getRightHandSide(),1);

      // before solving, fetch the residual
      resids[equation] = cdProblems_[equation]->get_discrete_residual();
      // solve equation with the new right hand side
      cdProblems_[equation]->solve();

      // Set back to uncoupled rhs.
      cdProblems_[equation]->get_rhs().copy(originalRightHandSides.at(equation).get_entries());
    }

    if(break_iteration(step, resids))
      break;

  }//endwhile bzw. endfor


  //descale the stiffness matrices of the problems, which also updates old_Au
  for (auto cd : cdProblems_)
  {
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    double theta_1 = TDatabase::TimeDB->THETA1;
    cd->descale_stiffness(tau, theta_1);
  }

  delete[] previousSolutions; //just delete the pointers array
}

void Coupled_Time_CDR_2D::assemble_initial_time(
    TFEFunction2D* velo1, TFEFunction2D* velo2)
{
  for (auto tcd : cdProblems_)
  {
    tcd->assemble_initial_time(velo1, velo2);
  }
}

void Coupled_Time_CDR_2D::assemble_uncoupled_part(
		TFEFunction2D* velo1, TFEFunction2D* velo2,
    std::vector<TFEFunction2D*> sources_and_sinks
    )
{
  if(sources_and_sinks.size() != nEquations_)
    ErrThrow("Number of given source-and-sink terms must match number of CDR equations!");
  int index=0;
  for (auto tcd : cdProblems_)
  {
    tcd->assemble(velo1, velo2, sources_and_sinks.at(index));
    ++index;
  }
}

std::vector<TFEFunction2D*> Coupled_Time_CDR_2D::get_fe_functions()
{
   std::vector<TFEFunction2D*> fe_fcts;
   for(auto eq : this->cdProblems_)
     fe_fcts.push_back(&eq->get_function());
   return fe_fcts;
}

std::vector<const TFEFunction2D*> Coupled_Time_CDR_2D::get_const_fe_functions() const
{
   std::vector<const TFEFunction2D*> fe_fcts;
   for(auto eq : this->cdProblems_)
     fe_fcts.push_back(&eq->get_function());
   return fe_fcts;
}

void Coupled_Time_CDR_2D::output(){
  //Let the work be done by the Time_CD objects for now.
  for (size_t equation = 0; equation<nEquations_;++equation){
    cdProblems_[equation]->output();
  }
}

/**
 * TODO
 * Es ist mir noch nicht klar, ob die Ermittlung des Residuals in Ordnung ist.
 * Bis jetzt wird in einer Gauss-Seidel-artigen Schleife geloest. Das Residuum
 * einer jeden gekoppelten Gleichung einzeln wird vor dem Loesen bestimmt, wobei
 * also die bis dahin bekannten anderen Loesungen eingehen. Ob die einzelnen
 * Residuen (d.h. ihr Maximum) klein genug sind, um die Loesungsschleife zu brechen,
 * wird am Ende einer vollen Loesungsschleife bestimmt. In die dort
 * beruecksichtigten Residuen gehen also nie die AKTUELLSTEN Loesungen ein -
 * nicht ueber die Kopplung in die rechte Seite, nicht wird die aktuelle Loesung
 * benutzt. Staendig zu updaten waere aber wahrscheinlich zu teuer.
 *
 */
bool Coupled_Time_CDR_2D::break_iteration(size_t step, std::vector<double> residuals)
{
  size_t max_it = db_["tcdre_system_solve_maxit"];
  double epsilon = db_["tcdre_system_solve_epsilon"];
  //find maximum residual of a single equation
  double max_res = *std::max_element(residuals.begin(), residuals.end());

  // Check whether target epsilon is hit.
  if(max_res < epsilon)
  {
    Output::stat("TCDRE SYSTEM SOLVE", "Target epsilon hit (", epsilon, ")");
    for(size_t p = 0 ; p < nEquations_; ++p)
      Output::dash("eq ", p ," discrete residual before solve: ", residuals[p]);
    return true;
  }

  // Check whether number of maximum iterations is hit.
  if(step >= max_it)
  {
    Output::stat("TCDRE SYSTEM SOLVE", "Number of maximum iterations hit (", step, ")");
    for(size_t p = 0 ; p < nEquations_; ++p)
      Output::dash("eq ", p ," discrete residual: ", residuals[p]);
    return true;
  }

  return false;

}














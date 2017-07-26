#include <Time_CD2D.h>
#include <Database.h>
#include <DirectSolver.h>
#include <MainUtilities.h>

#include <AlgebraicFluxCorrection.h>
#include <FEFunctionInterpolator.h>
#include <LocalAssembling2D.h>
#include <Assemble2D.h>
#include <LocalProjection.h>

#include <NSE2D_FixPo.h> //we need that include for a single "in-out function"
                         // (NSParamsVelo), maybe hard code the same thing here


/**************************************************************************** */
ParameterDatabase get_default_TCD2D_parameters()
{
  Output::print<5>("creating a default TCD2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default TCD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TCD2D parameter database");

  //TCD2D requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default afc database
  ParameterDatabase afc_db = AlgebraicFluxCorrection::default_afc_database();
  db.merge(afc_db, true);

  return db;
}

/**************************************************************************** */
Time_CD2D::System_per_grid::System_per_grid(const Example_TimeCD2D& example,
                                            TCollection& coll)
: fe_space(&coll, "time_cd2d space", "time_cd2d space", example.get_bc(0),
           TDatabase::ParamDB->ANSATZ_ORDER, nullptr),
           // TODO CB: Building the matrix here and rebuilding later is due to the
           // highly non-functional class TFEFunction2D, which does
           // neither provide default constructors nor working copy assignments
           stiff_matrix({&fe_space}),
           mass_matrix({&fe_space}),
           rhs(stiff_matrix, true),
           solution(stiff_matrix, false),
           old_Au(stiff_matrix, true),
           fe_function(&fe_space, "c", "a function subject to a time-dependent convection-diffusion equation",
                       solution.get_entries(), solution.length())
{
  stiff_matrix = BlockFEMatrix::CD2D(fe_space);
  mass_matrix = BlockFEMatrix::CD2D(fe_space);
}


/**************************************************************************** */

void Time_CD2D::System_per_grid::descale_stiff_matrix(double tau, double theta_1)
{
  if (tau==0 || theta_1 == 0)
  {
    ErrThrow("I will not divide by zero in descaling! tau=", tau, "theta_1=", theta_1);
  }
  //subtract the mass matrix (TODO smells like cancellation)...
  const FEMatrix& mass_block = *mass_matrix.get_blocks().at(0).get();
  //subtract the mass matrix...
  stiff_matrix.add_matrix_actives(mass_block, -1.0, {{0,0}}, {false});
  //...and descale the stiffness matrix...
  const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
  stiff_matrix.scale_blocks_actives(1./(tau*theta_1), cell_positions);
}

/**************************************************************************** */

void Time_CD2D::System_per_grid::update_old_Au()
{
  // the stiffness matrix must have been descaled (pure transport operator)
  stiff_matrix.apply(solution,old_Au);
}

/**************************************************************************** */
Time_CD2D::Time_CD2D(const TDomain& domain, const ParameterDatabase& param_db,
		int reference_id)
 : Time_CD2D(domain, param_db, Example_TimeCD2D(param_db), reference_id)
{
  
}

/**************************************************************************** */
Time_CD2D::Time_CD2D(const TDomain& domain, const ParameterDatabase& param_db,
		const Example_TimeCD2D& ex, int reference_id)
 : db(get_default_TCD2D_parameters()), solver(param_db), systems(), example(ex),
   errors(5, 0.0), timeDependentOutput(param_db), axisymmetric_(ex.is_axisymmetric())
{
  db.merge(param_db, false);
  set_parameters();
  // this is the L^inf(L^2(Omega)) error, initialize with some large number
  errors[4] = 1.e10; 
  
  bool usingMultigrid = solver.is_using_multigrid();
  if(!usingMultigrid)
  {
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
    // create finite element space and function, a matrix, rhs, and solution
    systems.emplace_back(example, *coll);
    
    TFESpace2D& space = systems.front().fe_space;
    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);
    Output::print<1>("N_Cells    : ", setw(12), coll->GetN_Cells());
    Output::print<1>("h(min,max) : ", setw(12), hmin, " ", setw(12), hmax);
    Output::print<1>("dof        : ", setw(12), space.GetN_DegreesOfFreedom());
    Output::print<1>("active dof : ", setw(12), space.GetN_ActiveDegrees());
    
    // old right hand side
    old_rhs.copy_structure(systems[0].rhs);
    
    // interpolate the initial data
    TFEFunction2D & fe_function = systems.front().fe_function;
    fe_function.Interpolate(example.get_initial_cond(0));
    // add the fe function to the output object. 
    timeDependentOutput.add_fe_function(&fe_function);
  }
  else
  {
    auto multigrid = solver.get_multigrid();
    if(multigrid->is_using_mdml())
      ErrThrow("mdml for TCD2D not yet implemented");
    
    size_t nMGLevel = multigrid->get_n_geometric_levels();
    int finest = domain.get_ref_level();
    int coarsest = finest-nMGLevel+1;
    
    std::list<BlockFEMatrix*> matrices;
    for(int grid_no = finest; grid_no >= coarsest; --grid_no)
    {
      TCollection *coll = domain.GetCollection(It_EQ, grid_no, -4711);
      systems.emplace_back(example, *coll);
      
      systems.front().fe_function.Interpolate(example.get_initial_cond(0));
      matrices.push_front(&systems.back().stiff_matrix);
    }
    multigrid->initialize(matrices);
  }
}

/**************************************************************************** */
void Time_CD2D::set_parameters()
{

  if(TDatabase::TimeDB->TIME_DISC == 0)
  {
    ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC 
          << " is not supported");
    throw("TIME_DISC: 0 is not supported. Chose 1 (backward Euler) or 2 (Crank-Nicolson).");
  }

  //set problem_type to Time_CD if not yet set
  if(!db["problem_type"].is(2))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 2;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to Time_CD."
          "It is now reset to the correct value for Time_CD (=2).");
      db["problem_type"] = 2;
    }
  }

  // an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }

  //////////////// Algebraic flux correction ////////////
  if(!db["algebraic_flux_correction"].is("none"))
  {//some kind of afc enabled
    if(!db["algebraic_flux_correction"].is("fem-fct-cn"))
    {
      db["algebraic_flux_correction"].set("fem-fct-cn");
      Output::print("Only kind of algebraic flux correction"
          " for TCD problems is Crank-Nicolson FEM-FCT (fem-fct-cn).");
    }

    if(solver.is_using_multigrid())
      ErrThrow("Multgrid and FEM-FCT are not yet enabled in Time_CD2D"); //TODO

    if(TDatabase::TimeDB->TIME_DISC !=2)
    {
      ErrThrow("Algebraic flux correction with FEM-FCT is "
          "implemented for Crank-Nicolson time stepping scheme only. "
          "Set TDatabase::TimeDB->TIME_DISC to 2.")
    }

    if(TDatabase::ParamDB->ANSATZ_ORDER != 1)
    {
      ErrThrow("Algebraic flux correction with FEM-FCT does only work for"
          "linear elements. Change ANSATZ_ORDER to 1!");
    }

    //make sure that galerkin discretization is used
    if (!db["space_discretization_type"].is("galerkin"))
    {//some other disctype than galerkin
      db["space_discretization_type"] = "galerkin";
      Output::warn<1>("Parameter 'space_discretization_type' changed to "
          " 'galerkin' because Algebraic Flux Correction is enabled.");
    }

    // when using afc, create system matrices as if all dofs were active
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
  }
}

/**************************************************************************** */
void Time_CD2D::assemble_initial_time(TFEFunction2D* velo1, TFEFunction2D* velo2)
{

  // chose the types for the local assembling objects
  LocalAssembling2D_type mass = LocalAssembling2D_type::TCD2D_Mass;
  LocalAssembling2D_type stiff_rhs = LocalAssembling2D_type::TCD2D;

  for(auto &s : this->systems)
  {
    // assemble mass matrix, stiffness matrix and rhs
    TFEFunction2D* fe_funct[1];
    fe_funct[0] = &s.fe_function;

    LocalAssembling2D la_mass(mass, fe_funct,
                              this->example.get_coeffs());
    LocalAssembling2D la_a_rhs(stiff_rhs, fe_funct,
                               this->example.get_coeffs());

    if (velo1 && velo2)
    { // modify the local assembling object for the stiffness matrix,
      // should a velocity field be given
      modify_and_call_assembling_routine(s, la_a_rhs, la_mass , true, velo1, velo2);
    }
    else
    {//no modifications necessary
      call_assembling_routine(s, la_a_rhs, la_mass , true);
    }

    // apply local projection stabilization method on stiffness matrix only!
    if(db["space_discretization_type"].is("local_projection")
        && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
    {
  	  if(axisymmetric_)
  		  ErrThrow("local_projection is untested for axisymmetric problems");

      if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
      {
        //fetch stiffness matrix as block
        std::vector<std::shared_ptr<FEMatrix>> stiff_blocks =
            s.stiff_matrix.get_blocks_uniquely();
        auto stiff_block = stiff_blocks.at(0).get();
        //call ultra local projection
        UltraLocalProjection((void *)stiff_block, false);
      }
      else
      {
        ErrThrow("LP_FULL_GRADIENT needs to be 1 to use LOCAL_PROJECTION");
      }

    }

    s.stiff_matrix.apply(s.solution,s.old_Au); //put up initial old_Au
  }

  System_per_grid& s = systems.front();
  old_rhs = s.rhs; //put up initial old_rhs
}

/**************************************************************************** */
void Time_CD2D::assemble(
	TFEFunction2D* velo1, TFEFunction2D* velo2,
    TFEFunction2D* sources_and_sinks)
{

  LocalAssembling2D_type stiff_rhs = LocalAssembling2D_type::TCD2D;  
  
  for(auto &s : this->systems)
  {

    TFEFunction2D * pointer_to_function = &s.fe_function;
    // create two local assembling object (second one will only be needed in SUPG case)
    LocalAssembling2D la_a_rhs(stiff_rhs, &pointer_to_function,
                               this->example.get_coeffs());

    if(db["space_discretization_type"].is("supg"))
    {
      ErrThrow("SUPG IS UNCHECKED IN THIS BRANCH (COUPLED CDR) AND NOT ADAPTED TO AXISYMMETRIC!");
      // In the SUPG case:
      // M = (u,v) + \tau (u,b.grad v)
      LocalAssembling2D_type mass_supg = LocalAssembling2D_type::TCD2D_Mass;;

      LocalAssembling2D la_m_supg(mass_supg, &pointer_to_function,
                               this->example.get_coeffs());

      if(velo1 && velo2 && sources_and_sinks)
      {//use both
        modify_and_call_assembling_routine(s, la_a_rhs, la_m_supg ,
                                           true, velo1, velo2, sources_and_sinks);
      }
      else
      {//or none
        //call assembling, including mass matrix (SUPG!)
        call_assembling_routine(s, la_a_rhs, la_m_supg , true);
      }
    }
    else //non SUPG
    {//call assembling, ignoring mass matrix part (third argument is not relevant)
      if(velo1 && velo2 && sources_and_sinks)
      {//use both
        modify_and_call_assembling_routine(s, la_a_rhs, la_a_rhs ,
                                           false, velo1, velo2, sources_and_sinks);
      }
      else
      {//or none: nothing external given, no modifications necessary
        call_assembling_routine(s, la_a_rhs, la_a_rhs , false);
      }
    }
  }
  
  // here the modifications due to time discretization begin
  if (!db["algebraic_flux_correction"].is("none") )
  {
    do_algebraic_flux_correction();
    return; // modifications due to time discretization are per-
            // formed inside the afc scheme, so step out here!
  }

  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  // preparing rhs 
  System_per_grid& s = this->systems.front();
  
  if(TDatabase::TimeDB->THETA4)
  {
    // scale by time step length and theta4 (only active dofs)
    s.rhs.scaleActive(tau*TDatabase::TimeDB->THETA4);    
    // add old right hand side scaled by time step length and theta3 (only 
    // active dofs)
    if(TDatabase::TimeDB->THETA3 != 0.)
      s.rhs.addScaledActive((this->old_rhs), tau*TDatabase::TimeDB->THETA3);	

    // save old right hand side (only if THETA3 != 0)
    if(TDatabase::TimeDB->THETA3)
    {
      this->old_rhs.addScaledActive(s.rhs, -1./(tau*TDatabase::TimeDB->THETA3));
      this->old_rhs.scaleActive(-TDatabase::TimeDB->THETA3/TDatabase::TimeDB->THETA4);
    }
  }
  else
  {
    if(TDatabase::TimeDB->TIME_DISC == 0)
    {
      ErrThrow("Forward Euler method is not supported. "
          "Choose TDatabase::TimeDB->TIME_DISC as 1 (bw Euler)"
          " or 2 (Crank-Nicoloson)");
    }
  }
  
  for(auto &s : this->systems)
  {
    // rhs += M*uold
    s.mass_matrix.apply_scaled_add_actives(s.solution, s.rhs, 1.0);
    // rhs -= tau*theta2*A_old*uold
    s.rhs.addScaledActive(s.old_Au, -tau*TDatabase::TimeDB->THETA2);
    
    // preparing the left hand side, i.e., the system matrix
    // stiffness matrix is scaled by tau*THETA1, after solving 
    // the matrix needs to be descaled if the coeffs does not depends
    // on time

    //scale  stiffness matrix...
    const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
    s.stiff_matrix.scale_blocks_actives(tau*TDatabase::TimeDB->THETA1, cell_positions);
    // ...and add the mass matrix
    const FEMatrix& mass_block = *s.mass_matrix.get_blocks().at(0).get();
    s.stiff_matrix.add_matrix_actives(mass_block, 1.0, {{0,0}}, {false});
  }  
  systems[0].solution.copy_nonactive(systems[0].rhs);
}

/**************************************************************************** */

void Time_CD2D::descale_stiffness(double tau, double theta_1)
{
    // restore stiffness matrix and store old_Au on all grids
    if(db["algebraic_flux_correction"].is("none"))
    {
      for(auto &s : this->systems)
      {
        s.descale_stiff_matrix(tau, theta_1);
        s.update_old_Au();
      }
    }
    else
    {
      //in AFC case the stiffness matrix is "ruined" by now -
      // TODO Keep this in mind - but do not print it all the time.
      // Output::print("AFC does not yet reset the stiffness"
      //  "matrix and old_Au correctly!");
    }
}

/**************************************************************************** */
void Time_CD2D::solve()
{
  double t = GetTime();
  System_per_grid& s = this->systems.front();
  
  this->solver.solve(s.stiff_matrix, s.rhs, s.solution);
  
  t = GetTime() - t;
  Output::print<1>("time for solving: ",  t);
  Output::print<2>("solution ", sqrt(Ddot(s.solution.length(),
                                          s.solution.get_entries(),
                                          s.solution.get_entries())) );
}

/**************************************************************************** */
void Time_CD2D::output()
{
  bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
  if(no_output)
    return;
  TFEFunction2D & fe_function = this->systems.front().fe_function;
  fe_function.PrintMinMax();
  
  if(db["output_compute_errors"])
  {
    double loc_e[5];
    TAuxParam2D aux; // NOTE: this is hjust a default constructed "aux" object
                     // TODO If an actual aux object is involved (meaning: another fe function
                     // enters the computation of coefficients, it must get here somehow!)

    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D* space = fe_function.GetFESpace2D();
    
    fe_function.GetErrors(this->example.get_exact(0), 3, AllDerivatives, 4,
                          SDFEMErrors, this->example.get_coeffs(), &aux, 1, 
                          &space, loc_e);
    
    Output::print<1>("time: ", TDatabase::TimeDB->CURRENTTIME);
    Output::print<1>("  L2: ", loc_e[0]);
    Output::print<1>("  H1-semi: ", loc_e[1]);
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    errors[0] += (loc_e[0]*loc_e[0] + errors[1])*tau*0.5;
    errors[1] = loc_e[0]*loc_e[0];
    Output::print<1>("  L2(0,T;L2) ", sqrt(errors[0]));
    errors[2] += (loc_e[1]*loc_e[1] + errors[3])*tau*0.5;
    errors[3] = loc_e[1]*loc_e[1];
    Output::print<1>("  L2(0,T;H1) ", sqrt(errors[2]));
    
    if(errors[4] < loc_e[0])
      errors[4] = loc_e[0];
    Output::print<1>("  Linfty(0,T;L2) ", errors[4]);
  }

  // write output
  timeDependentOutput.write(TDatabase::TimeDB->CURRENTTIME);

}

/* *************************************************************************** */
double Time_CD2D::get_discrete_residual() const
{
  const System_per_grid& s = systems.front();

  BlockVector defect = s.rhs;
  s.stiff_matrix.apply_scaled_add(s.solution, defect, -1.0);

  return defect.norm();
}

/* *************************************************************************** */
std::array< double, int(3) > Time_CD2D::get_errors() const
{
  std::array<double, int(3)> error_at_time_points;
  error_at_time_points.at(0) = sqrt(errors.at(0));
  error_at_time_points.at(1) = sqrt(errors.at(2));
  error_at_time_points.at(2) = sqrt(errors.at(4));
  return error_at_time_points;

    
    
}

/**************************************************************************** */

void Time_CD2D::do_algebraic_flux_correction()
{
  //determine which kind of afc to use
  if(db["algebraic_flux_correction"].is("default") ||
      db["algebraic_flux_correction"].is("fem-fct-cn"))
  {
    //TODO implement for multigrid!
    System_per_grid& s = systems.front();

    //get references to the relevant objects

    const TFESpace2D& feSpace = s.fe_space; //space
    const FEMatrix& mass = *s.mass_matrix.get_blocks().at(0).get(); // mass
    FEMatrix& stiff = *s.stiff_matrix.get_blocks_uniquely().at(0).get(); //stiffness
    //vector entry arrays
    const std::vector<double>& solEntries = s.solution.get_entries_vector();
    std::vector<double>& rhsEntries = s.rhs.get_entries_vector();
    std::vector<double>& oldRhsEntries = old_rhs.get_entries_vector();

    // fill a vector "neumannToDirichlet" with those rows that got
    // internally treated as Neumann although they are Dirichlet
    int firstDiriDof = feSpace.GetActiveBound();
    int nDiri = feSpace.GetN_Dirichlet();
    std::vector<int> neumToDiri(nDiri, 0);
    std::iota(std::begin(neumToDiri), std::end(neumToDiri), firstDiriDof);

    //determine prelimiter from Database
    AlgebraicFluxCorrection::Prelimiter prelim;
    switch((int) db["afc_prelimiter"])
    {
      case 1:
        prelim = AlgebraicFluxCorrection::Prelimiter::MIN_MOD;
        Output::print<2>("FEM-FCT: MinMod prelimiter chosen.");
        break;
      case 2:
        prelim = AlgebraicFluxCorrection::Prelimiter::GRAD_DIRECTION;
        Output::print<2>("FEM-FCT: Gradient direcetion prelimiter chosen.");
        break;
      case 3:
        prelim = AlgebraicFluxCorrection::Prelimiter::BOTH;
        Output::print<2>("FEM-FCT: Double prelimiting (MinMod & Gradient) chosen.");
        break;
      default:
        prelim = AlgebraicFluxCorrection::Prelimiter::NONE;
        Output::print<2>("FEM-FCT: No prelimiting chosen.");
    }

    //get current timestep length from Database
    double delta_t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

    // apply C-N FEM-FCT
    AlgebraicFluxCorrection::crank_nicolson_fct(
        mass, stiff,
        solEntries,
        rhsEntries, oldRhsEntries,
        delta_t,
        neumToDiri,
        prelim );

    //...and finally correct the entries in the Dirichlet rows
    AlgebraicFluxCorrection::correct_dirichlet_rows(stiff);
    //...and in the right hand side, too
    s.rhs.copy_nonactive(old_rhs);
  }
  else
  {
    ErrThrow("The chosen algebraic flux correction scheme is unknown "
        "to class Time_CD2D.");
  }
}

void Time_CD2D::call_assembling_routine(
    System_per_grid& s,
    LocalAssembling2D& la_stiff, LocalAssembling2D& la_mass,
    bool assemble_both)
{

  // Assemble mass matrix, stiffness matrix and rhs
  //...variables which are the same for both
  std::vector<const TFESpace2D* > fe_spaces = {&s.fe_space};
  int n_fe_spaces = 1;
  if((!axisymmetric_ && la_stiff.GetN_Parameters() == 2) || (axisymmetric_ && la_stiff.GetN_Parameters() == 3))
  {// this is an awful hack - yet it is my only sign, that
   // modify_and_call-assembling_routine was called WITHOUT source-and-sink term
    n_fe_spaces = 2;
    fe_spaces = {&s.fe_space, la_stiff.get_fe_function(1)->GetFESpace2D()};
  }
  else if((!axisymmetric_ && la_stiff.GetN_Parameters() == 3) || (axisymmetric_ && la_stiff.GetN_Parameters() == 4))
  {// this is an awful hack - yet it is my sign, that
   // modify_and_call-assembling_routine was called WITH source-and-sink term
    n_fe_spaces = 3;
    fe_spaces = {&s.fe_space,
                 la_stiff.get_fe_function(1)->GetFESpace2D(),
                 la_stiff.get_fe_function(3)->GetFESpace2D()};
  }

  // I checked: we need as many boundary_conditions as we have right hand sides (1),
  // those are the boundary conditions of the 0th fe space (the ansatz space)
  BoundCondFunct2D * boundary_conditions = fe_spaces[0]->GetBoundCondition();
  int N_Matrices = 1;
  double * rhs_entries = s.rhs.get_entries();

  BoundValueFunct2D * non_const_bound_value[1] {example.get_bd()[0]};

  //fetch stiffness matrix as block
  std::vector<std::shared_ptr<FEMatrix>> stiff_blocks = s.stiff_matrix.get_blocks_uniquely();
  TSquareMatrix2D * stiff_block[1]{reinterpret_cast<TSquareMatrix2D*>(stiff_blocks.at(0).get())};

  // Do the Assembling!

  //...stiffness matrix is second
  // reset right hand side and matrix to zero
  s.rhs.reset();
  stiff_block[0]->reset();
  Assemble2D(n_fe_spaces, &fe_spaces.at(0), N_Matrices, stiff_block, 0, NULL, 1, &rhs_entries,
             &fe_spaces.at(0), &boundary_conditions, non_const_bound_value, la_stiff);

  // If assemble_both is true, we also (re)assemble the mass matrix.
  if (assemble_both)
  {
    fe_spaces = {&s.fe_space};
    n_fe_spaces = 1;

    std::vector<std::shared_ptr<FEMatrix>> mass_blocks = s.mass_matrix.get_blocks_uniquely();
    TSquareMatrix2D * mass_block[1]{reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get())};

    mass_block[0]->reset();
    Assemble2D(n_fe_spaces, &fe_spaces.at(0), N_Matrices, mass_block, 0, NULL, 0, NULL,
               NULL, &boundary_conditions, non_const_bound_value, la_mass);
  }
}

// Used as a "ParamFunction" (a MooNMD specific oddity
// of the assembling process) in the following.
// ...shears away the first two parameters (x,y) and passes on the two next parameters
void TwoFEParametersFunction(double *in, double *out)
{
  out[0] = in[2];
  out[1] = in[3];
}
// ...shears away the first two parameters (x,y) and passes on the three next parameters
void ThreeFEParametersFunction(double *in, double *out)
{
  out[0] = in[2];
  out[1] = in[3];
  out[2] = in[4];
}
void AxisymmetricTwoFEParametersFunction(double *in, double *out)
{
	out[0] = in[1];
	out[1] = in[2];
	out[2] = in[3];
}
void AxisymmetricThreeFEParametersFunction(double *in, double *out)
{
	out[0] = in[1];
	out[1] = in[2];
	out[2] = in[3];
	out[3] = in[4];
}
void RadiusOnlyParametersFunction(double *in, double *out)
{
	out[0] = in[1];
}

void AxisymmetricLocalMatrixM(double Mult, double *coeff, double *param,
                           double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *MatrixRow;
  double ansatz00;
  double test00;
  double *Orig0;
  int i,j, N_;
  double r;

  r = param[0]; //y=r coordinate

  Matrix = LocMatrices[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    test00 = Orig0[i];

    for(j=0;j<N_;j++)
    {
      ansatz00 = Orig0[j];

      MatrixRow[j] += r * Mult*ansatz00*test00;
    } // endfor j
  } // endfor i
}

void AxisymmetricLocalMatrixARhs(double Mult, double *coeff, double *param,
                            double hK,
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, *Rhs, val, *MatrixRowA;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double r, eps, uz, ur, reac, f;// h;

  MatrixA = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  eps = coeff[0]; // eps
  uz = coeff[1]; // uz
  ur = coeff[2]; // ur
  reac = coeff[3]; // reaction coefficient
  f = coeff[4]; // f
  r  = coeff[5]; // y=r coordinate

  for(i=0;i<N_;i++)
  {
    MatrixRowA = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs[i] += Mult*test00*f*r;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = r * eps*(test10*ansatz10+test01*ansatz01);
      val += r * (uz*ansatz10+ur*ansatz01)*test00;
      val += r * reac*ansatz00*test00;

      MatrixRowA[j] += Mult * val;

    } // endfor j
  } // endfor i
}


void Time_CD2D::modify_and_call_assembling_routine(
    System_per_grid& s,
    LocalAssembling2D& la_stiff, LocalAssembling2D& la_mass,
    bool assemble_both,
	TFEFunction2D* velo1, TFEFunction2D* velo2,
    TFEFunction2D* sources_and_sinks)
{

  // NOTE: velo1, velo2 and sources_and_sinks must be defined on the right grid
  // - but currently there is no way to check that...

  // set up the input...
  std::vector<int> beginParameter = {0};

  TFEFunction2D* fe_funct[4]; //fill up the new fe function array (4th entry is optional, see below)
  fe_funct[0] = &s.fe_function;
  fe_funct[1] = velo1;
  fe_funct[2] = velo2;

  std::vector<int> feValueFctIndex = {1,2}; // to produce first fe value use fe function 1,
                                             // for second fe value use function 2
  std::vector<MultiIndex2D> feValueMultiIndex = {D00, D00}; // to produce first fe value use 0th derivative,
                                                            // for second fe value as well
  int N_parameters = 2; // two parameters (AFTER application of parameterFct...)
  int N_feValues = 2;   //..both of them stem from the evaluation of fe fcts
  int N_paramFct = 1;   // dealing with them is performed by 1 ParamFct

  // chose the parameter function ("in-out function") which shears away
  // the first to "in" values (x,y) and passes only u_x and u_y
  std::vector<ParamFct*> parameterFct = {TwoFEParametersFunction};

  if(axisymmetric_)
  {//parameter handling is different, y=r has to be passed on
	  N_parameters = 3;
	  parameterFct = {AxisymmetricTwoFEParametersFunction};
	  //this is Galerkin assembling for axisymmetric problem
	  la_stiff.setAssembleFctParam2D(AxisymmetricLocalMatrixARhs);
	  if(la_mass.get_type() == LocalAssembling2D_type::TCD2D_Mass) //yes, it's an actual mass matrix assembling method
	  {
		  la_mass.setAssembleFctParam2D(AxisymmetricLocalMatrixM);
		  la_mass.setBeginParameter({0});
		  la_mass.setN_Parameters(1);
		  la_mass.setN_FeValues(0);
		  la_mass.setN_ParamFct(1);
		  std::vector<ParamFct*> parameterFct = {RadiusOnlyParametersFunction};
		  la_mass.setParameterFct(parameterFct);
	  }
  }

  if(sources_and_sinks) // rhs source and sink terms are given
  {
    // NOTE: sources and sinks must be defined on the right grid
    // - but currently there is no way to check that...

    fe_funct[3] = sources_and_sinks;

    feValueFctIndex = {1,2,3}; // to produce first fe value use fe function 1,
                               // for second fe value use function 2,
                               // for third fe value use function 3

    feValueMultiIndex = {D00,D00,D00}; // to produce first fe value use 0th derivative,
                                       // for second and third fe value as well
    N_parameters = 3; // three parameters (AFTER application of parameterFct...)
    N_feValues = 3;   // ...all three of them stem from the evaluation of fe fcts
    N_paramFct = 1;   // dealing with them is performed by 1 ParamFct
    // chose the parameter function ("in-out function") which shears away
    // the first to "in" values (x,y) and passes u_x, u_y and f
    parameterFct = {ThreeFEParametersFunction};

    if(axisymmetric_)
    {//parameter handling is different, y=r has to be passed on
  	  N_parameters = 4;
  	  parameterFct = {AxisymmetricThreeFEParametersFunction};
    }

    //THIS IS DOUBLE CODE, but this stuff must be performed before
    // interpolated_sources_and_sinks and entries_source_and_sinks
    // go out of scope...
    // ...and call the corresponding setters
    la_stiff.setBeginParameter(beginParameter);
    la_stiff.setFeFunctions2D(fe_funct); //reset - now velo comp included
    la_stiff.setFeValueFctIndex(feValueFctIndex);
    la_stiff.setFeValueMultiIndex(feValueMultiIndex);
    la_stiff.setN_Parameters(N_parameters);
    la_stiff.setN_FeValues(N_feValues);
    la_stiff.setN_ParamFct(N_paramFct);
    la_stiff.setParameterFct(parameterFct);
    //...I expect that to do the trick.

    // step 4 - the assembling must be done before the velo functions
    // run out of scope
    call_assembling_routine(s, la_stiff, la_mass , assemble_both);
  }
  else
  {

  // ...and call the corresponding setters
  la_stiff.setBeginParameter(beginParameter);
  la_stiff.setFeFunctions2D(fe_funct); //reset - now velo comp included
  la_stiff.setFeValueFctIndex(feValueFctIndex);
  la_stiff.setFeValueMultiIndex(feValueMultiIndex);
  la_stiff.setN_Parameters(N_parameters);
  la_stiff.setN_FeValues(N_feValues);
  la_stiff.setN_ParamFct(N_paramFct);
  la_stiff.setParameterFct(parameterFct);
  //...I expect that to do the trick.

  // step 4 - the assembling must be done before the velo functions
  // run out of scope
  call_assembling_routine(s, la_stiff, la_mass , assemble_both);
  }
}

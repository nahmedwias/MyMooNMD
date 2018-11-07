#include <Time_CD2D.h>
#include <Database.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <AlgebraicFluxCorrection.h>
#include <Assemble2D.h>
#include <LocalProjection.h>
#include <AuxParam2D.h>
#include <TimeDiscretizations.h>


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
  
  // a default local assembling database
  db.merge(LocalAssembling2D::default_local_assembling_database(), true);
  
  db.merge(TimeDiscretization::default_TimeDiscretization_database(), true);

  return db;
}

/**************************************************************************** */
Time_CD2D::System_per_grid::System_per_grid(const Example_TimeCD2D& example,
                                            TCollection& coll)
: fe_space(&coll, "space", "time_cd2d space", example.get_bc(0),
           TDatabase::ParamDB->ANSATZ_ORDER)
{
  stiff_matrix = BlockFEMatrix::CD2D(fe_space);
  mass_matrix = BlockFEMatrix::CD2D(fe_space);

  rhs = BlockVector(stiff_matrix, true);
  solution = BlockVector(stiff_matrix, false);
  old_Au = BlockVector(stiff_matrix, true);
  fe_function = TFEFunction2D(&fe_space, "c", "c",
                              solution.get_entries(), solution.length());
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
   errors(5, 0.0), timeDependentOutput(param_db)
{
  db.merge(param_db, false);
  this->set_parameters();
  // this is the L^inf(L^2(Omega)) error, initialize with some large number
  errors[4] = 1.e10; 
  
  bool usingMultigrid = this->solver.is_using_multigrid();
  if(!usingMultigrid)
  {
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
    // create finite element space and function, a matrix, rhs, and solution
    this->systems.emplace_back(this->example, *coll);
    
    TFESpace2D& space = this->systems.front().fe_space;
    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);
    Output::print<1>("N_Cells    : ", setw(12), coll->GetN_Cells());
    Output::print<1>("h(min,max) : ", setw(12), hmin, " ", setw(12), hmax);
    Output::print<1>("dof        : ", setw(12), space.GetN_DegreesOfFreedom());
    Output::print<1>("active dof : ", setw(12), space.GetN_ActiveDegrees());
    
    // old right hand side
    old_rhs.copy_structure(this->systems[0].rhs);
    
    // interpolate the initial data
    TFEFunction2D & fe_function = this->systems.front().fe_function;
    fe_function.Interpolate(example.get_initial_cond(0));
    // add the fe function to the output object. 
    timeDependentOutput.add_fe_function(&fe_function);
  }
  else
  {
    auto multigrid = this->solver.get_multigrid();
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
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1; //FIXME THIS IS DANGEROUS, does not go well with BlockFEMatrix!
  }
}

/**************************************************************************** */
void Time_CD2D::assemble_initial_time()
{
  LocalAssembling_type mass = LocalAssembling_type::TCD2D_Mass;
  LocalAssembling_type stiff_rhs = LocalAssembling_type::TCD2D;

  for(auto &s : this->systems)
  {
    // assemble mass matrix, stiffness matrix and rhs
    TFEFunction2D* fe_funct[1] = {&s.fe_function}; //wrap up as '**'
    LocalAssembling2D la_mass(this->db, mass, fe_funct,
                              this->example.get_coeffs());
    LocalAssembling2D la_a_rhs(this->db, stiff_rhs, fe_funct,
                               this->example.get_coeffs());

    call_assembling_routine(s, la_a_rhs, la_mass , true);

    // apply local projection stabilization method on stiffness matrix only!
    if(db["space_discretization_type"].is("local_projection")
        && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
    {
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
void Time_CD2D::assemble()
{
  LocalAssembling_type stiff_rhs = LocalAssembling_type::TCD2D;  
  
  for(auto &s : this->systems)
  {

    TFEFunction2D * pointer_to_function = &s.fe_function;
    // create two local assembling object (second one will only be needed in SUPG case)
    LocalAssembling2D la_a_rhs(this->db, stiff_rhs, &pointer_to_function,
                               this->example.get_coeffs());

    if(db["space_discretization_type"].is("supg"))
    {
      // In the SUPG case:
      // M = (u,v) + \tau (u,b.grad v)
      LocalAssembling_type mass_supg = LocalAssembling_type::TCD2D_Mass;;

      LocalAssembling2D la_m_supg(this->db, mass_supg, &pointer_to_function,
                                  this->example.get_coeffs());

      //call assembling, including mass matrix (SUPG!)
      call_assembling_routine(s, la_a_rhs, la_m_supg , true);
    }
    else
    {//call assembling, ignoring mass matrix part (third argument is not relevant)
      call_assembling_routine(s, la_a_rhs, la_a_rhs , false);
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
  //this->systems[0].rhs.copy_nonactive(this->systems[0].solution);
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
      Output::print("AFC does not yet reset the stiffness"
          "matrix and old_Au correctly!");
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
    TAuxParam2D aux;
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

/**************************************************************************** */
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
  const TFESpace2D * fe_space = &s.fe_space;
  BoundCondFunct2D * boundary_conditions = fe_space->get_boundary_condition();
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
  Assemble2D(1, &fe_space, N_Matrices, stiff_block, 0, nullptr, 1, &rhs_entries,
             &fe_space, &boundary_conditions, non_const_bound_value, la_stiff);

  // If assemble_both is true, we also (re)assemble the mass matrix.
  if (assemble_both)
  {
    std::vector<std::shared_ptr<FEMatrix>> mass_blocks = s.mass_matrix.get_blocks_uniquely();
    TSquareMatrix2D * mass_block[1]{reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get())};

    mass_block[0]->reset();
    Assemble2D(1, &fe_space, N_Matrices, mass_block, 0, nullptr, 0, nullptr,
               nullptr, &boundary_conditions, non_const_bound_value, la_mass);
  }
}

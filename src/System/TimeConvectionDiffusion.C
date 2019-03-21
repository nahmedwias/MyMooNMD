#include "TimeConvectionDiffusion.h"
#include "AlgebraicFluxCorrection.h"
#include "Database.h"
#include "ConvDiff.h"
#include "Domain.h"
#include "Multigrid.h"
#ifdef __2D__
 #include "Assemble2D.h"
 #include "SquareMatrix2D.h"
 #include "AuxParam2D.h"
#else
 #include "Assemble3D.h"
 #include "SquareMatrix3D.h"
 #include "AuxParam3D.h"
#endif
#include "LocalProjection.h"
#ifdef _MPI
 #include "ParFECommunicator3D.h"
#endif


/* ************************************************************************* */
template <int d>
ParameterDatabase TimeConvectionDiffusion<d>::default_tcd_database()
{
  Output::print<5>("creating a default TCD2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default TCD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TimeConvectionDiffusion parameter database");

  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(AlgebraicFluxCorrection::default_afc_database());
  db.merge(LocalAssembling<d>::default_local_assembling_database(), true);
  db.merge(TimeDiscretization::default_TimeDiscretization_database(), true);
  return db;
}

/* ************************************************************************* */
template <int d>
TimeConvectionDiffusion<d>::SystemPerGrid::SystemPerGrid(
  const Example_TimeCD& example, TCollection& coll, int ansatz_order)
: fe_space(new FESpace(&coll, "space", "time_cd space", example.get_bc(0),
                       ansatz_order))
{
#ifdef _MPI
  fe_space->get_communicator().print_info();
#endif
#ifdef __3D__
  stiffness_matrix = BlockFEMatrix::CD3D(fe_space);
  mass_matrix = BlockFEMatrix::CD3D(fe_space);
#else
  stiffness_matrix = BlockFEMatrix::CD2D(fe_space);
  mass_matrix = BlockFEMatrix::CD2D(fe_space);
#endif

  rhs = BlockVector(stiffness_matrix, true);
  solution = BlockVector(stiffness_matrix, false);
  fe_function = FEFunction(fe_space, "c", "c", solution.get_entries(),
                           solution.length());
  solution_m1 = BlockVector(stiffness_matrix, false);
  u_m1 = FEFunction(fe_space,"um1","um1", solution_m1.get_entries(),
                    solution_m1.length());
  solution_m2 = BlockVector(stiffness_matrix, false);
  u_m2 = FEFunction(fe_space,"um2","um2", solution_m2.get_entries(),
                    solution_m2.length());
}

/* ************************************************************************* */
template <int d>
TimeConvectionDiffusion<d>::TimeConvectionDiffusion(
  const TDomain& domain, const ParameterDatabase& param_db)
 : TimeConvectionDiffusion<d>(domain, param_db, Example_TimeCD(param_db))
{
}

/* ************************************************************************* */
template <int d>
TimeConvectionDiffusion<d>::TimeConvectionDiffusion(
  const TDomain& domain, const ParameterDatabase &param_db,
  const Example_TimeCD& ex)
 : db(default_tcd_database()), solver(param_db), systems(), example(ex),
   old_rhs(), errors({}), outputWriter(param_db),
   time_stepping_scheme(param_db), rhs_from_time_disc()
{
  this->db.merge(param_db, false); // update this database with given values
  this->check_and_set_parameters();
  
  bool usingMultigrid = this->solver.is_using_multigrid();
  auto collections = domain.get_grid_collections();
  int ansatz_order = TDatabase::ParamDB->ANSATZ_ORDER;
  if(!usingMultigrid)
  {
    // the given collection for particular cell
    TCollection& cellCollection = *collections.front();
    systems.emplace_back(example, cellCollection, ansatz_order);
    // initial condition on the solution
    this->systems.front().fe_function.Interpolate(example.get_initial_cond(0));
  }
  else
  {
    auto multigrid = this->solver.get_multigrid();
    size_t nMgLevels = multigrid->get_n_geometric_levels();
    size_t n_grids = collections.size();
    if(n_grids < nMgLevels)
    {
      ErrThrow("Multigrid: expected ", nMgLevels, " collections, ", n_grids,
               " provided.");
    }
    // remove not needed coarser grid from list of collections
    for(size_t i = nMgLevels; i < n_grids; ++i)
    {
      collections.pop_back();
    }
    std::list<BlockFEMatrix*> matrices;
    // construct all SystemPerGrid and store them
    for(auto it : collections)
    {
      systems.emplace_back(example, *it, ansatz_order);
      systems.front().fe_function.Interpolate(example.get_initial_cond(0));
      matrices.push_front(&systems.back().stiffness_matrix);
    }
    multigrid->initialize(matrices);
  }// multigrid case
  // print useful information
  this->output_problem_size_info();
  outputWriter.add_fe_function(&this->systems.front().fe_function);
  this->rhs_from_time_disc.copy_structure(this->systems.front().rhs);
  errors.fill(0.); // initialize the array
   // initialize L_inf error to some negative number
  errors[TimeConvectionDiffusion::n_errors-1] = -1.;
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::output_problem_size_info() const
{
  // print some useful information
  const FESpace& space = *this->systems.front().fe_space;
  TCollection *coll = space.GetCollection();
#ifndef _MPI
  double hMin, hMax;
  coll->GetHminHmax(&hMin, &hMax);
  int n_cells = coll->GetN_Cells();
  int n_dof = space.GetN_DegreesOfFreedom();
#else // _MPI
  int root = 0; // root process number
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  int n_local_master_cells = coll->GetN_OwnCells();
  int n_cells;
  MPI_Reduce(&n_local_master_cells, &n_cells, 1, MPI_DOUBLE, MPI_SUM, root,
             MPI_COMM_WORLD);
  
  double local_hmin, local_hmax;
  coll->GetHminHmax(&local_hmin, &local_hmax);
  double hMin, hMax;
  MPI_Reduce(&local_hmin, &hMin, 1, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
  MPI_Reduce(&local_hmax, &hMax, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
  
  auto par_comm = systems.front().fe_space->get_communicator();
  int n_dof  = par_comm.get_n_global_dof();
  if(my_rank == root)
#endif
  {
    Output::stat("TimeConvectionDiffusion",
                 "Mesh data and problem size on finest grid");
    Output::dash("n cells     : ", setw(13), n_cells);
    Output::dash("h(min, max) : ", setw(13), hMin, " ", setw(13), hMax);
    Output::dash("n dofs      : ", setw(13), n_dof);
#ifndef _MPI
    Output::dash("n dof active: ", setw(13), space.GetActiveBound());
#endif
  }
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::check_and_set_parameters()
{
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

    if(solver.is_using_multigrid())
      ErrThrow("Multgrid and FEM-FCT are not yet enabled in "
               "TimeConvectionDiffusion"); /// @todo

    if(!db["algebraic_flux_correction"].is("fem-fct-cn"))
    {
      db["algebraic_flux_correction"].set("fem-fct-cn");
      Output::warn("TimeConvectionDiffusion::check_and_set_parameters()",
                   "Parameter 'algebraic_flux_correction' changed to "
                   "'fem-fct-cn', which is the only implemented kind of "
                   "algebraic flux correction for TCD problems.");
    }

    if(TDatabase::TimeDB->TIME_DISC != 2) /// @todo are we still using this?
    {
      ErrThrow("Algebraic flux correction with FEM-FCT is implemented for "
               "Crank-Nicolson time stepping scheme only. Set "
               "TDatabase::TimeDB->TIME_DISC to 2.")
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
      Output::warn("TimeConvectionDiffusion::check_and_set_parameters()",
                   "Parameter 'space_discretization_type' changed to 'galerkin'"
                   " because Algebraic Flux Correction is enabled.");
    }

    // when using afc, create system matrices as if all dofs were active
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1; //FIXME THIS IS DANGEROUS, does not go well with BlockFEMatrix!
  }
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::assemble_initial_time()
{
  LocalAssembling_type allMatrices = LocalAssembling_type::TCDStiffMassRhs;
  for(auto &s : this->systems)
  {
    FEFunction *feFunction = &s.fe_function;
    LocalAssembling<d> la(this->db, allMatrices, &feFunction,
                          example.get_coeffs());
    // Assemble stiffness, mass matrices and the rhs
    call_assembling_routine(s, la, true);
        // apply local projection stabilization method on stiffness matrix only!
    if(db["space_discretization_type"].is("local_projection")
        && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
    {
      if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
      {
        if(d == 3)
          ErrThrow("local_projection for TimeConvectionDiffusion<3> is not "
                   "implemented"); /// @todo
#ifdef __2D__
        //fetch stiffness matrix as block
        std::vector<std::shared_ptr<FEMatrix>> stiff_blocks =
            s.stiffness_matrix.get_blocks_uniquely();
        auto stiff_block = stiff_blocks.at(0).get();
        //call ultra local projection
        UltraLocalProjection((void *)stiff_block, false);
#endif // 2D
      }
      else
      {
        ErrThrow("LP_FULL_GRADIENT needs to be 1 to use LOCAL_PROJECTION");
      }

    }
  }
  SystemPerGrid& s = this->systems.front();
  old_rhs = s.rhs;
  s.solution_m1 = s.solution;
  s.solution_m2 = s.solution_m1;  
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::assemble()
{
  // In the case of SUPG: local assemble function itself take care of the 
  // number of matrices. One have to assemble also the weighted mass matrix 
  // which comes from the time discretization of the SUPG method. We have 
  // to assemble the Mass matrix also for each time step due to the convection 
  // field which might also depend on time as well.
  LocalAssembling_type la_type = LocalAssembling_type::TCDStiffRhs;
  bool assemble_both = false;
  for(auto &s : this->systems)
  {
    // call assembling routine 
    if(db["space_discretization_type"].is("supg"))
    {
      la_type = LocalAssembling_type::TCDStiffMassRhs;
      assemble_both = true;
    }
    FEFunction *feFunction = &s.fe_function;
    LocalAssembling<d> la(this->db, la_type, &feFunction, example.get_coeffs());
    call_assembling_routine(s, la, assemble_both);
  }
  
  // here the modifications due to time discretization begin
  if (!db["algebraic_flux_correction"].is("none") )
  {
    do_algebraic_flux_correction();
    rhs_from_time_disc = this->systems.front().rhs;
    return; // modifications due to time discretization are per-
            // formed inside the afc scheme, so step out here!
  }

  // preparing the right hand side discretized by the used time
  // stepping scheme
  SystemPerGrid& s = this->systems.front();
  rhs_from_time_disc.reset();
  rhs_from_time_disc = s.rhs;
  // all matrices are available
  unsigned int n_sols = time_stepping_scheme.n_old_solutions();
  std::vector<BlockVector> old_sols(n_sols);
  old_sols[0] = s.solution_m1;
  if(old_sols.size() == 2)
    old_sols[1] = s.solution_m2;
  std::vector<BlockVector> rhs(2);
  rhs[0] = rhs_from_time_disc;
  rhs[1] = old_rhs;
  // prepare the right hand side from the previous time step
  time_stepping_scheme.prepare_rhs_from_time_disc(s.stiffness_matrix, 
                                                  s.mass_matrix, rhs, old_sols);
  rhs_from_time_disc = rhs[0];
  old_rhs = s.rhs;
  rhs_from_time_disc.copy_nonactive(s.rhs);
  
  for(auto &s : this->systems)
    time_stepping_scheme.prepare_system_matrix(s.stiffness_matrix,
                                               s.mass_matrix);
  s.solution.copy_nonactive(s.rhs);
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::solve()
{
  SystemPerGrid& s = systems.front();
#ifndef _MPI
  solver.solve(s.stiffness_matrix, rhs_from_time_disc, s.solution);
#else
  if(solver.get_db()["solver_type"].is("direct"))
  {
    MumpsWrapper mumps_wrapper(s.stiffness_matrix);
    mumps_wrapper.solve(rhs_from_time_disc, s.solution);
  }
  else
    // same as sequential
    solver.solve(s.stiffness_matrix, rhs_from_time_disc, s.solution);
#endif
  // restore stiffness matrix
  if(db["algebraic_flux_correction"].is("none"))
  {
    for(auto &s : this->systems)
      time_stepping_scheme.reset_linear_matrices(s.stiffness_matrix,
                                                 s.mass_matrix);
  }
  else
  {
    //in AFC case the stiffness matrix is "ruined" by now -
    Output::print("AFC does not yet reset the stiffness matrix and old_Au "
                  "correctly!");
  }
  s.solution_m2 = s.solution_m1;
  s.solution_m1 = s.solution;
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::output()
{
  SystemPerGrid &s = this->systems.front();
  s.fe_function.PrintMinMax();
  bool i_am_root = true;
#ifdef _MPI
  int my_rank;
  int root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == root);
  // computing errors as well as writing vtk files requires a minimum 
  // consistency level of 1
  s.fe_space->get_communicator().consistency_update(s.solution.get_entries(), 1);
#endif // _MPI
  
  //write solution for visualization 
  outputWriter.write(time_stepping_scheme.current_time_);
  
  // compute errors 
  if(db["output_compute_errors"])
  {
    const int n_errors = 4;
    std::array<double, n_errors+1> locError;
#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D all_derivatives[4] = { D000, D100, D010, D001 };
#else
    TAuxParam2D aux;
    MultiIndex2D all_derivatives[3] = { D00, D10, D01 };
#endif
    const FESpace* space = s.fe_space.get();
    
    s.fe_function.GetErrors(example.get_exact(0), d+1, all_derivatives,
                            n_errors, conv_diff_l2_h1_linf_error<d>,
                            example.get_coeffs(), &aux, 1, &space,
                            locError.data());
#ifdef _MPI
    /// @todo the GetErrors method in TFEFunction3D should already do the 
    /// communication, it's surprising that in the mpi case the errors are 
    /// squared, while the square root has been taken already in the sequential 
    /// case.
    // global (across all processes) error, L2, H1, Linf, SD-error (for SUPG)
    std::vector<double> errorsReduced(n_errors);
    MPI_Reduce(locError.data(), errorsReduced.data(), n_errors-1, MPI_DOUBLE,
               MPI_SUM, root, MPI_COMM_WORLD);
    MPI_Reduce(&locError[n_errors-1], &errorsReduced[n_errors-1], 1,
               MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
    // correct values only on the root process!
    for(int i = 0; i < n_errors-1; ++i)
      locError[i] = std::sqrt(errorsReduced[i]);
    locError[n_errors-1] = errorsReduced[n_errors-1];
#endif // _MPI
    
    if(i_am_root)
    {
      Output::print<1>("time   : ", time_stepping_scheme.current_time_);
      Output::print<1>("  L2   : ", std::setprecision(14), locError[0]);
      Output::print<1>("  H1   : ", std::setprecision(14), locError[1]);
      Output::print<1>("  L_inf: ", std::setprecision(14), locError.at(3));
    }
    double tau = time_stepping_scheme.get_step_length();
    
    // errors[1] is the previous L2 error squared
    errors[0] += (locError[0] * locError[0] + errors[1] * errors[1]) * tau*0.5;
    errors[1] = locError[0];
    // errors[3] is the previous H1 error squared
    errors[2] += (locError[1] * locError[1] + errors[3] * errors[3]) * tau*0.5;
    errors[3] = locError[1];
    
    if(errors.at(4) < locError[3])
      errors[4] = locError[3];
    
    if(i_am_root)
    {
      Output::print<1>("  L2(0,T;L2)      : ", sqrt(errors[0]));
      Output::print<1>("  L2(0,T;H1)      : ", sqrt(errors[2]));
      Output::print<1>("  L_inf(0,T,L_inf): " , errors[4]);
    }
  }
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::call_assembling_routine(
  TimeConvectionDiffusion<d>::SystemPerGrid& system, LocalAssembling<d>& la,
  bool assemble_both)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  
  int nFESpaces = 1;
  const FESpace * fe_space = system.fe_space.get();
  
  auto stiffBlock = system.stiffness_matrix.get_blocks_uniquely()[0].get();
  auto massBlock = system.mass_matrix.get_blocks_uniquely()[0].get();
  
  int nsqMatrices = 2; // maximum number of matrices
  SquareMatrixD *sqMatrices[nsqMatrices];
  if(assemble_both)
  {
    sqMatrices[0] = reinterpret_cast<SquareMatrixD*>(stiffBlock);
    sqMatrices[0]->reset();
    sqMatrices[1] = reinterpret_cast<SquareMatrixD*>(massBlock);
    sqMatrices[1]->reset();
  }
  else
  {
    nsqMatrices = 1; // Mass and Stiffness Matrices at once
    sqMatrices[0] = reinterpret_cast<SquareMatrixD*>(stiffBlock);
    sqMatrices[0]->reset();
  }
  
  int nRectMatrices = 0;
  int nRhs = 1;
  double *rhsEntries = system.rhs.get_entries();
  system.rhs.reset(); // reset to zeros 
  
  auto * boundary_conditions = fe_space->get_boundary_condition();
  auto * boundary_value = example.get_bd(0);
  
#ifdef __3D__
    Assemble3D(
#else      
    Assemble2D(
#endif
               nFESpaces, &fe_space, nsqMatrices, sqMatrices,  nRectMatrices,
               nullptr, nRhs, &rhsEntries, &fe_space, &boundary_conditions,
               &boundary_value, la);
}

/* ************************************************************************* */
template <int d>
std::array<double, 3> TimeConvectionDiffusion<d>::get_errors() const
{
  std::array<double, 3> error_at_time_points;
  error_at_time_points.at(0) = sqrt(errors.at(0));
  error_at_time_points.at(1) = sqrt(errors.at(2));
  error_at_time_points.at(2) = sqrt(errors.at(4));
  return error_at_time_points;
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::do_algebraic_flux_correction()
{
  //determine which kind of afc to use
  if(db["algebraic_flux_correction"].is("default") ||
      db["algebraic_flux_correction"].is("fem-fct-cn"))
  {
    //TODO implement for multigrid!
    SystemPerGrid& s = systems.front();

    //get references to the relevant objects

    const FESpace& feSpace = *s.fe_space;
    const FEMatrix& mass = *s.mass_matrix.get_blocks().at(0).get();
    FEMatrix& stiff = *s.stiffness_matrix.get_blocks_uniquely().at(0).get();
    //vector entry arrays
#ifdef _MPI
    const TParFECommunicator3D& comm = s.fe_space->get_communicator();
    // solution is needed in consistency 2 for afc.
    comm.consistency_update(s.solution.block(0), 2);
#endif
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
    double delta_t = time_stepping_scheme.get_step_length();

    // apply C-N FEM-FCT
    AlgebraicFluxCorrection::crank_nicolson_fct(mass, stiff, solEntries,
                                                rhsEntries, oldRhsEntries,
                                                delta_t, neumToDiri, prelim);

    //...and finally correct the entries in the Dirichlet rows
    AlgebraicFluxCorrection::correct_dirichlet_rows(stiff);
    //...and in the right hand side, too
    s.rhs.copy_nonactive(old_rhs);

#ifdef _MPI
    // this update is necessary here, because copy_nonactive()
    // can disturb the consistency of a vector,
    // see Time_NSE3D::assemble_rhs() for more details
    comm.consistency_update(s.rhs.block(0), 2);
#endif
  }
  else
  {
    ErrThrow("The chosen algebraic flux correction scheme ",
             db["algebraic_flux_correction"],
             " is unknown to class TimeConvectionDiffusion.");
  }
}

#ifdef __3D__
template class TimeConvectionDiffusion<3>;
#else
template class TimeConvectionDiffusion<2>;
#endif

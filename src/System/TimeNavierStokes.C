#include "TimeNavierStokes.h"
#include "Database.h"

#include "LocalProjection.h"
#include "Hotfixglobal_AssembleNSE.h"
#include "GridTransfer.h"
#include "Multigrid.h"
#ifdef __2D__
#include "Upwind.h"
#include "Matrix2D.h"
#include "SquareMatrix2D.h"
#include "Assemble2D.h"
#include "AuxParam2D.h"
#else
#include "Upwind3D.h"
#include "Matrix3D.h"
#include "SquareMatrix3D.h"
#include "Assemble3D.h"
#include "AuxParam3D.h"
#endif
#ifdef _MPI
#include "ParFECommunicator3D.h"
#endif

/* ************************************************************************** */
template <int d>
ParameterDatabase TimeNavierStokes<d>::default_tnse_database()
{
  Output::print<5>("creating a default TNSE parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default tnse database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TimeNavierStokes parameter database");
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(ParameterDatabase::default_solution_in_out_database());
  db.merge(LocalAssembling<d>::default_local_assembling_database());
  return db;
}

/* ************************************************************************** */
template <int d>
TimeNavierStokes<d>::System_per_grid::System_per_grid(
  const Example_TimeNSE& example, TCollection& coll, std::pair<int, int> order)
 : velocity_space(new FESpace(&coll, "u", "velocity space", example.get_bc(0),
                              order.first)),
   pressure_space(new FESpace(&coll, "p", "pressure space", example.get_bc(d),
                              order.second))
{
  switch(TDatabase::ParamDB->NSTYPE)
  {
#ifdef __2D__
    case 1:
      matrix = BlockFEMatrix::NSE2D_Type1(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type1(*velocity_space, *pressure_space);
      break;
    case 2:
      matrix = BlockFEMatrix::NSE2D_Type2(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type2(*velocity_space, *pressure_space);
      break;
    case 3:
      matrix = BlockFEMatrix::NSE2D_Type3(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type3(*velocity_space, *pressure_space);
      break;
    case 4:
      matrix = BlockFEMatrix::NSE2D_Type4(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type4(*velocity_space, *pressure_space);
      break;
    case 14:
      matrix = BlockFEMatrix::NSE2D_Type14(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type4(*velocity_space, *pressure_space);
      break;
#else
    case 1:
      matrix = BlockFEMatrix::NSE3D_Type1(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type1(*velocity_space, *pressure_space);
      break;
    case 2:
      matrix = BlockFEMatrix::NSE3D_Type2(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type2(*velocity_space, *pressure_space);
      break;
    case 3:
      matrix = BlockFEMatrix::NSE3D_Type3(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type3(*velocity_space, *pressure_space);
      break;
    case 4:
      matrix = BlockFEMatrix::NSE3D_Type4(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type4(*velocity_space, *pressure_space);
      break;
    case 14:
      matrix = BlockFEMatrix::NSE3D_Type14(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type4(*velocity_space, *pressure_space);
      break;
#endif
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }
  rhs = BlockVector(matrix, true);
  solution = BlockVector(matrix, false);
  u = FEVectFunct(velocity_space.get(), "u", "u", solution.block(0),
                  solution.length(0), d);
  p = FEFunction(pressure_space.get(), "p", "p", this->solution.block(d),
                 solution.length(d));
  solution_m1 = BlockVector(matrix, false);
  u_m1 = FEVectFunct(velocity_space.get(), "u", "u", solution_m1.block(0),
                     solution_m1.length(0), d);
  p_m1 = FEFunction(pressure_space.get(), "p", "p", this->solution_m1.block(d),
                    solution_m1.length(d));
  solution_m2 = BlockVector(matrix, false);
  u_m2 = FEVectFunct(velocity_space.get(),"u","u", solution_m2.block(0),
                     solution_m2.length(0), d);
  p_m2 = FEFunction(pressure_space.get(), "p", "p", this->solution_m2.block(d),
                    solution_m2.length(d));

  time_avg_sol = BlockVector(matrix, false);
  u_time_avg = FEVectFunct(velocity_space.get(), "u_t_avg", "u time averaged",
                           time_avg_sol.block(0), time_avg_sol.length(0), d);
  p_time_avg = FEFunction(pressure_space.get(), "p_t_avg", "p time averaged",
                          this->time_avg_sol.block(d), time_avg_sol.length(d));

  combined_old_sols = BlockVector(matrix, false);
  comb_old_u = FEVectFunct(velocity_space.get(), "u", "u",
                           combined_old_sols.block(0),
                           combined_old_sols.length(0), d);
  extrapolate_sol = BlockVector(matrix, false);
  extrapolate_u = FEVectFunct(velocity_space.get(), "u", "u",
                              extrapolate_sol.block(0),
                              extrapolate_sol.length(0), d);
#ifdef _MPI
  //print some information
  velocity_space->get_communicator().print_info();
  pressure_space->get_communicator().print_info();
#endif
}

/* ************************************************************************** */
template <int d>
TimeNavierStokes<d>::System_per_grid::System_per_grid(
  const System_per_grid& other)
 : velocity_space(other.velocity_space), pressure_space(other.pressure_space),
   matrix(other.matrix), rhs(other.rhs), solution(other.solution)
{
  // the fe functions must be newly created, because copying would mean
  // referencing the BlockVectors in 'other'.
  u = FEVectFunct(velocity_space.get(), "u", "u", solution.block(0),
                  solution.length(0), d);
  p = FEFunction(pressure_space.get(), "p", "p", solution.block(d),
                 solution.length(d));
}

/* ************************************************************************** */
template <int d>
TimeNavierStokes<d>::TimeNavierStokes(const TDomain& domain,
                                      const ParameterDatabase& param_db)
: TimeNavierStokes<d>(domain, param_db, Example_TimeNSE(param_db))
{
}

/* ************************************************************************** */
template <int d>
TimeNavierStokes<d>::TimeNavierStokes(const TDomain& domain,
                                      const ParameterDatabase& param_db,
                                      const Example_TimeNSE& ex)
 : db(default_tnse_database()), systems(), outputWriter(param_db), example(ex),
   solver(param_db), defect(), old_residuals(), initial_residual(1e10),
   time_stepping_scheme(param_db), is_rhs_and_mass_matrix_nonlinear(false)
#ifdef __3D__
   , Lines()
#endif
{
  db.merge(param_db);
  this->check_and_set_parameters();

  std::pair<int,int> velo_pres_order(TDatabase::ParamDB->VELOCITY_SPACE,
                                     TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and pressure spaces
  // this function returns a pair which consists of
  // velocity and pressure order
  this->get_velocity_pressure_orders(velo_pres_order);
  
  bool usingMultigrid = this->solver.is_using_multigrid();
  auto collections = domain.get_grid_collections();
  TCollection *coll = collections.front(); // finest grid collection
  // create finite element space, functions, matrices, rhs and solution
  // at the finest grid
  this->systems.emplace_back(example, *coll, velo_pres_order);

  if(usingMultigrid)
  {
    // Construct multigrid object
    auto mg = solver.get_multigrid();
    size_t n_multigrid_levels = mg->get_n_geometric_levels();
    size_t n_grids=collections.size();
    if(n_multigrid_levels > n_grids)
    {
      ErrThrow("Wrong number of grids for multigrid! expecting ",
               n_multigrid_levels, " geometric grids but got", n_grids,".");
    }
    // remove not needed coarser grid from list of collections
    for(size_t i = n_multigrid_levels; i < n_grids; ++i)
    {
      collections.pop_back();
    }
    
    if(mg->is_using_mdml())
    {
      // change the discretization on the coarse grids to lowest order
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velo_pres_order = {-1, -4711};
      this->get_velocity_pressure_orders(velo_pres_order);
    }
    else
    {
      // for standard multigrid, pop the finest collection - it was already
      // used to construct a space before the "if(usingMultigrid)" clause
      // and will not (as in mdml) be used a second time with a different discretization
      collections.pop_front();
    }
    
    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    // matrix on finest grid is already constructed
    matrices.push_back(&systems.back().matrix);
    for(auto coll : collections)
    {
      systems.emplace_back(example, *coll, velo_pres_order);
      // prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    // initialize the multigrid object with all the matrices on all levels
    mg->initialize(matrices);
  }
  
  // initial solution on finest grid - read-in or interpolation
  if(db["read_initial_solution"].is(true))
  {
    if(!this->time_stepping_scheme.get_start_time())
    {
      Output::warn<1>("Initial Solution",
        "Restarting from existing solution but initial time is 0! This is "
        "probably not what you want! If for example your BC or RHS are "
        "time-dependent, you will apply initial BC or RHS to an already "
        "developed flow, instead of continuing your simulation! Set your "
        "time_start to the time of the binary file you are re-starting from "
        "(and don't forget to set continue_output_after_restart to true, and "
        "also read_metis parameters if needed). Or ignore this warning if you "
        "know what you are doing and are aware of the consequences (rhs and bc "
        "reset to t = 0, output restart from 0 and probably overwrites outputs "
        "from old simulation).");
    }
    std::string file = db["initial_solution_file"];
    Output::root_info("Initial Solution", "Reading initial solution from file ", file);
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    file += ".proc" + std::to_string(my_rank);
    Output::root_info("Initial Solution", "Appending .proc<RANK> to the "
        "expected initial solution file name.");
#endif
    systems.front().solution.read_from_file(file);
  }
  else
  {
    Output::info("Initial Solution", "Interpolating initial solution from example.");
    for(System_per_grid& s : this->systems)
    {
      for(int i = 0; i < d; ++i)
      {
        FEFunction * ui = s.u.GetComponent(i);
        ui->Interpolate(example.get_initial_cond(i));
        delete ui;
      }
    }
  }
  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);
  this->rhs_from_time_disc.copy_structure(this->systems.front().rhs);

  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());

  // print out the information (cells, dofs, etc)
  this->output_problem_size_info();
  this->errors.fill(0.);

#ifdef __3D__  
  if( db["output_along_line"] )
  {
    Lines = LinesEval<d>(domain, param_db);
  }
#endif
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::check_and_set_parameters()
{
  if(!db["problem_type"].is(6))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 6;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to "
                      "TimeNavierStokes. It is now reset to the correct value "
                      "for TimeNavierStokes (=6).");
      db["problem_type"] = 6;
    }
  }
  if(db["time_discretization"].is("forward_euler"))
  {
    ErrThrow("time discretization: ", db["time_discretization"],
             " is not supported");
  }

  if(db["imex_scheme_"])
  {
    // first two steps are performed with the full nonlinear 
    db["extrapolate_velocity"] = true;
  }
  // set the discretization parameters
  // standard Galerkin
  if(db["space_discretization_type"].is("galerkin"))
  {
    space_disc_global = 1;
    /// set scaling factor for B, BT's block
    time_stepping_scheme.n_scale_block = 2*d;
    time_stepping_scheme.b_bt_linear_nl = "linear";
  }
  
  if(db["space_discretization_type"].is("supg"))
  {
    // supg: NOTE: only tested with BDF2 so far
    if(!db["time_discretization"].is("bdf_two"))
    {
      ErrThrow("supg method is only implemented for BDF2 time stepping scheme");
    }
    space_disc_global = 2;
    /// set scaling factor for B, BT's block
    // depends on how to deal the nonlinearity in the 
    // test function: fully implicit case
    time_stepping_scheme.b_bt_linear_nl = "nonlinear";
    time_stepping_scheme.n_scale_block = d;
    if(TDatabase::ParamDB->NSTYPE == 14)
      time_stepping_scheme.n_scale_block = 2*d+1;
  }
  
  // Smagorinsky
  if(db["space_discretization_type"].is("smagorinsky"))
  {
    space_disc_global = 4;
    time_stepping_scheme.n_scale_block = 6;
    time_stepping_scheme.b_bt_linear_nl = "linear";
    // This is a hot fix: only the Smagorinsky routines with Laplace Type=1
    // (D(u):D(v) formulation) can be used properly. In the other case, the
    // block matrices passed to the routine are not correct, throwing a segfault.
    if(TDatabase::ParamDB->LAPLACETYPE==0)
    {
      ErrThrow("Smagorinsky works only with LAPLACETYPE = 1 so far.");
    }
  }

  if(db["space_discretization_type"].is("local_projection"))
  {
    if(d == 3)
      ErrThrow("local_projection in 3D not implemented");
    space_disc_global = 14;
  }
  // the only case where one have to re-assemble the right hand side
  if(db["space_discretization_type"].is("supg")
     && db["time_discretization"].is("bdf_two"))
  {
    is_rhs_and_mass_matrix_nonlinear = true;
  }  
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::get_velocity_pressure_orders(
  std::pair<int, int> &velo_pres_order) 
{
  int velocity_order = velo_pres_order.first;
  int pressure_order = velo_pres_order.second;
  int order = 0;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5:
    case 12: case 13: case 14: case 15:
      if(velocity_order > 10)
        order = velocity_order-10;
      else
        order = velocity_order;
      break;
    case -1: case -2: case -3: case -4: case -5: case -101:
      order = velocity_order;
      break;
    case 100: case 201: case 302: case 403: case 504:
      if(d == 3)
        ErrThrow("velocity_order ", velocity_order, " not supported in 3D");
      order = velocity_order;
      break;
    // conforming fe spaces with bubbles on triangles
    case 22: case 23: case 24:
      order = velocity_order;
      break;
      // discontinuous spaces
    case -11: case -12: case -13:
      order = velocity_order*10;
      break;
  }
  TDatabase::ParamDB->VELOCITY_SPACE = order;
  velo_pres_order.first = order;
  switch(pressure_order)
  {
    case -4711:
    {
      switch(velocity_order)
      {
        case -1: case -2: case -3: case -4:
          // nonconforming pw (bi)linear velo/ pw constant pressure
          // conforming pw (bi)linear velo/ pw constant pressure (not stable !!!)
          pressure_order = -velocity_order-1;
          break;
        case 1: // discontinuous space
          pressure_order = 0;
          break;
        case 2: case 3: case 4: case 5:
        // standard conforming velo and continuous pressure
          pressure_order = velocity_order-1;
          break;
          // discontinuous pressure spaces
          // standard conforming velo and discontinuous pressure
          // this is not stable on triangles !!!
        case 12: case 13: case 14: case 15:
        case -11: case -12: case -13: case -14:
          pressure_order = -(velocity_order-1)*10;
          break;
        case 22: case 23: case 24:
          pressure_order = -(velocity_order-11)*10;
          break;
        case 100: case 201: case 302: case 403: case 504:
          pressure_order = -(velocity_order%100 + 10)*10;
          break; 
      }
      break;
    }
    case 1: case 2: case 3: case 4: case 5:
      // pressure order is chosen correctly
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
    case 100: case 201: case 302: case 403: case 504:
      if(d == 3)
        ErrThrow("pressure_order ", pressure_order, " not supported in 3D");
      // pressure order is chosen correctly
      break;
    default:
      ErrThrow("pressure space is not chosen properly ", pressure_order);
  }
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  velo_pres_order.second = pressure_order;

  Output::print("velocity space ", setw(6), velo_pres_order.first);
  Output::print("pressure space ", setw(6), velo_pres_order.second);
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::assemble_initial_time()
{
  if(systems.size() > 1) // using multigrid
  {
    this->restrict_function();
  }
  for(auto &s : this->systems)
  {
    call_assembling_routine(s, LocalAssembling_type::TNSE3D_LinGAL);
    //update matrices for local projection stabilization
    if(db["space_discretization_type"].is("local_projection"))
      update_matrices_lps(s);
    // copy nonactives
    s.solution.copy_nonactive(s.rhs);

  /** After copy_nonactive, the solution vectors needs to be Comm-updated in 
   * MPI-case in order to be consistently saved. It is necessary that the vector
   * is consistently saved because it is the only way to ensure that its
   * multiplication with an inconsistently saved matrix (multiplication which
   * appears in the defect and rhs computations) give the correct results. When
   * we call copy_nonactive in MPI-case, we have to remember the following: it
   * can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active. That's why we have to update the values so that the vector
   * becomes consistent again.
   */
#ifdef _MPI
  for(int i = 0; i < d; ++i)
  {
    double *ui = s.solution.block(i);
    s.velocity_space->get_communicator().consistency_update(ui, 3);
  }
  double *p = s.solution.block(d);
  s.pressure_space->get_communicator().consistency_update(p, 3);
#endif
    
    s.solution_m1 = s.solution;
    s.solution_m2 = s.solution;
  }
  
  // the call to assembleslip is necessary here, in order to get
  // the correct old_rhs, i.e., zeros on the slip dofs
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1)
   this->modify_slip_bc(true, true);
  
  auto s = this->systems.front();
 
  // copy the current right hand side vector to the old_rhs
  this->old_rhs = s.rhs;
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::assemble_matrices_rhs(unsigned int it_counter)
{
  if(it_counter == 0)
  {
    // initialize the rhs from the time discretization
    System_per_grid& s = this->systems.front();
    rhs_from_time_disc = s.rhs;
    rhs_from_time_disc.reset();
    // only assemble the right-hand side
    call_assembling_routine(s, LocalAssembling_type::TNSE3D_Rhs);
    // copy the non active to the solution vector
    // since the rhs vector will be passed to the solver
    // and is modified with matrix vector multiplication
    // which also uses the non-actives
    s.solution.copy_nonactive(s.rhs);
    // copy the right hand side to the "rhs_from_time_disc"
    rhs_from_time_disc = s.rhs;
    // all matrices from the previous time step are available
    unsigned int n_sols = time_stepping_scheme.n_old_solutions();
    std::vector<BlockVector> oldsolutions(n_sols);
    oldsolutions[0] = s.solution_m1;
    if(oldsolutions.size() == 2)
      oldsolutions[1] = s.solution_m2;
    
    // modification of the matrices due to slip b.c
    if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1)
        this->modify_slip_bc(true, true);
    
    // one needs two right hand sides only for the crank-Nicolson
    // and fractional step theta schemes
    std::vector<BlockVector> all_rhs(2);
    all_rhs[0] = s.rhs; // current rhs
    all_rhs[1] = old_rhs;
   
    //NOTE: scale the B blocks only at the first iteration
    for(System_per_grid& sys : this->systems)
      time_stepping_scheme.scale_descale_all_b_blocks(sys.matrix, "scale");
   
    // prepare the right hand side for the solver
    time_stepping_scheme.prepare_rhs_from_time_disc(s.matrix, s.mass_matrix,
                                                    all_rhs, oldsolutions);
    rhs_from_time_disc = all_rhs[0];
    old_rhs = s.rhs;
    // copy the non-actives
    rhs_from_time_disc.copy_nonactive(s.solution);
    old_rhs.copy_nonactive(s.solution);
  
    /** After copy_nonactive, the solution vectors needs to be Comm-updated in 
     * MPI-case in order to be consistently saved. It is necessary that the vector
     * is consistently saved because it is the only way to ensure that its
     * multiplication with an inconsistently saved matrix (multiplication which
     * appears in the defect and rhs computations) give the correct results. When
     * we call copy_nonactive in MPI-case, we have to remember the following: it
     * can happen that some slave ACTTIVE DoFs are placed in the block of
     * NON-ACTIVE DoFs (because they are at the interface between processors).
     * Doing copy_nonactive changes then the value of these DOFs,although they are
     * actually active. That's why we have to update the values so that the vector
     * becomes consistent again.
     */
#ifdef _MPI
    for(int i = 0; i < d; ++i)
    {
      double *ui = s.solution.block(i);
      s.velocity_space->get_communicator().consistency_update(ui, 3);
    }
    double *p = s.solution.block(d);
    s.pressure_space->get_communicator().consistency_update(p, 3);
#endif // _MPI
  }
  // assemble the nonlinear matrices
  if(systems.size() > 1)
    this->restrict_function();
  for(System_per_grid & s : systems)
  {
    call_assembling_routine(s, LocalAssembling_type::TNSE3D_NLGAL);
    // update matrices with local projection term
    if(db["space_discretization_type"].is("local_projection"))
      update_matrices_lps(s);
  }

  // slip boundary modification of matrices
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1 )
  {
    if(is_rhs_and_mass_matrix_nonlinear)
      this->modify_slip_bc(true, true);
    else
      this->modify_slip_bc();
  }
  // prepare the right-hand side if it is nonlinear, assembling already 
  // done together with the matrices
  if(is_rhs_and_mass_matrix_nonlinear)
  {
    this->assemble_rhs_nonlinear();
  }
  //
  // also prepare the system matrices for the solver
  for(System_per_grid& s : this->systems)
  {
    // call the preparing method
    time_stepping_scheme.prepare_system_matrix(s.matrix, s.mass_matrix);
    if(db["space_discretization_type"].is("supg"))
      time_stepping_scheme.scale_nl_b_blocks(s.matrix);
  }
  Output::print<5>("Assembling of matrices and right hand side is done");
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::assemble_rhs_nonlinear()
{
  System_per_grid& s = this->systems.front();
  
  // copy the right hand side to the "rhs_from_time_disc"
  rhs_from_time_disc = s.rhs;
  // all matrices from the previous time step are available
  unsigned int n_sols = time_stepping_scheme.n_old_solutions();
  std::vector<BlockVector> oldsolutions(n_sols);
  oldsolutions[0] = s.solution_m1;
  if(oldsolutions.size() == 2)
    oldsolutions[1] = s.solution_m2;
  // one needs two right hand sides only for the crank-Nicolson
  // and fractional step theta schemes
  std::vector<BlockVector> rhs_(2);
  rhs_[0] = rhs_from_time_disc; // current rhs
  rhs_[1] = old_rhs; // old right hand side is needed for the Crank-Nicolson time stepping

  // prepare the right hand side for the solver
  time_stepping_scheme.prepare_rhs_from_time_disc(s.matrix, s.mass_matrix,
                                                  rhs_, oldsolutions);
  rhs_from_time_disc = rhs_[0];
  old_rhs = s.rhs;
  // copy the non-actives
  rhs_from_time_disc.copy_nonactive(s.rhs);
  s.solution.copy_nonactive(s.rhs);
  
  /** After copy_nonactive, the solution vectors needs to be Comm-updated in 
   * MPI-case in order to be consistently saved. It is necessary that the vector
   * is consistently saved because it is the only way to ensure that its
   * multiplication with an inconsistently saved matrix (multiplication which
   * appears in the defect and rhs computations) give the correct results. When
   * we call copy_nonactive in MPI-case, we have to remember the following: it
   * can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active. That's why we have to update the values so that the vector
   * becomes consistent again.
   */
#ifdef _MPI
  for(int i = 0; i < d; ++i)
  {
    double *ui = s.solution.block(i);
    s.velocity_space->get_communicator().consistency_update(ui, 3);
  }
  double *p = s.solution.block(d);
  s.pressure_space->get_communicator().consistency_update(p, 3);
#endif // _MPI
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::call_assembling_routine(
  TimeNavierStokes<d>::System_per_grid& s, LocalAssembling_type type)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  // set arrays of spaces for matrices and rhs
  std::vector<const FESpace*> spaces_mat;
  std::vector<const FESpace*> spaces_rhs;
  std::vector<FEFunction*> fefunctions;
  // call to routine to set arrays
  set_arrays(s, spaces_mat, spaces_rhs, fefunctions);

  // prepare matrices and rhs for assembling
  std::vector<SquareMatrixD*> sqMatrices;
  std::vector<MatrixD*> rectMatrices;
  std::vector<double*> rhs_array;
  // call the routine to prepare the matrices
  set_matrices_rhs(s, type, sqMatrices, rectMatrices, rhs_array);
  std::array<BoundaryConditionFunction*, d+1> boundCondition;
  std::array<BoundaryValuesFunction*, d+1> boundValues;
  for(int i = 0; i < d; ++i)
    boundCondition[i] = s.velocity_space->get_boundary_condition();
  boundCondition[d] = s.pressure_space->get_boundary_condition();
  for(int i = 0; i < d+1; ++i)
    boundValues[i] = example.get_bd(i);

  // local assembling settings
  LocalAssembling<d> la(this->db, type, fefunctions.data(),
                       this->example.get_coeffs(), space_disc_global);

  // default value
  assemble_nse = Hotfixglobal_AssembleNSE::WITH_CONVECTION;
  // find out if we have to do upwinding
  bool do_upwinding = false;
  if(type != LocalAssembling_type::TNSE3D_Rhs)
  {
    bool mdml = this->solver.is_using_multigrid()
                && this->solver.get_multigrid()->is_using_mdml();
    bool on_finest_grid = &systems.front() == &s;
    do_upwinding = (db["space_discretization_type"].is("upwind")
                    || (mdml && !on_finest_grid));
    /// @todo does upwinding work in tnse?
    if(do_upwinding)  // HOTFIX: Check the documentation!
      assemble_nse = Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION;
    else
      assemble_nse = Hotfixglobal_AssembleNSE::WITH_CONVECTION;
  }

#ifdef __3D__
    Assemble3D(
#else
    Assemble2D(
#endif
               spaces_mat.size(), spaces_mat.data(), sqMatrices.size(),
               sqMatrices.data(), rectMatrices.size(), rectMatrices.data(),
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
               boundCondition.data(), boundValues.data(), la);

  if(do_upwinding)
  {
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
      {
        // do upwinding with one matrix
#ifdef __2D__
        UpwindForNavierStokes(la.GetCoeffFct(), sqMatrices[0],
                              la.get_fe_function(0),
                              la.get_fe_function(1));
#else
        // the inverse of the example's diffusion coefficient
        double one_over_nu = 1./example.get_nu();
        UpwindForNavierStokes3D(sqMatrices[0], la.get_fe_function(0),
                                la.get_fe_function(1), la.get_fe_function(2),
                                one_over_nu);
#endif
        Output::print<3>("UPWINDING DONE");
        break;
      }
      case 3:
      case 4:
      case 14:
      {
        // do upwinding with d matrices
#ifdef __2D__
        UpwindForNavierStokes(la.GetCoeffFct(), sqMatrices[0],
                              la.get_fe_function(0),
                              la.get_fe_function(1));
        UpwindForNavierStokes(la.GetCoeffFct(), sqMatrices[1],
                              la.get_fe_function(0),
                              la.get_fe_function(1));
#else
        // the inverse of the example's diffusion coefficient
        double one_over_nu = 1./example.get_nu();
        UpwindForNavierStokes3D(sqMatrices[0], la.get_fe_function(0),
                                la.get_fe_function(1), la.get_fe_function(2),
                                one_over_nu);
        UpwindForNavierStokes3D(sqMatrices[1], la.get_fe_function(0),
                                la.get_fe_function(1), la.get_fe_function(2),
                                one_over_nu);
        UpwindForNavierStokes3D(sqMatrices[2], la.get_fe_function(0),
                                la.get_fe_function(1), la.get_fe_function(2),
                                one_over_nu);
#endif
        Output::print<3>("UPWINDING DONE");
        break;
      }
    } // endswitch
  }
  for (int i = 0; i < d; i++)
    delete fefunctions[i];
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::set_arrays(TimeNavierStokes<d>::System_per_grid& s,
                                     std::vector<const FESpace*>& spaces,
                                     std::vector<const FESpace*>& spaces_rhs,
                                     std::vector<FEFunction*>& functions)
{
  spaces.resize(2);
  spaces[0] = s.velocity_space.get();
  spaces[1] = s.pressure_space.get();
  
  spaces_rhs.resize(d+1);
  for(int i = 0; i < d; ++i)
    spaces_rhs[i] = s.velocity_space.get();
  spaces_rhs[d] = s.pressure_space.get();
  // standard for all methods.
  functions.resize(d+1);  
  for(int i = 0; i < d; ++i)
    functions[i] = s.u.GetComponent(i);
  functions[d] = &s.p;
  
  bool is_imex = imex_scheme();
  if(is_imex && db["extrapolate_velocity"])
  {
    if(db["space_discretization_type"].is("galerkin")
      || db["space_discretization_type"].is("local_projection"))
    {
      if(db["extrapolation_type"].is("constant_extrapolate"))
      {
        s.extrapolate_sol.reset();
        s.extrapolate_sol = s.solution_m1;
        
        for(int i = 0; i < d; ++i)
        {
          delete functions[i];
          functions[i] = s.extrapolate_u.GetComponent(i);
        }
      }
      else if(db["extrapolation_type"].is("linear_extrapolate"))
      {
        s.extrapolate_sol.reset();
        s.extrapolate_sol = s.solution_m1;
        s.extrapolate_sol.scale(2.);
        s.extrapolate_sol.add_scaled(s.solution_m2, -1.);
        
        for(int i = 0; i < d; ++i)
        {
          delete functions[i];
          functions[i] = s.extrapolate_u.GetComponent(i);
        }
      }
      else
      {
        ErrThrow("Only constant or linear extrapolation of velocity are used ",
                 db["extrapolation_type"]);
      }
    }
    // supg: NOTE: only tested with BDF2 so far
    if(db["space_discretization_type"].is("supg")
      && !db["time_discretization"].is("bdf_two"))
    {
      ErrThrow("supg method is only implemented for BDF2 time stepping scheme");
    }
    
    if(db["space_discretization_type"].is("supg"))
    {
      if(TDatabase::ParamDB->NSTYPE < 4)
      {
        ErrThrow("SUPG is implemented for equal-order only so far !!!");        
      }
      functions.resize(2*d);
      
      s.extrapolate_sol.reset();
      s.extrapolate_sol = s.solution_m1;
      s.extrapolate_sol.scale(2.);
      s.extrapolate_sol.add_scaled(s.solution_m2, -1.);
      for(int i = 0; i < d; ++i)
      {
        delete functions[i];
        functions[i] = s.extrapolate_u.GetComponent(i);
      }
      // combination of previous time solutions for assembling the right-hand
      // side, this is used for the pressure part.
      s.combined_old_sols.reset();
      // copy and scale the solution at previous time step with factor 2
      s.combined_old_sols = s.solution_m1;
      s.combined_old_sols.scale(2.);
      // subtract with right factor the solution at pre-previous solution
      s.combined_old_sols.add_scaled(s.solution_m2, -1./2.);
      
      for(int i = 0; i < d; ++i)
      {
        functions[i+d] = s.comb_old_u.GetComponent(i);
      }
    }
  }
  else 
  {
    // For the standard methods or symmetric stabilization schemes
    // there is no time derivative involved in combination with the 
    // pressure or velocity as in the supg case. So nothing to do 
    // for those schemes.
    
    // 
    if(db["space_discretization_type"].is("supg"))
    {
      functions.resize(2*d);
      if(time_stepping_scheme.pre_stage_bdf)
      {
        for(int i = 0; i < d; ++i)
          functions[i+d] = s.u_m1.GetComponent(i);
      }
      else
      {
        s.combined_old_sols.reset();
        // copy and scale the solution at previous time step with factor 2
        s.combined_old_sols = s.solution_m1;
        s.combined_old_sols.scale(2.);
        // subtract with right factor the solution at pre-previous solution
        s.combined_old_sols.add_scaled(s.solution_m2, -1./2.);
        
        for(int i = 0; i < d; ++i)
          functions[i+d] = s.comb_old_u.GetComponent(i);
      }
    }
  }
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::set_matrices_rhs(
  TimeNavierStokes<d>::System_per_grid& s, LocalAssembling_type type,
  std::vector<SquareMatrixD*>& sqMat, std::vector<MatrixD*>& reMat,
  std::vector<double*>& rhs_array)
{
  sqMat.resize(0);
  reMat.resize(0);
  // right hand side: for NSTYPE: 1,2 and 3, size is 2
  rhs_array.resize(d+1, nullptr);
  

  auto blocks = s.matrix.get_blocks_uniquely();
  auto mass_blocks = s.mass_matrix.get_blocks_uniquely(true);
  switch(type)
  {
    case LocalAssembling_type::TNSE3D_LinGAL:
    {
      // right hand sides are assembled for the initial time step
      // for the remaining time steps, they are assembled in
      // another function. so reset to zero here
      s.rhs.reset();
      for(int i = 0; i < d+1; ++i)
        rhs_array[i] = s.rhs.block(i);
      
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          if(blocks.size() != d+1)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size());
          }
          sqMat.resize(2);
          sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<SquareMatrixD*>(mass_blocks.at(0).get());
          // rectangular matrices
          reMat.resize(d);
          for(int i = 0; i < d; ++i)
            reMat[i] = reinterpret_cast<MatrixD*>(blocks.at(i+1).get());
          break;
        case 2:
          if(blocks.size() != 1+2*d)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size());
          }
          sqMat.resize(2);
          sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<SquareMatrixD*>(mass_blocks.at(0).get());
          // rectangular matrices
          reMat.resize(2*d);
          for(int i = 0; i < d; ++i) // first the lying B blocks
            reMat[i] = reinterpret_cast<MatrixD*>(blocks.at(d+1+i).get());
          for(int i = 0; i < d; ++i) // than the standing B blocks
            reMat[d+i] = reinterpret_cast<MatrixD*>(blocks.at(i+1).get());
          break;
        case 3:
          if(blocks.size() != d*d+d)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size());
          }
          sqMat.resize(d*d+d);
          for(int i = 0, j = 0; i < d*d; ++i, ++j)
          {
            if(i%d == 0 && i > 0)
              j++;
            sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
          }
          // mass matrix
          for(int i = 0; i < d; ++i)
            sqMat[d*d+i]
              = reinterpret_cast<SquareMatrixD*>(mass_blocks.at(i*(d+2)).get());
          // rectangular matrices
          reMat.resize(d);
          for(int i = 0; i < d; ++i)
            reMat[i] = reinterpret_cast<MatrixD*>(blocks.at(i*(d+1)+d).get());
          break;
        case 4:
        case 14:
          if( (TDatabase::ParamDB->NSTYPE == 14 && blocks.size() != (d+1)*(d+1))
           || (TDatabase::ParamDB->NSTYPE == 4  && blocks.size() != d*d+2*d))
          {
            ErrThrow("Wrong blocks.size() ", blocks.size());
          }
          sqMat.resize(d*d+d);
          for(int i = 0, j = 0; i < d*d; ++i, ++j)
          {
            if(i%d == 0 && i > 0)
              j++;
            sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
          }
          // mass matrix
          for(int i = 0; i < d; ++i)
            sqMat[d*d+i]
              = reinterpret_cast<SquareMatrixD*>(mass_blocks.at(i*(d+2)).get());
          if(TDatabase::ParamDB->NSTYPE == 14)
          {
            // C block
            sqMat.resize(d*d+d+1);
            sqMat[d*d+d] = reinterpret_cast<SquareMatrixD*>(blocks.at(d*d+2*d).get());
          }
          
          // rectangular matrices
          reMat.resize(2*d);
          for(int i = 0; i < d; ++i) // first the lying B blocks
            reMat[i] = reinterpret_cast<MatrixD*>(blocks.at(d*d+d+i).get());
          for(int i = 0; i < d; ++i) // than the standing B blocks
            reMat[d+i] = reinterpret_cast<MatrixD*>(blocks.at(i*(d+1)+d).get());
          break;
      }
      break;
    }
    case LocalAssembling_type::TNSE3D_NLGAL:
    {
      // no right-hand side needs to be assembled here (with a few exceptions)
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          sqMat.resize(1);
          sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
          reMat.resize(0);
          break;
        case 3:
        case 4:
          sqMat.resize(d);
          for(int i = 0; i < d; ++i)
            sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks.at(i*(d+2)).get());

          reMat.resize(0);
          if(db["space_discretization_type"].is("smagorinsky"))
          {
            sqMat.resize(d*d);
            rhs_array.resize(0);
            for(int i = 0, j = 0; i < d*d; ++i, ++j)
            {
              if(i%d == 0 && i > 0)
                j++;
              sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
            }
          }  
          // In the case of SUPG: together with the other contributions to 
          // the viscous and nonlinear terms, additional Mass matrix, BT-block, 
          // and right-had-side needs to be assembled during the nonlinear 
          // iteration due to the weighted test function.
          if(db["space_discretization_type"].is("supg"))
          {
            sqMat.resize(d+1);
            sqMat[d] = reinterpret_cast<SquareMatrixD*>(mass_blocks[0].get());
            reMat.resize(d);
            for(int i = 0; i < d; ++i)
              reMat[i] = reinterpret_cast<MatrixD*>(blocks.at(i*(d+1)+d).get());
            for(int i = 0; i < d; ++i)
              rhs_array[i] = s.rhs.block(i);
            s.rhs.reset(); // reset to zero
          }
          break;
        case 14:
          if(!db["space_discretization_type"].is("supg")
            || !db["space_discretization_type"].is("local_projection"))
          {
            ErrThrow("NSTYPE 14 only supports SUPG and Local Projection ");
          }
          // we need to re-assemble all the matrices due to the solution
          // dependency of the stabilization parameters
          sqMat.resize(d*d);
          for(int i = 0, j = 0; i < d*d; ++i, ++j)
          {
            if(i%d == 0 && i > 0)
              j++;
            sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
          }
          
          if(db["space_discretization_type"].is("supg") )
          {
            sqMat.resize(d*d+2);
            sqMat[d*d] = reinterpret_cast<SquareMatrixD*>(mass_blocks.at(0).get());
            // pressure-pressure block
            sqMat[d*d+1] = reinterpret_cast<SquareMatrixD*>(blocks.at(d*d+2*d).get());
          }
          
          reMat.resize(2*d);
          for(int i = 0; i < d; ++i) // first the lying B blocks
            reMat[i] = reinterpret_cast<MatrixD*>(blocks.at(d*d+d+i).get());
          for(int i = 0; i < d; ++i) // than the standing B blocks
            reMat[d+i] = reinterpret_cast<MatrixD*>(blocks.at(i*(d+1)+d).get());
          
          for(int i = 0; i < d; ++i)
            rhs_array[i] = s.rhs.block(i);
          s.rhs.reset();
          break;
      }// endswitch NSTYPE
      break;
    }// endswitch TNSE2D_NL
    break;
//---------------------------    
    case LocalAssembling_type::TNSE3D_Rhs:
    {
      // no matrices to be assembled
      sqMat.resize(0);
      reMat.resize(0);
      for(int i = 0; i < d+1; ++i)
        rhs_array[i] = s.rhs.block(i);
      s.rhs.reset();
      break;
    }
    default:
      ErrThrow("The assembling type ", type, 
               " is unknown to TimeNavierStokes.");
  }
  // reset matrices
  for(auto sm : sqMat)
    sm->reset();
  for(auto rm : reMat)
    rm->reset();
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::restrict_function()
{
  // assembling requires an approximate velocity solution on every grid
  for( int block = 0; block < d; ++block)
  {
    std::vector<const FESpace*> spaces;
    std::vector<double*> u_entries;
    std::vector<size_t> u_ns_dofs;
    for(auto &s : systems)
    {
      spaces.push_back(s.velocity_space.get());
      u_entries.push_back(s.solution.block(block));
      u_ns_dofs.push_back(s.solution.length(block));
    }
    GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
  }
}

/* ************************************************************************** */
template <int d>
bool TimeNavierStokes<d>::stop_it(unsigned int it_counter)
{
  const double old_norm_of_residual = this->get_full_residual();
  compute_residuals();
  
  bool i_am_root = true;
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == 0);
#endif
  const double norm_of_residual = this->get_full_residual();
  const double impulse_residual = this->get_impuls_residual();
  const double mass_residual = this->get_mass_residual();
  
  // some output:
  if(i_am_root)
  {
    Output::print<3>("nonlinear step  : " , setw(3), it_counter);
    Output::print<3>("impulse_residual: " , setw(12), impulse_residual);
    Output::print<3>("mass_residual   : " , setw(12), mass_residual);
    Output::print<3>("full residual   : " , setw(12), norm_of_residual);
  }

  if(it_counter == 0)
    initial_residual = norm_of_residual;
  else
    Output::print<3>("rate:           : " , setw(12),
                     norm_of_residual/old_norm_of_residual);

  // check if minimum number of iterations was performed already
  size_t min_it = db["nonlinloop_minit"];
  if(it_counter < min_it)
    return false;
  
  // hold the residual from 10 iterations ago
  const double very_old_norm_of_residual = old_residuals.front().fullResidual;
  size_t max_it = db["nonlinloop_maxit"];
  double convergence_speed = db["nonlinloop_slowfactor"];
  double limit = db["nonlinloop_epsilon"];
  bool slow_convergence = false;
  
  // check various stopping criteria
  
  if(norm_of_residual >= convergence_speed * very_old_norm_of_residual)
  {
    slow_convergence = true;
    Output::print("SLOW!! ", norm_of_residual);
  }
  
  if(db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= std::sqrt(this->get_size());
    if(i_am_root)
      Output::print("stopping tolerance for nonlinear iteration ", limit);
  }
  if(norm_of_residual <= limit || it_counter == max_it || slow_convergence)
  {
    reset_residuals();
    for(System_per_grid& s: this->systems)
    {
      s.solution_m2 = s.solution_m1;
      s.solution_m1 = s.solution;
    }
    // descale the matrices, since only the diagonal A block will
    // be reassembled in the next time step
    for(System_per_grid & s : this->systems)
    {
      time_stepping_scheme.reset_linear_matrices(s.matrix, s.mass_matrix);
      // descale if it's rescaled at the next time step for bdf schemes
      time_stepping_scheme.scale_descale_all_b_blocks(s.matrix, "descale");
    }
    return true;
  }
  else
    return false;
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::compute_residuals()
{
  System_per_grid& s = this->systems.front();
  unsigned int number_u_Dof = s.solution.length(0);
  unsigned int number_p_Dof = s.solution.length(d);

#ifdef _MPI
    //MPI: put solution in consistency level 3
    auto comms = s.matrix.get_communicators();
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(s.solution.block(bl), 3);
    }
#endif

  // copy rhs to defect and compute defect
  this->defect = rhs_from_time_disc;
  s.matrix.apply_scaled_add(s.solution, defect, -1.);

  if(s.matrix.pressure_projection_enabled())
  {
    FEFunction defect_fctn(s.pressure_space.get(), "p_defect",
                           "pressure defect function", &defect[d*number_u_Dof],
                           number_p_Dof);
    defect_fctn.project_into_L20();
  }

  // This is the calculation of the residual, given the defect.
  std::vector<unsigned int> velocity_blocks(d, 0);
  std::iota(std::begin(velocity_blocks), std::end(velocity_blocks), 0);
#ifdef _MPI
  std::vector<const TParFECommunicator3D*> velocity_comms(d, nullptr);
  for(int i = 0; i < d; ++i)
    velocity_comms[i] = comms[i];
  double impuls_residual_square = defect.norm(velocity_blocks, velocity_comms);
  double mass_residual_square = defect.norm({d}, {comms[d]});
#else
  double impuls_residual_square = defect.norm(velocity_blocks);
  double mass_residual_square = defect.norm({d});
#endif

  // the struct 'Residuals' takes squares:
  impuls_residual_square *= impuls_residual_square;
  mass_residual_square *= mass_residual_square;

  Residuals current_residuals(impuls_residual_square, mass_residual_square);
  old_residuals.add(current_residuals);
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::solve()
{
  System_per_grid& s = this->systems.front();
  // store previous solution for damping, it is a pointer so that we can avoid
  // the copy in case of no damping
  double damping = this->db["nonlinloop_damping_factor"];
  std::shared_ptr<BlockVector> old_solution(nullptr);
  if(damping != 1.0)
    old_solution = std::make_shared<BlockVector>(s.solution);
  // DEBUG CODE
  // #ifdef _MPI
  //   std::vector<const TParFECommunicator3D*> comms = s.matrix.get_communicators();
  //   Output::print("norm of rhs ", rhs_from_time_disc.norm(comms));
  //   Output::print("norm of sol before ", s.solution.norm(comms));
  //   Output::print("norm of matrix ", s.matrix.get_combined_matrix()->GetNorm(-2));
  // #else
  //   Output::print("norm of rhs ", rhs_from_time_disc.norm());
  //   Output::print("norm of sol before ", s.solution.norm());
  //   Output::print("norm of matrix ", s.matrix.get_combined_matrix()->GetNorm(-2));
  // #endif
  // END DEBUG CODE
  
#ifdef _MPI
  if(solver.get_db()["solver_type"].is("direct"))
  {
    MumpsWrapper mumps_wrapper(s.matrix);
    mumps_wrapper.solve(rhs_from_time_disc, s.solution);
  }
  else
#endif
  {
    solver.solve(s.matrix, rhs_from_time_disc, s.solution);
  }
  // DEBUG CODE
  // #ifdef _MPI
  //   Output::print("norm of sol after ", s.solution.norm(comms));
  // #else
  //   Output::print("norm of sol after ", s.solution.norm());
  // #endif
  // END DEBUG CODE
  
  if(damping != 1.0)
  {
    s.solution.scale(damping);
    s.solution.add_scaled(*old_solution, 1. - damping);
  }

  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11 and A22 matrices are
  // reset and assembled again but the A12 and A21 are scaled, so
  // for the next iteration we have to descale
  for(System_per_grid & s : this->systems)
    time_stepping_scheme.reset_linear_matrices(s.matrix, s.mass_matrix);
  
  if(s.matrix.pressure_projection_enabled())
       s.p.project_into_L20();
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::output()
{
  bool i_am_root = true;
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == 0);
#endif
  System_per_grid& s = this->systems.front();
  std::array<FEFunction*, d> velocity_components;
  for(int i = 0; i < d; ++i)
    velocity_components[i] = s.u.GetComponent(i);

  if((size_t)db["verbosity"] > 1)
  {
    for(int i = 0; i < d; ++i)
      velocity_components[i]->PrintMinMax();
    s.p.PrintMinMax();
  }
  
  double t = time_stepping_scheme.current_time_;

  if(db["output_compute_errors"])
  {
    std::vector<std::array<double, 5>> computed_errors(d+1);
#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D nsAllDerivs[d+1] = {D000, D100, D010, D001};
#else
    TAuxParam2D aux;
    MultiIndex2D nsAllDerivs[d+1] = {D00, D10, D01};
#endif
    const FESpace *velocity_space = &this->get_velocity_space();
    const FESpace *pressure_space = &this->get_pressure_space();
    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    
    // errors in the velocity components
    for(int i = 0; i < d; ++i)
    {
      auto ui = velocity_components[i];
      ui->GetErrors(example.get_exact(i), d+1, nsAllDerivs, d, L2H1Errors,
                    nullptr, &aux, 1, &velocity_space,
                    computed_errors[i].data());
    }
    // error in divergence
    double div = s.u.GetL2NormDivergenceError(example.get_exact(0),
                                              example.get_exact(1)
#ifdef __3D__
                                              , example.get_exact(2)
#endif
                                             );
    // errors in pressure
    s.p.GetErrors(example.get_exact(d), d+1, nsAllDerivs, d, L2H1Errors,
                  nullptr, &aux, 1, &pressure_space, computed_errors[d].data());
    
#ifdef _MPI
    int n_send = 2*(d+1)+1;
    double err_red[n_send]; //memory for global (across all processes) error
    double err_send[n_send]; //fill send buffer
    err_send[0] = computed_errors[0][0];
    err_send[1] = computed_errors[0][1];
    err_send[2] = computed_errors[1][0];
    err_send[3] = computed_errors[1][1];
    err_send[4] = computed_errors[2][0];
    err_send[5] = computed_errors[2][1];
    err_send[6] = div*div;
    err_send[7] = computed_errors[d][0];
    err_send[8] = computed_errors[d][1];

    MPI_Allreduce(err_send, err_red, n_send, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    for(int i = 0; i < n_send; i++)
    { //MPI: sqrt was skipped in GetErrors function - do it here globally!
      err_red[i] = std::sqrt(err_red[i]);
    }
    //fill the reduced errors back where they belong
    computed_errors[0][0] = err_red[0];
    computed_errors[0][1] = err_red[1];
    computed_errors[1][0] = err_red[2];
    computed_errors[1][1] = err_red[3];
    computed_errors[2][0] = err_red[4];
    computed_errors[2][1] = err_red[5];
    div = err_red[6];
    computed_errors[d][0] = err_red[7];
    computed_errors[d][1] = err_red[8];
#endif
    
    double l2_u = 0;
    double h1_u = 0;
    for(int i = 0; i < d; ++i)
    {
      l2_u += computed_errors[i][0] * computed_errors[i][0];
      h1_u += computed_errors[i][1] * computed_errors[i][1];
    }
    l2_u = std::sqrt(l2_u);
    h1_u = std::sqrt(h1_u);
    double l2_p = computed_errors[d][0];
    double h1_p = computed_errors[d][1];
    
    double old_l2_u = errors[0];
    double old_div = errors[1];
    double old_h1_u = errors[2];
    double old_l2_p = errors[3];
    double old_h1_p = errors[4];
    
    errors[0] = l2_u;
    errors[1] = div;
    errors[2] = h1_u;
    errors[3] = l2_p;
    errors[4] = h1_p;
    
    errors[5] += (l2_u*l2_u + old_l2_u*old_l2_u) * tau*0.5;
    errors[6] += (div*div + old_div*old_div) * tau*0.5;
    errors[7] += (h1_u*h1_u + old_h1_u*old_h1_u) * tau*0.5;
    errors[8] += (l2_p*l2_p + old_l2_p*old_l2_p) * tau*0.5;
    errors[9] += (h1_p*h1_p + old_h1_p*old_h1_p) * tau*0.5;
    
    if(i_am_root)
    {
      using namespace Output;
      stat("TimeNavierStokes", "Measured errors");
      dash(t, " L2(u)              : ", setprecision(14), errors[0]);
      dash(t, " L2(div(u))         : ", setprecision(14), errors[1]);
      dash(t, " H1-semi(u)         : ", setprecision(14), errors[2]);
      dash(t, " L2(p)              : ", setprecision(14), errors[3]);
      dash(t, " H1-semi(p)         : ", setprecision(14), errors[4]);
      dash(t, " L2(0,t,L2(u))      : ", setprecision(14), std::sqrt(errors[5]));
      dash(t, " L2(0,t,L2(div(u))) : ", setprecision(14), std::sqrt(errors[6]));
      dash(t, " L2(0,t,H1-semi(u)) : ", setprecision(14), std::sqrt(errors[7]));
      dash(t, " L2(0,t,L2(p))      : ", setprecision(14), std::sqrt(errors[8]));
      dash(t, " L2(0,t,H1-semi(p)) : ", setprecision(14), std::sqrt(errors[9]));
    }
  }

  if(db["output_compute_time_average"])
  {
    this->time_averaging();
    this->t_avg_output();
  }
#ifdef __3D__
  if(db["output_along_line"])
  {
    // fill a vector with all fe functions to be evaluated using this->Lines
    std::vector<const FEFunction*> fe_functions(velocity_components.begin(),
                                                velocity_components.end());
    fe_functions.push_back(&s.p);
    Lines.write_fe_values(fe_functions, t, "solution");
  }
#endif

  example.do_post_processing(*this);

  for(int i = 0; i < d; ++i)
    delete velocity_components[i];
  
  outputWriter.write(t);

  if(db["write_solution_binary"].is(true))
  {
    size_t interval = db["write_solution_binary_all_n_steps"];
    if(time_stepping_scheme.current_step_ % interval == 0)
    {
      //write solution to a binary file
      std::string file = db["write_solution_binary_file"];
      if(!db["overwrite_solution_binary"]) // create a new file every time
      {
        file += ".";
        file += std::to_string(t);
#ifdef _MPI
        file += ".proc" + std::to_string(my_rank);
#endif
      }
      Output::info("output", "Writing current solution to file ", file);
      systems.front().solution.write_to_file(file);
    }
  }
}

/* ************************************************************************** */
template<int d>
std::array<double, int(10)> TimeNavierStokes<d>::get_errors() const
{
  return this->errors;
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::output_problem_size_info() const
{
  int i_am_root = true;
#ifndef _MPI
  auto & velocity_space = *this->systems.front().velocity_space;
  auto & pressure_space = *this->systems.front().pressure_space;

  size_t nDofu = d*velocity_space.GetN_DegreesOfFreedom();
  size_t nDofp = pressure_space.GetN_DegreesOfFreedom();
  size_t nTotal = nDofu + nDofp;
  size_t nActive = d*velocity_space.GetActiveBound();

  TCollection* coll = velocity_space.GetCollection();
  double hmin, hmax;
  coll->GetHminHmax(&hmin, &hmax);
  int n_cells = coll->GetN_Cells();
#else
  int my_rank;
  int root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == root);
  auto velocity_comm = systems.front().velocity_space->get_communicator();
  auto pressure_comm = systems.front().pressure_space->get_communicator();
  int nDofu  = d*velocity_comm.get_n_global_dof();
  int nDofp  = pressure_comm.get_n_global_dof();
  int nTotal = nDofu + nDofp;
  
  TCollection* coll = systems.front().velocity_space->GetCollection();
  int n_local_master_cells = coll->GetN_OwnCells();
  int n_cells;
  MPI_Reduce(&n_local_master_cells, &n_cells, 1, MPI_DOUBLE, MPI_SUM, root,
             MPI_COMM_WORLD);
  double local_hmin, local_hmax;
  coll->GetHminHmax(&local_hmin, &local_hmax);
  double hmin, hmax;
  MPI_Reduce(&local_hmin, &hmin, 1, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
  MPI_Reduce(&local_hmax, &hmax, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
#endif

  if(i_am_root)
  {
    Output::print("N_Cells     : ", setw(10), n_cells);
    Output::print("h (min,max) : ", setw(12), hmin ," ", setw(12), hmax);
    Output::print("dof Velocity: ", setw(10), nDofu);
    Output::print("dof Pressure: ", setw(10), nDofp);
    Output::print("dof all     : ", setw(10), nTotal);
#ifndef _MPI
    Output::print("active dof  : ", setw(10), nActive);
#endif
  }
}


/* ************************************************************************** */
template <int d>
bool TimeNavierStokes<d>::imex_scheme()
{
  if((db["imex_scheme_"] && time_stepping_scheme.current_step_ >= 3))
  {
     // condition is here just to print it once
    if((size_t)db["nonlinloop_maxit"] > 1)
    {
      Output::info<1>("Nonlinear Loop MaxIteration",
                      "The parameter 'nonlinloop_maxit' was changed to 1."
                      " Only one non-linear iteration is done, because the "
                      "IMEX scheme was chosen.");
    }
    db["nonlinloop_maxit"] = 1;
    db["extrapolate_velocity"] = true;
    
    // this is typical for the problem with slip type boundary condition
    // NOTE: only tested for the mixing layer problem for the moment
    //if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
    //  is_rhs_and_mass_matrix_nonlinear = true;
    if(db["space_discretization_type"].is("supg"))
      space_disc_global = -2;
    return true;
  }
  else
    return false;
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::modify_slip_bc(bool BT_Mass, bool slip_A_nl)
{
  // modification of the matrices due to the
  // slip type boundary conditions: If the mass matrices,
  // the off-diagonal A-blocks , and the BT's block,
  // are unchanged during the time iteration, then this modification
  // is done only once in the time loop. However, in the SUPG
  // and residual based VMS method these matrices are also
  // updated during the time steps, so modification of all
  // of them including the right-hand side is necessary.The
  // modification of the diagonal A-blocks are necessary
  // in any case.
  if(TDatabase::ParamDB->NSTYPE < 4)
  {
    ErrThrow("Slip with friction b.c. is only implemented for NSTYPE 4 and 14");
  }
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  std::vector<const FESpace*> spaces_mat(1);
  std::vector<double*> rhs_array(d);
  std::vector<const FESpace*> rhs_space(d);

  for(System_per_grid& s: this->systems)
  {
    spaces_mat[0] = s.velocity_space.get();
    for(int i = 0; i < d; ++i)
      rhs_space[i] = spaces_mat[0];

    for(int i = 0; i < d; ++i)
      rhs_array[i] = s.rhs.block(i);

    auto blocks = s.matrix.get_blocks_uniquely();

    auto mass_blocks = s.mass_matrix.get_blocks_uniquely(true);
    
    
    std::array<BoundaryConditionFunction*, d+1> boundCondition;
    std::array<BoundaryValuesFunction*, d+1> boundValues;
    for(int i = 0; i < d; ++i)
      boundCondition[i] = s.velocity_space->get_boundary_condition();
    boundCondition[d] = s.pressure_space->get_boundary_condition();
    for(int i = 0; i < d+1; ++i)
      boundValues[i] = example.get_bd(i);

    std::vector<SquareMatrixD*> sqMat;
    std::vector<MatrixD*> reMat;
    sqMat.resize(d);
    
   /* all d*d A blocks at the first time step
    * ------------------
    * a11 a12 a12 b1t
    * a21 a22 a23 b2t
    * a31 a32 a33 b3t
    * b1  b2  b3  c
    * and only the first 2 within the nonlinear loop */
    for(int i = 0; i < d; ++i)
      sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[i*(d+2)].get()); // aii

    // if the off-diagonal are not changing within the non-linear loop
    // then dont need to assemble them again
    if(slip_A_nl)
    {
      sqMat.resize(d*d);
#ifdef __2D__
      sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());//a12
      sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());//a21
#else // __3D__
      sqMat[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());//a12
      sqMat[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());//a13
      sqMat[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());//a21
      sqMat[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());//a23
      sqMat[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());//a31
      sqMat[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());//a32
#endif   
    }

    // either at the first time step
    // or every time step if M and B's are changing
    reMat.resize(0);
    if(BT_Mass)
    {
#ifdef __2D__
      sqMat.resize(8);
      sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());//m11
      sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());//m22
      sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(1).get());//m12
      sqMat[7] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(3).get());//m21
      reMat.resize(2);
      reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
      reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
#else // __3D__
      sqMat.resize(18);
      sqMat[9]  = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get()); //m11
      sqMat[10] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(5).get()); //m22
      sqMat[11] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(10).get());//m33
  
      sqMat[12] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(1).get());//m12
      sqMat[13] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(2).get());//m13
      sqMat[14] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(4).get());//m21
      sqMat[15] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(6).get());//m23
      sqMat[16] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(8).get());//m31
      sqMat[17] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(9).get());//m32
      reMat.resize(3);
      reMat[0] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get()); // b1t
      reMat[1] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get()); // b2t
      reMat[2] = reinterpret_cast<TMatrix3D*>(blocks.at(11).get());// b3t
#endif  
    }

    auto u1 = s.u.GetComponent(0);
    auto u2 = s.u.GetComponent(1);
    // update the matrices and right hand side
#ifdef __2D__
    Assemble2DSlipBC(
#else // __3D__
    Assemble3DSlipBC(
#endif      
      spaces_mat.size(), spaces_mat.data(), sqMat.size(), sqMat.data(), 
      reMat.size(), reMat.data(), rhs_array.size(), rhs_array.data(), 
      rhs_space.data(), boundCondition.data(), boundValues.data()
#ifdef __2D__
      , u1, u2);
#else
    );
#endif
    
    delete u1;
    delete u2;
  }
}

/* ************************************************************************* */
template <int d>
typename TimeNavierStokes<d>::FEFunction*
  TimeNavierStokes<d>::get_velocity_component(int i)
{
  if(i == 0)
    return this->systems.front().u.GetComponent(0);
  else if(i == 1)
    return this->systems.front().u.GetComponent(1);
  else if(i == 2 && d == 3)
    return this->systems.front().u.GetComponent(2);
  else
    ErrThrow("There are only ", d, " velocity components!");
}


/* ************************************************************************* */
template <int d>
const Residuals& TimeNavierStokes<d>::get_residuals() const
{
  return old_residuals.back();
}

/* ************************************************************************* */
template <int d>
double TimeNavierStokes<d>::get_impuls_residual() const
{
  return old_residuals.back().impulsResidual;
}

/* ************************************************************************* */
template <int d>
double TimeNavierStokes<d>::get_mass_residual() const
{
  return old_residuals.back().massResidual;
}

/* ************************************************************************* */
template <int d>
double TimeNavierStokes<d>::get_full_residual() const
{
  return old_residuals.back().fullResidual;
}

/* ************************************************************************* */
template <int d>
void TimeNavierStokes<d>::reset_residuals()
{
  this->old_residuals = FixedSizeQueue<10, Residuals>();
}

/* ************************************************************************* */
template <int d>
void TimeNavierStokes<d>::update_matrices_lps(System_per_grid &s)
{
#ifdef __2D__
  std::vector<std::shared_ptr<FEMatrix>> blocks;
  blocks = s.matrix.get_blocks_uniquely();
  if(TDatabase::ParamDB->NSTYPE == 3 || TDatabase::ParamDB->NSTYPE == 4)
  {
    //update matrices for local projection stabilization
    std::vector<SquareMatrixD*> sqMat(2);
    sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
    sqMat[1] = reinterpret_cast<SquareMatrixD*>(blocks.at(4).get());
    UltraLocalProjection(sqMat[0], false);
    UltraLocalProjection(sqMat[1], false);
  }
  else
  {
    std::vector<SquareMatrixD*> sqMat(1);
    sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
    UltraLocalProjection(sqMat[0], false);
  }
#else // 2D->3D
  ErrThrow("TimeNavierStokes<d>::update_matrices_lps is only implemented in "
           "2D.");
#endif
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::time_averaging()
{
#ifdef __3D__
  double t           = time_stepping_scheme.current_time_;
  double tau         = TDatabase::TimeDB->TIMESTEPLENGTH;
  double t0          = db["time_start"];
  double t0_avg      = db["start_time_averaging_at"];
  System_per_grid& s = this->systems.front();

  if( t == t0 )
  {
    return; // in case of restart (i.e. continue_output_after_restart)
  }
  else if( t-tau >= t0_avg )
  {
    s.time_avg_sol.scale(t - t0_avg - tau);
    s.time_avg_sol.add_scaled(s.solution_m2, tau/2.);
    s.time_avg_sol.add_scaled(s.solution, tau/2.);
    s.time_avg_sol.scale(1./(t - t0_avg));
  }
#endif
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::t_avg_output()
{  
#ifdef __3D__
  // write the time averaged solution only for the last time step
  if(time_stepping_scheme.reached_final_time_step())
  {
    double t           = time_stepping_scheme.current_time_;
    System_per_grid& s = this->systems.front();
    if(db["output_along_line"])
    {
      std::vector<const FEFunction*> fe_functions;
      for( int i = 0 ; i < d ; ++i )
      {
        fe_functions.push_back(s.u_time_avg.GetComponent(i));
      }
      fe_functions.push_back(&s.p_time_avg);
      Lines.write_fe_values(fe_functions, t, "time_average");
      for( int i = 0 ; i < d ; ++i )
      {
        delete fe_functions[i];
      }
    }
    
    // create separate DataWrite object for writing the time average values, 
    // because it would not work with the case format in this->outputWriter.
    ParameterDatabase tmp(this->db);
    std::string basename = tmp["output_basename"];
    tmp["output_basename"].set(basename + "_time_average", false);
    tmp["output_write_case"] = false;
    if(tmp["output_write_vtk"].is(false) && tmp["output_write_vtu"].is(false))
    {
      tmp["output_write_vtu"] = true;
    }
    DataWriter<d> avg_writer(tmp);
    avg_writer.add_fe_vector_function(&s.u_time_avg);
    avg_writer.add_fe_function(&s.p_time_avg);
    avg_writer.write();

    //write solution to a binary file
    if(tmp["write_solution_binary"].is(true))
    {
      std::string file = tmp["write_solution_binary_file"];
      file += "_time_average.";
      file += std::to_string(t);
#ifdef _MPI
      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      file += ".proc" + std::to_string(my_rank);
#endif
      Output::info("output", "Writing time averaged solution to file ", file);
      s.time_avg_sol.write_to_file(file);
    }
  }
#endif /** #ifdef __3D__ */
} 

#ifdef __3D__
template class TimeNavierStokes<3>;
#else
template class TimeNavierStokes<2>;
#endif

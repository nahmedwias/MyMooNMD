#include "MultiPhase.h"
#include "Database.h"
#include "Example_TnseTcd2D.h"
#include "ParameterDatabase.h"
#include "Solver.h"
#include "TimeDiscretizations.h"
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

template <int d>
ParameterDatabase MultiPhase<d>::default_coupled_database()
{
  Output::print<5>("creating a default TNSE parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default tnse database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Coupled NS-CD parameter database");
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(ParameterDatabase::default_solution_in_out_database());
  db.merge(LocalAssembling<d>::default_local_assembling_database());
  return db;
}


template<int d> 
MultiPhase<d>::System_per_grid::System_per_grid(const Example_TnseTcd2D& example, TCollection& coll, std::tuple<int, int, int> order) 
 : velocity_space(new FESpace(&coll, "u", "u", example.get_bc(0), std::get<0>(order))),
   pressure_space(new FESpace(&coll, "p", "p", example.get_bc(d), std::get<1>(order))),
   temp_space(new FESpace(&coll, "c", "c", example.get_bc(0), std::get<2>(order)))
{
  switch(TDatabase::ParamDB->NSTYPE)
  {
 #ifdef __2D__
    case 4:
      matrix_NS = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      mass_matrixNS = BlockFEMatrix::Mass_NSE2D_Type4(velocity_space, pressure_space);
      break;
    case 14:
      matrix_NS = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
      mass_matrixNS = BlockFEMatrix::Mass_NSE2D_Type4(velocity_space, pressure_space);
      break;
#else
    case 4:
      matrix_NS = BlockFEMatrix::NSE3D_Type4(velocity_space, pressure_space);
      mass_matrixNS = BlockFEMatrix::Mass_NSE3D_Type4(velocity_space, pressure_space);
      break;
    case 14:
      matrix_NS = BlockFEMatrix::NSE3D_Type14(velocity_space, pressure_space);
      mass_matrixNS = BlockFEMatrix::Mass_NSE3D_Type4(velocity_space, pressure_space);
      break;
#endif
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }
  
  #ifdef __3D__
  stiffness_matrixtemp = BlockFEMatrix::CD3D(temp_space);
  mass_matrixtemp = BlockFEMatrix::CD3D(temp_space);
#else
  stiffness_matrixtemp = BlockFEMatrix::CD2D(temp_space);
  mass_matrixtemp = BlockFEMatrix::CD2D(temp_space);
#endif

  rhs_temp = BlockVector(stiffness_matrixtemp, true);
  solution_temp = BlockVector(stiffness_matrixtemp, false);
  temp = FEFunction(temp_space, (char*)"c", (char*)"c", solution_temp.block(0), solution_temp.length(0)); 
  
  solution_temp_m1 = BlockVector(stiffness_matrixtemp, false);
  temp_u_m1 = FEFunction(temp_space,"um1","um1", solution_temp_m1.get_entries(),
                    solution_temp_m1.length());
  solution_temp_m2 = BlockVector(stiffness_matrixtemp, false);
  temp_u_m2 = FEFunction(temp_space,"um2","um2", solution_temp_m2.get_entries(),
                    solution_temp_m2.length());
 
  rhs_NS = BlockVector(matrix_NS, true);
  solution_NS = BlockVector(matrix_NS, false);
  time_avg_sol = BlockVector(matrix_NS, false);
  combined_old_sols_NS = BlockVector(matrix_NS, false);
  extrapolate_sol_NS = BlockVector(matrix_NS, false);
  
  solution_NS_m1 = BlockVector(matrix_NS, false);
  solution_NS_m2 = BlockVector(matrix_NS, false);

  u = FEVectFunct(velocity_space, "u", "u", solution_NS.block(0),
                  solution_NS.length(0), d);
  p = FEFunction(pressure_space, "p", "p", this->solution_NS.block(d),
                 solution_NS.length(d));
  u_time_avg = FEVectFunct(velocity_space, "u_t_avg", "u time averaged",
                           time_avg_sol.block(0), time_avg_sol.length(0), d);
  p_time_avg = FEFunction(pressure_space, "p_t_avg", "p time averaged",
                          this->time_avg_sol.block(d), time_avg_sol.length(d));
  
  u_m1 = FEVectFunct(velocity_space, "u", "u", solution_NS_m1.block(0),
                     solution_NS_m1.length(0), d);
  p_m1 = FEFunction(pressure_space, "p", "p", this->solution_NS_m1.block(d),
                    solution_NS_m1.length(d));

  u_m2 = FEVectFunct(velocity_space,"u","u", solution_NS_m2.block(0),
                     solution_NS_m2.length(0), d);
  p_m2 = FEFunction(pressure_space, "p", "p", this->solution_NS_m2.block(d),
                    solution_NS_m2.length(d));
  
  comb_old_u = FEVectFunct(velocity_space, "u", "u",
                           combined_old_sols_NS.block(0),
                           combined_old_sols_NS.length(0), d);
  extrapolate_u = FEVectFunct(velocity_space, "u", "u",
                              extrapolate_sol_NS.block(0),
                              extrapolate_sol_NS.length(0), d);
    
#ifdef _MPI
  //print some information
  velocity_space->get_communicator().print_info();
  pressure_space->get_communicator().print_info();
  temp_space->get_communicator().print_info();
#endif
}


template <int d>
MultiPhase<d>::System_per_grid::System_per_grid(
  const System_per_grid& other)
 : velocity_space(other.velocity_space), pressure_space(other.pressure_space),
   matrix_NS(other.matrix_NS), rhs_NS(other.rhs_NS), solution_NS(other.solution_NS)
{
  u = FEVectFunct(velocity_space, "u", "u", solution_NS.block(0), solution_NS.length(0), d);
  p = FEFunction(pressure_space, "p", "p", this->solution_NS.block(d), solution_NS.length(d));
}


template<int d> 
MultiPhase<d>::MultiPhase(const TDomain& domain, const ParameterDatabase& param_db)
 : MultiPhase<d>(domain, param_db, Example_TnseTcd2D(param_db))
{
}




template<int d> 
MultiPhase<d>::MultiPhase(const TDomain& domain, const ParameterDatabase& param_db, 
				const Example_TnseTcd2D& ex)
: db(default_coupled_database()), solver(param_db), systems(), example(ex),
   old_rhs_temp(), errors_temp({}), outputWriter(param_db),
   time_stepping_scheme(param_db), rhs_temp_from_time_disc(),
   defect(), old_residuals(), initial_residual(1e10),
    is_rhs_and_mass_matrix_nonlinear(false),
   Lines()
{

  this->db.merge(param_db, false); // update this database with given values
  this->check_and_set_parameters();
  
  std::tuple<int,int,int> velo_pres_order(TDatabase::ParamDB->ANSATZ_ORDER,
                                          TDatabase::ParamDB->PRESSURE_SPACE,
                                          TDatabase::ParamDB->VELOCITY_SPACE);
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
      velo_pres_order = {-1, -4711, 0};
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
    matrices.push_back(&systems.back().matrix_NS);
    for(auto coll : collections)
    {
      systems.emplace_back(example, *coll, velo_pres_order);
      // prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix_NS);
    }
    // initialize the multigrid object with all the matrices on all levels
    mg->initialize(matrices);
  }
  else
    
  {
    // the given collection for particular cell
    TCollection& cellCollection = *collections.front();
    int ansatz_order = TDatabase::ParamDB->ANSATZ_ORDER;
    systems.emplace_back(example, cellCollection, velo_pres_order);
    // initial condition on the solution
    this->systems.front().temp.Interpolate(example.get_initial_cond(0));
  }
  
  // initial solution on finest grid - read-in or interpolation
  if(db["read_initial_solution"].is(true))
  {
    if(this->time_stepping_scheme.get_start_time() == 0.)
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
    Output::root_info("Initial Solution", "Reading initial solution from file ",
                      file);
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    file += ".proc" + std::to_string(my_rank);
    Output::root_info("Initial Solution", "Appending .proc<RANK> to the "
                      "expected initial solution file name.");
#endif
    systems.front().solution_NS.read_from_file(file);
    
    
    // read time average solution if it exists
    if(db["output_compute_time_average"])
    {
      Output::root_info("Time averaged solution",
                        "Appending _time_average.proc<RANK> to the expected "
                        "initial solution file name.");
      std::string file_average = db["initial_solution_file"];
      file_average += "_time_average";
#ifdef _MPI
      file_average += ".proc" + std::to_string(my_rank);
#endif
      try
      {
        systems.front().time_avg_sol.read_from_file(file_average);
      }
      catch(...)
      {
        Output::warn("Reading time averaged solution", " Could not open a "
                     "file ", file_average, ". Now I start with a zero time "
                     "average.");
      }
    }
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

  this->defect.copy_structure(this->systems.front().rhs_NS);
  this->rhs_NS_from_time_disc.copy_structure(this->systems.front().rhs_NS);
  this->rhs_temp_from_time_disc.copy_structure(this->systems.front().rhs_temp);
  
  outputWriter.add_fe_function(&this->systems.front().temp);
  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  if(db["output_compute_time_average"])
  {
    outputWriter.add_fe_vector_function(&systems.front().u_time_avg);
    outputWriter.add_fe_function(&systems.front().p_time_avg);
  }

  // print out the information (cells, dofs, etc)
  this->output_problem_size_info();
  this->errors_NS.fill(0.);
    errors_temp.fill(0.); // initialize the array
   // initialize L_inf error to some negative number
  errors_temp[MultiPhase::n_errors-1] = -1.;

  if( db["output_along_line"] )
  {
    Lines = LinesEval<d>(domain, param_db);
  }  
}
// 
// template<int d> 
// MultiPhase<d>::MultiPhase(const TDomain& domain, const ParameterDatabase& param_db, const Example_TnseTcd2D& ex)
// {
// }

/* ************************************************************************** */
template<int d>
void MultiPhase<d>::get_velocity_pressure_orders(
  std::tuple<int,int,int> &velocity_pressure_orders)
{
  int velocity_order = std::get<0>(velocity_pressure_orders);
  int pressure_order = std::get<1>(velocity_pressure_orders);
  int ansatz_order = std::get<2>(velocity_pressure_orders);
  int order;
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
  TDatabase::ParamDB->VELOCITY_SPACE = order;
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  TDatabase::ParamDB->ANSATZ_ORDER = ansatz_order;

  std::get<0>(velocity_pressure_orders) = velocity_order;
  std::get<1>(velocity_pressure_orders) = pressure_order;
  std::get<2>(velocity_pressure_orders) = ansatz_order;

  Output::print("velocity space ", setw(6), std::get<0>(velocity_pressure_orders));
  Output::print("pressure space ", setw(6), std::get<1>(velocity_pressure_orders));
  Output::print("temp   space ", setw(6), std::get<2>(velocity_pressure_orders));
}

/* ************************************************************************* */

template<int d>
void MultiPhase<d>::assemble_linear_terms()
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
    std::vector<FEFunction*> fefunctions;
      fefunctions.resize(d);  
      for(int i = 0; i < d; ++i)
        fefunctions[i] = s.u.GetComponent(i);
  
    LocalAssembling<d> la(this->db, la_type, fefunctions.data(),
                          example.get_coeffs());
    call_assembling_routine_ConvDiff(s, la, assemble_both);
  }
  
  // here the modifications due to time discretization begin
  if (!db["algebraic_flux_correction"].is("none") )
  {
    do_algebraic_flux_correction();
    rhs_temp_from_time_disc = this->systems.front().rhs_temp;
    return; // modifications due to time discretization are per-
            // formed inside the afc scheme, so step out here!
  }

  // preparing the right hand side discretized by the used time
  // stepping scheme
  System_per_grid& s = this->systems.front();
  rhs_temp_from_time_disc.reset();
  rhs_temp_from_time_disc = s.rhs_temp;
  // all matrices are available
  unsigned int n_sols = time_stepping_scheme.n_old_solutions();
  std::vector<BlockVector> old_sols_temp(n_sols);
  old_sols_temp[0] = s.solution_temp_m1;
  if(old_sols_temp.size() == 2)
    old_sols_temp[1] = s.solution_temp_m2;
  std::vector<BlockVector> rhs_temp(2);
  rhs_temp[0] = rhs_temp_from_time_disc;
  rhs_temp[1] = old_rhs_temp;
  // prepare the right hand side from the previous time step
  time_stepping_scheme.prepare_rhs_from_time_disc(s.stiffness_matrixtemp, 
                                                  s.mass_matrixtemp, rhs_temp, old_sols_temp);
  rhs_temp_from_time_disc = rhs_temp[0];
  old_rhs_temp = s.rhs_temp;
  rhs_temp_from_time_disc.copy_nonactive(s.rhs_temp);
  
  for(auto &s : this->systems)
    time_stepping_scheme.prepare_system_matrix(s.stiffness_matrixtemp,
                                               s.mass_matrixtemp);
  s.solution_temp.copy_nonactive(s.rhs_temp);


  for(System_per_grid & s : systems)
  {
    // assembling all matrices from NavierStokes model
    this->call_assembling_routine_NS(s, LocalAssembling_type::NSE3D_Linear);
  }
  
}

#ifdef __3D__
template class MultiPhase<3>;
#else
template class MultiPhase<2>;
#endif

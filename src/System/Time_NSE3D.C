#include <Time_NSE3D.h>
#include <Database.h>
#include <Assemble3D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <GridTransfer.h>
#include <Upwind3D.h>
#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!

#include <MainUtilities.h>
#include <Multigrid.h>

#include <BoundaryAssembling3D.h>

#include <sys/stat.h>

#ifdef _MPI
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

/* *************************************************************************** */
  //TODO  So far of this object only the nonlin it stuff is used - switch entirely!
ParameterDatabase get_default_TNSE3D_parameters()
{
  Output::print<5>("creating a default TNSE3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TNSE3D parameter database");

  //NSE3D requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default time database
//   ParameterDatabase time_db = ParameterDatabase::default_time_database();
//   db.merge(time_db,true);
 
  // a default solution in out database
  ParameterDatabase in_out_db = ParameterDatabase::default_solution_in_out_database();
  db.merge(in_out_db,true);
  
  // a default local assembling database
  db.merge(LocalAssembling3D::default_local_assembling_database());

  return db;
}

/* *************************************************************************** */
Time_NSE3D::System_per_grid::System_per_grid(const Example_TimeNSE3D& example,
                  TCollection& coll, std::pair< int, int > order, 
                  Time_NSE3D::Matrix type
#ifdef _MPI
                  , int maxSubDomainPerDof
#endif
)
 : velocitySpace_(&coll, "u", "velocity space",  example.get_bc(0),
                  order.first),
   pressureSpace_(&coll, "p", "pressure space", example.get_bc(3),
                  order.second)
{
  switch(type)
  {
    case Time_NSE3D::Matrix::Type1:
      matrix_ = BlockFEMatrix::NSE3D_Type1(velocitySpace_, pressureSpace_);
      massMatrix_ = BlockFEMatrix::Mass_NSE3D_Type1(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type2:
      matrix_ = BlockFEMatrix::NSE3D_Type2(velocitySpace_, pressureSpace_);
      massMatrix_ = BlockFEMatrix::Mass_NSE3D_Type2(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type3:
      matrix_ = BlockFEMatrix::NSE3D_Type3(velocitySpace_, pressureSpace_);
      massMatrix_ = BlockFEMatrix::Mass_NSE3D_Type3(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type4:
      matrix_ = BlockFEMatrix::NSE3D_Type4(velocitySpace_, pressureSpace_);
      massMatrix_ = BlockFEMatrix::Mass_NSE3D_Type4(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type14:
      matrix_ = BlockFEMatrix::NSE3D_Type14(velocitySpace_, pressureSpace_);
      massMatrix_ = BlockFEMatrix::Mass_NSE3D_Type4(velocitySpace_,pressureSpace_);
      break;
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }

  rhs_ = BlockVector(matrix_, true);
  solution_ = BlockVector(matrix_, false);

  u_ = TFEVectFunct3D(&velocitySpace_, "u", "u", solution_.block(0),
                     solution_.length(0), 3);
  p_ = TFEFunction3D(&pressureSpace_, "p", "p", solution_.block(3),
                    solution_.length(3));

#ifdef _MPI
  velocitySpace_.initialize_parallel(maxSubDomainPerDof);
  pressureSpace_.initialize_parallel(maxSubDomainPerDof);
  //print some information on the parallel infrastructure
  velocitySpace_.get_communicator().print_info();
  pressureSpace_.get_communicator().print_info();
#endif
}

/**************************************************************************** */
Time_NSE3D::Time_NSE3D(std::list< TCollection* > collections_, const ParameterDatabase& param_db, 
                       const Example_TimeNSE3D& ex
#ifdef _MPI
, int maxSubDomainPerDof
#endif  
)
: db_(get_default_TNSE3D_parameters()), outputWriter(param_db),systems_(), example_(ex),
   solver_(param_db), defect_(), old_residual_(), 
   initial_residual_(1e10), errors_(), oldtau_(), 
   time_stepping_scheme(param_db), is_rhs_and_mass_matrix_nonlinear(false),
   current_step_(0)
{
  db_.merge(param_db);
  this->check_and_set_parameters();
  std::pair <int,int>
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // get the velocity and pressure orders
  this->get_velocity_pressure_orders(velocity_pressure_orders);
  
  Time_NSE3D::Matrix type;
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case  1: type = Matrix::Type1;  break;
    case  2: type = Matrix::Type2;  break;
    case  3: type = Matrix::Type3;  break;
    case  4: type = Matrix::Type4;  break;
    case 14: type = Matrix::Type14; break;
    default:
      ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
               " That NSE Block Matrix Type is unknown to class Time_NSE3D.");
  }
  bool usingMultigrid = solver_.is_using_multigrid();
  TCollection *coll = collections_.front(); // finest grid collection
  // create finite element space and function, a matrix, rhs, and solution
#ifdef _MPI
  systems_.emplace_back(example_, *coll, velocity_pressure_orders, type,
                        maxSubDomainPerDof);
#else
  systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);
#endif
  if(usingMultigrid)
  {
    auto mg = this->solver_.get_multigrid();
    size_t n_multigrid_levels = mg->get_n_geometric_levels();
    size_t n_grids=collections_.size();
    if(n_multigrid_levels != n_grids)
    {
      ErrThrow("Wrong number of grids for multigrid! expecting ",
               n_multigrid_levels, " geometric grids but got", n_grids,".");
    }
    
    if(mg->is_using_mdml())
    {
      // change the discretization on the coarse grids to lowest order 
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velocity_pressure_orders={-1, -4711};
      this->get_velocity_pressure_orders(velocity_pressure_orders);
    }
    else
    {
      // for standard multigrid, pop the finest collection - it was already
      // used to construct a space before the "if(usingMultigrid)" clause
      // and will not (as in mdml) be used a second time with a different discretization
      collections_.pop_front();
    }
    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    // matrix on finest grid is already constructed
    matrices.push_back(&systems_.back().matrix_);
    
    for(auto coll : collections_)
    {
#ifdef _MPI
      systems_.emplace_back(example_, *coll, velocity_pressure_orders, type, 
                            maxSubDomainPerDof);
#else
      systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);
#endif
      // prepare input argument for multigrid 
      matrices.push_front(&systems_.back().matrix_);
    }
    // initialize multigrid object with matrices on all levels
    mg->initialize(matrices);
  }
  
  // initial solution on finest grid - read-in or interpolation
  if(db_["read_initial_solution"].is(true))
  {//initial solution is given
    std::string file = db_["initial_solution_file"];
    Output::root_info("Initial Solution", "Reading initial solution from file ", file);
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    file += ".proc" + std::to_string(my_rank);
    Output::root_info("Initial Solution", "Appending .proc<RANK> to the "
        "expected initial solution file name.");
#endif
    systems_.front().solution_.read_from_file(file);
  }
  else
  {//interpolate initial condition from the example
    Output::info("Initial Solution", "Interpolating initial solution from example.");
    for(System_per_grid& s : this->systems_)
    {
      s.u_.GetComponent(0)->Interpolate(example_.get_initial_cond(0));
      s.u_.GetComponent(1)->Interpolate(example_.get_initial_cond(1));
      s.u_.GetComponent(2)->Interpolate(example_.get_initial_cond(2));
    }
  }

  this->output_problem_size_info();
  // initialize the defect of the system. It has the same structure as
  // the rhs (and as the solution)
  this->defect_.copy_structure(this->systems_.front().rhs_);
}

///**************************************************************************** */
void Time_NSE3D::check_and_set_parameters()
{
 // Check problem_type
 if(!db_["problem_type"].is(6))
 {
   if (db_["problem_type"].is(0))
   {
     db_["problem_type"] = 6;
   }
   else
   {
     Output::warn<2>("The parameter problem_type doesn't correspond to Time_NSE."
         "It is now reset to the correct value for Time_NSE (=6).");
     db_["problem_type"] = 6;
   }
 }

 // Tell the user he is using IMEX
//  if(db_["time_discretization"].is(4))
//  {
//    if(solver_.is_using_multigrid())
//    {
//      ErrThrow("Multigrid with IMEX-scheme is not implemented yet");
//    }
//    else
//    {
//      Output::info<1>("check_and_set_parameters",
//                      "The IMEX scheme has been chosen as a time discretization scheme!\n");
//    }
//  }

 if(TDatabase::TimeDB->TIME_DISC == 0)
 {
   ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC
         << " does not supported");
   throw("TIME_DISC: 0 is not supported");
 }
 // standard method
 if(db_["space_discretization_type"].is("galerkin"))
 {
   space_disc_global = 1;
   // set scaling number of blocks that will be scaled by the 
   // factor times the time step length for each time discretization
   // shceme, e.g., 0.5*tau in Crank-Nicolson scheme
   time_stepping_scheme.n_scale_block = 6;
   time_stepping_scheme.b_bt_linear_nl = "linear";
 }
 // the supg case 
 if(db_["space_discretization_type"].is("supg"))
 {
   if(db_["time_discretization"].is("bdf_two"))
   {
     ErrThrow("supg method is only implemented for BDF2 time stepping scheme");
   }
   space_disc_global = 2;
   
   /// set scaling factor for B, BT's block
    // depends on how to deal the nonlinearity in the 
    // test function: fully implicit case
    time_stepping_scheme.b_bt_linear_nl = "nonlinear";
    // inf-sup stable case
    time_stepping_scheme.n_scale_block = 3;
    // equal order case
    if(TDatabase::ParamDB->NSTYPE==14)
      time_stepping_scheme.n_scale_block = 7;
 }
 // Smagorinsky
 if(db_["space_discretization_type"].is("smagorinsky"))
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
 
 // the only case where one have to re-assemble the right hand side
  if(db_["space_discretization_type"].is("supg") && db_["time_discretization"].is("bdf_two"))
  {
    is_rhs_and_mass_matrix_nonlinear = true;
  }
}

/**************************************************************************** */
void Time_NSE3D::get_velocity_pressure_orders(std::pair< int, int > &velocity_pressure_orders)
{
  int velocity_order = velocity_pressure_orders.first;
  int pressure_order = velocity_pressure_orders.second;
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
    case -1: case -2: case -3: case -4: case -5:
    case -101:
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
  velocity_pressure_orders.first = order;
  switch(pressure_order)
  {
    case -4711:
      switch(velocity_order)
      {
        case -1:
        case -2:
        case -3:
        case -4:
          // nonconforming pw (bi)linear velo/ pw constant pressure
          // conforming pw (bi)linear velo/ pw constant pressure (not stable !!!)
          pressure_order = -velocity_order-1;
          break;
        case 1: // discontinuous space
          pressure_order = 0;
          Output::print<1>("Warning: The P1/P0 element pair (Q1/Q0 on hexa) is "
              " not stable. Make sure to use stabilization!");
          break;
        case 2: case 3: case 4: case 5:
        // standard conforming velocity and continuous pressure
          pressure_order = velocity_order-1;
          break;
          // discontinuous pressure spaces
          // standard conforming velo and discontinuous pressure
          // this is not stable on triangles !!!
        case 12: case 13: case 14: case 15:
          pressure_order = -(velocity_order-1)*10;
          break;
        case 22: case 23: case 24:
          pressure_order = -(velocity_order-11)*10;
          break;
      }
      break;
    // discontinuous spaces
    case 1:case 2: case 3: case 4: case 5:
      // TODO CHECK IF THIS IS CORRECT!!! IN NSE3D,pressure_order=1  ?!!
      pressure_order = -(velocity_order-1)*10;
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
  }
  TDatabase::ParamDB->PRESSURE_SPACE  = pressure_order;
  velocity_pressure_orders.second = pressure_order;

  Output::print("velocity space", setw(10), velocity_pressure_orders.first);
  Output::print("pressure space", setw(10), velocity_pressure_orders.second);
}

/**************************************************************************** */
void Time_NSE3D::assemble_initial_time()
{
  if(systems_.size() > 1)
  {
    this->restrict_function();
  }
  
  for(auto &s : this->systems_)
  {
    // assemble the initial matrices and right hand side
    call_assembling_routine(s, LocalAssembling_type::TNSE3D_LinGAL);
    // manage dirichlet condition by copying non-actives DoFsfrom rhs to solution
    s.solution_.copy_nonactive(s.rhs_);
  }
  
  /** After copy_nonactive, the solution vectors needs to be Comm-updated
   * in MPI-case in order to be consistently saved. It is necessary that
   * the vector is consistently saved because it is the only way to
   * ensure that its multiplication with an inconsistently saved matrix
   * (multiplication which appears in the defect and rhs computations)
   * give the correct results.
   * When we call copy_nonactive in MPI-case, we have to remember the following:
   * it can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active.
   * That's why we have to update the values so that the vector becomes consistent again.
   */
  #ifdef _MPI
    double *u1  = this->systems_.front().solution_.block(0);
    double *u2  = this->systems_.front().solution_.block(1);
    double *u3  = this->systems_.front().solution_.block(2);
    double *p   = this->systems_.front().solution_.block(3);
    this->systems_.front().velocitySpace_.get_communicator().consistency_update(u1, 3);
    this->systems_.front().velocitySpace_.get_communicator().consistency_update(u2, 3);
    this->systems_.front().velocitySpace_.get_communicator().consistency_update(u3, 3);
    this->systems_.front().pressureSpace_.get_communicator().consistency_update(p, 3);
  #endif
  // copy the last right hand side and solution vectors to the old ones
  this->old_rhs_      = this->systems_.front().rhs_;
  this->solution_m1 = this->systems_.front().solution_;
}
/**************************************************************************** */
void Time_NSE3D::restrict_function()
{
  for( int block = 0; block < 3 ;++block)
  {
    std::vector<const TFESpace3D*> spaces;
    std::vector<double*> u_entries;
    std::vector<size_t> u_ns_dofs;
    for(auto &s : systems_ )
    {
      spaces.push_back(&s.velocitySpace_);
      u_entries.push_back(s.solution_.block(block));
      u_ns_dofs.push_back(s.solution_.length(block));
    }
    GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
  }
}

/**************************************************************************** */
void Time_NSE3D::assemble_rhs()
{
  System_per_grid& s = this->systems_.front();
  // reset the right hand side of the grid of interest (finest)
  s.rhs_.reset();
  rhs_from_time_disc.copy_structure(s.rhs_);
  // preparation of the right hand side for the system solve
  rhs_from_time_disc.reset();
  // calling assembling routine
  call_assembling_routine(s, LocalAssembling_type::TNSE3D_Rhs);
  BlockVector temp = s.rhs_;
  // copy right hand side 
  rhs_from_time_disc = s.rhs_;
  unsigned int n_sols = time_stepping_scheme.n_old_solutions();
  std::vector<BlockVector> old_sols (n_sols);
  old_sols[0] = s.solution_;
  // this is needed for the BDF2 method
  if(old_sols.size() == 2)
    old_sols[1] = solution_m2;
  // the right hand side vectors: old rhs is needed for the 
  // Crank-Nicolson and Fractional Step schemes
  std::vector<BlockVector> all_rhs(2);
  all_rhs[0] = s.rhs_;
  all_rhs[1] = old_rhs_;
  
  // prepare the right hand side vector
  time_stepping_scheme.prepare_rhs_from_time_disc(s.matrix_, s.massMatrix_, 
						  all_rhs, old_sols);
  // on the way back, the system rhs is stored in all_rhs[0], copy to the 
  // rhs_from_time_disc vector
  rhs_from_time_disc = all_rhs[0];
  // update the rhs for the next time step
  old_rhs_ = s.rhs_;
  
  // scaling of all B-block is Important:
  // If all B-blocks are linear, we have to scale them only once. So it's done only
  // at the first time. In the case of nonlinearity, one have to do this in 
  // each nonlinear iteration. For example, the SUPG case.
  if(time_stepping_scheme.current_step_==1 && (time_stepping_scheme.b_bt_linear_nl.compare("linear")==0))
  {
    for(System_per_grid &s : this->systems_)
      time_stepping_scheme.scale_descale_all_b_blocks(s.matrix_, "scale");
    if(time_stepping_scheme.current_step_ >1)
    {
      ErrThrow("wrong time stepping");
    }
  }
  
  // retrieve the non active from "temporary" into rhs vector
  rhs_from_time_disc.copy_nonactive(temp);
  old_rhs_.copy_nonactive(temp);
  // copy the non active to the solution vector
  s.solution_.copy_nonactive(temp);

  /** After copy_nonactive, the solution vectors needs to be Comm-updated
   * in MPI-case in order to be consistently saved. It is necessary that
   * the vector is consistently saved because it is the only way to
   * ensure that its multiplication with an inconsistently saved matrix
   * (multiplication which appears in the defect and rhs computations)
   * give the correct results.
   * When we call copy_nonactive in MPI-case, we have to remember the following:
   * it can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active.
   * That's why we have to update the values so that the vector becomes consistent again.
   */
  #ifdef _MPI
    double *u1 = this->systems_.front().solution_.block(0);
    double *u2 = this->systems_.front().solution_.block(1);
    double *u3 = this->systems_.front().solution_.block(2);
    double *p  = this->systems_.front().solution_.block(3);
    this->systems_.front().velocitySpace_.get_communicator().consistency_update(u1, 3);
    this->systems_.front().velocitySpace_.get_communicator().consistency_update(u2, 3);
    this->systems_.front().velocitySpace_.get_communicator().consistency_update(u3, 3);
    this->systems_.front().pressureSpace_.get_communicator().consistency_update(p, 3);
  #endif
}
/**************************************************************************** */
void Time_NSE3D::assemble_nonlinear_term()
{
  if(systems_.size()>1)
    this->restrict_function();
  for(System_per_grid &s : this->systems_)
  {
    call_assembling_routine(s, LocalAssembling_type::TNSE3D_NLGAL);
  }
}

/**************************************************************************** */
void Time_NSE3D::assemble_system()
{
  for(System_per_grid &s : this->systems_)
    time_stepping_scheme.prepare_system_matrix(s.matrix_, s.massMatrix_);

  Output::info<5>("Assemble System", "Assembled the system matrix which"
      " will be passed to the solver");
}

/**************************************************************************** */
bool Time_NSE3D::stop_it(unsigned int iteration_counter)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif

  // compute, update and display defect and residuals
  //  compute_residuals();

  // stores current norm of the residual. They are normed per default in
  // the class Residuals
  const double normOfResidual        = this->get_full_residual();
  //  const double normOfImpulseResidual = this->get_impulse_residual();
  //  const double normOfMassResidual    = this->get_mass_residual();
  //  const double oldNormOfResidual     = this->old_residual_[1].fullR;
  // hold the residual from up to 10 iterations ago
  const double veryOldNormOfResidual  = this->old_residual_.front().fullResidual;

  //  Output::print("nonlinear step  :  " , setw(3), iteration_counter);
  //  Output::print("impulse_residual:  " , setw(3), normOfImpulseResidual);
  //  Output::print("mass_residual   :  " , setw(3), normOfMassResidual);
  //  Output::print("residual        :  " , setw(3), normOfResidual);
  System_per_grid& s = this->systems_.front();
  size_t nu=s.solution_.length(0);
  size_t np=s.solution_.length(3);
  
  Output::print<5>("B " , Ddot(3*nu+np,s.solution_.get_entries(),s.solution_.get_entries()), " ",
                Ddot(3*nu,rhs_from_time_disc.get_entries(),rhs_from_time_disc.get_entries()) , " "  , 
                Ddot(np,rhs_from_time_disc.get_entries()+3*nu,rhs_from_time_disc.get_entries()+3*nu)," ",
                Ddot(3*nu+np,rhs_from_time_disc.get_entries(),rhs_from_time_disc.get_entries()));
  
  // this is the convergence ratio between actual step and last step
  // TODO : correct oldNormOfResidual to be the residual of last step
  // and print it out
  //  if ( iteration_counter > 0 )
  //  {
  //  Output::print("convergence rate:  " , setw(3), normOfResidual/oldNormOfResidual);
  //  }

  /* Some parameters have to be set up at the first nonlinear iteration */
  if(iteration_counter == 0)
  {
    initial_residual_ = normOfResidual;

    // saves the solution from previous time step with nonActive of current step
    this->solution_m1 = this->systems_.front().solution_;
  }
  // check if minimum number of iterations was performed already
  size_t min_it = db_["nonlinloop_minit"];
  if(iteration_counter < min_it)
	  return false;

  /** Parameters for stopping criteria (desired precision epsilon, max number
   *  of iteration, convergence rate, IMEX_scheme : if IMEX is used, we only
   *  need to solve the nonlinear iterations for the first time step in order
   *  to obtain both initial solution and the one after. Then, the extrapolated
   *  scheme is solved just once, as a linear system. ) */
  double epsilon    = db_["nonlinloop_epsilon"];
  size_t max_It     = db_["nonlinloop_maxit"];
  double conv_speed = db_["nonlinloop_slowfactor"];
  bool slow_conv    = false;
  
  if ( db_["nonlinloop_scale_epsilon_with_size"] )
  {
    epsilon *= sqrt(this->get_size());
    if (my_rank==0)
      Output::print("stopping tolerance for nonlinear iteration ", epsilon);
  }

//   if ( normOfResidual >= conv_speed*veryOldNormOfResidual )
//   {
//     slow_conv = true;
//   }
  // Stopping criteria
  if ( (normOfResidual <= epsilon) || (iteration_counter == max_It)
      || (slow_conv) )
  {
    if(slow_conv && my_rank==0)
      Output::print<1>(" SLOW !!! ", normOfResidual/veryOldNormOfResidual);

    if(my_rank==0)
    {
      Output::print("Last nonlinear iteration : ", setw(3), iteration_counter,
                    "\t\t", "Residual: ", normOfResidual,
                    "\t\t", "Reduction: ", normOfResidual/initial_residual_);
    }
    // copy solution for the next time step: BDF2 needs two previous solutions 
    solution_m2 = solution_m1;
    solution_m1 = s.solution_;
    // descale the matrices, since only the diagonal A block will
    // be reassembled in the next time step
    if (this->imex_scheme(0) && iteration_counter>0)
      return true; // in these conditions, the matrix are already descaled
    else
    {
      for(System_per_grid &s : this->systems_)
      {
	time_stepping_scheme.reset_linear_matrices(s.matrix_, s.massMatrix_);
        // descale if it's rescaled at the next time step for bdf schemes
	time_stepping_scheme.scale_descale_all_b_blocks(s.matrix_, "descale");
      }
      return true;
    }
  }
  else
    return false;
}

/**************************************************************************** */
void Time_NSE3D::compute_residuals()
{
  System_per_grid& s = this->systems_.front();
  unsigned int number_u_Dof = s.solution_.length(0);
  unsigned int number_p_Dof = s.solution_.length(3);

#ifdef _MPI
    //MPI: put solution in consistency level 3
    auto comms = s.matrix_.get_communicators();
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(s.solution_.block(bl), 3);
    }
#endif

  // copy rhs to defect and compute defect
  this->defect_ = rhs_from_time_disc;
  s.matrix_.apply_scaled_add(s.solution_, defect_,-1.);

  if(s.matrix_.pressure_projection_enabled())
  {
    TFEFunction3D defect_fctn(&s.pressureSpace_,
                              "p_def","pressure defect function",
                              &defect_[3*number_u_Dof], number_p_Dof);
    defect_fctn.project_into_L20();
  }

  // This is the calculation of the residual, given the defect.
  BlockVector defect_impuls({number_u_Dof,number_u_Dof,number_u_Dof});
  BlockVector defect_mass({number_p_Dof});
  //copy the entries (BlockVector offers no functionality to do this more nicely)
  for(size_t i = 0; i<3*number_u_Dof ;++i)
    defect_impuls.get_entries()[i] = defect_.get_entries()[i];
  for(size_t i =0 ; i<number_p_Dof ; ++i)
    defect_mass.get_entries()[i] = defect_.get_entries()[3*number_u_Dof + i];

#ifdef _MPI
  double impuls_residual_square = defect_impuls.norm_global({comms[0],comms[1],comms[2]});
  impuls_residual_square *= impuls_residual_square;
  double mass_residual_square = defect_mass.norm_global({comms[3]});
  mass_residual_square *= mass_residual_square;
#else
  double impuls_residual_square = defect_impuls.norm();
  impuls_residual_square *= impuls_residual_square;
  double mass_residual_square = defect_mass.norm();
  mass_residual_square *= mass_residual_square;
#endif

  Residuals current_residual(impuls_residual_square, mass_residual_square);
  old_residual_.add(current_residual);
}

/**************************************************************************** */
void Time_NSE3D::solve()
{
  System_per_grid& s = systems_.front();
  // store previous solution for damping, it is a pointer so that we can avoid
  // the copy in case of no damping
  double damping = this->db_["nonlinloop_damping_factor"];
  std::shared_ptr<BlockVector> old_solution(nullptr);
  if(damping != 1.0)
    old_solution = std::make_shared<BlockVector>(s.solution_);
#ifndef _MPI
  solver_.solve(s.matrix_, rhs_from_time_disc, s.solution_);
#elif defined(_MPI)
  if(solver_.get_db()["solver_type"].is("direct"))
  {
    if(damping != 1.0)
      Output::warn("Time_NSE3D::solve", "damping in an MPI context is not tested");

    //set up a MUMPS wrapper
    MumpsWrapper mumps_wrapper(s.matrix_);
    //kick off the solving process
    mumps_wrapper.solve(rhs_from_time_disc, s.solution_);
  }
  else
    solver_.solve(s.matrix_, rhs_from_time_disc, s.solution_); // same as sequential
#endif

  // apply damping if prescribed
  if(damping != 1.0)
  {
    s.solution_.scale(damping);
    s.solution_.add_scaled(*old_solution, 1-damping);
  }

  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11, A22 and A33 matrices are
  // reset and assembled again but the non-diagonal blocks are scaled, so
  // for the next iteration we have to descale, see assemble_system()
  for(System_per_grid &s : this->systems_)
    time_stepping_scheme.reset_linear_matrices(s.matrix_, s.massMatrix_);

  if(s.matrix_.pressure_projection_enabled())
     s.p_.project_into_L20();
}

/**************************************************************************** */
void Time_NSE3D::output(int m, int &image)
{
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  	bool no_output = !db_["output_write_vtk"] &&
 	                 !db_["output_compute_errors"] &&
 	                 !db_["write_solution_binary"];
	if(no_output)
		return;

  System_per_grid& s = this->systems_.front();
  TFEFunction3D* u1 = s.u_.GetComponent(0);
  TFEFunction3D* u2 = s.u_.GetComponent(1);
  TFEFunction3D* u3 = s.u_.GetComponent(2);

  Output::print("** COMPUTING FLUX **");
  double flux;
  s.u_.compute_flux(3, flux);
  Output::print("** FLUX ON 3 = ", flux,
		" ** Expected: ", 2*3.1415*TDatabase::TimeDB->CURRENTTIME);
  
  if((size_t)db_["verbosity"]> 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    u3->PrintMinMax();
    s.p_.PrintMinMax();
  }

  // write solution to a vtk file
  outputWriter.add_fe_function(&s.p_);
  outputWriter.add_fe_vector_function(&s.u_);
  outputWriter.write(image);

  // Measure errors to known solution
  // if an exact solution is not known, it is usually set to be zero, so that
  // in such a case, here only integrals of the solution are computed.
  if(db_["output_compute_errors"])
  {
    double err_u1[4];  // FIXME? Of these arrays only the 2 first entries are
    double err_u2[4];  // used. But the evil GetErrors() will corrupt memory if
    double err_u3[4];  // these have not at least size 4.
    double err_p[4];

    TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    MultiIndex3D allderiv[4]= {D000, D100, D010, D001};
    const TFESpace3D *v_space = &this->get_velocity_space();
    const TFESpace3D *p_space = &this->get_pressure_space();

   double tau = TDatabase::TimeDB->TIMESTEPLENGTH;

    // Errors in velocity components and pressure
    u1 ->GetErrors(example_.get_exact(0), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &v_space, err_u1);
    u2 ->GetErrors(example_.get_exact(1), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &v_space, err_u2);
    u3 ->GetErrors(example_.get_exact(2), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &v_space, err_u3);
    s.p_.GetErrors(example_.get_exact(3), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &p_space, err_p);

#ifdef _MPI
    double err_red[8]; //memory for global (across all processes) error
    double err_send[8]; //fill send buffer
    err_send[0]=err_u1[0];
    err_send[1]=err_u1[1];
    err_send[2]=err_u2[0];
    err_send[3]=err_u2[1];
    err_send[4]=err_u3[0];
    err_send[5]=err_u3[1];
    err_send[6]=err_p[0];
    err_send[7]=err_p[1];

    MPI_Allreduce(err_send, err_red, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    int j;
    for(j=0;j<8;j++)
    {//MPI: sqrt was skipped in GetErrors function - do it here globally!
      err_red[j] = sqrt(err_red[j]);
    }
    //fill the reduced errors back where they belong
    err_u1[0] = err_red[0];
    err_u1[1] = err_red[1];
    err_u2[0] = err_red[2];
    err_u2[1] = err_red[3];
    err_u3[0] = err_red[4];
    err_u3[1] = err_red[5];
    err_p[0] = err_red[6];
    err_p[1] = err_red[7];
#else
    int my_rank = 0;
#endif

    errors_.at(0) = err_u1[0]*err_u1[0] + err_u2[0]*err_u2[0] +
                    err_u3[0]*err_u3[0];  // (L2-norm)^2 for u
    errors_.at(1) = err_u1[1]*err_u1[1] + err_u2[1]*err_u2[1] +
                    err_u3[1]*err_u3[1];  // (H1-semi)^2 for u
    errors_.at(2) = err_p[0]*err_p[0];  // (L2-norm)^2 for p
    errors_.at(3) = err_p[1]*err_p[1];  // (H1-norm)^2 for p

    errors_.at(4) += (errors_[0] + this->errors_.at(5))*tau*0.5; //l2(0,t,l2)(u)
    errors_.at(5) = errors_.at(0);
    errors_.at(6) += (errors_.at(1) + errors_.at(7))*tau*0.5; // l2(0,t,h1)(u)
    errors_.at(7) = errors_.at(1);
    
    errors_.at(8) += (errors_.at(2)+errors_.at(9))*tau*0.5; //l2(0,t,l2) (p)
    errors_.at(9) = errors_.at(2);
    errors_.at(10) += (errors_.at(3) + errors_.at(11))*tau*0.5;//l2(0,t,h1)(p)
    errors_.at(11) = errors_.at(3);
    
    // print errors
    if (my_rank == 0 )
    {
      Output::print<1>("L2(u)         : ", setprecision(10), sqrt(errors_.at(0)));
      Output::print<1>("H1-semi(u)    : ", setprecision(10), sqrt(errors_.at(1)));
      
      Output::print<1>("L2(0,t,L2(u)) : ", setprecision(10), sqrt(errors_.at(4)));
      Output::print<1>("L2(0,t,H1-semi(u)) : ", 
                       setprecision(10), sqrt(errors_.at(6)));
      Output::print<1>("L2(p)      : ", setprecision(10), sqrt(errors_.at(2)));
      Output::print<1>("H1-semi(p)): ", setprecision(10), sqrt(errors_.at(3)));
      
      Output::print<1>("L2(0,t,L2(p)) : ", setprecision(10), sqrt(errors_.at(8)) );
      Output::print<1>("L2(0,t,H1-semi(p)) : ", setprecision(10), sqrt(errors_.at(10)) );
    }
  }
   delete u1;
   delete u2;
   delete u3;

   // do post-processing step depending on what the example implements, if needed
   example_.do_post_processing(*this);
   
   if(db_["write_solution_binary"].is(true))
   { size_t interval = db_["write_solution_binary_all_n_steps"];
    if(m % interval == 0)
    {//write solution to a binary file
      std::string file = db_["write_solution_binary_file"];
      if(!db_["overwrite_solution_binary"]) //create a new file every time
      {
    	  file += ".";
    	  file += std::to_string(TDatabase::TimeDB->CURRENTTIME);
#ifdef _MPI
    	  file += ".proc" + std::to_string(my_rank);
#endif
      }
      Output::info("output", "Writing current solution to file ", file);
      systems_.front().solution_.write_to_file(file);
    }
  }
}

/**************************************************************************** */
const Residuals& Time_NSE3D::get_residuals() const
{
  return old_residual_.back();
}

/**************************************************************************** */
double Time_NSE3D::get_impulse_residual() const
{
  return old_residual_.back().impulsResidual;
}

/**************************************************************************** */
double Time_NSE3D::get_mass_residual() const
{
  return old_residual_.back().massResidual;
}

/**************************************************************************** */
double Time_NSE3D::get_full_residual() const
{
  return old_residual_.back().fullResidual;
}

/**************************************************************************** */
std::array< double, int(6) > Time_NSE3D::get_errors() const
{
  std::array<double, int(6)> error_at_time_points;
  error_at_time_points[0] = sqrt(this->errors_[0]); // L2 velocity error
  error_at_time_points[1] = sqrt(this->errors_[1]); // H1 velocity error
  error_at_time_points[2] = sqrt(this->errors_[2]); // L2 pressure error
  error_at_time_points[3] = sqrt(this->errors_[3]); // H1 pressure error

  return error_at_time_points;
}

/**************************************************************************** */
void Time_NSE3D::output_problem_size_info() const
{
  // print out some information about number of DoFs and mesh size
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 3 * n_u + n_p; // total number of degrees of freedom
  int nActive = this->get_velocity_space().GetN_ActiveDegrees();
  double h_min, h_max;
  TCollection * coll = this->get_velocity_space().GetCollection();
  coll->GetHminHmax(&h_min, &h_max);

  Output::print("N_Cells     : ", setw(10), coll->GetN_Cells());
  Output::print("h (min,max) : ", setw(10), h_min ," ", setw(12), h_max);
  Output::print("dof Velocity: ", setw(10), 3* n_u);
  Output::print("dof Pressure: ", setw(10), n_p   );
  Output::print("dof all     : ", setw(10), n_dof );
  Output::print("active dof  : ", setw(10), 3*nActive);
}

/**************************************************************************** */
void Time_NSE3D::construct_extrapolated_solution()
{
  this->extrapolated_solution_.reset();
  this->extrapolated_solution_ = this->solution_m1;
  this->extrapolated_solution_.scale(-1.);
  this->extrapolated_solution_.add_scaled(this->systems_.front().solution_,2.);
  this->extrapolated_solution_.copy_nonactive(this->systems_.front().rhs_);
  // Now extrapolated_solution_ = 2*u(t-1)-u(t-2), only on the finest mesh
}

/**************************************************************************** */
TFEFunction3D* Time_NSE3D::get_velocity_component(int i)
{
  if(i==0)
    return this->systems_.front().u_.GetComponent(0);
  else if(i==1)
    return this->systems_.front().u_.GetComponent(1);
  else  if(i==2)
    return this->systems_.front().u_.GetComponent(2);
  else
    throw std::runtime_error("There are only three velocity components!");
}

/**************************************************************************** */
bool Time_NSE3D::imex_scheme(bool print_info)
{
  // change maximum number of nonlin_iterations to 1 in IMEX case
  if((db_["imex_scheme_"] && time_stepping_scheme.current_step_ >= 3))
  {
    db_["nonlinloop_maxit"] = 1;
    if(print_info) // condition is here just to print it once
      Output::info<1>("Nonlinear Loop MaxIteration",
                    "The parameter 'nonlinloop_maxit' was changed to 1."
                    " Only one non-linear iteration is done, because the IMEX scheme was chosen.\n");
    return true;
  }
  else
    return false;
}

/**************************************************************************** */
void Time_NSE3D::call_assembling_routine(Time_NSE3D::System_per_grid& s, 
                          LocalAssembling_type type)
{
  std::vector<const TFESpace3D*> space_mat;
  std::vector<const TFESpace3D*> space_rhs;
  std::vector<TFEFunction3D*> fefunctions;
  // set the memory 
  set_arrays(s, space_mat, space_rhs, fefunctions);
  
  // prepare matrices and rhs
  std::vector<TSquareMatrix3D*> sqMat;
  std::vector<TMatrix3D*> reMat;
  std::vector<double*> rhs_array;
  
  set_matrices_rhs(s, type, sqMat, reMat, rhs_array);
  // find out if we have to do upwinding
  bool do_upwinding = false;
  if(type != LocalAssembling_type::TNSE3D_Rhs)
  {
      bool mdml =  solver_.is_using_multigrid()
                  && solver_.get_multigrid()->is_using_mdml();
      bool on_finest_grid = &systems_.front() == &s;
      do_upwinding = (db_["space_discretization_type"].is("upwind")
                     || (mdml && !on_finest_grid));
    
    if(do_upwinding)  //HOTFIX: Check the documentation!
      assemble_nse = Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION;
    else
      assemble_nse = Hotfixglobal_AssembleNSE::WITH_CONVECTION;
  }
  // Boundary conditions and value
  BoundCondFunct3D * boundary_conditions[4] = {
    s.velocitySpace_.get_boundary_condition(),
    s.velocitySpace_.get_boundary_condition(),
    s.velocitySpace_.get_boundary_condition(),
    s.pressureSpace_.get_boundary_condition() };

  std::array<BoundValueFunct3D*, 4> boundary_values;
  boundary_values[0] = example_.get_bd(0);
  boundary_values[1] = example_.get_bd(1);
  boundary_values[2] = example_.get_bd(2);
  boundary_values[3] = example_.get_bd(3);
  
  LocalAssembling3D localAssembling(this->db_, type, fefunctions.data(),
                                    this->example_.get_coeffs(),
                                    this->get_space_disc_global());

  if (type == LocalAssembling_type::TNSE3D_Rhs) {
    Output::print(" ** START ASSEMBLE PRESSURE BC ON RHS **");

    const TFESpace3D *v_space = &s.velocitySpace_;
    const TFESpace3D *p_space = &s.pressureSpace_;

    std::vector<TBoundFace*> boundaryFaceList;

    BoundaryAssembling3D bi;
    for (int k = 0; k < TDatabase::ParamDB->n_neumann_boundary; k++)
    {
      double t=TDatabase::TimeDB->CURRENTTIME;
      double PI = acos(-1.0);
      double pressure_of_t = TDatabase::ParamDB->neumann_boundary_value[k]*sin(2*PI*t);

      Output::print(" ** set value ", pressure_of_t,
		    " on boundary ",TDatabase::ParamDB->neumann_boundary_id[k]);
      
      bi.rhs_g_v_n(s.rhs_, v_space, nullptr, boundaryFaceList,
		   TDatabase::ParamDB->neumann_boundary_id[k], pressure_of_t);
    }
  }
  
  // assemble all the matrices and right hand side 
  Assemble3D(space_mat.size(), space_mat.data(),
	     sqMat.size(), sqMat.data(), reMat.size(), reMat.data(),
             rhs_array.size(), rhs_array.data(), space_rhs.data(),
             boundary_conditions, boundary_values.data(), localAssembling);
  
  
  if(do_upwinding && type != LocalAssembling_type::TNSE3D_Rhs)
  {
    double one_over_nu = 1/example_.get_nu(); //the inverse of the example's diffusion coefficient
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        UpwindForNavierStokes3D(sqMat[0], fefunctions[0], fefunctions[1],
                                fefunctions[2], one_over_nu);
        Output::print<3>("UPWINDING DONE with 1 square matrices.");
      break;
      case 3:
      case 4:
      case 14:
        UpwindForNavierStokes3D(sqMat[0], fefunctions[0], fefunctions[1],
                                fefunctions[2], one_over_nu);
        UpwindForNavierStokes3D(sqMat[1], fefunctions[0], fefunctions[1],
                                fefunctions[2], one_over_nu);
        UpwindForNavierStokes3D(sqMat[2], fefunctions[0], fefunctions[1],
                                fefunctions[2], one_over_nu);
        Output::print<3>("UPWINDING DONE with 3 square matrices.");
        break;
    }
  }
  
}
/**************************************************************************** */
void Time_NSE3D::set_arrays(Time_NSE3D::System_per_grid& s, 
        std::vector<const TFESpace3D*> &spaces, std::vector< const TFESpace3D* >& spaces_rhs,
        std::vector< TFEFunction3D*> &functions)
{
  spaces.resize(2);
  spaces[0] = &s.velocitySpace_;
  spaces[1] = &s.pressureSpace_;
  spaces_rhs.resize(3);
  spaces_rhs[0] = &s.velocitySpace_;
  spaces_rhs[1] = &s.velocitySpace_;
  spaces_rhs[2] = &s.velocitySpace_;
  
  if(TDatabase::ParamDB->NSTYPE==14)
  {
    spaces_rhs.resize(4);
    spaces_rhs[3] = &s.pressureSpace_;
  }
  functions.resize(4);  
  functions[0] = s.u_.GetComponent(0);
  functions[1] = s.u_.GetComponent(1);
  functions[2] = s.u_.GetComponent(2);
  functions[3] = &s.p_;
}
/**************************************************************************** */
void Time_NSE3D::set_matrices_rhs(Time_NSE3D::System_per_grid& s, LocalAssembling_type type,
        std::vector<TSquareMatrix3D*> &sqMat, std::vector<TMatrix3D*> &reMat,
        std::vector<double*> &rhs_array)
{
  rhs_array.resize(0);
  sqMat.resize(0);
  reMat.resize(0);
  
  std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix_.get_blocks_uniquely();
  std::vector<std::shared_ptr<FEMatrix>> mass_blocks
         = s.massMatrix_.get_blocks_uniquely(true);
  
  switch(type)
  {
    case LocalAssembling_type::TNSE3D_LinGAL:
    {
      rhs_array.resize(3);
      rhs_array[0] = s.rhs_.block(0);
      rhs_array[1] = s.rhs_.block(1);
      rhs_array[2] = s.rhs_.block(2);
      s.rhs_.reset();
      switch(TDatabase::ParamDB->NSTYPE)
      {
	case 1:
	  if(blocks.size() != 4)
	  {
	    ErrThrow("Wrong blocks.size() ", blocks.size(), " instead of 4.");
	  }
	  sqMat.resize(2);
          sqMat[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
          // rectangular matrices
          reMat.resize(3);
          reMat[0] = reinterpret_cast<TMatrix3D*>(blocks.at(1).get());
          reMat[1] = reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
	  reMat[2] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());
	  break;
	case 2:
	  if(blocks.size() != 7)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size(), " instead of 7.");
          }
          sqMat.resize(2);
          sqMat[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          // mass matrix
          sqMat[1] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
          
          
          // rectangular matrices
          reMat.resize(6);
          reMat[0] = reinterpret_cast<TMatrix3D*>(blocks.at(4).get()); //first the lying B blocks
          reMat[1] = reinterpret_cast<TMatrix3D*>(blocks.at(5).get());
          reMat[2] = reinterpret_cast<TMatrix3D*>(blocks.at(6).get());
          reMat[3] = reinterpret_cast<TMatrix3D*>(blocks.at(1).get()); //then the standing B blocks
          reMat[4] = reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
          reMat[5] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());
	  break;
	case 3:
	  if(blocks.size() != 12)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size(), " instead of 12.");
          }
          sqMat.resize(12);
          sqMat[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
          sqMat[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
          sqMat[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
          sqMat[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
          sqMat[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
          sqMat[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
          sqMat[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
          sqMat[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
          // mass matrices
          sqMat[9] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
	  sqMat[10] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(5).get());
	  sqMat[11] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(10).get());
          
          // rectangular matrices
          reMat.resize(3);
          reMat[0] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());  // standing B blocks
          reMat[1] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
          reMat[2] = reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
	  break;
	case 4:
	  if(blocks.size() != 15)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size(), " instead of 15.");
          }
          sqMat.resize(12);
          sqMat[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
          sqMat[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
          sqMat[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
          sqMat[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
          sqMat[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
          sqMat[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
          sqMat[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
          sqMat[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
          // mass matrices
          sqMat[9] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
	  sqMat[10] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(5).get());
	  sqMat[11] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(10).get());
          
          // rectangular matrices
          reMat.resize(6);
          reMat[0] = reinterpret_cast<TMatrix3D*>(blocks.at(12).get()); //first the lying B blocks
          reMat[1] = reinterpret_cast<TMatrix3D*>(blocks.at(13).get());
          reMat[2] = reinterpret_cast<TMatrix3D*>(blocks.at(14).get());
          reMat[3] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());  //than the standing B blocks
          reMat[4] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
          reMat[5] = reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
	  break;
      }// endswitch NSTYPE
    }
    break;// case LocalAssembling3D_type::TNSE3D_LinGAL
    //===============================================
    case LocalAssembling_type::TNSE3D_NLGAL:
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
	case 1:
	case 2:
	  blocks = s.matrix_.get_blocks_uniquely({{0,0},{1,1},{2,2}});
	  sqMat.resize(1);
          sqMat[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
	  break;
	case 3:
	case 4:
	  if(db_["space_discretization_type"].is("smagorinsky"))
          {
            sqMat.resize(9);
            blocks = s.matrix_.get_blocks_uniquely({{0,0}, {0,1}, {0,2}, 
                                                    {1,0}, {1,1}, {1,2},
                                                    {2,0}, {2,1}, {2,2}});
            sqMat[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
            sqMat[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(3).get());
            sqMat[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
            sqMat[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
            sqMat[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
            sqMat[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(7).get());
            sqMat[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
          }
          else
          {
            sqMat.resize(3); 
            blocks = s.matrix_.get_blocks_uniquely({{0,0},{1,1},{2,2}});
            sqMat[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
          }
	  break;
      }// endswitch NSTYPE
    }
    break;
    //===============================================
    case LocalAssembling_type::TNSE3D_Rhs:
    {
      rhs_array.resize(3);
      rhs_array[0] = s.rhs_.block(0);
      rhs_array[1] = s.rhs_.block(1);
      rhs_array[2] = s.rhs_.block(2);
      if(TDatabase::ParamDB->NSTYPE==14)
      {
	rhs_array.resize(4);
	rhs_array[3] = s.rhs_.block(3);
      }
      s.rhs_.reset();
    }
    break;
    //===============================================
    default:
      ErrThrow("Local Assembling type is not supported");
  }
  
  // reset matrices
  for(auto sm : sqMat)
    sm->reset();
  for(auto rm : reMat)
    rm->reset();
}
/**************************************************************************** */
/** ************************************************************************ */

#include <Time_NSE3D.h>
#include <Database.h>
#include <Assemble3D.h>
#include <LocalAssembling3D.h>
#include <LinAlg.h>
#include <Output3D.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <Multigrid.h>
#include <Variational_MultiScale3D.h>

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
  ParameterDatabase time_db = ParameterDatabase::default_time_database();
  db.merge(time_db,true);
  
  ParameterDatabase turb_mod_db = ParameterDatabase::default_turbulence_model_database();
  db.merge(turb_mod_db,true);

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
 : velocitySpace_(&coll, (char*)"u", (char*)"velocity space",  example.get_bc(0),
                  order.first),
   pressureSpace_(&coll, (char*)"p", (char*)"pressure space", example.get_bc(3),
                  order.second),
   matrix_({&velocitySpace_, &velocitySpace_, &velocitySpace_, &pressureSpace_}),
   massMatrix_({&velocitySpace_, &velocitySpace_, &velocitySpace_}),
   rhs_(matrix_, true),
   solution_(matrix_, false),
   u_(&velocitySpace_, (char*)"u", (char*)"u", solution_.block(0),
     solution_.length(0), 3),
   p_(&pressureSpace_, (char*)"p", (char*)"p", solution_.block(3),
     solution_.length(3)), 
   MatrixK({&velocitySpace_, &velocitySpace_, &velocitySpace_, &pressureSpace_}),
   Old_Sol(matrix_, false),
   u_old(&velocitySpace_, (char*) "uold", (char*)"uold", Old_Sol.block(0), 
         Old_Sol.length(0), 3)
{
  // Mass Matrix
  // Output::increaseVerbosity(5);

  massMatrix_ = BlockFEMatrix::Mass_NSE3D(velocitySpace_);
      
  switch(type)
  {
    case Time_NSE3D::Matrix::Type1:
      matrix_ = BlockFEMatrix::NSE3D_Type1(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type2:
      matrix_ = BlockFEMatrix::NSE3D_Type2(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type3:
      matrix_ = BlockFEMatrix::NSE3D_Type3(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type4:
      matrix_ = BlockFEMatrix::NSE3D_Type4(velocitySpace_, pressureSpace_);
      if(TDatabase::ParamDB->DISCTYPE == RESIDUAL_VMS)
      {
        MatrixK = BlockFEMatrix::NSE3D_Type4(velocitySpace_, pressureSpace_);
      }
      break;
    case Time_NSE3D::Matrix::Type14:
      matrix_ = BlockFEMatrix::NSE3D_Type14(velocitySpace_, pressureSpace_);
      if(TDatabase::ParamDB->DISCTYPE == RESIDUAL_VMS)
      {
        MatrixK = BlockFEMatrix::NSE3D_Type4(velocitySpace_, pressureSpace_);
      }
      break;
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }
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
: db_(get_default_TNSE3D_parameters()), systems_(), example_(ex),
   solver_(param_db), defect_(), old_residual_(), 
   initial_residual_(1e10), errors_(), oldtau_()
{
  db_.merge(param_db, false);
  this->check_parameters();
  
  // get the velocity and pressure orders
  std::pair <int,int>
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                               TDatabase::ParamDB->PRESSURE_SPACE);
  this->get_velocity_pressure_orders(velocity_pressure_orders);
  
  // get matrix type
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
  
  if((TDatabase::ParamDB->DISCTYPE == SUPG || 
    TDatabase::ParamDB->DISCTYPE == RESIDUAL_VMS )  && 
    (type != Matrix::Type14 && type !=Matrix::Type4))
  {
    ErrThrow("The SUPG method is only implemented for NSTYPE 4 and 14");
  }
  
  // multigrid method 
  bool usingMultigrid = solver_.is_using_multigrid();
  TCollection *coll = collections_.front(); // finest grid collection
  // create finite element space and function, a matrix, rhs, and solution
#ifdef _MPI
  systems_.emplace_back(example_, *coll, velocity_pressure_orders, type,
                        maxSubDomainPerDof);
#else
  systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);
#endif  
  // projection-based VMS, only on the finest level 
  if(TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
  {
    // VMS projection order 
    int projection_order = db_["vms_projection_space_order"];
    if(projection_order < 0)
      projection_order = 0;
    
    // projection space
    projection_space_ = 
       std::make_shared<TFESpace3D>(collections_.front(), (char*)"L",  
                           (char*)"vms projection space", example_.get_bc(2), 
                           DiscP_PSpace, projection_order);
    int ndofP = projection_space_->GetN_DegreesOfFreedom();
   
    // create vector for vms projection, used if small resolved scales are needed
    this->vms_small_resolved_scales.resize(6*ndofP);
    // finite element vector function for vms projection
    this->vms_small_resolved_scales_fefct = 
       std::make_shared<TFEVectFunct3D>(projection_space_.get(), (char*)"v", (char*)"v", 
                                        &vms_small_resolved_scales[0], ndofP, 6);

 
    // create the label space, used in the adaptive method   
    int n_cells = collections_.front()->GetN_Cells();
    // initialize the piecewise constant vector
    if(projection_order == 0)
      this->label_for_local_projection.resize(n_cells, 0);
    else if(projection_order == 1)
      this->label_for_local_projection.resize(n_cells, 1);
    else
      ErrThrow("local projection space is not defined");
    
    // create fefunction for the labels such that they can be passed to the assembling routines
    label_for_local_projection_space_ = 
      std::make_shared<TFESpace3D>(collections_.front(), (char*)"label_for_local_projection", 
                                   (char*)"label_for_local_projection", example_.get_bc(2), DiscP_PSpace, 0);
    // finite element function for local projection space
    this->label_for_local_projection_fefct = 
      std::make_shared<TFEFunction3D>(label_for_local_projection_space_.get(), (char*)"vms_local_projection_space_fefct", 
                                      (char*)"vms_local_projection_space_fefct", &label_for_local_projection[0], 
                                      n_cells);
   
    // matrices for the vms method needs to assemble only on the 
    // finest grid. On the coarsest grids, the SMAGORINSKY model 
    // is used and for that matrices are available on all grids:
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1: case 2:
        ErrThrow("VMS projection cannot be supported for NSTYPE  ", 
                 TDatabase::ParamDB->NSTYPE);
        break;
      case 3: case 4: case 14:
        const TFESpace3D& velocity_space = this->get_velocity_space();
        // matrices G_tilde
        matrices_for_turb_mod.at(0) = std::make_shared<FEMatrix>(&velocity_space, projection_space_.get());
        matrices_for_turb_mod.at(1) = std::make_shared<FEMatrix>(&velocity_space, projection_space_.get());
        matrices_for_turb_mod.at(2) = std::make_shared<FEMatrix>(&velocity_space, projection_space_.get());
        // matrices G
        matrices_for_turb_mod.at(3) = std::make_shared<FEMatrix>(projection_space_.get(), &velocity_space);
        matrices_for_turb_mod.at(4) = std::make_shared<FEMatrix>(projection_space_.get(), &velocity_space);
        matrices_for_turb_mod.at(5) = std::make_shared<FEMatrix>(projection_space_.get(), &velocity_space);
        // mass matrix 
        matrices_for_turb_mod.at(6) = std::make_shared<FEMatrix>(projection_space_.get(), projection_space_.get());
        break;
    }
  }
  // prepartaion of the multigrid 
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
  
  this->output_problem_size_info();
  // initialize the defect of the system. It has the same structure as
  // the rhs (and as the solution)
  this->defect_.copy_structure(this->systems_.front().rhs_);
  for(System_per_grid& s : this->systems_)
  {
    s.u_.GetComponent(0)->Interpolate(example_.get_initial_cond(0));
    s.u_.GetComponent(1)->Interpolate(example_.get_initial_cond(1));
    s.u_.GetComponent(2)->Interpolate(example_.get_initial_cond(2));
  }
}

///**************************************************************************** */
void Time_NSE3D::interpolate()
{
//   TFEFunction3D *u1 = this->systems_.front().u_.GetComponent(0);
//   TFEFunction3D *u2 = this->systems_.front().u_.GetComponent(1);
//   TFEFunction3D *u3 = this->systems_.front().u_.GetComponent(2);
//   
//   u1->Interpolate(example_.get_initial_cond(0));
//   u2->Interpolate(example_.get_initial_cond(1));
//   u3->Interpolate(example_.get_initial_cond(2));  
}

///**************************************************************************** */
void Time_NSE3D::check_parameters()
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
 if(db_["time_discretization"].is(4))
 {
   if(solver_.is_using_multigrid())
   {
     ErrThrow("Multigrid with IMEX-scheme is not implemented yet");
   }
   else
   {
     Output::info<1>("check_parameters",
                     "The IMEX scheme has been chosen as a time discretization scheme!\n");
   }
 }

 if(TDatabase::TimeDB->TIME_DISC == 0)
 {
   ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC
         << " does not supported");
   throw("TIME_DISC: 0 is not supported");
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
      // pressure order is chosen correctly
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
  std::vector<TSquareMatrix3D*> sqMatrices;
  std::vector<TMatrix3D*>       rectMatrices;
  std::vector<double*> rhsArray(3);
  // from front to back (fine to coarse grid)
  for(System_per_grid& s : this->systems_) 
  {
    call_assembling_routine(s, LocalAssembling3D_type::TNSE3D_LinGAL);
    
    /** manage dirichlet condition by copying non-actives DoFs
     * from rhs to solution of front grid (=finest grid)
     * Note: this operation can also be done inside the loop, so that
     * the s.solution is corrected on every grid. This is the case in
     * TNSE2D.
     * TODO: CHECK WHAT IS THE DIFFERENCE BETWEEN doing this on every grid
     * and doing it only on the finest grid!
     **/
    s.solution_.copy_nonactive(s.rhs_);
    // copy the solution to the old solution for the residual computations
    // used in the RESIDUAL_VMS method
    s.Old_Sol = s.solution_;

    if(TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
    {
      std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix_.get_blocks_uniquely();
      // update mass matrix of projection
      LumpMassMatrixToDiagonalMatrix3D(matrices_for_turb_mod.at(6));
      // update stiffness matrix
      VMS_ProjectionUpdateMatrices3D(blocks, matrices_for_turb_mod);
      // reset flag for projection-based VMS method such that Smagorinsky LES method
      // is used on coarser grids 
      TDatabase::ParamDB->DISCTYPE = SMAGORINSKY_COARSE;
    }
  }// end for system per grid - the last system is the finer one (front)
  // reset   DISCTYPE to VMS_PROJECTION to be correct in the next assembling
  if(TDatabase::ParamDB->DISCTYPE == SMAGORINSKY_COARSE)
    TDatabase::ParamDB->DISCTYPE = VMS_PROJECTION;

  /** manage dirichlet condition by copying non-actives DoFs
  * from rhs to solution of front grid (=finest grid)
  * Note: this operation can also be done inside the loop, so that
  * the s.solution is corrected on every grid. This is the case in
  * TNSE2D.
  * TODO: CHECK WHAT IS THE DIFFERENCE BETWEEN doing this on every grid
  * and doing it only on the finest grid!
  * **/
  // this->systems_.front().solution_.copy_nonactive(this->systems_.front().rhs_);

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
  this->old_solution_ = this->systems_.front().solution_;
}

/**************************************************************************** */
void Time_NSE3D::assemble_rhs()
{
  System_per_grid& s = this->systems_.front();

  // TODO Should it be timesteplength or currenttimesteplength
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
//  const double theta1 = TDatabase::TimeDB->THETA1;
  const double theta2 = TDatabase::TimeDB->THETA2;
  const double theta3 = TDatabase::TimeDB->THETA3;
  const double theta4 = TDatabase::TimeDB->THETA4;
  // reset the right hand side of the grid of interest (finest)
  s.rhs_.reset();
  call_assembling_routine(s, LocalAssembling3D_type::TNSE3D_Rhs);
  /* just a temporary vector which is going to be used at the end to
   * retrieve the nonActive of rhs_. If we dont do it, the nonActive of rhs
   * will be changed during the following matrix.vector operations and we'll
   * loose the values of the Dirichlet nodes. */
  BlockVector temporary = s.rhs_;

  // now it is this->systems[i].rhs = f^k
  // scale by time step length and theta4 (only active dofs)
  s.rhs_.scaleActive(tau*theta4);

  // add rhs from previous time step
  if(theta3 != 0)
  {
    s.rhs_.addScaledActive((this->old_rhs_), tau*theta3);

    // now it is this->systems[i].rhs = tau*theta3*f^{k-1} + tau*theta4*f^k
    // next we want to set old_rhs to f^k (to be used in the next time step)
    this->old_rhs_.addScaledActive(s.rhs_, -1./(tau*theta3));
    this->old_rhs_.scaleActive(-theta3/theta4);
    this->old_rhs_.copy_nonactive(s.rhs_);
  }

  // FIXME Find other solution than this submatrix method.
  // M u^{k-1} NOTE : here s.solution_ is exactly u^{k-1}
  if(TDatabase::ParamDB->DISCTYPE == RESIDUAL_VMS)
    s.MatrixK.apply_scaled_submatrix(old_solution_, s.rhs_, 3, 3, 1.0);
  else
    s.massMatrix_.apply_scaled_submatrix(s.solution_, s.rhs_, 3, 3, 1.0);
  // -tau*theta2 * A u^{k-1} NOTE : here, s.solution_ is u^{k-1}
  double factor = -tau*theta2;
  s.matrix_.apply_scaled_submatrix(s.solution_, s.rhs_, 3, 3, factor);

  // scale the BT blocks with time step length
  for(System_per_grid& s : this->systems_)
  {
    if(tau != oldtau_)
    {
      // TODO: change the factor to be THETA1*tau; why??
      factor = /*theta1*/tau;
      if(this->oldtau_ != 0.0)
      {
        factor /= this->oldtau_;
        Output::print<1>("change in tau", this->oldtau_, "->", tau);
      }
      // scale the BT transposed blocks with the current time step
      const std::vector<std::vector<size_t>> cell_positions = {{0,3},
                                                               {1,3},
                                                               {2,3}};
      s.matrix_.scale_blocks(factor, cell_positions);
      if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
      {
        const std::vector<std::vector<size_t>> cell_positions_t = {{3,0},
                                                                   {3,1},
                                                                   {3,2}};
        s.matrix_.scale_blocks(factor, cell_positions_t);
      }
    }
  }

  this->oldtau_ = tau;

  // retrieve the non active from "temporary" into rhs vector
  s.rhs_.copy_nonactive(temporary);

  // copy the non active to the solution vector
  s.solution_.copy_nonactive(s.rhs_);

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

  /* Reset old_residual_ for this time step iteration
  * otherwise, ones compares with the old_residual_ from
  * the previous time iteration, which is not correct. */
  this->old_residual_ = FixedSizeQueue<10,Residuals>();

  Output::info<5>("Assemble_rhs()", "End of the assembling of right hand side.");
}

/**************************************************************************** */
void Time_NSE3D::assemble_nonlinear_term()
{
  // subtract the right hand which comes from the SUPG contribution, assemble and 
  // add it to the rhs vector for nonlinear iteration
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double t1 = TDatabase::TimeDB->THETA1;
  double t2 = TDatabase::TimeDB->THETA2;
  double t3 = TDatabase::TimeDB->THETA3;
  double t4 = TDatabase::TimeDB->THETA4;
  int disctype = TDatabase::ParamDB->DISCTYPE;
  for(System_per_grid& s : this->systems_)
  {
    call_assembling_routine(s, LocalAssembling3D_type::TNSE3D_NLGAL);
    // scale the standing blocks due to nonlinearity
    if(disctype == SUPG || disctype == RESIDUAL_VMS)
    {
      const std::vector<std::vector<size_t>> cell_posi = {{0,3}, {1,3}, {2,3}};
      s.matrix_.scale_blocks(t1*tau, cell_posi);
    }
    
    if(disctype == VMS_PROJECTION)
    {
      std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix_.get_blocks_uniquely();
      // update mass matrix of projection
      LumpMassMatrixToDiagonalMatrix3D(matrices_for_turb_mod.at(6));
      // update stiffness matrix
      VMS_ProjectionUpdateMatrices3D(blocks, matrices_for_turb_mod);
      // reset flag for projection-based VMS method such that Smagorinsky LES method
      // is used on coarser grids 
      TDatabase::ParamDB->DISCTYPE = SMAGORINSKY_COARSE;
    }
  }
  // update the right hand side for the next iteration: All A blocks, standing B blocks 
  // and the right hand sides needs to be re-assemble for the SUPG and RESIDUAL_VMS
  if(disctype == SUPG || disctype == RESIDUAL_VMS)
  {
    System_per_grid& s = this->systems_.front();
    s.rhs_.scaleActive(tau*t4);
    if(TDatabase::ParamDB->DISCTYPE == RESIDUAL_VMS)
      s.MatrixK.apply_scaled_submatrix(old_solution_, s.rhs_, 3, 3, 1.0);
    else 
      s.massMatrix_.apply_scaled_submatrix(old_solution_, s.rhs_, 3, 3, 1.0);
  }
  // reset   DISCTYPE to VMS_PROJECTION to be correct in the next assembling
  if(TDatabase::ParamDB->DISCTYPE == SMAGORINSKY_COARSE)
    TDatabase::ParamDB->DISCTYPE = VMS_PROJECTION;
  Output::info<5>("Assemble non linear terms", "End of the assembling of the nonlinear matrix.");
}

/**************************************************************************** */
void Time_NSE3D::assemble_system()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;

  for(System_per_grid& s : this->systems_)
  {
    const std::vector<std::vector<size_t>>
      cell_positions = {{0,0}, {0,1}, {0,2},
                        {1,0}, {1,1}, {1,2},
                        {2,0}, {2,1}, {2,2}};

    // note: declaring the auxiliary cell_positions is needed by the compiler
    // to sort out the overriding of the function scale_blocks_actives(...,...)
    s.matrix_.scale_blocks_actives(factor, cell_positions);
    if(TDatabase::ParamDB->DISCTYPE == RESIDUAL_VMS)
    {
      const FEMatrix& m00 = *s.MatrixK.get_blocks().at(0).get();
      s.matrix_.add_matrix_actives(m00, 1.0, {{0,0}}, {false});
      
      const FEMatrix& m01 = *s.MatrixK.get_blocks().at(1).get();
      s.matrix_.add_matrix_actives(m01, 1.0, {{0,1}}, {false});
      
      const FEMatrix& m02 = *s.MatrixK.get_blocks().at(2).get();
      s.matrix_.add_matrix_actives(m02, 1.0, {{0,2}}, {false});
      
      const FEMatrix& m10 = *s.MatrixK.get_blocks().at(4).get();
      s.matrix_.add_matrix_actives(m10, 1.0, {{1,0}}, {false});
      
      const FEMatrix& m11 = *s.MatrixK.get_blocks().at(5).get();
      s.matrix_.add_matrix_actives(m11, 1.0, {{1,1}}, {false});
      
      const FEMatrix& m12 = *s.MatrixK.get_blocks().at(6).get();
      s.matrix_.add_matrix_actives(m12, 1.0, {{1,2}}, {false});
      
      const FEMatrix& m20 = *s.MatrixK.get_blocks().at(8).get();
      s.matrix_.add_matrix_actives(m20, 1.0, {{2,0}}, {false});
      
      const FEMatrix& m21 = *s.MatrixK.get_blocks().at(9).get();
      s.matrix_.add_matrix_actives(m21, 1.0, {{2,1}}, {false});
      
      const FEMatrix& m22 = *s.MatrixK.get_blocks().at(10).get();
      s.matrix_.add_matrix_actives(m22, 1.0, {{2,2}}, {false});
    }
    else
    {
      const FEMatrix& mass_blocks =
         *s.massMatrix_.get_blocks().at(0).get();
         
      s.matrix_.add_matrix_actives(mass_blocks, 1.0,
                                 {{0,0}, {1,1}, {2,2}},
                                 {false, false, false});
    }
  }

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
  // hold the residual from up to 10 iterations ago
  const double veryOldNormOfResidual  = this->old_residual_.front().fullResidual;
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
    //this->old_solution_ = this->systems_.front().solution_;
  }

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

  if ( normOfResidual >= conv_speed*veryOldNormOfResidual )
  {
    slow_conv = true;
  }

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
    // descale the matrices, since only the diagonal A block will
    // be reassembled in the next time step
    if (this->imex_scheme(0) && iteration_counter>0)
      return true; // in these conditions, the matrix are already descaled
    else
    {
      this->descale_matrices();
      this->old_solution_ = this->systems_.front().solution_;
      this->systems_.front().Old_Sol = this->systems_.front().solution_;
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
    //MPI: solution in consistency level 3 (TODO: maybe this is superfluous here
    // (because solution might be in level 3 consistency already)!)
    auto comms = s.matrix_.get_communicators();
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(s.solution_.block(bl), 3);
    }
#endif

  // copy rhs to defect and compute defect
  this->defect_ = s.rhs_;
  s.matrix_.apply_scaled_add(s.solution_, defect_,-1.);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    IntoL20Vector3D(&defect_[3*number_u_Dof], number_p_Dof,
                      TDatabase::ParamDB->PRESSURE_SPACE);
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
  
#ifndef _MPI
  solver_.solve(s.matrix_, s.rhs_, s.solution_);
#endif
  
#ifdef _MPI
  if(solver_.get_db()["solver_type"].is("direct"))
  {
    //set up a MUMPS wrapper
    MumpsWrapper mumps_wrapper(s.matrix_);
    //kick off the solving process
    mumps_wrapper.solve(s.rhs_, s.solution_);
  }
  else
    this->solver_.solve(s.matrix_, s.rhs_, s.solution_);
#endif
  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11, A22 and A33 matrices are
  // reset and assembled again but the non-diagonal blocks are scaled, so
  // for the next iteration we have to descale, see assemble_system()
  this->descale_matrices();

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       s.p_.project_into_L20();
}

/**************************************************************************** */
void Time_NSE3D::descale_matrices()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;
  for(System_per_grid& s : this->systems_)
  {
    if(TDatabase::ParamDB->DISCTYPE != SUPG || 
      TDatabase::ParamDB->DISCTYPE != RESIDUAL_VMS)
    {
      const FEMatrix& mass_blocks = *s.massMatrix_.get_blocks().at(0).get();
      s.matrix_.add_matrix_actives(mass_blocks, -1.0,
                                   {{0,0}, {1,1}, {2,2}},
                                   {false, false, false});
      const std::vector<std::vector<size_t>>
      cell_positions = {{0,0}, {0,1}, {0,2},
                        {1,0}, {1,1}, {1,2},
                        {2,0}, {2,1}, {2,2}};
      // note: declaring the auxiliary cell_positions is needed by the compiler
      // to sort out the overriding of the function scale_blocks_actives(...,...)
      s.matrix_.scale_blocks_actives(1./factor, cell_positions);
    }
    // else::Descaling of the matrices are not needed 
    // for the SUPG and RESIDUAL_VMS methods since the matrices will 
    // be reassemble during the nonlinear iteration
  }
}

/**************************************************************************** */
void Time_NSE3D::output(int m, int &image)
{
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    bool no_output = !db_["output_write_vtk"] && !db_["output_compute_errors"];
    if(no_output)
      return;

  System_per_grid& s = this->systems_.front();

  TFEFunction3D* u1 = s.u_.GetComponent(0);
  TFEFunction3D* u2 = s.u_.GetComponent(1);
  TFEFunction3D* u3 = s.u_.GetComponent(2);

  if((size_t)db_["verbosity"]> 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    u3->PrintMinMax();
    s.p_.PrintMinMax();
  }

  if((m==0) || (m/TDatabase::TimeDB->STEPS_PER_IMAGE) )
  {
    if(db_["output_write_vtk"])
    {
      // last argument in the following is domain but is never used in this class
      TOutput3D output(5, 5, 2, 1, NULL);
      output.AddFEFunction(&s.p_);
      output.AddFEVectFunct(&s.u_);
#ifdef _MPI
      char SubID[] = "";
      if(my_rank == 0)
        mkdir(db_["output_vtk_directory"], 0777);
      std::string dir = db_["output_vtk_directory"];
      std::string base = db_["output_basename"];
      output.Write_ParVTK(MPI_COMM_WORLD, image, SubID, dir, base);
      image++;
#else
    // Create output directory, if not already existing.
    mkdir(db_["output_vtk_directory"], 0777);
    std::string filename = db_["output_vtk_directory"];
    filename += "/" + db_["output_basename"].value_as_string();

      if(image<10) filename += ".0000";
      else if(image<100) filename += ".000";
      else if(image<1000) filename += ".00";
      else if(image<10000) filename += ".0";
      else filename += ".";
      filename += std::to_string(image) + ".vtk";
      output.WriteVtk(filename.c_str());
      image++;
#endif
    }
  }

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
    double t = TDatabase::TimeDB->CURRENTTIME;
    // print errors
    if (my_rank == 0 )
    {
      Output::print<1>("time: ", t, " L2(u)         : ", setprecision(10), sqrt(errors_.at(0)));
      Output::print<1>("time: ", t, " H1-semi(u)    : ", setprecision(10), sqrt(errors_.at(1)));
      
      Output::print<1>("time: ", t, " L2(0,t,L2(u)) : ", setprecision(10), sqrt(errors_.at(4)));
      Output::print<1>("time: ", t, " L2(0,t,H1-semi(u)) : ", 
                       setprecision(10), sqrt(errors_.at(6)));
      Output::print<1>("time: ", t, " L2(p)      : ", setprecision(10), sqrt(errors_.at(2)));
      Output::print<1>("time: ", t, " H1-semi(p)): ", setprecision(10), sqrt(errors_.at(3)));
      
      Output::print<1>("time: ", t, " L2(0,t,L2(p)) : ", setprecision(10), sqrt(errors_.at(8)) );
      Output::print<1>("time: ", t, " L2(0,t,H1-semi(p)) : ", setprecision(10), sqrt(errors_.at(10)) );
    }
  }
   delete u1;
   delete u2;
   delete u3;

   // do post-processing step depending on what the example implements, if needed
   example_.do_post_processing(*this);

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
  double h_min, h_max;
  TCollection * coll = this->get_velocity_space().GetCollection();
  coll->GetHminHmax(&h_min, &h_max);
#ifndef _MPI
  // print out some information about number of DoFs and mesh size
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 3 * n_u + n_p; // total number of degrees of freedom
  int nActive = this->get_velocity_space().GetN_ActiveDegrees();

  Output::print("N_Cells     : ", setw(10), coll->GetN_Cells());
  Output::print("h (min,max) : ", setw(10), h_min ," ", setw(12), h_max);
  Output::print("dof Velocity: ", setw(10), 3* n_u);
  Output::print("dof Pressure: ", setw(10), n_p   );
  Output::print("dof all     : ", setw(10), n_dof );
  Output::print("active dof  : ", setw(10), 3*nActive);
#else
  int size, my_rank;
  auto comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &size);
  
  double global_hmin = 0, global_hmax=0;
  MPI_Reduce(&h_min,&global_hmin, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&h_max,&global_hmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  
  std::vector<int> n_hal_cells(size,0);// hallo cells
  std::vector<int> n_own_cells(size,0);
  int l_hal_cells = coll->GetN_HaloCells();
  int l_own_cells = coll->GetN_OwnCells();
  // gather the halo cells to root
  MPI_Gather(&l_hal_cells, 1, MPI_INT, // send
             &n_hal_cells.at(0), 1, MPI_INT,// receive
             0, comm); // control
  // gather the own cells to root
  MPI_Gather(&l_own_cells, 1, MPI_INT, // send
             &n_own_cells.at(0), 1, MPI_INT, // receive
             0, comm); // control

  auto velocity_comm = systems_.front().velocitySpace_.get_communicator();
  auto pressure_comm = systems_.front().pressureSpace_.get_communicator();
  int n_loc_master_dof_velo = velocity_comm.GetN_Master();
  std::vector<int> n_master_dof_velo(size,0);
  MPI_Gather(&n_loc_master_dof_velo, 1, MPI_INT, // send
             &n_master_dof_velo.at(0), 1, MPI_INT, // receive
             0, comm); // control
  int n_loc_master_dof_pres = pressure_comm.GetN_Master();
  std::vector<int> n_master_dof_pres(size,0);
  MPI_Gather(&n_loc_master_dof_pres, 1, MPI_INT, // send
             &n_master_dof_pres.at(0), 1, MPI_INT, // receive
             0, comm); // control
  
  if(my_rank == 0)
  {
    Output::stat("Time_NSE3D", "information on the FE space");
    size_t n_total_cells=0; size_t n_dofs_velo = 0; size_t n_dofs_pres = 0;
    for(int i=0; i<size; ++i)
    {
      Output::dash("Process ", i, "\t n_own_cells ", n_own_cells.at(i), 
                   "\t n_halo_cells ", n_hal_cells.at(i), 
                   "\t h_min/h_max ", h_min, "/", h_max,
                   "\t n_master_dof_velo ", n_master_dof_velo.at(i), 
                   "\t n_master_dof_pres ", n_master_dof_pres.at(i));
      n_total_cells += n_own_cells.at(i);
      n_dofs_velo += n_master_dof_velo.at(i);
      n_dofs_pres += n_master_dof_pres.at(i);
    }
    Output::dash("Total number of cells:                       ", n_total_cells);
    Output::dash("Total number of velcoity degrees of freedom: ", 3*n_dofs_velo);
    Output::dash("Total number of pressure degrees of freedom: ", n_dofs_pres);
    Output::dash("Total number of degrees of freedom:          ", 3*n_dofs_velo+n_dofs_pres);
  }
#endif
  
}

/**************************************************************************** */
void Time_NSE3D::construct_extrapolated_solution()
{
  this->extrapolated_solution_.reset();
  this->extrapolated_solution_ = this->old_solution_;
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
  //IMEX-scheme needs to get out of the iteration directly after the 1st solve()
  bool interruption_condition  = (db_["time_discretization"].is(4))*
                      (this->current_step_>=3);

  // change maximum number of nonlin_iterations to 1 in IMEX case
  if (interruption_condition)
  {
    db_["nonlinloop_maxit"] = 1;
    if(print_info) // condition is here just to print it once
      Output::info<1>("Nonlinear Loop MaxIteration",
                    "The parameter 'nonlinloop_maxit' was changed to 1."
                    " Only one non-linear iteration is done, because the IMEX scheme was chosen.\n");
  }
  return interruption_condition;
}

void Time_NSE3D::call_assembling_routine(Time_NSE3D::System_per_grid& s, 
                                         LocalAssembling3D_type type)
{
  // set arrays of spaces used to 
  // assemble matrices and right hand side 
  std::vector<const TFESpace3D*> spaces_mat(2);
  spaces_mat[0] = &s.velocitySpace_;
  spaces_mat[1] = &s.pressureSpace_;
  // spaces for right hand side
  std::vector<const TFESpace3D*> spaces_rhs(3);
  spaces_rhs[0] = &s.velocitySpace_;
  spaces_rhs[1] = &s.velocitySpace_;
  spaces_rhs[2] = &s.velocitySpace_;
  // finite element functions for nonlinear convective term
  // for initial assembling, they correspond to the initial condition    
  std::vector<TFEFunction3D*> fefunctions(3);
  int disctype = TDatabase::ParamDB->DISCTYPE;
  // The assembly with the extrapolated velocity of IMEX_scheme begins
  // at step 3
  bool is_imex = this->imex_scheme(0);
  if(!is_imex)
  {
    fefunctions[0]=s.u_.GetComponent(0); fefunctions[1]=s.u_.GetComponent(1);
    fefunctions[2]=s.u_.GetComponent(2); 

    if(disctype == VMS_PROJECTION)
    {
      spaces_mat.resize(4);
      // projection space
      spaces_mat[2]=projection_space_.get();
      // label space that indicates local projection space
      spaces_mat[3]=label_for_local_projection_space_.get();
      // append function for the labels of the local projection
      fefunctions.resize(4);
      fefunctions[3]=label_for_local_projection_fefct.get();
    }
    // 
    if(disctype == RESIDUAL_VMS)
    {
      fefunctions.resize(7);
      fefunctions[3]=&s.p_;
      fefunctions[4]=s.u_old.GetComponent(0);
      fefunctions[5]=s.u_old.GetComponent(1);
      fefunctions[6]=s.u_old.GetComponent(2);
    }
  }
  else
  {
    // construct the extrapolated solution 2*u(t-1)-u(t-2) in case of IMEX-scheme
    // Note : in this function, the non active Dofs are taken care of.
    // Namely, extrapolated_solution takes the nonActive of the current Rhs.
    this->construct_extrapolated_solution();
    TFEVectFunct3D extrapolated_velocity_vector(&this->systems_.front().velocitySpace_,
                                                (char*)"", (char*)"",
                                                extrapolated_solution_.block(0),
                                                extrapolated_solution_.length(0), 3);

    // Construct now the corresponding fe_functions for local assembling
    fefunctions[0]=extrapolated_velocity_vector.GetComponent(0);
    fefunctions[1]=extrapolated_velocity_vector.GetComponent(1);
    fefunctions[2]=extrapolated_velocity_vector.GetComponent(2);
  }
  // assign boundary conditions and boundary values
  std::vector<BoundCondFunct3D*>  bound_cond;
  std::vector<BoundValueFunct3D*> bound_valu;
  boundary_data(s,bound_cond,bound_valu);
  // prepare matrices and rhs for assembling 
  std::vector<TSquareMatrix3D*> sqMatrices;
  std::vector<TMatrix3D*> rectMatrices;
  std::vector<double*> rhs_array;
  prepare_matrices_rhs(s, type, sqMatrices, rectMatrices, rhs_array);
      // local assembling object - used in Assemble3D
  const LocalAssembling3D la(type, fefunctions.data(),this->example_.get_coeffs());

  // assemble all the matrices and right hand side
  Assemble3D(spaces_mat.size(), spaces_mat.data(),
             sqMatrices.size(), sqMatrices.data(),
             rectMatrices.size(), rectMatrices.data(),
             rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
             bound_cond.data(), bound_valu.data(), la);
}


void Time_NSE3D::boundary_data(System_per_grid& s, 
                               std::vector<BoundCondFunct3D*> &bc, 
                               std::vector<BoundValueFunct3D*> &bv)
{
  // assign boundary conditions
  bc.resize(3); 
  bc.at(0) = s.velocitySpace_.getBoundCondition();
  bc.at(1) = s.velocitySpace_.getBoundCondition();
  bc.at(2) = s.velocitySpace_.getBoundCondition();
  // bc.at(3) = s.pressureSpace_.getBoundCondition();
  // assign boundary values  
  bv.resize(3);
  bv.at(0) = example_.get_bd(0);
  bv.at(1) = example_.get_bd(1);
  bv.at(2) = example_.get_bd(2);
  //bv.at(3) = example_.get_bd(3);
}

void Time_NSE3D::prepare_matrices_rhs(System_per_grid& s,
                                  LocalAssembling3D_type type,
                                  std::vector<TSquareMatrix3D*> &sqMatrices,
                                  std::vector<TMatrix3D*> &rectMatrices, 
                                  std::vector<double*> &rhs_array)
{
  rectMatrices.resize(0);
  sqMatrices.resize(0);
  rhs_array.resize(0);
  
  std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix_.get_blocks_uniquely();
  int disctype = TDatabase::ParamDB->DISCTYPE;
  int nstype = TDatabase::ParamDB->NSTYPE;
  switch(type)
  {
    case LocalAssembling3D_type::TNSE3D_LinGAL:
    {
      // right had side which is common for all nstypes
      rhs_array.resize(3);
      rhs_array[0]=s.rhs_.block(0);
      rhs_array[1]=s.rhs_.block(1);
      rhs_array[2]=s.rhs_.block(2);
      // get the blocks of the mass matrix 
      std::vector<std::shared_ptr<FEMatrix>> mass_blocks
         = s.massMatrix_.get_blocks_uniquely();
  
      switch(nstype)
      {
        case 1:
          if(blocks.size() != 4)
          {
            ErrThrow("NSTYPE 1: Wrong blocks.size() ", blocks.size(), " instead of 4.");
          }
          sqMatrices.resize(2);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          // mass matrix
          sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
          
          rectMatrices.resize(3);
          rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(1).get());
          rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
          rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get());
          break;
        case 2:
          if(blocks.size() != 7)
          {
            ErrThrow("NSTYPE 2: Wrong blocks.size() ", blocks.size(), " instead of 7.");
          }
          sqMatrices.resize(2);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          // mass matrix
          sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
          
          rectMatrices.resize(6);
          rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(4).get()); //first the lying B blocks
          rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(5).get());
          rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(6).get());
          rectMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks.at(1).get()); //then the standing B blocks
          rectMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
          rectMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get());
          break;
        case 3:
          if(blocks.size() != 12)
          {
            ErrThrow("NSTYPE 3: Wrong blocks.size() ", blocks.size(), " instead of 12.");
          }
          sqMatrices.resize(10);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
          sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
          sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
          sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
          sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
          sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
          sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
          sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
          // mass matrix 
          sqMatrices[9]=reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
          
          // rectangular matrices 
          rectMatrices.resize(3);
          rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get());  // standing B blocks
          rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
          rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
          break;
        case 4:
        case 14:
          if(nstype==4 && blocks.size() != 15)
          {
            ErrThrow("NSTYPE 4: Wrong blocks.size() ", blocks.size(), " instead of 15.");
          }
          if(nstype==14 && blocks.size() != 16)
          {
            ErrThrow("NSTYPE 14: Wrong blocks.size() ", blocks.size(), " instead of 16.");
          }
          // square matrices
          sqMatrices.resize(10);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
          sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
          sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
          sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
          sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
          sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
          sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
          sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
          // mass matrix 
          sqMatrices[9] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());

          // additional matrices of the VMS method
          // mass matrix L of projection space 
          if(disctype == VMS_PROJECTION)
          {
            sqMatrices.resize(11);
            sqMatrices[10]=reinterpret_cast<TSquareMatrix3D*>(
                                            matrices_for_turb_mod.at(6).get());
          }
          // rectangular matrices 
          rectMatrices.resize(6);
          rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(12).get()); // standing B blocks
          rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(13).get());
          rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(14).get());
          rectMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get()); //than the standing B blocks
          rectMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
          rectMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
          if(disctype == VMS_PROJECTION)
          {
            rectMatrices.resize(12);
            // matrices  \tilde G_11, \tilde G_24, \tilde G_36
            rectMatrices[6]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(0).get());
            rectMatrices[7]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(1).get());
            rectMatrices[8]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(2).get());
            // matrices G_11, G_42, G_63
            rectMatrices[9]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(3).get());
            rectMatrices[10]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(4).get());
            rectMatrices[11]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(5).get());
          }
          if(disctype == RESIDUAL_VMS)
          {
            sqMatrices.resize(18);
            std::vector<std::shared_ptr<FEMatrix>> m_kblocks
               = s.MatrixK.get_blocks_uniquely();
            sqMatrices[9] =reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(0).get());
            sqMatrices[10]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(1).get());
            sqMatrices[11]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(2).get());
            sqMatrices[12]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(4).get());
            sqMatrices[13]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(5).get());
            sqMatrices[14]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(6).get());
            sqMatrices[15]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(8).get());
            sqMatrices[16]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(9).get());
            sqMatrices[17]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(10).get());
          }
          if(nstype ==14 )
          {
            rhs_array.resize(4);
            rhs_array[3] = s.rhs_.block(3);
          }
          break;
        default:
          ErrThrow("TDatabase::ParamDB->NSTYPE = ", nstype ,
               " That NSE Block Matrix Type is unknown to class Time_NSE3D.");
      }
      s.rhs_.reset();
      break; 
    } //LocalAssembling3D_type::TNSE3D_LinGAL
    case LocalAssembling3D_type::TNSE3D_NLGAL: // nonlinear
    {
      std::vector<std::shared_ptr<FEMatrix>> mass_blocks
             = s.massMatrix_.get_blocks_uniquely();
      switch(nstype)
      {
        case 1: case 2:
          blocks = s.matrix_.get_blocks_uniquely({{0,0},{1,1},{2,2}});
          sqMatrices.resize(1);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          break;
        case 3: case 4:
          switch(disctype)
          {
            case GALERKIN:
              sqMatrices.resize(3);
              blocks = s.matrix_.get_blocks_uniquely({{0,0},{1,1},{2,2}});
              sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
              sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
              sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
              break;
            case SMAGORINSKY:
            case VMS_PROJECTION:
            case SMAGORINSKY_COARSE:
            case SUPG:
              sqMatrices.resize(9);
              sqMatrices[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
              sqMatrices[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
              sqMatrices[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
              sqMatrices[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
              sqMatrices[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
              sqMatrices[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
              sqMatrices[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
              sqMatrices[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
              sqMatrices[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
              // rectangular matrices for vms projection 
              if(disctype==VMS_PROJECTION)
              {
                rectMatrices.resize(3);
                rectMatrices[0]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(0).get());
                rectMatrices[1]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(1).get());
                rectMatrices[2]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(2).get());
              }
              if(disctype == SUPG)
              {
                // weighted mass matrix
                sqMatrices.resize(10);
                sqMatrices[9]=reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
                // rectangular matrices:: standing B-blocks
                rectMatrices.resize(3);
                rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get()); 
                rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
                rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
                // right hand side
                rhs_array.resize(3);
                rhs_array[0]=s.rhs_.block(0);
                rhs_array[1]=s.rhs_.block(1);
                rhs_array[2]=s.rhs_.block(2);
                s.rhs_.reset();
              }
              break;
            case RESIDUAL_VMS:
              std::vector<std::shared_ptr<FEMatrix>> m_kblocks
               = s.MatrixK.get_blocks_uniquely();
              sqMatrices.resize(18);
              sqMatrices[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
              sqMatrices[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
              sqMatrices[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
              sqMatrices[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
              sqMatrices[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
              sqMatrices[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
              sqMatrices[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
              sqMatrices[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
              sqMatrices[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
              
              sqMatrices[9] =reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(0).get());
              sqMatrices[10]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(1).get());
              sqMatrices[11]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(2).get());
              sqMatrices[12]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(4).get());
              sqMatrices[13]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(5).get());
              sqMatrices[14]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(6).get());
              sqMatrices[15]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(8).get());
              sqMatrices[16]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(9).get());
              sqMatrices[17]=reinterpret_cast<TSquareMatrix3D*>(m_kblocks.at(10).get());
              
              rectMatrices.resize(3);
              rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get()); 
              rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
              rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
              // right hand side
              rhs_array.resize(3);
              rhs_array[0]=s.rhs_.block(0);
              rhs_array[1]=s.rhs_.block(1);
              rhs_array[2]=s.rhs_.block(2);
              s.rhs_.reset();
              break;// RESIDUAL_VMS
          }
          break;
      } // switch over NSTYPE
      break; // nonlinear
    }//LocalAssembling3D_type::TNSE3D_NLGAL:
    case LocalAssembling3D_type::TNSE3D_Rhs:
    {
      // no matrices to be assembled
      sqMatrices.resize(0);
      rectMatrices.resize(0);
      // right hand side 
      rhs_array.resize(3);
      rhs_array[0]= s.rhs_.block(0);
      rhs_array[1]= s.rhs_.block(1);
      rhs_array[2]= s.rhs_.block(2);
      if(nstype == 14) // TODO remove the case 4: no need the pressure block 
      {
        rhs_array.resize(4);
        rhs_array[3]=s.rhs_.block(3);
      }
      s.rhs_.reset();
      break;
    } //LocalAssembling3D_type::TNSE3D_Rhs:
  }  
  // reset matrices
  for(auto mat : sqMatrices)
    mat->reset();
  for(auto remat : rectMatrices)
    remat->reset();
}
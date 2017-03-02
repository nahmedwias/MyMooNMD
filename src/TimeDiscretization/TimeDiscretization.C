#include <TimeDiscretization.h>
#include <Database.h>
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

/* ****************************************************************************/
ParameterDatabase TimeDiscretization::default_TimeDiscretization_database()
{
  ParameterDatabase db("default ParMooN time parameters database");

  db.add("time_start", 0.,
         "This is the start time of the simulation. The total simulated time "
         "interval is [start_time, end_time].", -1000.,1000.);

  db.add("time_end", 1.,
         "This is the end time of the simulation. The total simulated time "
         "interval is [start_time, end_time].", -1000.,1000.);

  db.add("time_step_length", 0.05,
         "This is the time step length. Without time adaptivity, it "
         "is the same constant value for all time steps, whereas for "
         "time adaptivity, it only corresponds to the initial value.",
         0., 0.5);

  db.add("extrapolate_velocity", "constant_extrapolate",
         "This is the extrapolate velcotity parameter"
         "These values can be used with BDF's schemes or Crank_Nicholson",
         {"constant_extrapolate", "linear_extrapolate", "quadratic_extrapolate"});

  db.add("ansatz_test_extrapolate", "no_extrapolation",
         "This is the parameter for deciding the extrapolation in test"
         " or ansatz functions. This is due to the non-linear terms appearing"
        "in the SUPG or Residual Based VMS methods",
         {"only_velocity_test", "no_extrapolation"});

  db.add("time_discretization", "backward_euler",
         "Determine the time stepping scheme: currently supported are",
         {"backward_euler", "crank_nicolson", "bdf_two", "fractional_step"});

  return db;
}

/* ****************************************************************************/
TimeDiscretization::TimeDiscretization(const ParameterDatabase & param_db)
: db(default_TimeDiscretization_database())
{
  db.merge(param_db, false);

  if(db["time_discretization"].is("bdf_two"))
  {
    bdf_coefficients = {4./3., -1./3., 2./3.};
    // pre-stage computation is done either by the
    // backwar-Euler or Crank-Nicolson scheme
    // currently it's done by the backward Euler scheme.
    // This parameter is set to False in the main code
    // for the second time step
    pre_stage_bdf = true;
  }
  else
  {
    rk = std::make_shared<RungeKuttaTable>(db["time_discretization"]);
  }
  current_time_step_length = db["time_step_length"];
  TDatabase::TimeDB->TIMESTEPLENGTH = current_time_step_length;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = current_time_step_length;
  
  n_scale_block = 0;
  
  b_bt_linear_nl= "linear";

  Output::print<5>("Time discretization is initialized with ",
                db["time_discretization"]);
}

/* ****************************************************************************/
void TimeDiscretization::prepare_timestep(BlockFEMatrix& system_matrix,
  const BlockFEMatrix& mass_matrix, std::vector<BlockVector> & rhs,
  const std::vector<BlockVector> old_solutions)
{
  //TODO: here the idea can be implemented of multisteps
}

/* ****************************************************************************/
void TimeDiscretization::prepare_rhs_from_time_disc(
     BlockFEMatrix& system_matrix, const BlockFEMatrix& mass_matrix,
     std::vector<BlockVector> & rhs, const std::vector<BlockVector> old_solutions)
{
  // at the moment it's implemented in the way it was done in MooNMD
  // but have to be adopted to be fit according to the Butscher table
  if(db["time_discretization"].is("bdf_two") && !pre_stage_bdf)
  {
    // right hand side with tau*..
    rhs[0].scaleActive(bdf_coefficients[2]*current_time_step_length);

    mass_matrix.apply_scaled_add_actives(old_solutions[0], rhs[0],
                                           bdf_coefficients[0]);
    mass_matrix.apply_scaled_add_actives(old_solutions[1], rhs[0],
                                           bdf_coefficients[1]);
    Output::print<5>(RED, "bdf2_stage ", pre_stage_bdf, BLACK);
  }
  else if(db["time_discretization"].is("backward_euler") || pre_stage_bdf)
  {
    rhs[0].scaleActive(current_time_step_length);
    // mass matrix times old solution goes to right hand side
    mass_matrix.apply_scaled_add_actives(old_solutions[0], rhs[0], 1.);

    if(db["time_discretization"].is("bdf_two") && pre_stage_bdf)
    {
      Output::print<5>("First step in BDF2 scheme is performed by the BDF1 scheme");
    }
  }
  else if(db["time_discretization"].is("crank_nicolson"))
  {
    ErrThrow("not yet implemented");
  }
}

/* ****************************************************************************/
void TimeDiscretization::prepare_system_matrix(
  BlockFEMatrix& system_matrix, const BlockFEMatrix& mass_matrix, 
  unsigned int it_counter)
{
  if(db["time_discretization"].is("backward_euler") || pre_stage_bdf)
  {
    const std::vector<std::vector<size_t>> cells = {{0,0},{0,1},{1,0},{1,1}};
    system_matrix.scale_blocks_actives(current_time_step_length, cells);
  }
  else if(db["time_discretization"].is("bdf_two") && !pre_stage_bdf)
  {
    const std::vector<std::vector<size_t>> cells = {{0,0},{0,1},{1,0},{1,1}};
    system_matrix.scale_blocks_actives(current_time_step_length*bdf_coefficients[2], cells);
  }
  else if(db["time_discretization"].is("crank_nicolson"))
  {
    const std::vector<std::vector<size_t>> cells = {{0,0},{0,1},{1,0},{1,1}};
    system_matrix.scale_blocks_actives(current_time_step_length*0.5, cells);
  }
  else
    ErrThrow("Time stepping scheme ", db["time_discretization"], " is not supported");
  // add the scaled block matrices
  system_matrix.add_blockfe_matrix(mass_matrix);

  /// @DETAILS: scaling of the B and BT's blocks (also "C")
  double factor;
  if(db["time_discretization"].is("backward_euler") || pre_stage_bdf)
    factor = current_time_step_length;
  else if(db["time_discretization"].is("bdf_two") && !pre_stage_bdf)
    factor = current_time_step_length*bdf_coefficients[2];
  else if(db["time_discretization"].is("crank_nicolson"))
    factor = current_time_step_length*0.5;//this have to be adopted according to the stages
    else
      ErrThrow("Time stepping scheme ", db["time_discretization"], " is not supported");
  //BEGIN DEBUG
  // cout<< "pre_stage_bdf  " << pre_stage_bdf << " time_discretization " 
  // << db["time_discretization"] << "  " << factor<< endl;
  // system_matrix.get_blocks().at(6)->Print("B1T");exit(0);
  //END DEBUG
      
  // standard case when all B-blocks have to be scaled only once for all 
  // time steps if equi-disctinct time stepping is used
  if( (b_bt_linear_nl.compare("linear")==0) && (current_step_ == 1) 
      && (it_counter==0))
  {// scale the B, BT blocks once 
    const std::vector<std::vector<size_t>> cells = {{0,2},{1,2},{2,0},{2,1}};
    system_matrix.scale_blocks(factor, cells);
  }
  else if(b_bt_linear_nl.compare("nonlinear")==0)
  {
    // inf-sup stable case 
    if(n_scale_block==2)
    {// scale all four blocks at first nonlinear-iteration step 
      if(it_counter==0){
        const std::vector<std::vector<size_t>> cells = {{0,2},{1,2},{2,0},{2,1}};
        system_matrix.scale_blocks(factor, cells);
      }
      else{//scale only the BT blocks due to the re-assembling 
        const std::vector<std::vector<size_t>> cells = {{0,2},{1,2}};
        system_matrix.scale_blocks(factor, cells);
      }
    }
    if(n_scale_block==5)
    {//
      const std::vector<std::vector<size_t>> cells = {{0,2},{1,2},{2,0},{2,1}, {2,2}};
      system_matrix.scale_blocks(factor, cells);
    }
  }
  else if(b_bt_linear_nl.compare("solution_dependent") ==0)
  {
    // inf-sup case
    if(n_scale_block==4){
      if(it_counter==0){
        const std::vector<std::vector<size_t>> cells = {{0,2},{1,2},{2,0},{2,1}};
        system_matrix.scale_blocks(factor, cells);
      }
      else{
        const std::vector<std::vector<size_t>> cells = {{2,0},{2,1}};
        system_matrix.scale_blocks(factor, cells);
      }
    }
    else{
      if(it_counter==0){
        const std::vector<std::vector<size_t>> cells = {{0,2},{1,2},{2,0},{2,1}, {2,2}};
        system_matrix.scale_blocks(factor, cells);
      }
      else{
        const std::vector<std::vector<size_t>> cells = {{2,0},{2,1}};
        system_matrix.scale_blocks(factor, cells);
      }
    }
  }
  else{
    ErrThrow("Please check the scaling of the B, BT blocks");
  }
  
  Output::print<5>("addition of the mass and stiffness matrix done ");
}
/* ****************************************************************************/
void TimeDiscretization::reset_linear_matrices(BlockFEMatrix& matrix,
                                               const BlockFEMatrix& mass)
{
  double step_length = db["time_step_length"];
  // first subtract the mass matrix from the matrix
  matrix.add_blockfe_matrix(mass, -1.);
  double factor;
  if(db["time_discretization"].is("backward_euler") || pre_stage_bdf)
    factor = step_length;
  else if(db["time_discretization"].is("bdf_two") && !pre_stage_bdf)
    factor = bdf_coefficients[2]*step_length;
  else if(db["time_discretization"].is("crank_nicolson"))
    factor = step_length*0.5;
  else
    ErrThrow("Time stepping scheme ", db["time_discretization"], " is not supported");

  // reset the matrices
  const std::vector<std::vector<size_t>> cells = {{0,0},{0,1},{1,0},{1,1}};
  matrix.scale_blocks_actives(1./factor, cells);
  // No need to reset the B, BT and C blocks because they are only scaled once during
  // the time step or re-assembled and scaled in the nonlinear case.
}

/* ****************************************************************************/
void TimeDiscretization::set_time_disc_parameters()
{
  // the first step needs to be done by the backward Euler
  // or Crank-Nicolson scheme
  if(db["time_discretization"].is("backward_euler"))
  {
    TDatabase::TimeDB->THETA1 = 1.0;
    TDatabase::TimeDB->THETA2 = 0.0;
    TDatabase::TimeDB->THETA3 = 0.0;
    TDatabase::TimeDB->THETA4 = 1.0;
  }
  if(db["time_discretization"].is("bdf_two") && (current_step_<=1) )
  {
    TDatabase::TimeDB->THETA1 = 1.0;
    TDatabase::TimeDB->THETA2 = 0.0;
    TDatabase::TimeDB->THETA3 = 0.0;
    TDatabase::TimeDB->THETA4 = 1.0;
  }
  if(db["time_discretization"].is("bdf_two") &&
          (this->current_step_ >= 2))
  {
    this->pre_stage_bdf = false;
    TDatabase::TimeDB->THETA1 = 2./3;
    TDatabase::TimeDB->THETA2 = 0.0;
    TDatabase::TimeDB->THETA3 = 0.0;
    TDatabase::TimeDB->THETA4 = 2./3.;
  }
  // set the global parameters used in the local assembling routines 
  // and at some other places
  TDatabase::TimeDB->TIMESTEPLENGTH = current_time_step_length;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = current_time_step_length;
  Output::print<5>("BDF2 scheme with the steps perform ", pre_stage_bdf);
}

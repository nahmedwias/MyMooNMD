#include "ROM_TCDR2D_Rep.h"
#include "Database.h"
#include "Assemble2D.h"
#include "AuxParam2D.h"
#include "MainUtilities.h"

/**************************************************************************** */
ParameterDatabase get_default_tcd2d_rom_parameters()
{
  Output::print<5>("creating a default tcd2d_rom parameter database");
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("tcd2d_rom parameter database");
  db.merge(LocalAssembling2D::default_local_assembling_database(), true);
  db.merge(TimeDiscretization::default_TimeDiscretization_database(), true);
  return db;
}

/**************************************************************************** */
ROM_TCDR2D_Rep::System_per_grid::System_per_grid(const Example_TimeCD2D& example,
                                            TCollection& coll)
: space_(&coll, "space", "time_cd2d space", example.get_bc(0),
           TDatabase::ParamDB->ANSATZ_ORDER)
{
  gramian_matrix_=BlockFEMatrix::CD2D(space_);
  rhs_ = BlockVector(gramian_matrix_, true);
  solution_ = BlockVector(gramian_matrix_, false);
  fe_function_ = TFEFunction2D(&space_, "c", "c",
                               solution_.get_entries(), solution_.length());  
}
/**************************************************************************** */
ROM_TCDR2D_Rep::ROM_TCDR2D_Rep(const TDomain& domain, 
                               const ParameterDatabase& param_db, 
                               const Example_TimeCD2D& ex)
: ROM(param_db),
  db(get_default_tcd2d_rom_parameters()), solver(param_db),systems(), 
  example(ex), time_stepping_scheme(param_db), errors(5, 0.), outputWriter(param_db)
{
  db.merge(param_db, false);
  this->set_parameters();
  
  bool usingMultigrid = this->solver.is_using_multigrid();
  if(!usingMultigrid)
  {
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, -4711);
    // create finite element space and function, a matrix, rhs, and solution
    this->systems.emplace_back(this->example, *coll);
    
    TFESpace2D& space = this->systems.front().space_;
    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);
    Output::print<1>("N_Cells    : ", setw(12), coll->GetN_Cells());
    Output::print<1>("h(min,max) : ", setw(12), hmin, " ", setw(12), hmax);
    Output::print<1>("dof        : ", setw(12), space.GetN_DegreesOfFreedom());
    Output::print<1>("active dof : ", setw(12), space.GetN_ActiveDegrees());
    
    // interpolate the initial data
    TFEFunction2D & fe_function = this->systems.front().fe_function_;
    fe_function.Interpolate(example.get_initial_cond(0));
    
    // add the fe function to the output object
    outputWriter.add_fe_function(&fe_function);
  }
  else
  {
    ErrThrow(" multigrid is not supported yet ");
  }
}

/**************************************************************************** */
void ROM_TCDR2D_Rep::set_parameters()
{
  Output::print("set me");
}

/**************************************************************************** */
void ROM_TCDR2D_Rep::compute_initial_solution()
{
  if(rom_db["rom_init_regularized"])
  {
    /* see S.Giere, PhD Thesis 2016 */
    Output::print<1>("Type of ROM initial condition: regularized");
    /* assemble and reduce the stiffness matrix (without diffusion parameter) */
    LocalAssembling_type la_type = LocalAssembling_type::TCDGradGradOnly;
    assemble_and_reduce(la_type);
    // reduce the stiffness matrix
    ublas::matrix<double> reduced_h1_mat = ROM::reduce(temp_matrix.get_combined_matrix());
    ublas::vector<double> reduced_h1_mat_mean = ROM::reduce_mat_mean(temp_matrix.get_combined_matrix());
    double mu = rom_db["differential_filter_width"];
    /* Set system matrix and system rhs for Helmholtz equation */
    ublas::matrix<double> initial_system_matrix = 
    mass_mat_ + mu*mu * reduced_h1_mat;
    // set initial right hand side
    ublas::vector<double> initial_system_rhs = prod(mass_mat_, red_sol_);
    initial_system_rhs -= mu * mu * reduced_h1_mat_mean;
    
    // solve Helmholtz equation
    red_sol_ = ROM::solve(initial_system_matrix, initial_system_rhs);
  }
  else
  {
    Output::print<1>("Reducing initial condition...");    
    /* reduce initial solution */
    ROM::reduce_solution(systems.front().solution_.get_entries(),red_sol_);
  }
}

/**************************************************************************** */
void ROM_TCDR2D_Rep::assemble_matrices_rhs(bool mat_depend_time)
{
  if(mat_depend_time)
  {
    // assemble the stiffness matrix and right hand side
    assemble_and_reduce(LocalAssembling_type::TCDStiffRhs);
    // assemble the mass matrix 
    assemble_and_reduce(LocalAssembling_type::TCDMassOnly);
  }
  else
  {
    // assemble the right hand side only
    assemble_and_reduce(LocalAssembling_type::TCDRhsOnly);
  }
  
}

/**************************************************************************** */
void ROM_TCDR2D_Rep::set_system_matrix()
{
  sys_mat_.clear();
  double tau = time_stepping_scheme.get_step_length();
  if(db["time_discretization"].is("backward_euler"))
    sys_mat_ = tau*cdr_mat_;
  else if(db["time_discretization"].is("crank_nicolson"))
    sys_mat_ = tau*0.5*cdr_mat_;
  else
  {
    ErrThrow("time discretization: ", db["time_discretization"], " is not supported");
  }
  sys_mat_ += mass_mat_;
  /* NOTE prepare_solve() makes only sense for time-independent problem coefficients */
  ROM::prepare_solve(sys_mat_);
}

/**************************************************************************** */
void ROM_TCDR2D_Rep::set_system_rhs(bool reassemble)
{
  double tau = time_stepping_scheme.get_step_length();
  sys_rhs_.clear();
  sys_rhs_ = prod(mass_mat_, red_sol_);

  if(db["time_discretization"].is("backward_euler"))
    sys_rhs_ -= tau * prod(cdr_mat_, red_sol_);
  else if(db["time_discretization"].is("crank_nicolson"))
  {
    sys_rhs_ -= tau * 0.5 * prod(cdr_mat_, red_sol_);
    sys_rhs_ += tau * 0.5 * reduced_rhs_;
  }
  else
  {
    ErrThrow("time discretization: ", db["time_discretization"], " is not supported");
  }
  if(rom_db["pod_fluct"])
    sys_rhs_ -= tau * cdr_mat_mean_;
  /* for time-independent source term, the next line assemble_reduce(true) can be commented out
   * NOTE: Online (= within the time loop) assembling involving FE dimension is not efficient
   * in the ROM context. To avoid it, one could alternatively store all FE coefficients of the
   * source term and reduce them offline, and use here the corresponding reduced-order vectors.
   * For source terms that are formulated in the separated time-space form one can pre-assemble
   * and reduce the space part before the time loop and multiply it in every time step with the
   * appropriate temporal factor.
   */
  if(reassemble)
    assemble_and_reduce(LocalAssembling_type::TCDRhsOnly);
  if(db["time_discretization"].is("backward_euler"))
    sys_rhs_ += tau * reduced_rhs_;
  else if(db["time_discretization"].is("crank_nicolson"))
    sys_rhs_ += tau * 0.5 * reduced_rhs_;
}

/**************************************************************************** */
void ROM_TCDR2D_Rep::solve()
{
  red_sol_ = ROM::solve(sys_rhs_);
}

/**************************************************************************** */
void ROM_TCDR2D_Rep::assemble_and_reduce(LocalAssembling_type type)
{
  temp_matrix = BlockFEMatrix::CD2D(systems.front().space_);
  int n_fespaces = 1;
  const TFESpace2D * _fe_space_ = &systems.front().space_;  
  std::vector<TSquareMatrix2D*> sqMatrices;
  std::vector<double*> rhs_array(0);
  std::vector<std::shared_ptr<FEMatrix>> block;
  switch(type)
  {
    case LocalAssembling_type::TCDMassOnly:
      sqMatrices.resize(1);
      block = temp_matrix.get_blocks_uniquely();
      sqMatrices[0] = {reinterpret_cast<TSquareMatrix2D*>(block.at(0).get())};
      break;
    case LocalAssembling_type::TCDGradGradOnly:
      if (rom_db["pod_inner_product"].get<std::string>() == "l2")
      {
        sqMatrices.resize(1);
        block = temp_matrix.get_blocks_uniquely();
        sqMatrices[0] = {reinterpret_cast<TSquareMatrix2D*>(block.at(0).get())};
      }
      else if (rom_db["pod_inner_product"].get<std::string>() == "eucl")
      {
        Output::print<1>("For given parameter 'pod_inner_product' no assembling needed.");
        return;
      }
      break;
    case LocalAssembling_type::TCDRhsOnly:
      rhs_array.resize(1);
      rhs_array[0] = systems.front().rhs_.get_entries();
      break;
    case LocalAssembling_type::TCDStiffRhs:
      sqMatrices.resize(1);
      block = temp_matrix.get_blocks_uniquely();
      sqMatrices[0] = {reinterpret_cast<TSquareMatrix2D*>(block.at(0).get())};
      rhs_array.resize(1);
      rhs_array[0] = systems.front().rhs_.get_entries();      
      break;
    default:
      ErrThrow("LocalAssembling_type ", type, " not supported");
  }
  for(auto sm : sqMatrices)
    sm->reset();
  // local assembling 
  TFEFunction2D * pointer_to_function = &systems.front().fe_function_;
  LocalAssembling2D la(db, type, &pointer_to_function, example.get_coeffs());
  // boundary conditions and boundary values
  BoundCondFunct2D * bc = _fe_space_->get_boundary_condition();
  BoundValueFunct2D * bv[1] {example.get_bd()[0]};
  
  Assemble2D(n_fespaces, &_fe_space_, sqMatrices.size(), sqMatrices.data(), 
             0, nullptr, rhs_array.size(), rhs_array.data(),
             &_fe_space_, &bc, bv, la);
  
  switch(type)
  {
    case LocalAssembling_type::TCDMassOnly:
      Output::print<5>("Reducing finite element mass matrix...");
      mass_mat_ = ROM::reduce(temp_matrix.get_combined_matrix());
      ROM::set_gramian(temp_matrix.get_combined_matrix());
      break;
    case LocalAssembling_type::TCDRhsOnly:
      Output::print<5>("Assembling finite element source term...");
      this->reduced_rhs_ = ROM::reduce(systems.front().rhs_.get_entries());
      break;
    case LocalAssembling_type::TCDStiffRhs:
      Output::print<5>("Reducing finite element convection-diffusion-reaction matrix and rhs...");
      cdr_mat_ = ROM::reduce(temp_matrix.get_combined_matrix());
      cdr_mat_mean_ = ROM::reduce_mat_mean(temp_matrix.get_combined_matrix());
      // reduce right hand side
      Output::print<5>("Reducing source term...");
      reduced_rhs_ = ROM::reduce(systems.front().rhs_.get_entries());
      break;
    default:
      ErrThrow("LocalAssembling_type ", type, " not supported");
  }
}

/**************************************************************************** */
void ROM_TCDR2D_Rep::output(int time_step)
{
  bool no_output = !rom_db["output_write_vtk"] && !rom_db["output_compute_errors"];
  if(no_output)
    return;
  System_per_grid& s = this->systems.front();

  ROM::get_full_solution(this->red_sol_, s.solution_.get_entries());
  TFEFunction2D & _fe_function = s.fe_function_;
  _fe_function.PrintMinMax();

  if( time_step % TDatabase::TimeDB->STEPS_PER_IMAGE == 0 )
  {
    Output::print<1>("Writing vtk output at time ",time_stepping_scheme.current_time_);
    // write output
    this->outputWriter.write(time_stepping_scheme.current_time_);
  }

  if(rom_db["output_compute_errors"])
  {
    double loc_e[5];
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D* space = _fe_function.GetFESpace2D();
    
    _fe_function.GetErrors(this->example.get_exact(0), 3, AllDerivatives, 4,
                           SDFEMErrors, this->example.get_coeffs(), &aux, 1,
                           &space, loc_e);
    
    Output::print<1>("time: ", time_stepping_scheme.current_time_);
    Output::print<1>("  L2: ", loc_e[0]);
    Output::print<1>("  H1-semi: ", loc_e[1]);
    double tau = time_stepping_scheme.get_step_length();
    errors[0] += (loc_e[0]*loc_e[0] + errors[1])*tau*0.5;
    errors[1] = loc_e[0]*loc_e[0];
    Output::print<1>("  L2(0,T;L2) ", sqrt(errors[0]));
    errors[2] += (loc_e[1]*loc_e[1] + errors[3])*tau*0.5;
    errors[3] = loc_e[1]*loc_e[1];
    Output::print<1>("  L2(0,T;H1) ", sqrt(errors[2]));
    
    if(errors[4] < loc_e[0])
      errors[4] = loc_e[0];
    Output::print<1>("  Linfty(0,T;L2) ", errors[4]);
    
    double min, max;
    s.fe_function_.MinMax(min,max);
    std::string error_file = rom_db["outfile"].get<std::string>() + ".errors";
    /** @todo Add additional functionality in namespace Output to write information
     * into seond (or generally several) output file (like here Output::redirect(error_file)) */
    
    //TODO: Check the redirect function for true and false case
    double start_time = db["time_start"];
    if (time_stepping_scheme.current_time_== start_time)
      Output::redirect(error_file);//, true);
    // else
    //   Output::redirect(error_file, false);
    // Output: "time" "L2 error" "H01 error" "min" "max"
    Output::print<1>(time_stepping_scheme.current_time_, " ",
                     loc_e[0], " ", loc_e[1], " ", min, " ", max);
    Output::resetOutfile();
  }
}

/**************************************************************************** */

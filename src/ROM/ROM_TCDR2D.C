#include <ROM_TCDR2D.h>
#include <Database.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <AlgebraicFluxCorrection.h>
#include <Assemble2D.h>
#include <LocalProjection.h>
#include <AuxParam2D.h>
#include <TimeDiscRout.h>
#include <TCD2D_POD_LocalAssemble.h>

/**************************************************************************** */
ROM_TCDR2D::ROM_TCDR2D(TCollection& coll, const ParameterDatabase& param_db,
		const Example_TimeCD2D& ex)
  : ROM(param_db),
    example(ex),
    fe_space(&coll, (char*)"space", (char*)"cd2d space", example.get_bc(0), 
             TDatabase::ParamDB->ANSATZ_ORDER),
    gramian_matrix({&fe_space}),
	full_sol(this->gramian_matrix, false),
    fe_function(&this->fe_space, (char*)"c", (char*)"c", this->full_sol.get_entries(),
    		     this->full_sol.length()),
    outputWriter(param_db),
	errors(5, 0.0)
{
  // TODO initialize all ublas matrices and vectors
  set_parameters();

  // this is the L^inf(L^2(Omega)) error, initialize with some large number
  // CHECK IF IT IS THE ERROR WE WANT TO MEASURE
  errors[4] = 1.e10;

  // Assemble gramian matrix and set it for ROM class
  assemble_gramian();
  ROM::set_gramian(this->gramian_matrix.get_combined_matrix());


  int ndof = this->fe_space.GetN_DegreesOfFreedom();
  double hmin, hmax;
  coll.GetHminHmax(&hmin, &hmax);
  Output::print<1>("N_Cells    : ", setw(12), coll.GetN_Cells());
  Output::print<1>("h(min,max) : ", setw(12), hmin, " ", setw(12), hmax);
  Output::print<1>("dof        : ", setw(12), ndof);
  Output::print<1>("active dof : ", setw(12), this->fe_space.GetN_ActiveDegrees());

  // Initialize the output object, add the fe function to it.
  outputWriter.add_fe_function(&this->fe_function);

  assemble_reduce();
}
/**************************************************************************** */
ROM_TCDR2D::~ROM_TCDR2D()
{
}

/**************************************************************************** */
void ROM_TCDR2D::set_parameters()
{

  if(TDatabase::TimeDB->TIME_DISC == 0)
  {
    ErrThrow("TIME_DISC: ", TDatabase::TimeDB->TIME_DISC, " is not supported.\n"
    		 "Choose 1 (backward Euler) or 2 (Crank-Nicolson).");
  }

  //set problem_type to Time_CD if not yet set
  if(!rom_db["problem_type"].is(2))
  {
    if (rom_db["problem_type"].is(0))
    {
      rom_db["problem_type"] = 2;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to Time_CD."
          "It is now reset to the correct value for Time_CD (=2).");
      rom_db["problem_type"] = 2;
    }
  }

  // an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    ErrThrow("Ansatz order 0 is no use in convection diffusion "
    		 "reaction problems! (Vanishing convection and diffusion term).");
  }

  // time discretization parameters
  if(TDatabase::TimeDB->TIME_DISC == 1)
    {
      TDatabase::TimeDB->THETA1 = 1.0;
      TDatabase::TimeDB->THETA2 = 0.0;
      TDatabase::TimeDB->THETA3 = 0.0;
      TDatabase::TimeDB->THETA4 = 1.0;
    }
    else if(TDatabase::TimeDB->TIME_DISC == 2)
    {
      TDatabase::TimeDB->THETA1 = 0.5;
      TDatabase::TimeDB->THETA2 = 0.5;
      TDatabase::TimeDB->THETA3 = 0.5;
      TDatabase::TimeDB->THETA4 = 0.5;
    }
}

/**************************************************************************** */
void ROM_TCDR2D::assemble_gramian()
{
  Output::print<1>("Assembling the gramian matrix for reduction of full-order solution...\n"
		           "(Note that gramian matrix has to be the same as used for POD computation)");

  this->gramian_matrix = BlockFEMatrix::CD2D(this->fe_space);

  if (rom_db["pod_inner_product"].get<std::string>() == "l2")
  {
	  this->assemble(this->gramian_matrix, mat_p_q);
  }
  else if (rom_db["pod_inner_product"].get<std::string>() == "eucl")
  {
	  Output::print<1>("For given parameter 'pod_inner_product' no assembling needed.");
	  return;
  }
}

/**************************************************************************** */
void ROM_TCDR2D::assemble_reduce(bool only_rhs)
{
  BlockFEMatrix fe_mass_mat = BlockFEMatrix::CD2D(fe_space);
  BlockFEMatrix fe_cdr_mat  = BlockFEMatrix::CD2D(fe_space);
  BlockVector fe_rhs(fe_cdr_mat, true);
  
  // bool supg = (TDatabase::ParamDB->DISCTYPE==2) ? true : false;
  bool supg = (rom_db["space_discretization_type"].is("supg")) ? true : false; 
  
  if (!only_rhs)
  {
    Output::print<1>("Assembling finite element mass matrix...");
    this->assemble(fe_mass_mat, (supg) ? mat_p_q_supg : mat_p_q);
    Output::print<1>("Reducing finite element mass matrix...");
    this->mass_mat = ROM::reduce(fe_mass_mat.get_combined_matrix());
    
    Output::print<1>("Assembling finite element convection-diffusion-reaction matrix...");
    this->assemble(fe_cdr_mat, (supg) ? mat_cdr_supg : mat_cdr);
    
    Output::print<1>("Reducing finite element convection-diffusion-reaction matrix...");
    this->cdr_mat      = ROM::reduce(fe_cdr_mat.get_combined_matrix());
    
    Output::print<1>("Reducing finite element convection-diffusion-reaction matrix times snaps mean...");
    this->cdr_mat_mean = ROM::reduce_mat_mean(fe_cdr_mat.get_combined_matrix());
  }
  Output::print<1>("Assembling finite element source term...");
  this->assemble(fe_rhs, (supg) ? rhs_f_q_supg : rhs_f_q);
  
  Output::print<1>("Reducing source term...");
  this->source = ROM::reduce(fe_rhs.get_entries());
}

/**************************************************************************** */
void ROM_TCDR2D::compute_initial_solution()
{
  Output::print<1>("Computing initial condition...");

  /** interpolate initial condition (if given in example file) */
  TFEFunction2D & fe_function = this->fe_function;
  fe_function.Interpolate(example.get_initial_cond(0));

  /* reduce initial solution */
  ROM::reduce_solution(this->full_sol.get_entries(),this->red_sol);

  if (rom_db["rom_init_regularized"])
  {
    /* see S.Giere, PhD Thesis 2016 */
    Output::print<1>("Type of ROM initial condition: regularized");
    
    /* assemble and reduce the stiffness matrix (without diffusion parameter) */
    BlockFEMatrix h1_mat = BlockFEMatrix::CD2D(fe_space);
    this->assemble(h1_mat, mat_gradp_gradq);
    ublas::matrix<double> red_h1_mat      = ROM::reduce(h1_mat.get_combined_matrix());
    ublas::vector<double> red_h1_mat_mean = ROM::reduce_mat_mean(h1_mat.get_combined_matrix());
    
    double mu =rom_db["differential_filter_width"];
    
    /* Set system matrix and system rhs for Helmholtz equation */
    ublas::matrix<double> init_sys_mat = this->mass_mat + mu*mu * red_h1_mat;
    
    ublas::vector<double> init_sys_rhs = prod(this->mass_mat,this->red_sol);
    init_sys_rhs -=	mu*mu * red_h1_mat_mean;

    /* Solve Helmholtz equation */
    this->red_sol = ROM::solve(init_sys_mat, init_sys_rhs);
  }
  else
    Output::print<1>("Type of ROM initial condition: standard");
}

/**************************************************************************** */
void ROM_TCDR2D::set_system_matrix()
{
  this->sys_mat.clear();
  this->sys_mat  = TDatabase::TimeDB->THETA1 * TDatabase::TimeDB->TIMESTEPLENGTH * this->cdr_mat;
  this->sys_mat += this->mass_mat;
  /* NOTE prepare_solve() makes only sense for time-independent problem coefficients */
  ROM::prepare_solve(sys_mat);
}

/**************************************************************************** */
void ROM_TCDR2D::set_system_rhs(bool reassemble)
{
  this->sys_rhs.clear();
  this->sys_rhs  = prod(this->mass_mat, this->red_sol);
  this->sys_rhs -= TDatabase::TimeDB->THETA2 *
		           TDatabase::TimeDB->TIMESTEPLENGTH *
				   prod( this->cdr_mat, this->red_sol);
  if(rom_db["pod_fluct"])
  {
	this->sys_rhs -= (TDatabase::TimeDB->THETA1 + TDatabase::TimeDB->THETA2) *
			          TDatabase::TimeDB->TIMESTEPLENGTH *
					  this->cdr_mat_mean;
  }
  this->sys_rhs +=  TDatabase::TimeDB->THETA3 * TDatabase::TimeDB->TIMESTEPLENGTH * this->source;

  /* for time-independent source term, the next line assemble_reduce(true) can be commented out
   * NOTE: Online (= within the time loop) assembling involving FE dimension is not efficient
   * in the ROM context. To avoid it, one could alternatively store all FE coefficients of the
   * source term and reduce them offline, and use here the corresponding reduced-order vectors.
   * For source terms that are formulated in the separated time-space form one can pre-assemble
   * and reduce the space part before the time loop and multiply it in every time step with the
   * appropriate temporal factor.
   */
  if (reassemble)
	assemble_reduce(true);
  this->sys_rhs +=  TDatabase::TimeDB->THETA4 * TDatabase::TimeDB->TIMESTEPLENGTH * this->source;
}

/**************************************************************************** */
void ROM_TCDR2D::solve()
{
  this->red_sol = ROM::solve(this->sys_rhs);
}

/**************************************************************************** */
void ROM_TCDR2D::output(int time_step)
{
  bool no_output = !rom_db["output_write_vtk"] && !rom_db["output_compute_errors"];
  if(no_output)
	return;

  ROM::get_full_solution(this->red_sol, this->full_sol.get_entries());
  TFEFunction2D & _fe_function = this->fe_function;
  _fe_function.PrintMinMax();

  if( time_step % TDatabase::TimeDB->STEPS_PER_IMAGE == 0 )
  {
	Output::print<1>("Writing vtk output at time ",TDatabase::TimeDB->CURRENTTIME);
	// write output
	this->outputWriter.write(TDatabase::TimeDB->CURRENTTIME);
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
    
    double min, max;
    fe_function.MinMax(min,max);
    std::string error_file = rom_db["outfile"].get<std::string>() + ".errors";
    /** @todo Add additional functionality in namespace Output to write information
     * into seond (or generally several) output file (like here Output::redirect(error_file)) */
    
    //TODO: Check the redirect function for true and false case
    if (TDatabase::TimeDB->CURRENTTIME==TDatabase::TimeDB->STARTTIME)
      Output::redirect(error_file);//, true);
    // else
    //   Output::redirect(error_file, false);
    // Output: "time" "L2 error" "H01 error" "min" "max"
    Output::print<1>(TDatabase::TimeDB->CURRENTTIME, " ",
                     loc_e[0], " ", loc_e[1], " ", min, " ", max);
    Output::resetOutfile();
  }
}

/**************************************************************************** */
void ROM_TCDR2D::assemble(BlockFEMatrix &mat, AssembleFctParam2D *local_assemble_param)
{
  int n_matrices = 1;
  std::vector<int> row_space = { 0 };
  std::vector<int> col_space = { 0 };
  int n_rhs = 0;
  std::vector<int> rhs_space = { 0 };
  ManipulateFct2D *manipulate = NULL;
  int n_terms = 5;
  std::vector<MultiIndex2D> derivatives = { D00, D10, D01, D02, D20 };
  bool *needs_second_derivatives;
  needs_second_derivatives = new bool[1];
  needs_second_derivatives[0] = true;
  std::vector<int> fe_space_number = { 0, 0, 0, 0, 0 };
  int n_param = 0;
  int n_param_fct = 0;
  std::vector<ParamFct*> param_fct = {};
  int n_fe_values = 0;
  std::vector<int> fe_value_fct_index = {};
  std::vector<MultiIndex2D> fe_value_multi_index = {};
  std::vector<int> begin_param = {};
  
  TFEFunction2D * pointer_to_function = &this->fe_function;
  
  LocalAssembling2D la_mat(n_terms, derivatives, fe_space_number, row_space, col_space,
                           rhs_space, this->example.get_coeffs(), local_assemble_param,
                           manipulate, n_matrices, n_rhs, n_param_fct, param_fct, begin_param,
                           n_param, &pointer_to_function, n_fe_values, fe_value_fct_index,
                           fe_value_multi_index);

  const TFESpace2D * _fe_space = &this->fe_space;
  BoundCondFunct2D * boundary_conditions = _fe_space->get_boundary_condition();
  BoundValueFunct2D * non_const_bound_value[1] {this->example.get_bd()[0]};
  //fetch matrix as block
  std::vector<std::shared_ptr<FEMatrix>> mat_blocks = mat.get_blocks_uniquely();
  TSquareMatrix2D * mat_block[1]{reinterpret_cast<TSquareMatrix2D*>(mat_blocks.at(0).get())};
  
  //do the assembling
  mat_block[0]->reset();
  Assemble2D(1, &_fe_space, n_matrices, mat_block, 0, NULL, 0, NULL,
             NULL, &boundary_conditions, non_const_bound_value, la_mat);
}


/**************************************************************************** */
void ROM_TCDR2D::assemble(BlockVector &rhs, AssembleFctParam2D *local_assemble_param)
{
  int n_matrices = 0;
  std::vector<int> row_space = { 0 };
  std::vector<int> col_space = { 0 };
  int n_rhs = 1;
  std::vector<int> rhs_space = { 0 };
  ManipulateFct2D *manipulate = NULL;
  int n_terms = 3;
  std::vector<MultiIndex2D> derivatives = { D00, D10, D01 };
  bool *needs_second_derivatives;
  needs_second_derivatives = new bool[1];
  needs_second_derivatives[0] = false;
  std::vector<int> fe_space_number = { 0, 0, 0 };
  int n_param = 0;
  int n_param_fct = 0;
  std::vector<ParamFct*> param_fct = {};
  int n_fe_values = 0;
  std::vector<int> fe_value_fct_index = {};
  std::vector<MultiIndex2D> fe_value_multi_index = {};
  std::vector<int> begin_param = {};
  
  TFEFunction2D * pointer_to_function = &this->fe_function;
  
  LocalAssembling2D la_rhs(n_terms, derivatives, fe_space_number, row_space, col_space,
                           rhs_space, this->example.get_coeffs(), local_assemble_param,
                           manipulate, n_matrices, n_rhs, n_param_fct, param_fct, begin_param,
                           n_param, &pointer_to_function, n_fe_values, fe_value_fct_index,
                           fe_value_multi_index);

  const TFESpace2D * _fe_space = &this->fe_space;
  BoundCondFunct2D * boundary_conditions = _fe_space->get_boundary_condition();
  BoundValueFunct2D * non_const_bound_value[1] {this->example.get_bd()[0]};
  
  double * rhs_entries = rhs.get_entries();
  rhs.reset();
  
  //do the assembling
  Assemble2D(1, &_fe_space, 0, NULL, 0, NULL, n_rhs, &rhs_entries,
             &_fe_space, &boundary_conditions, non_const_bound_value, la_rhs);
}

#include <CD2D.h>
#include <Database.h>
#include <Multigrid.h>
#include <MainUtilities.h> // L2H1Errors
#include <AlgebraicFluxCorrection.h>
#include <PostProcessing2D.h>
#include <LocalAssembling2D.h>
#include <Assemble2D.h>
#include <Upwind.h>
#include <LocalProjection.h>

ParameterDatabase get_default_CD2D_parameters()
{
  Output::print<5>("creating a default CD2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default CD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("CD2D parameter database");
  
  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default afc database
  ParameterDatabase afc_db = AlgebraicFluxCorrection::default_afc_database();
  db.merge(afc_db, true);

  return db;
}
/** ************************************************************************ */
CD2D::System_per_grid::System_per_grid(const Example_CD2D& example,
                                       TCollection& coll, int ansatz_order)
: fe_space(&coll, (char*)"space", (char*)"cd2d fe_space", example.get_bc(0),
           ansatz_order, nullptr),
           // TODO CB: Building the matrix here and rebuilding later is due to the
           // highly non-functional class TFEVectFunction2D (and TFEFunction2D,
           // which do neither provide default constructors nor working copy assignments.)
           matrix({&fe_space}),
           rhs(this->matrix, true),
           solution(this->matrix, false),
           fe_function(&this->fe_space, (char*)"c", (char*)"c",
                       this->solution.get_entries(), this->solution.length())
{
  
  matrix = BlockFEMatrix::CD2D(fe_space);
}

/** ************************************************************************ */
CD2D::CD2D(const TDomain& domain, const ParameterDatabase& param_db,
           int reference_id)
 : CD2D(domain, param_db, Example_CD2D(param_db), reference_id)
{
  eps=1.0,t=0, eps_old=1.0;
  omega_min=0.01, omega=1.0, omega_max=1.0;
  c1=1.001, c2=1.1, c3=1.001, c4=0.9;
  time_newton=0, time_rhs=0;
  rhs_flag=0, newton_flag=0;
  newton_iterate=0, rhs_iterate=0;
  up_param=1e-05;
}

/** ************************************************************************ */
CD2D::CD2D(const TDomain& domain, const ParameterDatabase& param_db,
           const Example_CD2D& example, int reference_id)
 : systems(), example(example), db(get_default_CD2D_parameters()),
   outputWriter(param_db), solver(param_db), errors()
{
  this->db.merge(param_db, false); // update this database with given values
  this->set_parameters();
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  // create finite element space and function, a matrix, rhs, and solution
  int ansatz_order = TDatabase::ParamDB->ANSATZ_ORDER;
  this->systems.emplace_back(this->example, *coll, ansatz_order);

  outputWriter.add_fe_function(&this->get_function());
  

  // print out some information
  TFESpace2D & space = this->systems.front().fe_space;
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  Output::print<1>("N_Cells    : ", setw(12), coll->GetN_Cells());
  Output::print<2>("h (min,max): ", setw(12), h_min, " ", setw(12), h_max);
  Output::print<1>("dof all    : ", setw(12), space.GetN_DegreesOfFreedom());
  Output::print<2>("dof active : ", setw(12), space.GetN_ActiveDegrees());

  
  // done with the constructor in case we're not using multigrid
  if(!this->solver.is_using_multigrid())
    return;
  // else multigrid
  
  auto mg = this->solver.get_multigrid();
  bool mdml = mg->is_using_mdml();
  if(mdml)
  {
    // change the discretization to lowest order
    /// @todo for mdml: is P1/Q1 the correct space on the other grids? Maybe 
    /// what we really need is say Q3/P3, Q2/P2, Q1/P1 on the finest grid and 
    /// Q1/P1 on all coarser grids.
    if(ansatz_order == -1 || ansatz_order == 1)
    {
      // - using non conforming P1 already, it makes no sense to use another 
      //   discretization on the finest grid. 
      // - using conforming P1, we don't do another multigrid level with non 
      //   conforming P1 elements on the finest grid, because this space is 
      //   typically larger that conforming P1.
      // Either way we just do regular multigrid
      mdml = false;
    }
    else
      ansatz_order = -1;
  }
  if(mdml)
    /// @todo mdml for CD2D: We need a special assembling function which 
    /// does not assemble the convection term. Instead one then calls an upwind 
    /// method.
    ErrThrow("mdml is currently not working.");
  
  // number of multigrid levels
  size_t n_levels = mg->get_n_geometric_levels();
  // index of finest grid
  int finest = domain.get_ref_level(); // -> there are finest+1 grids
  // index of the coarsest grid used in this multigrid environment
  int coarsest = finest - n_levels + 1;
  if(mdml)
  {
    coarsest++;
  }
  else
  {
    // only for mdml there is another matrix on the finest grid, otherwise
    // the next system to be created is on the next coarser grid
    finest--;
  }
  if(coarsest < 0 )
  {
    ErrThrow("the domain has not been refined often enough to do multigrid "
             "on ", n_levels, " levels. There are only ",
             domain.get_ref_level() + 1, " grid levels.");
  }
  
  // Construct systems per grid and store them, finest level first
  std::list<BlockFEMatrix*> matrices;
  // matrix on finest grid is already constructed
  matrices.push_back(&systems.back().matrix);
  for (int grid_no = finest; grid_no >= coarsest; --grid_no)
  {
    TCollection *coll = domain.GetCollection(It_EQ, grid_no, reference_id);
    systems.emplace_back(example, *coll, ansatz_order);
    //prepare input argument for multigrid object
    matrices.push_front(&systems.back().matrix);
  }
  mg->initialize(matrices);
  eps=1.0,t=0, eps_old=1.0;
  omega_min=0.01, omega=1.0, omega_max=1.0;
  c1=1.001, c2=1.1, c3=1.001, c4=0.9;
  time_newton=0, time_rhs=0;
  rhs_flag=0, newton_flag=0;
  newton_iterate=0, rhs_iterate=0;
  up_param=1e-05;
}

/** ************************************************************************ */
CD2D::~CD2D()
{
  // delete the collections created during the contructor
  for(auto & s : this->systems)
    delete s.fe_space.GetCollection();
}

/** ************************************************************************ */
void CD2D::set_parameters()
{
  //set problem_type to CD if not yet set
  if(!db["problem_type"].is(1))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 1;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to CD."
          "It is now reset to the correct value for CD (=1).");
      db["problem_type"] = 1;
    }
  }
  //////////////// Algebraic flux correction ////////////
  if(!db["algebraic_flux_correction"].is("none"))
  {//some kind of afc enabled
    if(!db["algebraic_flux_correction"].is("afc"))
    {
      db["algebraic_flux_correction"].set("afc");
      Output::print("Only kind of algebraic flux correction"
          " for CD problems is AFC (afc).");
    }
    //make sure that galerkin discretization is used
    if (!db["space_discretization_type"].is("galerkin"))
    {//some other disctype than galerkin
      db["space_discretization_type"] = "galerkin";
      Output::warn<1>("Parameter 'space_discretization_type' changed to 'galerkin' "
          "because Algebraic Flux Correction is enabled.");
    }
    // when using afc, create system matrices as if all dofs were active
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
  }
}

/** ************************************************************************ */
void CD2D::assemble()
{
  LocalAssembling2D_type t = LocalAssembling2D_type::ConvDiff;
  bool mdml = this->solver.is_using_multigrid()
             && this->solver.get_multigrid()->is_using_mdml();
  // in case of mdml, we need to change the local assembling, (not yet 
  // implemented)

  // this loop has more than one iteration only in case of multigrid
  for(auto & s : this->systems)
  {
    TFEFunction2D * pointer_to_function = &s.fe_function;
    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling2D la(t, &pointer_to_function, example.get_coeffs());

    // assemble the system matrix with given local assembling, solution and rhs
    const TFESpace2D * fe_space = &s.fe_space;
    BoundCondFunct2D * boundary_conditions = fe_space->GetBoundCondition();
    int N_Matrices = 1;
    double * rhs_entries = s.rhs.get_entries();

    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();
    TSquareMatrix2D * matrix = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

    BoundValueFunct2D * non_const_bound_value[1] {example.get_bd()[0]};

      // reset right hand side and matrix to zero (just in case)
      s.rhs.reset();
      matrix->reset();

      // assemble
      Assemble2D(1, &fe_space, N_Matrices, &matrix, 0, NULL, 1, &rhs_entries,
                 &fe_space, &boundary_conditions, non_const_bound_value, la);

      // apply local projection stabilization method
      if(db["space_discretization_type"].is("local_projection")
         && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
      {
        if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
        {
          UltraLocalProjection((void *)&matrix, false);
        }
        else
        {
          ErrThrow("LP_FULL_GRADIENT needs to be one to use LOCAL_PROJECTION");
        }
      }
      
      bool finest_grid = &systems.front() == &s;
      if(mdml && !finest_grid)
      {
        UpwindForConvDiff(la.GetCoeffFct(), matrix, rhs_entries, fe_space, 
                          nullptr, nullptr, false);
      }

      // copy Dirichlet values from rhs to solution vector (this is not really
      // necessary in case of a direct solver)
      s.solution.copy_nonactive(s.rhs);
  }

  // when using afc, do it now
  if(!db["algebraic_flux_correction"].is("none"))
  {
    do_algebraic_flux_correction();
  }

}

bool CD2D::solve(const int iteration)
{

  t = GetTime();
  System_per_grid& s = this->systems.front();
  // compute norm of residual
  BlockVector res = s.rhs;
  //The Dynamic damping factor is used for Zalesak Limiter because of slow convergence fin LINEAR_PRESERVE_BJK17
  //For iteration=0 we would have got an error which is highlighted below
  if(!db["algebraic_flux_correction"].is("none") && iteration>0 && !db["afc_limiter"].is("linearity_preserving_BJK17"))
  {
    Output::print<2>("nonlinear step ", iteration, " residual before solve: ", eps_old);
    if(eps_old>(double)db["afc_nonlinloop_epsilon"])
    {
      first_damp=1;
      abc:                                                                    // label to use for goto
      //The error comes in this line as for iteration=0, we don't have any value for old_solution
      AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, s.solution, omega);
      s.matrix.apply_scaled_add(s.solution , res ,-1.0);                      //Calculation of residue r_k+1
      eps = res.norm();
      if(eps<eps_old || omega<=c1*omega_min)                                  //eps_old is the residue r_k
      {
	if(eps<eps_old &&first_damp==1)
        {
	  omega_max=std::min(1.0,c3*omega_max);
       	  omega=std::min(omega_max,c2*omega);
        }
      }
      else
      {
	omega=std::max(omega_min,omega/2);
        if(first_damp==1)
        {
	  omega_max=std::max(omega_min,c4*omega_max);
	  first_damp=0;	  
	}
	//Output::print<4>("Epsilon Old :", eps_old, " and Epsilon new is ", eps," at iteration: ",iteration, " and omega :", omega);
	goto abc;
      }
      //Output::print<4>("Value of Omega ", iteration, " : ", omega);
    }
    //When the solution converges
    else
    {
      Output::print<1>("NONLINEAR iteration converged: iterations ", iteration, " norm of residual ", eps);
      return(TRUE);
    } 
  }  
  //For initial iteration and when the limiter is LINEAR_PRESERVE_BJK17
  if(iteration==0 || db["afc_limiter"].is("linearity_preserving_BJK17"))
  {
    s.matrix.apply_scaled_add(s.solution , res ,-1.0);
    eps = res.norm();
    Output::print<2>("nonlinear step ", iteration, " residual before solve: ", eps);
  }
  //Stopping criteria for LINEAR_PRESERVE_BJK17
  if(!db["algebraic_flux_correction"].is("none") && db["afc_limiter"].is("linearity_preserving_BJK17"))
  {
    // stop iteration if norm of residual is sufficiently small 
    if (eps <= (double)db["afc_nonlinloop_epsilon"])
    {
      Output::print<1>("NONLINEAR iteration converged: iterations ", iteration, " norm of residual ", eps);
      return(TRUE);
    }
  }
  //Storage of old solution
  if(!db["algebraic_flux_correction"].is("none"))
    old_solution = s.solution;
  
  this->solver.solve(s.matrix, s.rhs, s.solution);
  Output::print<2>("SOLVE DONE");
  t = GetTime() - t;
  Output::print<3>(" solving of a CD2D problem done in ", t, " seconds");
  //Damping of the solution using the damping factor in database
  if(!db["algebraic_flux_correction"].is("none") && db["afc_limiter"].is("linearity_preserving_BJK17"))
  {
    AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, s.solution, db["afc_nonlinloop_damping_factor"]);
  }
  Output::print<4>("solve done ", iteration);
  //Storage of old residue
  eps_old=eps;
  return(FALSE);   
}

/** ************************************************************************ */
void CD2D::output(int i)
{
  // print the value of the largest and smallest entry in the finite element 
  // vector
  TFEFunction2D & fe_function = this->systems.front().fe_function;
  fe_function.PrintMinMax();
  
  // write solution to a vtk file or in case-format
  outputWriter.write(i);

  /*
  // implementation with the old class TOutput2D
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(1, 1, 0, 0, NULL);
    Output.AddFEFunction(&fe_function);

    // Create output directory, if not already existing.
    mkdir(db["output_vtk_directory"], 0777);
    std::string filename = this->db["output_vtk_directory"];
    filename += "/" + this->db["output_basename"].value_as_string();

    if(i >= 0)
      filename += "_" + std::to_string(i);
    filename += ".vtk";
    Output.WriteVtk(filename.c_str());
  }
  */
  

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(this->db["output_compute_errors"])
  {
    // this should be a little longer than this->errors, because of a bug in
    // FEFunction::GetErrors. Otherwise we could use this->errors directly.
    // Note that we can not write 
    // 'constexpr size_t n_errors = errors.max_size();'. The reason is that the
    // method 'max_size' is not marked const in c++11, but it is in c++14. We
    // should switch to that.
    std::array<double, 5> errors;
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D* space = fe_function.GetFESpace2D();
    
    fe_function.GetErrors(this->example.get_exact(0), 3, AllDerivatives, 4,
                          SDFEMErrors, this->example.get_coeffs(), &aux, 1, 
                          &space, errors.data());
    
    Output::print<1>("L2     : ", errors[0]);
    Output::print<1>("H1-semi: ", errors[1]);
    Output::print<1>("SD     : ", errors[2]);
    Output::print<1>("L_inf  : ", errors[3]);
    // copy local variable to member variable
    std::copy(errors.begin(), errors.end()-1, this->errors.begin());
  } 
}

AlgebraicFluxCorrection::Limiter string_to_limiter(std::string afc_limiter);
AlgebraicFluxCorrection::Iteration_Scheme string_to_it_scheme(std::string afc_iteration_scheme);

/** ************************************************************************ */
void CD2D::do_algebraic_flux_correction()
{
  Output::print<4>("AFC: enter do_algebraic_flux_correction");
  for(auto & s : this->systems) // do it on all levels.
  {
    bool compute_D_and_gamma = FALSE;
    
    //determine which kind of afc to use
    if(db["algebraic_flux_correction"].is("default") ||
        db["algebraic_flux_correction"].is("afc"))
    {
      
      //
      AlgebraicFluxCorrection::Limiter limiter = string_to_limiter(db["afc_limiter"]);
      AlgebraicFluxCorrection::Iteration_Scheme it_scheme = string_to_it_scheme(db["afc_iteration_scheme"]);
      
      //USE OF NEWTON METHOD WHEN RESIDUE BECOMES LESS THAN E-05
      /*There is change of scheme when time taken by the scheme exceeds 2 secomds
       * When the time step is approximated to 2 seconds in Newton we bring back the time for FIXED_POINT_RHS to 0 seconds
       * For FIXED_POINT_RHS we used a different method because as soon as we have time_rhs>1 and if we used ceiling function
       * we would have moved back to newton method and hence we use floor function. As soon as we get time_rhs>2 
       * we reset time_newton to 0 and hence we move back to NEWTON method.
       */
      
      //CONDITION WHEN ONE ITERATE OF NEWTON AS WELL AS FIXED POINT RHS TAKES MORE THAN 2 SECONDS
      /*There can be cases when one iterate of Newton as well as Fixed point RHS takes more than 2 seconds. For that
       *case we are using newton_iterate and rhs_iterate which will count the number of iterations that happened while
       *running that scheme. If newton as well as rhs takes 2 seconds for one iteration then we make rhs_flag as 1, which will help
       * us to identify that both the iterations are slow and hence we move to Newton Method as it takes less number of steps. 
       */
      
      //CHANGE IN TOLERANCE WHEN NEWTON MOVES BACK TO FIXED_POINT_RHS
      /*When the iteration changes from newton to fixed_point_rhs we increase the tolerance limit by a factor of 1e-02 
       */
      if(eps_old<up_param && time_newton<=2.0 && rhs_flag==0)
      {
	Output::print<1>("NEWTON IMPLEMENTED");
	it_scheme=AlgebraicFluxCorrection::Iteration_Scheme::NEWTON;
	time_newton+=t;
	newton_iterate++;
	if(std::ceil(time_newton)==2.0)
	{
	  time_rhs=0.0;
	  rhs_iterate=0; 
	}
	newton_flag=1;
      }
      else if(eps_old<up_param && time_rhs<=2.5 && rhs_flag==0)
      {
	it_scheme=AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS;
	time_rhs+=t;
	rhs_iterate++;
	if(std::floor(time_rhs)==2.0 && rhs_iterate==1 && newton_iterate==1)
	  rhs_flag=1;
	if(std::floor(time_rhs)==2.0)
	{
	  newton_iterate=0;
	  time_newton=0.0;
	}
      }
      if(eps_old<up_param && rhs_flag==1)
	it_scheme=AlgebraicFluxCorrection::Iteration_Scheme::NEWTON;
      if(eps_old>=up_param && newton_flag==1)
      {
	up_param*=1e-02;
	newton_flag=0;
      }
      Output::print<4>("Iterate RHS: ", rhs_iterate);
      Output::print<4>("Time RHS: ", time_rhs);
      Output::print<4>("Iterate Newton: ", newton_iterate);
      Output::print<4>("Time NEWTON: ", time_newton);
      
      //get pointers/references to the relevant objects
      TFESpace2D& feSpace = s.fe_space;
      FEMatrix& one_block = *s.matrix.get_blocks_uniquely().at(0).get();
      const std::vector<double>& solEntries = s.solution.get_entries_vector();
      std::vector<double>& rhsEntries = s.rhs.get_entries_vector();

      // fill a vector "neumannToDirichlet" with those rows that got
      // internally treated as Neumann although they are Dirichlet
      int firstDiriDof = feSpace.GetActiveBound();
      int nDiri = feSpace.GetN_Dirichlet();

      std::vector<int> neumToDiri(nDiri, 0);
      std::iota(std::begin(neumToDiri), std::end(neumToDiri), firstDiriDof);
      
      // if necessary, set up vector gamma and matrix D
      if(s.afc_matrix_D_entries.empty())
      {
	Output::print<4>("AFC: allocate matrix D");
	s.afc_matrix_D_entries.resize(one_block.GetN_Entries(),0.0);
	compute_D_and_gamma = TRUE;
      }
      if ((limiter== AlgebraicFluxCorrection::Limiter::LINEAR_PRESERVE_BJK17) && (s.afc_gamma.empty()))
      {
	Output::print<4>("AFC: vector gamma");
	s.afc_gamma.resize(feSpace.GetN_DegreesOfFreedom(),0.0);
      }
      
      // apply FEM-TVD
      AlgebraicFluxCorrection::steady_state_algorithm(
          one_block,
          solEntries,rhsEntries,
          neumToDiri, 
	  s.afc_matrix_D_entries, s.afc_gamma, compute_D_and_gamma,
	  limiter, it_scheme);

      //...and finally correct the entries in the Dirchlet rows
      AlgebraicFluxCorrection::correct_dirichlet_rows(one_block);
      //...and in the right hand side, too, assum correct in solution vector
      s.rhs.copy_nonactive(s.solution);
    }
    else
    {
      ErrThrow("The chosen algebraic flux correction scheme is unknown "
          "to class CD2D.");
    }
  }
}

/** ************************************************************************ */
AlgebraicFluxCorrection::Limiter string_to_limiter(std::string afc_limiter)
{
  if (afc_limiter == std::string("zalesak"))
    return AlgebraicFluxCorrection::Limiter::ZALESAK;
  else if (afc_limiter == std::string("linearity_preserving_BJK17"))
    return AlgebraicFluxCorrection::Limiter::LINEAR_PRESERVE_BJK17;
  else 
  {
    ErrThrow("afc_limiter ", afc_limiter, " not implemented!!!"); 
  }
}

/** ************************************************************************ */
AlgebraicFluxCorrection::Iteration_Scheme string_to_it_scheme(std::string afc_iteration_scheme)
{
  if (afc_iteration_scheme == std::string("fixed_point_rhs"))
    return AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS;
  else if (afc_iteration_scheme == std::string("fixed_point_matrix"))
    return AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_MATRIX;
  else if (afc_iteration_scheme == std::string("newton"))
    return AlgebraicFluxCorrection::Iteration_Scheme::NEWTON;
  else 
  {
    ErrThrow("afc_iteration_scheme ", afc_iteration_scheme, " not implemented!!!"); 
  }
}

/** ************************************************************************ */
double CD2D::get_L2_error() const
{
  return this->errors[0];
}

/** ************************************************************************ */
double CD2D::get_H1_semi_error() const
{
  return this->errors[1];
}

/** ************************************************************************ */
double CD2D::get_SD_error() const
{
  return this->errors[2];
}

/** ************************************************************************ */
double CD2D::get_L_inf_error() const
{
  return this->errors[3];
}


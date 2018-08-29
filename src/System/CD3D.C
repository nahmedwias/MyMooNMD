#include <CD3D.h>
#include <Example_CD3D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <AlgebraicFluxCorrection.h>
#include <LinAlg.h>
#include <LocalAssembling3D.h>
#include <Assemble3D.h>
#include <LocalProjection.h>
#include <anderson.h>

#include <PostProcessing3D.h>
#include <Upwind3D.h>
#include <Multigrid.h>
#include <MainUtilities.h>

#include <DirectSolver.h>

#include <MainUtilities.h> // L2H1Errors

#include <Multigrid.h>

#include <sys/stat.h>

#ifdef _MPI
#include <MumpsWrapper.h>
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

ParameterDatabase get_default_CD3D_parameters()
{
  Output::print<5>("creating a default CD3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default CD3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("CD3D parameter database");
  
  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default afc database
  ParameterDatabase afc_db = AlgebraicFluxCorrection::default_afc_database();
  db.merge(afc_db, true);
  
  return db;
}

#ifdef _MPI
  CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example,
                                     TCollection& coll, int maxSubDomainPerDof)
   : feSpace_(new TFESpace3D(&coll, "space", "cd3d fe_space", example.get_bc(0),
              TDatabase::ParamDB->ANSATZ_ORDER))
  {
    //inform the fe space about the maximum number of subdomains per dof
    feSpace_->initialize_parallel(maxSubDomainPerDof);
    feSpace_->get_communicator().print_info();

    // set the matrix with named constructor
    matrix_ = BlockFEMatrix::CD3D(*feSpace_);

    rhs_ = BlockVector(matrix_, true);
    solution_ = BlockVector(matrix_, false);

    feFunction_ = TFEFunction3D(feSpace_.get(), "c", "c",
                                solution_.get_entries(), solution_.length());

  }
#else
  /* ************************************************************************ */
  CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example,
                                     TCollection& coll)
   : feSpace_(new TFESpace3D(&coll, "space", "cd3d fe_space", example.get_bc(0),
              TDatabase::ParamDB->ANSATZ_ORDER))
  {
    // set the matrix with named constructor
    matrix_ = BlockFEMatrix::CD3D(*feSpace_);

    rhs_ = BlockVector(matrix_, true);
    solution_ = BlockVector(matrix_, false);

    feFunction_ = TFEFunction3D(feSpace_.get(), "c", "c", 
                                solution_.get_entries(), solution_.length());
  }
#endif

  /** ************************************************************************ */
  CD3D::CD3D(std::list<TCollection*> collections,
             const ParameterDatabase& param_db, const Example_CD3D& example
#ifdef _MPI
             ,int maxSubDomainPerDof
#endif
  )
  : systems_(), example_(example), db(get_default_CD3D_parameters()),
    outputWriter(param_db), solver(param_db), errors_(),  alphas_x_i(), old_solution()
  {
    this->db.merge(param_db, false); // update this database with given values
    this->checkParameters();
    // The construction of the members differ, depending on whether
    // a multigrid solver will be used or not.
    bool usingMultigrid = this->solver.is_using_multigrid();

    if (!usingMultigrid)
    {
      //Check at least if the collections list contains exactly one Collection.
      if(collections.size() != 1 )
      {
        ErrThrow("Non-multigrid: Expected exactly one collection!");
      }

      // Get the one given collection.
      TCollection& cellCollection = *collections.front();

#ifdef _MPI
      // create finite element space and function, a matrix, rhs, and solution
      systems_.emplace_back(example_, cellCollection, maxSubDomainPerDof);
#else
      // create finite element space and function, a matrix, rhs, and solution
      systems_.emplace_back(example_, cellCollection);      
#endif
    }
    else
    {// we are using multigrid
        size_t n_levels = collections.size();
        if(!param_db["multigrid_n_levels"].is(n_levels))
           ErrThrow("Number of collection does not equal number of multigrid levels!");

        auto mg = this->solver.get_multigrid();
        // Construct systems per grid and store them, finest level first
        std::list<BlockFEMatrix*> matrices;
        for (auto coll : collections)
        {
#ifdef _MPI
          systems_.emplace_back(example, *coll, maxSubDomainPerDof);
#else
          systems_.emplace_back(example, *coll);
#endif
          //prepare input argument for multigrid object
          matrices.push_front(&systems_.back().matrix_);
        }
        mg->initialize(matrices);
    }
    alphas_x_i=this->systems_.front().solution_;
    old_solution=this->systems_.front().solution_;
      // print useful information
      this->output_problem_size_info();
      time_newton=0, time_rhs=0;
      rhs_flag=0, newton_flag=0;
      newton_iterate=0, rhs_iterate=0;
      up_param=1e-25;
      rejected_steps = 0;
      is_not_afc_fixed_point_rhs=1;
  }

/** ************************************************************************ **/
//==============================================================================
void CD3D::output_problem_size_info() const
{
  // print some useful information
  auto& space = *this->systems_.front().feSpace_;
  double hMin, hMax;
  TCollection *coll = space.GetCollection();
  coll->GetHminHmax(&hMin, &hMax);
  Output::print<1>("N_Cells    : ", setw(13), coll->GetN_Cells());
  Output::print<1>("h(min, max): ", setw(13), hMin, " ", setw(13), hMax);
  Output::print<1>("dofs all   : ", setw(13), space.GetN_DegreesOfFreedom());
  Output::print<1>("dof active : ", setw(13), space.GetActiveBound());
}
/** ************************************************************************ **/
void CD3D::assemble(const int iteration)
{
  //determine the local assembling type to be CD3D
  LocalAssembling3D_type t = LocalAssembling3D_type::CD3D;
  bool mdml = this->solver.is_using_multigrid()
             && this->solver.get_multigrid()->is_using_mdml();
  // this loop has more than one iteration only in case of multigrid
  if(iteration==0)
    time_total=GetTime();
  
  for(auto & s : systems_)
  {
    TFEFunction3D * pointer_to_function = &s.feFunction_;
    int afc_ini_supg = 0;
 
    std::shared_ptr<LocalAssembling3D> la;
    
    if(!db["algebraic_flux_correction"].is("none") && iteration==0 
     && db["afc_initial_iterate"].is("supg"))
    {
      global_space_type = 2;
      afc_ini_supg = 1; 
      // create a local assembling object which is needed to assemble the matrix
      la = std::make_shared<LocalAssembling3D>(t, &pointer_to_function, example_.get_coeffs(), global_space_type);
    }
    else
    {
      // create a local assembling object which is needed to assemble the matrix
      //LocalAssembling3D la(t, &pointer_to_function, example.get_coeffs(),global_space_type);
      //global_space_type = (int) db["space_discretization_type"];
      Output::print<2>("assembling discretization ",global_space_type);
      la = std::make_shared<LocalAssembling3D>(t, &pointer_to_function, example_.get_coeffs(),global_space_type);
    }
    /*if(!db["algebraic_flux_correction"].is("none") && iteration==0 )
    {
      if (db["afc_limiter"].is("kuzmin"))
      {
	db["afc_nonlinloop_kuzmin_first"].set<>("no");
      }
      else
      {
        if (db["afc_nonlinloop_kuzmin_first"].is("yes"))
	{
          db["afc_limiter"].set<>("kuzmin");
	  Output::print<2>("limiter set to kuzmin");
	}
      }
    }*/
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_uniquely();
    TSquareMatrix3D * matrix = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
    
    //BoundValueFunct3D * non_const_bound_value[1] {example_.get_bd()[0]};
    
    if(is_not_afc_fixed_point_rhs==1)  
    {
      call_assembling_routine(s, *la);
      if(iteration==1)
        rhs_copy=s.rhs_;  
    }
    //for FIXED_POINT_RHS and iteration>1
    else
    {
      s.rhs_=rhs_copy;
    }
    // apply local projection stabilization method
    if(db["space_discretization_type"].is("local_projection")
      && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
    {
      if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
      {
        UltraLocalProjection3D((void *)&matrix, false);
      }
      else
      {
        ErrThrow("LP_FULL_GRADIENT needs to be one to use LOCAL_PROJECTION");
      }  
    }
    bool finest_grid = &systems_.front() == &s;
    if(is_not_afc_fixed_point_rhs==1)
    {
      if ((mdml && !finest_grid) || (db["space_discretization_type"].is("upwind")))
      {
        Output::print<2>("upwind for convection-diffusion equation");
        UpwindForConvDiff(matrix, s.rhs_.get_entries(), s.feSpace_.get(), 
                          *la);  
      }
      if (!db["algebraic_flux_correction"].is("none")&&(iteration==0)&&
          db["afc_initial_iterate"].is("upwind"))
      {
        Output::print<2>("upwind for convection-diffusion equation");
        UpwindForConvDiff(matrix, s.rhs_.get_entries(), s.feSpace_.get(), 
                          *la);  
      }
    }
    if (afc_ini_supg == 1)
      db["space_discretization_type"].set<>("galerkin");
    // copy Dirichlet values from rhs to solution vector (this is not really
    // necessary in case of a direct solver)
    s.solution_.copy_nonactive(s.rhs_);
    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling3D laObject(t, &pointer_to_function, example_.get_coeffs());

    // assemble the system matrix with given local assembling, solution and rhs
    //s.matrix_.assemble(laObject, s.solution_, s.rhs_);
    //call_assembling_routine(s, laObject);
  }

  if (iteration == 0)
  {
   db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    // the following flags are for the BAIL proceedings
    // switch to fixed point rhs
    if (((int)db["afc_nonlinloop_switch_to_newton_scheme"]==10)||
        ((int)db["afc_nonlinloop_switch_to_newton_scheme"]==20))
    {
      db["afc_fixed_point_matrix_weight"] = 0.0;
      db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }
    // switch to fixed point matrix 
    if (((int)db["afc_nonlinloop_switch_to_newton_scheme"]==11)||
       ((int)db["afc_nonlinloop_switch_to_newton_scheme"]==21)) 
    {
      db["afc_fixed_point_matrix_weight"] = 1.0;
      db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }
  }

  if(!db["algebraic_flux_correction"].is("none"))
    do_algebraic_flux_correction(iteration, is_not_afc_fixed_point_rhs);
  
  //Previous Implementation
  /*if(!db["algebraic_flux_correction"].is("none"))
    do_algebraic_flux_correction(iteration);*/
  
}
/** ***************************************************************************************************** **/
AlgebraicFluxCorrection::Iteration_Scheme string_to_it_scheme(std::string afc_iteration_scheme);
/** ***************************************************************************************************** **/

bool CD3D::solve(const int iteration)
{
  double t = GetTime();
  SystemPerGrid& s = this->systems_.front();
  BlockVector res = s.rhs_;
  BlockVector current_rhs = s.rhs_;
  BlockVector current_sol = s.solution_;
  AlgebraicFluxCorrection::Iteration_Scheme it_scheme = string_to_it_scheme(db["afc_iteration_scheme"]);
  is_not_afc_fixed_point_rhs=0;
  
  //To assemble for 1st iteration withj FIXED_POINT_RHS
  if(it_scheme!=AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS
    || iteration == 1)
    is_not_afc_fixed_point_rhs=1;

  // special treatment of linear discretization or first iteration
  if (iteration==0)
  {
    // compute residual vector
     s.matrix_.apply_scaled_add(s.solution_ , res ,-1.0); 
     
     //s.matrix.get_combined_matrix()->Print("mat");
     //s.rhs.print("rhs");
     //s.solution.print("sol");
     //exit(1);
     
     // compute the norm of the residual vector, residual_old is the norm of residual r_k
     residual_old = res.norm();
     old_solution = s.solution_;
     this->solver.solve(s.matrix_, s.rhs_, s.solution_);
     
     t = GetTime() - t;
     Output::print("  solving of a CD3D problem done in ", t, " seconds");
     
     if(db["algebraic_flux_correction"].is("none"))
     {
       // THIS STATEMENT ASSUMES THAT THE SOLVING PROCESS WAS SUFFICIENTLY ACCURATE
       return(true);
     }
     else
     {
       double h_min, h_max;
       const TFESpace3D & space = *this->systems_.front().feSpace_;
       TCollection *coll = space.GetCollection();
       coll->GetHminHmax(&h_min, &h_max);
       double sigma = (double)db["afc_newton_regu_sigma"];
       sigma *= h_max * h_max * h_max * h_max;
       db["afc_newton_regu_sigma"] = sigma;
       Output::print<2>("afc_newton_regu_sigma changed to ", db["afc_newton_regu_sigma"]);
       
       double dof = space.GetN_DegreesOfFreedom();
       double eps = (double)db["afc_nonlinloop_epsilon"];
       eps *= sqrt(dof);
       db["afc_nonlinloop_epsilon"] = eps;
       Output::print<2>("afc_nonlinloop_epsilon normalized to ", db["afc_nonlinloop_epsilon"]);
 
       Output::print<2>("nonlinear step ", iteration, " residual: ", residual_old);
       
       return(false);
     }
  }
  
  if(!db["algebraic_flux_correction"].is("none"))
  {
    if ((db["afc_nonlinloop_anderson_acc"].is("yes"))&&(iteration >= 
      (int)db["afc_nonlinloop_anderson_acc_start"]))
    {
      int N_Unknowns=res.length();
      anderson_acceleration_damping(N_Unknowns,iteration, solAnderson, deltaAnderson);
    }
    
    //Constant Damping
    if (db["afc_nonlinloop_damping_factor_constant"].is("yes"))
    {
      // compute proposal for the next solution 
      AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, s.solution_, db);
      assemble(iteration);
      // calculation of residual r_k+1
      res = s.rhs_;
      s.matrix_.apply_scaled_add(s.solution_ , res ,-1.0); 
      // compute the norm of the residual vector, residual_old is the norm of residual r_k
      residual = res.norm();
      Output::print<4>("  residual for proposed new iterate ", residual);
    }
    //Dynamic Damping from [JK08]
    else 
      dynamic_damping(iteration);
  }

  old_solution.add_scaled(s.solution_,-1.0); 
  Output::print<2>("nonlinear step ", iteration, " residual: ", residual, " reduction ", residual/residual_old ," change of sol ", old_solution.norm());
  // stopping criterion satisfied
  if ((residual < (double)db["afc_nonlinloop_epsilon"])&&(iteration >1))
  {
    if (db["afc_limiter"].is("kuzmin") && db["afc_nonlinloop_kuzmin_first"].is("yes"))
    {
      db["afc_limiter"].set<>("BJK17");
      Output::print<2>("ite ", iteration, " limiter set to BJK17");
      db["afc_nonlinloop_kuzmin_first"].set<>("no");
      return(false);
    }
    time_total = GetTime()-time_total;
    Output::print<1>("NONLINEAR ITERATION: ite ", iteration, " res ", residual, " rejections ", rejected_steps, " time ", time_total, " t/it ", time_total/(iteration+rejected_steps));
    return(true);
  }
  //maximal number of iterations
  if (iteration == (int)db["afc_nonlinloop_maxit"])
  {
    time_total = GetTime()-time_total;
    Output::print<1>("MAX_NONLINEAR ITERATION: ite ", iteration, " res ", residual, " rejections ", rejected_steps, " time ", time_total, " t/it ", time_total/(iteration+rejected_steps));
    return(true);    
  }
  // storage of old solution
  old_solution = s.solution_;
  
  // solve the linear system
  //matrix is factorised only once and then stored in the system
  if(is_not_afc_fixed_point_rhs==1)
    this->solver.update_matrix(s.matrix_);
  this->solver.solve(s.rhs_, s.solution_);
  //Previous Implementation
  //this->solver.solve(s.matrix_, s.rhs_, s.solution_);
  t = GetTime() - t;
  Output::print<4>("  iteration done in ", t, " seconds");
  
  
  if ((db["afc_iteration_scheme"].is("newton"))&&((int)db["afc_nonlinloop_switch_to_newton_scheme"]==1))
  {
    double threshold = (double)db["afc_change_method_threshold"];
    if (residual <= threshold)
    {
      double omega_derivative = (double)db["afc_fixed_point_derivative_weight_factor"];
      // first application
      if (omega_derivative < 1e-3)
        omega_derivative = 0.25;
      if (residual/residual_old > 0.99)
      {
        omega_derivative *=0.999;
        if (omega_derivative < 0.1)
          omega_derivative = 0.1;
      }
      if (residual/residual_old < 0.99)
      {
        omega_derivative *=1.001;
        if (omega_derivative >1.0)
          omega_derivative = 1.0;
      }
      db["afc_fixed_point_derivative_weight_factor"] = omega_derivative;
    }
    // go back to former method
    if (residual > 100*threshold)
    {
      db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }
    Output::print<2>("afc_fixed_point_derivative_weight_factor ",  (double)db["afc_fixed_point_derivative_weight_factor"]);
  }

  /*********************************************************************************************/
  // this scheme is for the BAIL proceedings 
  /*********************************************************************************************/
  if ((db["afc_iteration_scheme"].is("newton"))&&((int)db["afc_nonlinloop_switch_to_newton_scheme"]>=10))
  {
    double threshold = (double)db["afc_change_method_threshold"];
    // switch to formal Newton 
    if (residual <= threshold)
    {
        db["afc_fixed_point_matrix_weight"] = 1.0;
        db["afc_fixed_point_derivative_weight_factor"] = 1.0;
    }
    /*else
    {
        if (residual > 10*threshold)
    {
        db["afc_fixed_point_matrix_weight"] = 0.0;
        db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }
    }*/
    if (residual > 100*threshold)
    {
      if ((int)db["afc_nonlinloop_switch_to_newton_scheme"]==20)
      {
        db["afc_fixed_point_matrix_weight"] = 0.0;
        db["afc_fixed_point_derivative_weight_factor"] = 0.0;
      }
      if ((int)db["afc_nonlinloop_switch_to_newton_scheme"]==21)
      {
        db["afc_fixed_point_matrix_weight"] = 1.0;
        db["afc_fixed_point_derivative_weight_factor"] = 0.0;
      }
    }
    Output::print<2>("afc_fixed_point_matrix_weight ",  (double)db["afc_fixed_point_matrix_weight"]);
  }
  residual_old=residual;   
  return(false);   
}
/** ************************************************************************ */
void CD3D::dynamic_damping(const int iteration)
{
  SystemPerGrid& s = this->systems_.front();
  double omega = db["afc_nonlinloop_damping_factor"];
  int first_damp=1;
  BlockVector res = s.rhs_;
  BlockVector current_rhs = s.rhs_;
  BlockVector current_sol = s.solution_;
  int check_anderson=0;
  if(iteration<=(int)db["afc_nonlinloop_anderson_acc_vec"]+(int)db["afc_nonlinloop_anderson_acc_start"])
    check_anderson=1;
  while(1)
  {
    // damping from WN11
    if(db["afc_anderson_damping"].is("yes"))
    {
      if(db["afc_nonlinloop_anderson_acc"].is("yes"))
      {
	if(check_anderson==1)
	  AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, s.solution_, db);
	//starts anderson damping after the vectors are stored
	else
	{
	  s.solution_.scale(omega);
	  s.solution_.add_scaled(alphas_x_i, 1.-omega);
	}
      }
      else
	AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, s.solution_, db);
      assemble(iteration);
      // calculation of residual r_k+1
      res = s.rhs_;
      s.matrix_.apply_scaled_add(s.solution_ , res ,-1.0); 
      // compute the norm of the residual vector, residual_old is the norm of residual r_k
      residual = res.norm();
      Output::print<4>("  residual for proposed new iterate ", residual);
      // accept the first damping parameter
      if (iteration==1)
	break;
      if(db["afc_nonlinloop_anderson_acc"].is("yes") && check_anderson==1)
	break;
    }
    // dynamic damping from JK08
    else
    {
       AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, s.solution_, db);
       assemble(iteration);
       // calculation of residual r_k+1
       res = s.rhs_;
       s.matrix_.apply_scaled_add(s.solution_ , res ,-1.0); 
       // compute the norm of the residual vector, residual_old is the norm of residual r_k
       residual = res.norm();
       Output::print<4>("  residual for proposed new iterate ", residual);
       // accept the first damping parameter
       if (iteration==1)
	 break;
    }
    // if the norm of the residual vector decreases or if the damping parameter is already very small
    // then accept the new iterate 
    if (residual<residual_old || omega<=(double)db["afc_nonlinloop_damping_factor_min_tol"]*(double)db["afc_nonlinloop_damping_factor_min"])
    {
      // if the norm of the residual decreases without having decreased the damping 
      // parameter in this step, then increase the damping parameter if possible
      if(residual<residual_old &&first_damp==1)
      {
	db["afc_nonlinloop_damping_factor_max"]=std::min((double)db["afc_nonlinloop_damping_factor_max_global"],
							 (double)db["afc_nonlinloop_damping_factor_max_increase"]
							 *(double)db["afc_nonlinloop_damping_factor_max"]);
	omega=std::min((double)db["afc_nonlinloop_damping_factor_max"],(double)db["afc_nonlinloop_damping_factor_increase"]*omega);   
      }
      Output::print<2>("  iterate accepted, damping factor ", omega, " ", db["afc_nonlinloop_damping_factor_max"]);
      db["afc_nonlinloop_damping_factor"] = omega;
      break;  
    }
    else
    {
      if (db["afc_nonlinloop_anderson_acc"].is("no"))
      {
        rejected_steps++;
        // get starting situation back 
        s.solution_ = current_sol;
        s.rhs_ = current_rhs;
      }
      // reduce damping factor  
      omega=std::max((double)db["afc_nonlinloop_damping_factor_min"],
		     omega*(double)db["afc_nonlinloop_damping_factor_decrease"]);
      // reduce maximal value for damping factor
      if(first_damp==1)
      {
	db["afc_nonlinloop_damping_factor_max"]=std::max((double)db["afc_nonlinloop_damping_factor_min"],
							 (double)db["afc_nonlinloop_damping_factor_max_decrease"]*(double)db["afc_nonlinloop_damping_factor_max"]);
	first_damp=0;	  
      }
      Output::print<2>("  iterate rejected, res old ", residual_old, " res new: ", residual, " damping factor :", omega);
      db["afc_nonlinloop_damping_factor"] = omega;  
      if(newton_iterate==1 && omega<=(double)db["afc_damping_bound_newton"])
	newton_flag=1;
    } 
    // if Anderson acceleration, then only reduction of the damping parameter
    // but no computation of different update
    if ((db["afc_nonlinloop_anderson_acc"].is("yes"))&&(iteration >= 
      (int)db["afc_nonlinloop_anderson_acc_start"])) 
      break;
  }  
}
/** ************************************************************************ */
void CD3D::anderson_acceleration_damping(int N_Unknowns, const int iteration, 
     std::list<std::vector<double>> & solAnderson,
     std::list<std::vector<double>> & deltaAnderson)
{
  SystemPerGrid& s = this->systems_.front();
  int k;
  
  std::vector <double> newSol(N_Unknowns);
  std::vector <double> newDelta(N_Unknowns);
  
  int lwork, ierr;
  lwork=-1;
  double temp;
  
  int ndim, kdim;
  double *residuals, *fpast, *work;
  int kAnderson = (int)db["afc_nonlinloop_anderson_acc_vec"];
  BlockVector update, new_solution;
  update = s.solution_;
  
  // update =\hat u^{\nu+1}-\tilde{u^{\nu}}
  update.add_scaled(old_solution,-1.0);

  //collect data for previous iterations
  for(int i=0;i<N_Unknowns;i++)
  {
    newDelta[i]=s.solution_[i];
    newSol[i]=update[i];  
  }
  //update list deltaAnderson
  //add new element at the end
  deltaAnderson.push_back(newDelta);
  //if list has sufficiently many entries, remove first entry
  if((int)deltaAnderson.size()>kAnderson)
    deltaAnderson.pop_front();
  
  //update list solAnderson
  solAnderson.push_back(newSol);
  if((int)solAnderson.size()>kAnderson)
    solAnderson.pop_front();
  
  //Use Anderson acceleration if there are sufficiently many entries
  if((iteration>kAnderson+(int)db["afc_nonlinloop_anderson_acc_start"]) && ((int)solAnderson.size()==kAnderson))
  {
    ndim=N_Unknowns;
    kdim=kAnderson;
    //entries for solution
    residuals=new double[2*ndim*kdim];
    //entries for updates
    fpast=residuals+ndim*kdim;
    
    //copy arrays
    std::list <std::vector<double>>::iterator i_delta = deltaAnderson.begin();
    std::list <std::vector<double>>::iterator i_sol = solAnderson.begin();
    for(int i=0;i<kdim;i++)
    {
      //advance list iterators to next element
      for(int j=0;j<ndim;j++)
      {
	if(i<kdim) 
	{
	  k=i*ndim+j;
	  fpast[k]=i_delta->at(j);
	  residuals[k]=i_sol->at(j);
	}
      }
      std::advance(i_delta,1);
      std::advance(i_sol,1);  
    }     //end of loop i
    //set initial length of work array=-1
    if(lwork==-1)
    {
      //the first call sets the length of work
      //anderson_acceleration(ndim, kdim, &s.solution[0], fpast, residuals,
	//		    &temp, lwork, &ierr);  
      anderson_acceleration(ndim, kdim, &s.solution_[0], fpast, residuals,
			    &temp, lwork, &ierr, &alphas_x_i[0]);  
    }
    lwork=(int)temp;
    work=new double[lwork];
    Output::print<2>("Call anderson acceleration, iteration=", iteration-kAnderson-(int)db["afc_nonlinloop_anderson_acc_start"]);
    anderson_acceleration(ndim, kdim, &s.solution_[0], fpast, residuals,
			  work, lwork, &ierr, &alphas_x_i[0]);
    //s.solution = new_solution;
    delete[] work;
    delete[] residuals;  
  }
}
/** ************************************************************************ */

/** ************************************************************************ */
void CD3D::output(int i)
{
#ifdef _MPI
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

	bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
	if(no_output)
		return;

  SystemPerGrid& syst = systems_.front() ;

  // print the value of the largest and smallest entry in the FE vector
  syst.feFunction_.PrintMinMax();

#ifdef _MPI
  // computing errors as well as writing vtk files requires a minimum 
  // consistency level of 1
  syst.feSpace_->get_communicator().consistency_update(
    syst.solution_.get_entries(), 1);
#endif // _MPI
  
  // write solution to a vtk file
  outputWriter.add_fe_function(&syst.feFunction_);
  outputWriter.write();

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    double errors[5];
    TAuxParam3D aux(1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
    MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
    const TFESpace3D* space = syst.feFunction_.GetFESpace3D();

    syst.feFunction_.GetErrors(example_.get_exact(0), 4, AllDerivatives,
                               2, L2H1Errors, example_.get_coeffs(),
                               &aux, 1, &space, errors);
#ifdef _MPI
    double errorsReduced[4]; //memory for global (across all processes) error

    MPI_Allreduce(errors, errorsReduced, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(i=0;i<2;i++)
      errors[i] = sqrt(errorsReduced[i]);
#else
    int my_rank =0;
#endif

    //store errors
    errors_.at(0) = errors[0]; //L2
    errors_.at(1) = errors[1]; //H1-semi

    //print errors
    if(my_rank == 0)
    {
      Output::print("");
      Output::print( "L2: ", errors_.at(0));
      Output::print( "H1-semi: ", errors_.at(1));
    }
  } // if(this->db["compute_errors"]S)
}

/** ******************************************************************************** **/
AlgebraicFluxCorrection::Limiter string_to_limiter(std::string afc_limiter);
/** ******************************************************************************** **/

void CD3D::do_algebraic_flux_correction(const int iteration, const int is_not_afc_fixed_point_rhs)
{
  Output::print<4>("AFC: enter do_algebraic_flux_correction");
  for(auto & s : this->systems_) // do it on all levels.
  {
    bool compute_D_and_gamma = false;
    
    //determine which kind of afc to use
    if(db["algebraic_flux_correction"].is("default") ||
        db["algebraic_flux_correction"].is("afc"))
    {
      // determine which kind of afc to use
      AlgebraicFluxCorrection::Limiter limiter = string_to_limiter(db["afc_limiter"]);  
      // determine which kind of iteration scheme is active
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
      /*
      if(residual<up_param && time_newton<=2.0 && rhs_flag==0 && db["afc_iteration_scheme_automatic"]=="yes")
      {
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
      else if(residual<up_param && time_rhs<=2.5 && rhs_flag==0&& db["afc_iteration_scheme_automatic"]=="yes")
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
      if(residual<up_param && rhs_flag==1&& db["afc_iteration_scheme_automatic"]=="yes")
      {
	it_scheme=AlgebraicFluxCorrection::Iteration_Scheme::NEWTON;
      }
      if(residual>=up_param && newton_flag==1&& db["afc_iteration_scheme_automatic"]=="yes")
      {
	up_param*=1e-02;
	newton_flag=0;
      }
      Output::print<4>("Iterate RHS: ", rhs_iterate);
      Output::print<4>("Time RHS: ", time_rhs);
      Output::print<4>("Iterate Newton: ", newton_iterate);
      Output::print<4>("Time NEWTON: ", time_newton);
      */
      
     // automatic choice of the scheme for the next iterate
     if (db["afc_iteration_scheme_automatic"].is("yes"))
     {
       if (db["afc_iteration_scheme"].is("fixed_point_matrix") 
	 && (double)db["afc_nonlinloop_damping_factor"] <=  
	 (double)db["afc_nonlinloop_switch_newton_to_fprhs"] * (double)db["afc_nonlinloop_damping_factor_min"])
       {
	 db["afc_iteration_scheme"].set<>("fixed_point_rhs");
	 Output::print<2>("afc_iteration_scheme changed to fixed_point_rhs");
       }
       if (db["afc_iteration_scheme"].is("fixed_point_rhs") 
	 && (double)db["afc_nonlinloop_damping_factor"]  > (double)db["afc_nonlinloop_switch_fprhs_to_newton"] 
	 * (double)db["afc_nonlinloop_damping_factor_min"])
       {
	 db["afc_iteration_scheme"].set<>("fixed_point_matrix");
	 Output::print<2>("afc_iteration_scheme changed to fixed_point_matrix");
       }
      }
      
      //get pointers/references to the relevant objects
      TFESpace3D& feSpace = *s.feSpace_;
      FEMatrix& one_block = *s.matrix_.get_blocks_uniquely().at(0).get();
      const std::vector<double>& solEntries = s.solution_.get_entries_vector();
      std::vector<double>& rhsEntries = s.rhs_.get_entries_vector();

      // fill a vector "neumannToDirichlet" with those rows that got
      // internally treated as Neumann although they are Dirichlet
      int firstDiriDof = feSpace.GetActiveBound();
      int nDiri = feSpace.GetN_Dirichlet();

      std::vector<int> neumToDiri(nDiri, 0);
      std::iota(std::begin(neumToDiri), std::end(neumToDiri), firstDiriDof);
      

      
      if (iteration >0 || db["afc_initial_iterate"].is("afc_zero"))
      {
        // if necessary, set up vector gamma and matrix D
        if(s.afc_matrix_D_entries.empty())
        {
  	  Output::print<4>("AFC: allocate matrix D");
	  s.afc_matrix_D_entries.resize(one_block.GetN_Entries(),0.0);
	  compute_D_and_gamma = true;
        }
        if ((limiter== AlgebraicFluxCorrection::Limiter::BJK17) && (s.afc_gamma.empty()))
        {
	  Output::print<4>("AFC: vector gamma");
	  s.afc_gamma.resize(feSpace.GetN_DegreesOfFreedom(),0.0);
        }
      
        // apply AFC 
        if(is_not_afc_fixed_point_rhs==1)
	{
	  AlgebraicFluxCorrection::steady_state_algorithm(
	    one_block,
	    solEntries,rhsEntries,
	    neumToDiri, 
	    s.afc_matrix_D_entries, s.afc_gamma, compute_D_and_gamma, db, 
	    limiter, it_scheme, is_not_afc_fixed_point_rhs);
	  //performed only once in the whole iteration process
	  if(iteration==1 && it_scheme==AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS)
	  {
	    //matrix_copy=A+D 
	    matrix_copy=s.matrix_;
	  }	 
	}
	//case for fixed point rhs and iteration>1
	else
	{
	  /*
	   * the matrix that is used for the AFC scheme needs to have the correct Dirichlet entries
	   * and hence after the first iteration sending these values.
	   */
	  FEMatrix& one_block1 = *matrix_copy.get_blocks_uniquely().at(0).get();
	  AlgebraicFluxCorrection::steady_state_algorithm(
	    one_block1,
	    solEntries,rhsEntries,
	    neumToDiri, 
	    s.afc_matrix_D_entries, s.afc_gamma, compute_D_and_gamma, db, 
	    limiter, it_scheme, is_not_afc_fixed_point_rhs);
	}
        // Previous Implementation
         /*AlgebraicFluxCorrection::steady_state_algorithm(
          one_block,
          solEntries,rhsEntries,
          neumToDiri, 
	  s.afc_matrix_D_entries, s.afc_gamma, compute_D_and_gamma, db,
	  limiter, it_scheme, is_not_afc_fixed_point_rhs);
	Output::print<1>("Matrix Norm: ",s.matrix_.get_blocks_uniquely({{0,0}})[0]->GetNorm());*/
      }
      //Previous Implementation
      //AlgebraicFluxCorrection::correct_dirichlet_rows(one_block);
      if (is_not_afc_fixed_point_rhs)
	AlgebraicFluxCorrection::correct_dirichlet_rows(one_block);
      //...and in the right hand side, too, assum correct in solution vector
      s.rhs_.copy_nonactive(s.solution_);
    }
    else
    {
      ErrThrow("The chosen algebraic flux correction scheme is unknown "
          "to class CD3D.");
    }
  }
}

/** ************************************************************************ */
AlgebraicFluxCorrection::Limiter string_to_limiter(std::string afc_limiter)
{
  if (afc_limiter == std::string("kuzmin"))
    return AlgebraicFluxCorrection::Limiter::KUZMIN;
  else if (afc_limiter == std::string("BJK17"))
    return AlgebraicFluxCorrection::Limiter::BJK17;
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
  else if (afc_iteration_scheme == std::string("newton_regu"))
    return AlgebraicFluxCorrection::Iteration_Scheme::NEWTON_REGU;
  else 
  {
    ErrThrow("afc_iteration_scheme ", afc_iteration_scheme, " not implemented!!!"); 
  }
}
/* *************************************************************************** */

void CD3D::checkParameters()
{
  //check if the correct problem type is set, change eventually
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
  if(db["space_discretization_type"].is("galerkin"))
    global_space_type=1;
  else if(db["space_discretization_type"].is("supg"))
    global_space_type=2;
  else if(db["space_discretization_type"].is("upwind"))
    global_space_type=1;
  else
    ErrThrow("space_discretization_type ", db["space_discretization_type"].get_name(), " not implemented !");   
  

  //an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }
}

void CD3D::call_assembling_routine(SystemPerGrid& s, LocalAssembling3D& local_assem)
{//FIXME the body of this function was copy and paste

  const TFESpace3D * fe_space = s.feSpace_.get();
  BoundCondFunct3D * boundary_conditions = fe_space->getBoundCondition();
  int N_Matrices = 1;
  double * rhs_entries = s.rhs_.get_entries();

  BoundValueFunct3D * non_const_bound_value[1] {example_.get_bd()[0]};

  //fetch stiffness matrix as block
  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_uniquely();
  TSquareMatrix3D * block[1]{reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get())};

  // Do the Assembling!

  
  // reset right hand side and matrix to zero
  s.rhs_.reset();
  block[0]->reset();
  //and call the method
  Assemble3D(1, &fe_space, N_Matrices, block, 0, nullptr, 1, &rhs_entries,
             &fe_space, &boundary_conditions, non_const_bound_value, local_assem);

}

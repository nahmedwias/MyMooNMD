#include "ConvectionDiffusion_AFC.h"
#include "Database.h"
#include "AlgebraicFluxCorrection.h"
#include "LinAlg.h"
#include "anderson.h"


#ifdef __2D__
 #include "Assemble2D.h"
 #include "SquareMatrix2D.h"
 #include "Example_CD2D.h"
 #include "AuxParam2D.h"
 #include "Upwind.h"
#else
 #include "Assemble3D.h"
 #include "SquareMatrix3D.h"
 #include "AuxParam3D.h" 
 #include "Example_CD3D.h"
 #include "Upwind3D.h"
#endif

#ifdef _MPI
#include "ParFECommunicator3D.h"
#endif

/* ************************************************************************* */
/**
 * @brief: While using AFC, system matrices is created as if
 *         all DOFs were active.
 */
const ParameterDatabase& test(const ParameterDatabase& db)
{
  if(!db["algebraic_flux_correction"].is("none"))
  {                                               
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
  }
  return db;
}

/* ************************************************************************ */

template<int d>
ConvectionDiffusion_AFC<d>::ConvectionDiffusion_AFC(const TDomain& domain, 
						 const ParameterDatabase& param_db)
: ConvectionDiffusion<d>(domain, test(param_db), Example_CD(param_db)),
                         alphas_x_i(), old_solution()
{
  // A default AFC database. Use to merge AFC database with Parameter DB
  auto afc_db(AlgebraicFluxCorrection::default_afc_database());
  afc_db.merge(param_db);
  this->db.merge(afc_db, true);
  
  set_AFC_parameters();
  //Currently, AFC doesn't support Multigrid
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("AFC doesn't support Multigrid!");
  }
  alphas_x_i=this->ConvectionDiffusion<d>::systems.front().solution;
  old_solution=this->ConvectionDiffusion<d>::systems.front().solution;
  is_not_afc_fixed_point_rhs=1;
  rejected_steps=0;
  
  /*
   * Creates an array to store the interpolation of the exact solution. Useful
   * in finding AFC norm and d_h( ; , ) error.
   */
  
  auto& space = this->ConvectionDiffusion<d>::systems.front().fe_space;
  exact_interpolant.resize(space->GetN_DegreesOfFreedom());
  
  //creates a FEFUnction  
#ifdef __2D__
  TFEFunction2D exact_interpolation(space, "interpolant",
                                    "interpolant of exact solution",
                                    exact_interpolant.data(),
                                    space->GetN_DegreesOfFreedom());
#else
  TFEFunction3D exact_interpolation(space, "interpolant",
				    "interpolant of exact solution",
				    exact_interpolant.data(),   space->GetN_DegreesOfFreedom());
#endif
  //Interpolates and store the value in exact_interpolant
  exact_interpolation.Interpolate(this->example.get_exact(0));  
}

/** ************************************************************************ */
template<int d>
void ConvectionDiffusion_AFC<d>::set_AFC_parameters()
{
  if(!this->db["algebraic_flux_correction"].is("none"))
  {
    if(!this->db["algebraic_flux_correction"].is("afc"))
    {
      this->db["algebraic_flux_correction"].set("afc");
      Output::print("Only kind of algebraic flux correction"
        " for CD problems is AFC (afc).");
    }
    //make sure that galerkin discretization is used
    if (!this->db["space_discretization_type"].is("galerkin"))
    {                                            
      this->db["space_discretization_type"] = "galerkin";
      Output::warn<1>("Parameter 'space_discretization_type' changed to", 
       "Galerkin because Algebraic Flux Correction is enabled.");
    }
  }
}
/** ************************************************************************ */
template<int d>
void ConvectionDiffusion_AFC<d>::assemble(const int iteration)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  //determine the local assembling type to be ConvectionDiffusion
  LocalAssembling_type laTypet = LocalAssembling_type::ConvDiff;
  
  if(iteration==0)
    time_total=GetTime();
  
  auto & s = this->ConvectionDiffusion<d>::systems.front();
  FEFunction * feFunctionPtr = &s.fe_function;
  
  int afc_ini_supg = 0;
  if(!this->db["algebraic_flux_correction"].is("none") && iteration==0 
    && this->db["afc_initial_iterate"].is("supg"))
  {
    afc_ini_supg = 1;   
  }
  
  LocalAssembling<d> laObject(this->db, laTypet, &feFunctionPtr,
                                this->example.get_coeffs());
  
  //fetch stiffness matrix as block
  auto blocks = s.matrix.get_blocks_uniquely();
  SquareMatrixD* block = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
  
  /**
   * @brief: For Fixed Point RHS(FPR), the matrix and the RHS has to assembled 
   * only once. Becasue in FPR the rhs gets modified in 
   * do_algebraic_flux_correction, where the modification is only added and
   * hence we copy the old rhs.
   * 
   */
  if(is_not_afc_fixed_point_rhs==1)  
  {
    ConvectionDiffusion<d>::call_assembling_routine(s, laObject);
    if(iteration==1)
      rhs_copy=s.rhs;    
  }
  //for FIXED_POINT_RHS and iteration>1
  else
  {
    s.rhs=rhs_copy;  
  }
  if(is_not_afc_fixed_point_rhs==1)
  {
    if (!this->db["algebraic_flux_correction"].is("none")&&(iteration==0)&&
          this->db["afc_initial_iterate"].is("upwind"))
    {
      Output::print<2>("upwind for convection-diffusion equation");
#ifdef __2D__
      UpwindForConvDiff(laObject.GetCoeffFct(), block, s.rhs.get_entries(),
                        s.fe_space.get(),nullptr, nullptr, false);
#else
      UpwindForConvDiff(block, &s.rhs.get_entries_vector()[0], 
                          s.fe_space.get(),laObject);  
#endif  
    }  
  }
  if (afc_ini_supg == 1)
    this->db["space_discretization_type"].set("galerkin");
  
  if (iteration == 0)
  {
   this->db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    /*// the following flags are for the BAIL proceedings
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
    }*/
  }

  if(!this->db["algebraic_flux_correction"].is("none"))
    do_algebraic_flux_correction(iteration, is_not_afc_fixed_point_rhs);
  
}

/* ************************************************************************ */
AlgebraicFluxCorrection::Iteration_Scheme string_to_it_scheme(
                                          std::string afc_iteration_scheme);
/* ************************************************************************ */
template<int d>
bool ConvectionDiffusion_AFC<d>::solve(const int iteration)
{
  double t = GetTime();
  auto& s = this->systems.front();
  BlockVector res = s.rhs;
  BlockVector current_rhs = s.rhs;
  BlockVector current_sol = s.solution;
  AlgebraicFluxCorrection::Iteration_Scheme it_scheme =
                       string_to_it_scheme(this->db["afc_iteration_scheme"]);
  
  
  /*
   * In case of not FPR, after the first iteration this has to be zero as then
   * we have to assemble the matrix again and again.
   */
                       
  is_not_afc_fixed_point_rhs=0;

  // To assemble for 1st iteration with FPR
  if(it_scheme!=AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS
    || iteration == 1)
    is_not_afc_fixed_point_rhs=1;

  // Special treatment of first iteration
  if (iteration==0)
  {
    // Compute residual vector
    s.matrix.apply_scaled_add(s.solution , res ,-1.0);

    // Compute the norm of the residual vector, residual_old 
    //  is the norm of residual r_k
    residual_old = res.norm();
    old_solution = s.solution;
    this->ConvectionDiffusion<d>::solver.solve(s.matrix, s.rhs, s.solution);

    t = GetTime() - t;
    Output::print(" solving of a Convection Diffusion problem done in ",
                   t, " seconds");
    double h_min, h_max;
    
#ifdef __2D__
    const TFESpace2D & space = *this->systems.front().fe_space;
#else
    const TFESpace3D & space = *this->systems.front().fe_space;
#endif
    auto coll = space.GetCollection();
    coll->GetHminHmax(&h_min, &h_max);
    
    // For reqularized Newton method, the sigma is depending on h. Idea from
    // [BB17.CMAME]
    double sigma = (double)this->db["afc_newton_regu_sigma"];
    sigma *= h_max * h_max * h_max * h_max;
    this->db["afc_newton_regu_sigma"] = sigma;
    Output::print<2>("afc_newton_regu_sigma changed to ", 
                                    this->db["afc_newton_regu_sigma"]);   
    
    // The stopping tolerance is dependent on mesh refinment
    double dof = space.GetN_DegreesOfFreedom();
    double eps = (double)this->db["afc_nonlinloop_epsilon"];
    eps *= sqrt(dof);
    this->db["afc_nonlinloop_epsilon"] = eps;
    Output::print<2>("afc_nonlinloop_epsilon normalized to ", 
                     this->db["afc_nonlinloop_epsilon"]);
    Output::print<2>("nonlinear step ", iteration, " residual: ", residual_old);
    return(false);
  }

  if(!this->db["algebraic_flux_correction"].is("none"))
  {
    //Anderson acceleration
    if ((this->db["afc_nonlinloop_anderson_acc"].is("yes"))&&(iteration >=
      (int)this->db["afc_nonlinloop_anderson_acc_start"]))
    {
      int N_Unknowns=res.length();
      anderson_acceleration_damping(N_Unknowns,iteration, 
                                    solAnderson, deltaAnderson);
    }

    //Constant Damping
    if (this->db["afc_nonlinloop_damping_factor_constant"].is("yes"))
    {
      // compute proposal for the next solution
      AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, 
                                                       s.solution, this->db);
      assemble(iteration);
      // calculation of residual r_k+1
      res = s.rhs;
      s.matrix.apply_scaled_add(s.solution , res ,-1.0);
      // compute the norm of the residual vector, residual_old is the 
      // norm of residual r_k
      residual = res.norm();
      Output::print<4>("  residual for proposed new iterate ", residual);
    }
    
    //Dynamic Damping from [JK08.CMAME]
    else
      dynamic_damping(iteration);
  }

  old_solution.add_scaled(s.solution,-1.0);
  Output::print<2>("nonlinear step ", iteration, " residual: ", residual, 
                   " reduction ",residual/residual_old ," change of sol ", 
                                                      old_solution.norm());
  
  // stopping criterion satisfied
  if ((residual < (double)this->db["afc_nonlinloop_epsilon"])&&(iteration >1))
  {
    time_total = GetTime()-time_total;
    Output::print<1>("NONLINEAR ITERATION: ite ", iteration, " res ", 
                                                    residual, " rejections ",
		     rejected_steps, " time ", time_total, " t/it ", 
		     time_total/(iteration+rejected_steps));
    return(true);
  }
  //maximal number of iterations
  if (iteration == (int)this->db["afc_nonlinloop_maxit"])
  {
    time_total = GetTime()-time_total;
    Output::print<1>("MAX_NONLINEAR ITERATION: ite ", iteration, " res ", 
                     residual, " rejections ", rejected_steps, " time ",
                     time_total, " t/it ", time_total/(iteration+rejected_steps));
    return(true);
  }
  // storage of old solution
  old_solution =s.solution;

  // Solve the linear system
  //Matrix is factorised only once and then stored in the system
  if(is_not_afc_fixed_point_rhs==1)
    this->ConvectionDiffusion<d>::solver.update_matrix(s.matrix);
  this->ConvectionDiffusion<d>::solver.solve(s.rhs, s.solution);
  
  //Previous Implementation
  //this->solver.solve(s.matrix_, s.rhs_, s.solution_);
  
  t = GetTime() - t;
  Output::print<4>("  iteration done in ", t, " seconds");

  //Dynamic choice of Newton damping parameter
  if ((this->db["afc_iteration_scheme"].is("newton"))||
       (this->db["afc_iteration_scheme"].is("newton_regu")))
  {
    double threshold = (double)this->db["afc_change_method_threshold"];
    if (residual <= threshold)
    {
      double omega_derivative = 
                 (double)this->db["afc_fixed_point_derivative_weight_factor"];
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
      this->db["afc_fixed_point_derivative_weight_factor"] = omega_derivative;
    }
    // go back to former method
    if (residual > 100*threshold)
    {
      this->db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }
    Output::print<2>("afc_fixed_point_derivative_weight_factor ",
		     (double)this->db["afc_fixed_point_derivative_weight_factor"]);
  }

  /* ************************************************************************ */
  // this scheme is for the BAIL proceedings
  /* ************************************************************************ */
  /*if ((db["afc_iteration_scheme"].is("newton"))&&
     ((int)db["afc_nonlinloop_switch_to_newton_scheme"]>=10))
  {
    double threshold = (double)db["afc_change_method_threshold"];
    // switch to formal Newton
    if (residual <= threshold)
    {
      db["afc_fixed_point_matrix_weight"] = 1.0;
      db["afc_fixed_point_derivative_weight_factor"] = 1.0;
    }
    else
    {
        if (residual > 10*threshold)
    {
        db["afc_fixed_point_matrix_weight"] = 0.0;
        db["afc_fixed_point_derivative_weight_factor"] = 0.0;
    }
    }
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
    Output::print<2>("afc_fixed_point_matrix_weight ",  
                          (double)db["afc_fixed_point_matrix_weight"]);
  }*/
  
  //Storage of old residual
  residual_old=residual;
  return(false);
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion_AFC<d>::dynamic_damping(const int iteration)
{
  auto& s = this->systems.front();
  double omega = this->db["afc_nonlinloop_damping_factor"];
  int first_damp=1;
  BlockVector res = s.rhs;
  BlockVector current_rhs = s.rhs;
  BlockVector current_sol = s.solution;
  int check_anderson=0;
  if(iteration<=(int)this->db["afc_nonlinloop_anderson_acc_vec"]+
                 (int)this->db["afc_nonlinloop_anderson_acc_start"])
    check_anderson=1;
  while(1)
  {
    // Damping from [WN11.SINUM]
    if(this->db["afc_anderson_damping"].is("yes"))
    {
      if(this->db["afc_nonlinloop_anderson_acc"].is("yes"))
      {
        if(check_anderson==1)
          AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, 
                                                           s.solution, this->db);
        //starts anderson damping after the vectors are stored
        else
        {
          s.solution.scale(omega);
          s.solution.add_scaled(alphas_x_i, 1.-omega);
        }
      }
      else
        AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, 
                                                         s.solution, this->db);
      assemble(iteration);
      // calculation of residual r_k+1
      res = s.rhs;
      s.matrix.apply_scaled_add(s.solution , res ,-1.0);
      // compute the norm of the residual vector, 
      // residual_old is the norm of residual r_k
      residual = res.norm();
      Output::print<4>("  residual for proposed new iterate ", residual);
      // accept the first damping parameter
      if (iteration==1)
        break;
      if(this->db["afc_nonlinloop_anderson_acc"].is("yes") && check_anderson==1)
        break;
    }
    // Dynamic damping from [JK08.CMAME]
    else
    {
      AlgebraicFluxCorrection::AFC_Compute_New_Iterate(old_solution, s.solution, 
                                                       this->db);
      assemble(iteration);
      // calculation of residual r_k+1
      res = s.rhs;
      s.matrix.apply_scaled_add(s.solution , res ,-1.0);
      // compute the norm of the residual vector, 
      //    residual_old is the norm of residual r_k
      residual = res.norm();
      Output::print<4>("  residual for proposed new iterate ", residual);
      // accept the first damping parameter
      if (iteration==1)
        break;
    }
    /*
     * if the norm of the residual vector decreases or if the damping parameter 
     * is already very small then accept the new iterate
     */
    if (residual<residual_old ||
               omega<=(double)this->db["afc_nonlinloop_damping_factor_min_tol"]*
      (double)this->db["afc_nonlinloop_damping_factor_min"])
    {
      /*
       * if the norm of the residual decreases without having decreased the
       * damping parameter in this step, then increase the damping parameter 
       * if possible
       */
      if(residual<residual_old &&first_damp==1)
      {
        this->db["afc_nonlinloop_damping_factor_max"]=
        std::min((double)this->db["afc_nonlinloop_damping_factor_max_global"],
          (double)this->db["afc_nonlinloop_damping_factor_max_increase"]
          *(double)this->db["afc_nonlinloop_damping_factor_max"]);
        omega=std::min((double)this->db["afc_nonlinloop_damping_factor_max"],
        (double)this->db["afc_nonlinloop_damping_factor_increase"]*omega);
      }
      Output::print<2>("  iterate accepted, damping factor ", omega, " ",
          this->db["afc_nonlinloop_damping_factor_max"]);
      this->db["afc_nonlinloop_damping_factor"] = omega;
      break;
    }
    else
    {
      if (this->db["afc_nonlinloop_anderson_acc"].is("no"))
      {
        rejected_steps++;
        // get starting situation back
        s.solution = current_sol;
        s.rhs = current_rhs;
      }
      // reduce damping factor
      omega=std::max((double)this->db["afc_nonlinloop_damping_factor_min"],
        omega*(double)this->db["afc_nonlinloop_damping_factor_decrease"]);
      // reduce maximal value for damping factor
      if(first_damp==1)
      {
        this->db["afc_nonlinloop_damping_factor_max"]=
        std::max((double)this->db["afc_nonlinloop_damping_factor_min"],
          (double)this->db["afc_nonlinloop_damping_factor_max_decrease"]
              *(double)this->db["afc_nonlinloop_damping_factor_max"]);
        first_damp=0;
      }
      Output::print<2>("  iterate rejected, res old ", residual_old,
          " res new: ", residual, " damping factor :", omega);
      this->db["afc_nonlinloop_damping_factor"] = omega;
    }
    // if Anderson acceleration, then only reduction of the damping parameter
    // but no computation of different update
    if ((this->db["afc_nonlinloop_anderson_acc"].is("yes"))&&(iteration >=
      (int)this->db["afc_nonlinloop_anderson_acc_start"]))
      break;
  }
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion_AFC<d>::anderson_acceleration_damping(
  int N_Unknowns, const int iteration,
  std::list<std::vector<double>> & solAnderson,
  std::list<std::vector<double>> & deltaAnderson)
{
  auto& s = this->systems.front();
  int k;

  std::vector <double> newSol(N_Unknowns);
  std::vector <double> newDelta(N_Unknowns);

  int lwork, ierr;
  lwork=-1;
  double temp;

  int ndim, kdim;
  double *residuals, *fpast, *work;
  int kAnderson = (int)this->db["afc_nonlinloop_anderson_acc_vec"];
  BlockVector update, new_solution;
  update = s.solution;

  // update =\hat u^{\nu+1}-\tilde{u^{\nu}}
  update.add_scaled(old_solution,-1.0);

  //collect data for previous iterations
  for(int i=0;i<N_Unknowns;i++)
  {
    newDelta[i]=s.solution[i];
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
  if((iteration>kAnderson+(int)this->db["afc_nonlinloop_anderson_acc_start"]) &&
    ((int)solAnderson.size()==kAnderson))
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
    }                                             //end of loop i
    //set initial length of work array=-1
    if(lwork==-1)
    {
      //the first call sets the length of work
      //anderson_acceleration(ndim, kdim, &s.solution[0], fpast, residuals,
      //		    &temp, lwork, &ierr);
      anderson_acceleration(ndim, kdim, &s.solution[0], fpast, residuals,
        &temp, lwork, &ierr, &alphas_x_i[0]);
    }
    lwork=(int)temp;
    work=new double[lwork];
    Output::print<2>("Call anderson acceleration, iteration=", 
		     iteration-kAnderson-(int)this->db["afc_nonlinloop_anderson_acc_start"]);
    anderson_acceleration(ndim, kdim, &s.solution[0], fpast, residuals,
      work, lwork, &ierr, &alphas_x_i[0]);
    //s.solution = new_solution;
    delete[] work;
    delete[] residuals;
  }
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion_AFC<d>::output(int i)
{
  ConvectionDiffusion<d>::output(i);
  if(this->db["output_compute_errors"])
  {    
    //AFC Errors
    if(this->db["algebraic_flux_correction"].is("afc"))
    {
      double epsilon=this->ConvectionDiffusion<d>::example.get_nu();
      Output::print<1>("sigma*L2     : ", 
                       setprecision(14), ConvectionDiffusion<d>::errors[0]);
      Output::print<1>("epsilon*H1-semi: ", 
                       setprecision(14), std::sqrt(epsilon)*
                       ConvectionDiffusion<d>::errors[1]);     
      Output::print<1>("Energy Norm^2: ", setprecision(14), 
                       ConvectionDiffusion<d>::errors[0]*
                       ConvectionDiffusion<d>::errors[0]
                       +epsilon*ConvectionDiffusion<d>::errors[1]*
                       ConvectionDiffusion<d>::errors[1]);
    }
  }
}

/* ************************************************************************* */
AlgebraicFluxCorrection::Limiter string_to_limiter(std::string afc_limiter);
/* ************************************************************************* */

template<int d>
void ConvectionDiffusion_AFC<d>::do_algebraic_flux_correction(
  const int iteration, const int is_not_afc_fixed_point_rhs)
{
  Output::print<4>("AFC: enter do_algebraic_flux_correction");
  auto & s = this->ConvectionDiffusion<d>::systems.front();
  bool compute_D_and_gamma = false;
  
  //determine which kind of afc to use
  if(this->db["algebraic_flux_correction"].is("default") ||
      this->db["algebraic_flux_correction"].is("afc"))
  {
    // Determine which kind of limiter to use
    AlgebraicFluxCorrection::Limiter limiter = 
                                    string_to_limiter(this->db["afc_limiter"]);
    // Determine which kind of iteration scheme is active
    AlgebraicFluxCorrection::Iteration_Scheme it_scheme = 
                          string_to_it_scheme(this->db["afc_iteration_scheme"]);
    
    //get pointers/references to the relevant objects
    auto& feSpace = *s.fe_space;
    FEMatrix& one_block = *s.matrix.get_blocks_uniquely().at(0).get();
    const std::vector<double>& solEntries = s.solution.get_entries_vector();
    std::vector<double>& rhsEntries = s.rhs.get_entries_vector();
    
    // fill a vector "neumannToDirichlet" with those rows that got
    // internally treated as Neumann although they are Dirichlet
    int firstDiriDof = feSpace.GetActiveBound();
    int nDiri = feSpace.GetN_Dirichlet();
    
    std::vector<int> neumToDiri(nDiri, 0);
    std::iota(std::begin(neumToDiri), std::end(neumToDiri), firstDiriDof);
    
    //AFC only applied after the 0th iteration or the initial iterate is zero.
    if (iteration >0 || this->db["afc_initial_iterate"].is("afc_zero"))
    {
      // if necessary, set up vector gamma and matrix D
      if(afc_matrix_D_entries.empty())
      {
        Output::print<4>("AFC: allocate matrix D");
        afc_matrix_D_entries.resize(one_block.GetN_Entries(),0.0);
        compute_D_and_gamma = true;  
      }
      if ((limiter== AlgebraicFluxCorrection::Limiter::BJK17) 
                          && (afc_gamma.empty()))
      {
        Output::print<4>("AFC: vector gamma");
        afc_gamma.resize(feSpace.GetN_DegreesOfFreedom(),0.0);  
      }
      
      // apply AFC
      if(is_not_afc_fixed_point_rhs==1)
      {
        AlgebraicFluxCorrection::steady_state_algorithm(
          one_block,
          solEntries,rhsEntries,
          neumToDiri,
          afc_matrix_D_entries, afc_gamma, compute_D_and_gamma, this->db,
          limiter, it_scheme, is_not_afc_fixed_point_rhs);
        //performed only once in the whole iteration process
        if(iteration==1 && 
           it_scheme==AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS)
        {
          //matrix_copy=A+D
          //In AlgebraicFluxCorrection.C we don't need to add D again
          matrix_copy=s.matrix;  
        }  
      }
      //case for fixed point rhs and iteration>1
      else
      {
        /*
         * the matrix that is used for the AFC scheme needs to have the correct 
         * Dirichlet entries and hence after the first iteration sending these 
         * values.
         */
        FEMatrix& one_block1 = *matrix_copy.get_blocks_uniquely().at(0).get();
        AlgebraicFluxCorrection::steady_state_algorithm(
          one_block1,
          solEntries,rhsEntries,
          neumToDiri,
          afc_matrix_D_entries, afc_gamma, compute_D_and_gamma, this->db,
          limiter, it_scheme, is_not_afc_fixed_point_rhs);  
      }
      // Previous Implementation
      /*AlgebraicFluxCorrection::steady_state_algorithm(
       *         one_block,
       *         solEntries,rhsEntries,
       *         neumToDiri,
       *        s.afc_matrix_D_entries, s.afc_gamma, compute_D_and_gamma, db,
       *        limiter, it_scheme, is_not_afc_fixed_point_rhs);
       *        Output::print<1>("Matrix Norm: ",
       *                  s.matrix_.get_blocks_uniquely({{0,0}})[0]->GetNorm());
       */  
    }
    //Previous Implementation
    //AlgebraicFluxCorrection::correct_dirichlet_rows(one_block);
    if (is_not_afc_fixed_point_rhs)
      AlgebraicFluxCorrection::correct_dirichlet_rows(one_block);
    //...and in the right hand side, too, assum correct in solution vector
    s.rhs.copy_nonactive(s.solution);  
  }
  else
  {
    ErrThrow("The chosen algebraic flux correction scheme is unknown "
              "to class ConvectionDiffusion<",d,">.");  
  }
}


/* ************************************************************************* */

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

/* ************************************************************************* */
AlgebraicFluxCorrection::Iteration_Scheme string_to_it_scheme(
                                              std::string afc_iteration_scheme)
{
  if (afc_iteration_scheme == std::string("fixed_point_rhs"))
    return AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_RHS;
  else if (afc_iteration_scheme == std::string("fixed_point_matrix"))
    return AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_MATRIX;
  else if (afc_iteration_scheme == std::string("newton"))
    return AlgebraicFluxCorrection::Iteration_Scheme::NEWTON;
  else if (afc_iteration_scheme == std::string("newton_no_damp"))
    return AlgebraicFluxCorrection::Iteration_Scheme::NEWTON;
  else if (afc_iteration_scheme == std::string("newton_regu"))
    return AlgebraicFluxCorrection::Iteration_Scheme::NEWTON_REGU;
  else
  {
    ErrThrow("afc_iteration_scheme",afc_iteration_scheme, "not implemented!!!");
  }
}

#ifdef __3D__
template class ConvectionDiffusion_AFC<3>;
#else
template class ConvectionDiffusion_AFC<2>;
#endif
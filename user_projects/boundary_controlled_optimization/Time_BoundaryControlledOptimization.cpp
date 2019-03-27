#include "Time_BoundaryControlledOptimization.hpp"
#include "TimeDiscretizations.h"
#include "LoopInfo.h"
#include "Database.h"
#include "BaseCell.h"
#include "BoundEdge.h"
#include <algorithm>
#include <BoundEdge.h>
#include <BaseCell.h>

template<int d>
ParameterDatabase Time_BoundaryControlledOptimization<d>::default_BCO_database()
{
  // is only needed for the adjoint problem, but there is only one global 
  // parameter for this.
  TDatabase::ParamDB->NSTYPE = 4;
  // this is set here to have the full structure and assembled entries in the
  // Dirichlet rows. This is really only needed for the adjoint matrix to 
  // compute the gradient of the functional. So, for the assembling of the 
  // primal system this is tunred off again.
  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
  auto db = ParameterDatabase::parmoon_default_database();
  db.add("controlled_boundary_component", 1, 
         "The boundary component on which the control is defined. All degrees "
         "of freedom on this component are then part of the control space.");
  
  db.add("alpha_cost", 1., 
         "The scalar alpha in the functional J_hat, which is multiplied with "
         "the control term.", 0., 1.0e10);
  
  db.add("cost_functional", {1., 0., 0.},
         "Switch between a few possible cost functionals via weights. The "
         "order of the functionals is: L2_norm_of_curl, backward_facing_step, " 
         "L2_norm_diff_stokes", 0., 1.);

  db.add("restricted_curl_functional", false, 
         "The L2_norm_of_curl functional can be restricted such the the "
         "integral is taken only on a subdomain Omega intersected with x<10 ");
  
  db.add("control_depends_on_time", true,
         "Activate it if you want to allow the control dofs to vary in time. The "
         "control space size will be multiplied by the number of time steps.");

  db.add("control_depends_on_space", true,
         "Activate it if you want to allow the control dofs to vary in space. The "
         "control space size will be multiplied by the number of time steps.");

  db.add("control_in_x_direction", true,
         "Activate it if you want that the control space includes the x component"
         "at the boundary.");

  db.add("control_in_y_direction", true,
         "Activate it if you want that the control space includes the y component"
         "at the boundary.");

  return db;
}

template<int d>
bool Time_BoundaryControlledOptimization<d>::check_input_parameter_consistency(const ParameterDatabase& param_db)
{
  ParameterDatabase db = Time_BoundaryControlledOptimization::default_BCO_database();
  db.merge(param_db, false);

  bool consistent = true;

  if (!db["control_in_x_direction"] && !db["control_in_y_direction"])
  {
    Output::warn("check_input_parameter_consistency", "the control should have at "
        "least one component, either x or y!");
    consistent = false;
  }

  std::vector<double> functional_weights = db["cost_functional"];
  double total_weights = 0;
  for (unsigned int i = 0; i < functional_weights.size(); ++i)
    total_weights += functional_weights.at(i);
  if (total_weights == 0 )
  {
    Output::warn("check_input_parameter_consistency", "the cost functional is zero!");
    consistent = false;
  }

  if (db["control_depends_on_time"] && !db["control_depends_on_space"])
  {
    Output::info<5>("check_input_parameter_consistency", " The control is "
        "only time-dependent!");
  }

  if (!db["control_depends_on_time"] && db["control_depends_on_space"])
  {
    Output::info<5>("check_input_parameter_consistency", " The control is "
        "only space-dependent, although you are running optimization on time-dependent "
        "problems!");
  }

  if (!db["control_depends_on_time"] && !db["control_depends_on_space"])
  {
    Output::info<5>("check_input_parameter_consistency", " The control is "
        "independent from time and space!");
  }

  return consistent;
}

template<int d>
ParameterDatabase Time_BoundaryControlledOptimization<d>::get_primal_database(const ParameterDatabase& param_db)
{
  return param_db;
}

template<int d>
ParameterDatabase Time_BoundaryControlledOptimization<d>::get_adjoint_database(const ParameterDatabase& param_db)
{
  ParameterDatabase adjoint_db(ParameterDatabase::parmoon_default_database());
  adjoint_db.merge(ParameterDatabase::default_output_database(), true);
  adjoint_db.merge(param_db, true, true);
  std::string adjoint_basename = param_db["output_basename"];
  adjoint_basename += std::string("bco_adjoint");
  adjoint_db["output_basename"].set(adjoint_basename, false);
  adjoint_db["nonlinloop_damping_factor"] = 1.;
  return adjoint_db;
}

template<int d>
Time_BoundaryControlledOptimization<d>::tnse_primal_solution::tnse_primal_solution()
{

}

template<int d>
Time_BoundaryControlledOptimization<d>::tnse_primal_solution::tnse_primal_solution(
    const BlockVector& solution_vector, int step, const FEVectFunct& u_sol,
    const FEFunction& p_sol)
: vector_at_timestep_t_(solution_vector), timestep_t_(step)
{
  // the fe functions must be newly created, because copying would mean
  // referencing the BlockVectors in 'other'.
#ifdef __2D__
  u_at_timestep_t_ = TFEVectFunct2D(u_sol.GetFESpace2D(), "u", "u",
                                    vector_at_timestep_t_.block(0),
                                    vector_at_timestep_t_.length(0), 2);
  p_at_timestep_t_ = TFEFunction2D(p_sol.GetFESpace2D(), "p", "p",
                                   vector_at_timestep_t_.block(2),
                                   vector_at_timestep_t_ .length(2));
#else
  u_at_timestep_t_ = TFEVectFunct3D(u_sol.GetFESpace3D(), "u", "u",
                                    vector_at_timestep_t_.block(0),
                                    vector_at_timestep_t_.length(0), 3);
  p_at_timestep_t_ = TFEFunction3D(p_sol.GetFESpace3D(), "p", "p",
                                   vector_at_timestep_t_.block(3),
                                   vector_at_timestep_t_ .length(3)); 
#endif
}

template<int d>
Time_BoundaryControlledOptimization<d>::Time_BoundaryControlledOptimization(
  const TDomain& domain, const ParameterDatabase& param_db)
 : db(default_BCO_database()), n_control(0), n_components(0), 
   n_dof_per_component(0), n_control_per_time_step(0), 
   control_index_x(), control_index_y(), control_dofs(),
   tnse_primal(domain, get_primal_database (param_db)),
   tnse_adjoint(tnse_primal, get_adjoint_database(param_db)),
   stokes_fe_vector(new BlockVector(tnse_primal.get_solution())),
   stokes_sol(
     new FEVectFunct(tnse_primal.get_velocity_space(), "stokes_sol", "",
                     stokes_fe_vector->get_entries(),
                     tnse_primal.get_velocity_space()->GetN_DegreesOfFreedom(),
                     2)),
   n_time_steps_(0),
   tnse_primal_solutions_(),
   optimization_info("optimization", true, true, 1), // 1 -> full verbosity
   nonlinear_info("nonlinear", true, true, 3),
   current_J_hat(std::numeric_limits<double>::infinity()), 
   control_old(), n_calls(0),
   n_computation_derivative(0)
{
  Output::print<5>("Creating the Time_BoundaryControlledOptimization object");
  db.merge(param_db, false);
  if (!check_input_parameter_consistency(db))
    ErrThrow("Inconsistent input parameters!");
  
  // Set the number of time steps from database of tnse_primal member
  /* TODO: DONE! (N.A. 07/11/2018, needs to be checked before removing this todo) 
   * this step can be done much more nicely if the global parameter
   * is deleted and replaced by the local "end_time" stored in tnse_primal.
   * Unfortunately, the latter one is not used it. Once it is the case,
   * n_time_steps_ can be directly set in the initialization of the members,
   * just above (tnse_primal.get_time_stepping_scheme().get_end_time() instead
   * of 0)
   *
   */
  double dt = this->tnse_primal.get_time_stepping_scheme().get_step_length();
  n_time_steps_ = (this->tnse_primal.get_time_stepping_scheme().get_end_time()/ dt) + 1 ;

  // assemble linear (Stokes) parts
  // TODO: REMOVE THIS STOKES SOLUTION?
  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 0;
//  tnse_primal.assemble();
  // solve a Stokes system once, to store its solution (as 'desired state')
//  tnse_primal.solve();
  *stokes_fe_vector = tnse_primal.get_solution(); // copy
//  tnse_primal.add_to_output(stokes_sol.get());
  
  // the id of the component on which the control is to be applied
  int controlled_boundary_component = db["controlled_boundary_component"];
  // find the degrees of freedom (dof) which are to be controlled
  auto fe_space = tnse_primal.get_velocity_space();
  auto collection = fe_space->GetCollection();
  auto n_cells = collection->GetN_Cells();
  for(int i_cell = 0; i_cell < n_cells; ++i_cell)
  {
    auto cell = collection->GetCell(i_cell);
    auto n_joints = cell->GetN_Joints();
    for(int i_joint = 0; i_joint < n_joints; ++i_joint)
    {
      auto joint = cell->GetJoint(i_joint);
      if(joint->GetType() == BoundaryEdge)
      {
        auto *boundary_joint = dynamic_cast<const TBoundEdge*>(joint);
        int boundary_id = boundary_joint->GetBoundComp()->GetID();
        if(boundary_id == controlled_boundary_component)
        {
          // found a joint which has is on the boundary where the control is at
#ifdef __2D__
          auto fe_descriptor = fe_space->get_fe(i_cell).GetFEDesc2D();
#else
          auto fe_descriptor = fe_space->get_fe(i_cell).GetFEDesc3D();
#endif
          int n_joint_dof = fe_descriptor->GetN_JointDOF();
          auto local_joint_dofs = fe_descriptor->GetJointDOF(i_joint);
          auto local_to_global_dof = fe_space->GetGlobalDOF(i_cell);
          for(int i_joint_dof = 0; i_joint_dof < n_joint_dof; ++i_joint_dof)
          {
            int global_dof = local_to_global_dof[local_joint_dofs[i_joint_dof]];
            control_dofs.push_back(global_dof);
          }
        }
      }
    }
  }
  // second loop of this kind to remove dofs which are also on other boundaries
  for(int i_cell = 0; i_cell < n_cells; ++i_cell)
  {
    auto cell = collection->GetCell(i_cell);
    auto n_joints = cell->GetN_Joints();
    for(int i_joint = 0; i_joint < n_joints; ++i_joint)
    {
      auto joint = cell->GetJoint(i_joint);
      if(joint->GetType() == BoundaryEdge)
      {
        auto *boundary_joint = dynamic_cast<const TBoundEdge*>(joint);
        int boundary_id = boundary_joint->GetBoundComp()->GetID();
        if(boundary_id != controlled_boundary_component)
        {
          // found a joint which has is on an other boundary
#ifdef __2D__
          auto fe_descriptor = fe_space->get_fe(i_cell).GetFEDesc2D();
#else
          auto fe_descriptor = fe_space->get_fe(i_cell).GetFEDesc3D();
#endif
          int n_joint_dof = fe_descriptor->GetN_JointDOF();
          auto local_joint_dofs = fe_descriptor->GetJointDOF(i_joint);
          auto local_to_global_dof = fe_space->GetGlobalDOF(i_cell);
          for(int i_joint_dof = 0; i_joint_dof < n_joint_dof; ++i_joint_dof)
          {
            int global_dof = local_to_global_dof[local_joint_dofs[i_joint_dof]];
            // see https://stackoverflow.com/a/3385251
            auto it = std::remove(control_dofs.begin(), control_dofs.end(),
                                  global_dof);
            control_dofs.erase(it, control_dofs.end()); 
          }
        }
      }
    }
  }
  // sort and remove duplicates
  std::sort(control_dofs.begin(), control_dofs.end());
  auto it = std::unique(control_dofs.begin(), control_dofs.end());
  control_dofs.resize(std::distance(control_dofs.begin(), it));
  control_dofs.shrink_to_fit();
  
//   for(auto e : control_dofs)
//   {
//     Output::print("control_dof ", e);
//   }

  // initialize the control size depending on input parameters
  initialize_control();
   
  Output::print<3>("Created the Time_BoundaryControlledOptimization object, ",
                   "n_control = ", n_control);
}


template<int d>
void Time_BoundaryControlledOptimization<d>::initialize_control()
{
  bool depends_on_time  = db["control_depends_on_time"];
  bool depends_on_space = db["control_depends_on_space"];
  bool x_direction = db["control_in_x_direction"];
  bool y_direction = db["control_in_y_direction"];
  
  // Set the size of the control
  // minimum value: 1 (time- and space-INdependent)
  // maximum value: 2 * control_dofs.size()*n_time_steps (time- and space-dependent)
  n_control = 1;
  if (depends_on_time)
    n_control *= n_time_steps_;

  n_components = 0;
  if (x_direction)
    n_components++;
  if (y_direction)
    n_components++;

  n_dof_per_component = depends_on_space ? control_dofs.size() : 1;
  n_control_per_time_step = n_components * n_dof_per_component;

  n_control *= n_control_per_time_step;  
  control_old = std::vector<double>(n_control, 0.0);
  
  control_index_x.resize(control_dofs.size());
  control_index_y.resize(control_dofs.size());
}

template<int d>
double Time_BoundaryControlledOptimization<d>::compute_functional_and_derivative(
  unsigned int n, const double* x, double* grad)
{
  if(n != n_control)
  {
    ErrThrow("the given dimension ", n, " must match the dimension ", n_control,
             " of the control space");
  }
  ++n_calls;
  Output::print<1>("\n#########################################"
                   "\nCall to 'compute_functional_and_derivative' ", n_calls,
                   "\n#########################################"); 
  if(n_calls > 1)
  {
    double norm_diff_control = std::sqrt(
      std::inner_product(x, x+n, control_old.begin(), 0.0, std::plus<double>(),
                         [](double a, double b){ return (a-b) * (a-b); }));
    double norm_mean_control = std::sqrt(
      std::inner_product(x, x+n, control_old.begin(), 0.0, std::plus<double>(),
                         [](double a, double b){ return (a+b) * (a+b) / 4.; }));
    Output::print<1>("diff of control (in l2): ", norm_diff_control,
                  "  mean: ", norm_mean_control);
  }
  std::copy(x, x+n, control_old.begin());
  Output::print<2>("Control vector (", n, " dofs):");
  for(auto i = 0u; i < n_time_steps_; ++i)
  {
    for(auto j = 0u; j < n_dof_per_component; ++j)
    {
    Output::print<3>("Time step i=", i, " x[i]=", x[j+i*2*n_dof_per_component],
                  " y[i]=", x[j+i*2*n_dof_per_component+n_dof_per_component]);
    }
  }
  
  apply_control_and_solve(x);

  // Write all the solutions - for validation and debugging purposes
//  write_all_solutions();

  Output::print<5>("Resetting time to 0 after solve!");
  TDatabase::TimeDB->CURRENTTIME = 0;

  current_J_hat = compute_functional_in_time(0, n_time_steps_);
  if(grad != nullptr)
  {
    ++n_computation_derivative;
    solve_adjoint_equation();
    compute_derivative(x, grad);
    double norm_of_grad = std::sqrt(std::accumulate(
      grad, grad+n, 0., [](double a, double b){ return a + b*b; }));
    Output::print("norm of gradient: ", norm_of_grad);
//     for(auto i = 0u; i < n_control; ++i)
//     {
//       Output::print("computed_gradient[", i, "] = ", grad[i]);
//     }
  }
  
  // clear all primal solutions of current optimization iteration (VERY IMPORTANT LINE!)
  tnse_primal_solutions_.clear();
  
  if(n_calls == 1)
  {
    optimization_info.restart("optimization", current_J_hat);
  }
  optimization_info.print(n_calls, current_J_hat);
  return current_J_hat;
}

template<int d>
void Time_BoundaryControlledOptimization<d>::apply_control_and_solve(const double* x)
{
  // make sure Dirichlet rows are handled as usual (ones on the diagonal)
  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 0;
  
  Output::print<2>("primal solve");
  TimeDiscretization& tss = tnse_primal.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.current_time_ = 0.;
  TDatabase::TimeDB->CURRENTTIME = 0.;
  tss.set_time_disc_parameters();

  // it is necessary to re-set the initial conditions here to re-start correctly the
  // optimization iteration
  if(tnse_primal.get_db()["read_initial_solution"].is(true))
  {//initial solution is given
    std::string file = tnse_primal.get_db()["initial_solution_file"];
    Output::info("Initial Solution", "Reading initial solution from file ", file);
    BlockVector& sol_vector = tnse_primal.get_solution();
    sol_vector.read_from_file(file);
  }
  else
  {///interpolate initial condition from the example
    Output::info<5>("Initial Solution", "Interpolating initial solution from example.");
    for(int i = 0; i < d; ++i)
    {
      FEFunction * ui = tnse_primal.get_velocity_component(i);
      ui->Interpolate(tnse_primal.get_example().get_initial_cond(i));
      delete ui;
    }
  }

  // apply control x
  impose_control_in_rhs_and_sol(x, tss.current_step_);

  // save solution for later use (cost functional and adjoint problem)
  tnse_primal_solutions_.emplace_back(tnse_primal.get_solution(),tss.current_step_,
                                      tnse_primal.get_velocity(),tnse_primal.get_pressure());

  tnse_primal.assemble_initial_time();
  impose_control_in_rhs_and_sol(x, tss.current_step_); // to overwrite assembled rhs values
  tnse_primal.output();
  LoopInfo loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1;
  int linear_iteration=0;

  // solve the primal system, time loop
  while(!tss.reached_final_time_step())
  {
    tss.current_step_++;
    // set the time parameters
    tss.set_time_disc_parameters();
    double tau = tnse_primal.get_db()["time_step_length"];
    tss.current_time_ += tss.get_step_length();
    TDatabase::TimeDB->CURRENTTIME += tau;
    Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

    tnse_primal.assemble_matrices_rhs(0);
    impose_control_in_rhs_and_sol(x,tss.current_step_);

    LoopInfo loop_info("nonlinear");
    loop_info.print_time_every_step = true;
    loop_info.verbosity_threshold = 1;
    for(unsigned int i=0;; i++)
    {
      if(tnse_primal.stop_it(i))
      {
        loop_info.finish(i,tnse_primal.get_full_residual());
        linear_iteration +=i;
        loop_info_time.print(linear_iteration, tnse_primal.get_full_residual());
        break;
      }
      else
        loop_info.print(i, tnse_primal.get_full_residual());

      tnse_primal.solve();

      if(tnse_primal.imex_scheme())
        continue;

      tnse_primal.assemble_matrices_rhs(i+1);
      impose_control_in_rhs_and_sol(x,tss.current_step_);
    }
    // save solution for later use (cost functional and adjoint problem)
    tnse_primal_solutions_.emplace_back(tnse_primal.get_solution(),tss.current_step_,
                                        tnse_primal.get_velocity(),tnse_primal.get_pressure());

    tnse_primal.output(); // use this line to output(tss.current_step_) everything
  }
//  tnse_primal.output();  // use this line to output(n_calls) only the last step
}

template<int d>
void Time_BoundaryControlledOptimization<d>::impose_control_in_rhs_and_sol(const double* x,
                                                                        int current_time_step)
{
  auto& rhs = this->tnse_primal.get_rhs_from_time_disc();
  auto& sol = this->tnse_primal.get_solution();
  auto& old_sol = this->tnse_primal.get_old_solution();
  int length = rhs.length(0);

  //adjust control indices according to space- and time-dependency (see input parameters)
  set_control_indices(current_time_step); 

  for(auto i = 0u; i < control_dofs.size(); ++i)
  {
    rhs[control_dofs[i]] = db["control_in_x_direction"] ? x[control_index_x[i]] : 0.;
    sol[control_dofs[i]] = db["control_in_x_direction"] ? x[control_index_x[i]] : 0.;
    old_sol[control_dofs[i]] = db["control_in_x_direction"] ? x[control_index_x[i]] : 0.;

    rhs[control_dofs[i] + length] = db["control_in_y_direction"] ? x[control_index_y[i]]: 0.;
    sol[control_dofs[i] + length] = db["control_in_y_direction"] ? x[control_index_y[i]]: 0.;
    old_sol[control_dofs[i] + length] = db["control_in_y_direction"] ? x[control_index_y[i]]: 0.;

    Output::print("Control_dofs x_component: i=",i, " in time step=", current_time_step,
                  " equals:", rhs[control_dofs[i]]);
    Output::print("Control_dofs y_component: i=", i, " in time step=",
    current_time_step, " equals:", rhs[control_dofs[i] + length]);
   }
//    if(n_calls >=2)
//    {
//      Output::print("EXITING AT SECOND CALL TO FORWARD PROBLEM");
//      exit(0);
//    }
}

template<int d>
void Time_BoundaryControlledOptimization<d>::set_control_indices(int current_time_step) const
{
  int time_index = db["control_depends_on_time"] ? current_time_step*n_control_per_time_step : 0;
  
  for(auto i = 0u; i < control_dofs.size(); ++i)
  {
    control_index_x[i] = time_index;
    if(db["control_depends_on_space"])
      control_index_x[i] += i;

    if(db["control_in_x_direction"] && control_index_x[i] >= n_control)
          ErrThrow("ERROR", control_index_x[i], " ", n_control);

    if(db["control_in_y_direction"])
      control_index_y[i] = control_index_x[i] + n_dof_per_component;

    if(db["control_in_y_direction"] && control_index_y[i] >= n_control)
      ErrThrow("ERROR", control_index_y[i], " ", n_control);
   }
}

template<int d>
double Time_BoundaryControlledOptimization<d>::compute_functional_in_time(int t0, int t1) const
{
  double functional_value = 0;
  for (int i = t0; i < t1; ++i)
  {
    functional_value += compute_functional_at_t(i);
  }

  if(n_calls > 1)
  {
    Output::print("difference to previous functional ",
                  current_J_hat - functional_value);
  }
  return functional_value;
}

template<int d>
double Time_BoundaryControlledOptimization<d>::compute_functional_at_t(int time_step) const
{
  if (time_step > n_time_steps_ - 1)
    ErrThrow("The functional can not be computed at t = ", time_step, " > ", n_time_steps_-1);

  auto u = tnse_primal_solutions_.at(time_step).u_at_timestep_t_;
  auto u1 = u.GetComponent(0);
  auto u2 = u.GetComponent(1);
#ifdef __3D__
  auto u3 = u.GetComponent(2);
#endif
  //auto p = tnse_primal_solutions_.at(time_step).p_at_timestep_t_();

  // the id of the component on which the control is to be applied
  int comp = db["controlled_boundary_component"];
  double alpha = db["alpha_cost"];
  std::vector<double> cost_functional = db["cost_functional"];
  double functional_value = 0;
#ifdef __2D__
  auto l2_norm_on_boundary1 = u1->get_L2_norm_on_boundary(comp);
  auto l2_norm_on_boundary2 = u2->get_L2_norm_on_boundary(comp);
  auto l2_norm_on_boundary = l2_norm_on_boundary1*l2_norm_on_boundary1
   + l2_norm_on_boundary2*l2_norm_on_boundary2;

  auto compute_L2_norm_of_curl =
  [&](){
    if(db["restricted_curl_functional"])
    {
      std::vector<double> values(1);
      auto f = [](std::vector<double>& v, std::array<double, 8> e)
               {
                 double curl = (e[0] < 8. && e[0] > 4.) ? -e[6] + e[5] : 0.;
                 v[0] += curl;
               };
      u.get_functional_value(values, f);
      return values[0] * values[0];
    }
    else
    {
      auto div_curl = u.get_L2_norm_divergence_curl();
      return std::get<1>(div_curl)*std::get<1>(div_curl);
    }
  };
  auto compute_backward_facing_step =
  [&](){
    std::vector<double> values(1);
    auto f = [](std::vector<double>& v, std::array<double, 8> e)
             {
               double min_u1_0 = e[2] < 0. ? e[2]*e[2] : 0.; // min{u1,0}^2
               double max_u2_0 = e[3] > 0. ? e[3]*e[3] : 0.; // max{u2,0}^2
               v[0] += std::sqrt(min_u1_0 + max_u2_0);
             };
    u.get_functional_value(values, f);
    return values[0]*values[0];
  };
  auto compute_L2_norm_diff_stokes =
  [&](){
    const BlockVector& primal_sol = tnse_primal.get_solution();
    *this->stokes_fe_vector -= primal_sol;
    *this->stokes_fe_vector *= -1.;
    auto u1_diff = this->stokes_sol->GetComponent(0);
    auto u2_diff = this->stokes_sol->GetComponent(1);
    auto u1_diff_l2_norm = u1_diff->get_L2_norm();
    auto u2_diff_l2_norm = u2_diff->get_L2_norm();
    auto l2_norm_diff = u1_diff_l2_norm * u1_diff_l2_norm
                      + u2_diff_l2_norm * u2_diff_l2_norm;
    *this->stokes_fe_vector *= -1.;
    *this->stokes_fe_vector += primal_sol;
    delete u1_diff;
    delete u2_diff;
    return l2_norm_diff;
  };

  auto curl = compute_L2_norm_of_curl();
  auto min_max_bfs = compute_backward_facing_step();
  auto l2_norm_diff = compute_L2_norm_diff_stokes();
  functional_value = 0.5 * alpha * l2_norm_on_boundary;
  functional_value += 0.5 * curl * cost_functional[0];
  functional_value += 0.5 * min_max_bfs * cost_functional[1];
  functional_value += 0.5 * l2_norm_diff * cost_functional[2];
  Output::print("functional parts at time step i ", time_step, " : curl ", curl, ",   min_max ", min_max_bfs,
                ",   l2_norm_diff ", l2_norm_diff, ",   cost ",
                alpha * l2_norm_on_boundary);
  
  delete u1;
  delete u2;
#else
  ErrThrow("No 3D yet!");
#endif
  return functional_value;
}

template<int d>
void Time_BoundaryControlledOptimization<d>::solve_adjoint_equation()
{
  Output::print<2>("\n#########################################"
                   "\n#### STARTING SOLVE_ADJOINT_EQUATION ####"
                   "\n#########################################"); 
  // make sure Dirichlet rows are handled differently (1e30 on the diagonal)
  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
  
  Output::print<2>("adjoint solve");
  TimeDiscretization& tss = tnse_adjoint.get_time_stepping_scheme();
  tss.current_step_ = n_time_steps_; // last step first / backward in time
  tss.set_time_disc_parameters();

  {
    // the "initial" (t=T) condition is zero unless the objective has a term
    // involving u(T)
    Output::info<5>("Initial Solution", "Interpolating initial solution from example.");
    for(int i = 0; i < d; ++i)
    {
      FEFunction * ui = tnse_adjoint.get_velocity_component(i);
      ui->Interpolate(tnse_adjoint.get_example().get_initial_cond(i));
      delete ui;
    }
    FEFunction&  p = tnse_adjoint.get_pressure();
    p.Interpolate(tnse_adjoint.get_example().get_initial_cond(d));
  }

  auto final_forward_solution = tnse_primal_solutions_.at(
    tnse_primal.get_time_stepping_scheme().current_step_);
  tnse_adjoint.assemble_initial_time(final_forward_solution.u_at_timestep_t_,
                                     final_forward_solution.p_at_timestep_t_);
  tnse_adjoint.output();
  double end_time = 0.;
  TDatabase::TimeDB->CURRENTTIME = tnse_primal.get_time_stepping_scheme().get_end_time();
  LoopInfo loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1;
  int linear_iteration=0;

  // solve the primal system, time loop
  while(tss.current_time_ > end_time + 1e-10)
  {
    tss.current_step_--;
    // set the time parameters
    tss.set_time_disc_parameters();
    double tau = tnse_adjoint.get_db()["time_step_length"];
    TDatabase::TimeDB->CURRENTTIME -= tau;
    Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

    auto forward_solution = tnse_primal_solutions_.at(tss.current_step_);
    tnse_adjoint.assemble_matrices_rhs(forward_solution.u_at_timestep_t_,
                                       forward_solution.p_at_timestep_t_);

    tnse_adjoint.solve();
    // save solution for later use (cost functional and adjoint problem)
    //tnse_adjoint_solutions_.emplace_back(tnse_adjoint.get_solution(),tss.current_step_,
    //                                    tnse_adjoint.get_velocity(),tnse_adjoint.get_pressure());

    tnse_adjoint.output(); // use this line to output(tss.current_step_) everything
  }
}

template<int d>
void Time_BoundaryControlledOptimization<d>::compute_derivative(const double* x,
                                                        double* grad) const
{
//  double alpha = db["alpha_cost"];
//  auto & adjoint_solution = tnse_adjoint.get_solution();
//  auto adjoint_residual = tnse_adjoint.get_rhs(); // copy
//  auto & adjoint_mat = tnse_adjoint.get_matrix();
//  adjoint_mat.BlockMatrix::apply_scaled_add(adjoint_solution, adjoint_residual, -1.);
//  auto length = adjoint_solution.length(0);
//  auto n = control_dofs.size(); // number of dof per component
//  for(auto i = 0u; i < n; ++i)
//  {
//    grad[i] = alpha * x[i] - adjoint_residual[control_dofs[i]];
//    grad[i+n] = alpha * x[i+n] - adjoint_residual[control_dofs[i] + length];
////     Output::print("derivative x ", alpha * x[i], " ", alpha * x[i+n],
////                   " adjoint ", adjoint_residual[control_dofs[i]], " ",
////                   adjoint_residual[control_dofs[i] + length]);
//  }
//   for(auto i = 0u; i < adjoint_solution.length(); ++i)
//   {
//     Output::print("i ", i, " \t sol ", adjoint_solution[i], " \t rhs ", 
//                   tnse_adjoint.get_rhs()[i], " \t residual ",
//                   adjoint_residual[i]);
//   }
}

template<int d>
void Time_BoundaryControlledOptimization<d>::write_all_solutions()
{
  std::string name;
  if (tnse_primal_solutions_.size()!=n_time_steps_)
  {
    ErrThrow("The number of saved solutions must be equal to the number of "
        "time steps!");
  }
  else
  {
    for (unsigned int i = 0; i < tnse_primal_solutions_.size(); ++i)
    {
      Output::print(tnse_primal_solutions_.at(i).timestep_t_);
      name = "sol_vectors";
      name += std::to_string(i);
      tnse_primal_solutions_.at(i).vector_at_timestep_t_.write(name);
      tnse_primal_solutions_.at(i).u_at_timestep_t_.WriteSol(i,"./","sol_fefuncts");
    }
  }
//  exit(0);
}

#ifdef __3D__
template class Time_BoundaryControlledOptimization<3>;
#else
template class Time_BoundaryControlledOptimization<2>;
#endif

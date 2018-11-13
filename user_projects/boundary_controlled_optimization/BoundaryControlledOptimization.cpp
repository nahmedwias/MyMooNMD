#include "BoundaryControlledOptimization.hpp"
#include "LoopInfo.h"
#include "Database.h"
#include <algorithm>

ParameterDatabase BoundaryControlledOptimization::default_BCO_database()
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
  
  return db;
}


ParameterDatabase BoundaryControlledOptimization::get_primal_database(const ParameterDatabase& param_db)
{
  return param_db;
}

ParameterDatabase BoundaryControlledOptimization::get_adjoint_database(const ParameterDatabase& param_db)
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


BoundaryControlledOptimization::BoundaryControlledOptimization(
  const TDomain& domain, const ParameterDatabase& param_db)
 : db(default_BCO_database()), n_control(0), control_dofs(),
   nse_primal(domain, get_primal_database (param_db)), 
   nse_adjoint( nse_primal, get_adjoint_database(param_db)), 
   stokes_fe_vector(new BlockVector(nse_primal.get_solution())),
   stokes_sol(
     new TFEVectFunct2D(&nse_primal.get_velocity_space(), "stokes_sol", 
                        "", stokes_fe_vector->get_entries(), 
                        nse_primal.get_velocity_space().GetN_DegreesOfFreedom(), 
                        2)),
   optimization_info("optimization", true, true, 1), // 1 -> full verbosity
   nonlinear_info("nonlinear", true, true, 3),
   current_J_hat(std::numeric_limits<double>::infinity()), 
   control_old(), n_calls(0),
   n_computation_derivative(0)
{
  Output::print<5>("Creating the BoundaryControlledOptimization object");
  db.merge(param_db, false);
  
  // assemble linear (Stokes) parts
  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 0;
  nse_primal.assemble();
  // solve a Stokes system once, to store its solution (as 'desired state')
  nse_primal.solve();
  *stokes_fe_vector = nse_primal.get_solution(); // copy
  
  nse_primal.add_to_output(stokes_sol.get());
  
  // the id of the component on which the control is to be applied
  int controlled_boundary_component = db["controlled_boundary_component"];
  // find the degrees of freedom (dof) which are to be controlled
  auto& fe_space = nse_primal.get_velocity_space();
  auto collection = fe_space.GetCollection();
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
        auto *boundary_joint = dynamic_cast<TBoundEdge*>(joint);
        int boundary_id = boundary_joint->GetBoundComp()->GetID();
        if(boundary_id == controlled_boundary_component)
        {
          // found a joint which has is on the boundary where the control is at
          auto fe_descriptor = fe_space.get_fe(i_cell).GetFEDesc2D();
          int n_joint_dof = fe_descriptor->GetN_JointDOF();
          auto local_joint_dofs = fe_descriptor->GetJointDOF(i_joint);
          auto local_to_global_dof = fe_space.GetGlobalDOF(i_cell);
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
        auto *boundary_joint = dynamic_cast<TBoundEdge*>(joint);
        int boundary_id = boundary_joint->GetBoundComp()->GetID();
        if(boundary_id != controlled_boundary_component)
        {
          // found a joint which has is on an other boundary
          auto fe_descriptor = fe_space.get_fe(i_cell).GetFEDesc2D();
          int n_joint_dof = fe_descriptor->GetN_JointDOF();
          auto local_joint_dofs = fe_descriptor->GetJointDOF(i_joint);
          auto local_to_global_dof = fe_space.GetGlobalDOF(i_cell);
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
  // sort and remove dublicates
  std::sort(control_dofs.begin(), control_dofs.end());
  auto it = std::unique(control_dofs.begin(), control_dofs.end());
  control_dofs.resize(std::distance(control_dofs.begin(), it));
  control_dofs.shrink_to_fit();
  
//   for(auto e : control_dofs)
//   {
//     Output::print("control_dof ", e);
//   }
  n_control = 2 * control_dofs.size(); // space dimension 2
  
  control_old = std::vector<double>(n_control, 0.0);
  Output::print<3>("Created the BoundaryControlledOptimization object, ",
                   "n_control = ", n_control);
}


double BoundaryControlledOptimization::compute_functional_and_derivative(
  unsigned int n, const double* x, double* grad)
{
  if(n != n_control)
  {
    ErrThrow("the given dimension ", n, " must match the dimension ", n_control,
             " of the control space");
  }
  ++n_calls;
  Output::print<2>("Call to 'compute_functional_and_derivative' ", n_calls);
  if(n_calls > 1)
  {
    double norm_diff_control = std::sqrt(
      std::inner_product(x, x+n, control_old.begin(), 0.0, std::plus<double>(),
                         [](double a, double b){ return (a-b) * (a-b); }));
    double norm_mean_control = std::sqrt(
      std::inner_product(x, x+n, control_old.begin(), 0.0, std::plus<double>(),
                         [](double a, double b){ return (a+b) * (a+b) / 4.; }));
    Output::print("diff of control (in l2): ", norm_diff_control,
                  "  mean: ", norm_mean_control);
  }
  std::copy(x, x+n, control_old.begin());
  
  apply_control_and_solve(x);
  current_J_hat = compute_functional();
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
  
  
  if(n_calls == 1)
  {
    optimization_info.restart("optimization", current_J_hat);
  }
  optimization_info.print(n_calls, current_J_hat);
  return current_J_hat;
}


void BoundaryControlledOptimization::apply_control_and_solve(const double* x)
{
  // apply control x
  auto& rhs = nse_primal.get_rhs();
  int length = rhs.length(0);
  auto n_dof_per_component = control_dofs.size();
  for(auto i = 0u; i < n_dof_per_component; ++i)
  {
    rhs[control_dofs[i]] = x[i];
    rhs[control_dofs[i] + length] = x[i+n_dof_per_component];
  }
  // make sure Dirichlet rows are handled as usual (ones on the diagonal)
  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 0;
  
  Output::print<2>("primal solve");
  // solve the primal system, nonlinear loop
  nse_primal.reset_residuals(); // forget previous iterations
  nse_primal.assemble_nonlinear_term();
  // check if solution is already set, compute residuals
  if(!nse_primal.stopIt(0)) 
  {
    bool stokes = false;
    nonlinear_info.restart("nonlinear", nse_primal.getFullResidual());
    nonlinear_info.print(0, nse_primal.getFullResidual());
    for(unsigned int k = 1;; k++)
    {
      //Output::print(); // new line for a new nonlinear iteration
      nse_primal.solve();
      
      //no nonlinear iteration for Stokes problem
      if(stokes)
        break;
      
      nse_primal.assemble_nonlinear_term();
      
      if( nse_primal.stopIt(k))
      {
        nonlinear_info.finish(k, nse_primal.getFullResidual());
        break;
      }
      else
      {
        nonlinear_info.print(k, nse_primal.getFullResidual());
      }
    } // end for k
  }
  nse_primal.output(n_calls);
}


double BoundaryControlledOptimization::compute_functional() const
{
  auto u = nse_primal.get_velocity();
  auto u1 = u.GetComponent(0);
  auto u2 = u.GetComponent(1);
  //auto p = nse_primal.get_pressure();
  // the id of the component on which the control is to be applied
  int comp = db["controlled_boundary_component"];
  double alpha = db["alpha_cost"];
  std::vector<double> cost_functional = db["cost_functional"];
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
    const BlockVector& primal_sol = nse_primal.get_solution();
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
  double functional_value = 0.5 * alpha * l2_norm_on_boundary;
  functional_value += 0.5 * curl * cost_functional[0];
  functional_value += 0.5 * min_max_bfs * cost_functional[1];
  functional_value += 0.5 * l2_norm_diff * cost_functional[2];
  Output::print("functional parts: curl ", curl, ",   min_max ", min_max_bfs,
                ",   l2_norm_diff ", l2_norm_diff, ",   cost ",
                alpha * l2_norm_on_boundary);
  if(n_calls > 1)
  {
    Output::print("difference to previous functional ", 
                  current_J_hat - functional_value);
  }
  
  delete u1;
  delete u2;
  return functional_value;
}


void BoundaryControlledOptimization::solve_adjoint_equation()
{
  auto u = nse_primal.get_velocity();
  auto p = nse_primal.get_pressure();
  std::vector<double> cost_functional_weights = db["cost_functional"];
  bool restricted_curl_functional = db["restricted_curl_functional"];
  
  Output::print<2>("adjoint solve ", n_computation_derivative);
  nse_adjoint.assemble(u, p, *stokes_sol, cost_functional_weights,
                       restricted_curl_functional);
  nse_adjoint.solve();
  nse_adjoint.output(n_calls);
}


void BoundaryControlledOptimization::compute_derivative(const double* x, 
                                                        double* grad) const
{
  double alpha = db["alpha_cost"];
  auto & adjoint_solution = nse_adjoint.get_solution();
  auto adjoint_residual = nse_adjoint.get_rhs(); // copy
  auto & adjoint_mat = nse_adjoint.get_matrix();
  adjoint_mat.BlockMatrix::apply_scaled_add(adjoint_solution, adjoint_residual, -1.);
  auto length = adjoint_solution.length(0);
  auto n = control_dofs.size(); // number of dof per component
  for(auto i = 0u; i < n; ++i)
  {
    grad[i] = alpha * x[i] - adjoint_residual[control_dofs[i]];
    grad[i+n] = alpha * x[i+n] - adjoint_residual[control_dofs[i] + length];
//     Output::print("derivative x ", alpha * x[i], " ", alpha * x[i+n], 
//                   " adjoint ", adjoint_residual[control_dofs[i]], " ",
//                   adjoint_residual[control_dofs[i] + length]);
  }
//   for(auto i = 0u; i < adjoint_solution.length(); ++i)
//   {
//     Output::print("i ", i, " \t sol ", adjoint_solution[i], " \t rhs ", 
//                   nse_adjoint.get_rhs()[i], " \t residual ", 
//                   adjoint_residual[i]);
//   }
}


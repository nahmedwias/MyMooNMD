#include <Solver.h>
#include <IterativeMethod.h>
#include <Iteration_bicgstab.h>
#include <Iteration_cg.h>
#include <Iteration_cgs.h>
#include <Iteration_gmres.h>
#include <Iteration_jacobi.h>
#include <Iteration_multigrid.h>
#include <Iteration_richardson.h>
#include <Iteration_sor.h>
#include <Preconditioner.h>
#include <Saddle_point_preconditioner.h>
#include <DirectSolver.h>

template <class L, class V>
ParameterDatabase Solver<L, V>::default_solver_database()
{
  Output::print<3>("creating a default solver parameter database");
  ParameterDatabase db("default solver database");
  
  db.add("solver_type", std::string("direct"),
         "Determine which kind of solver should be used. This can be an "
         "iterative or a direct solver", {"direct", "iterative"});
  
  db.add("direct_solver_type", std::string("umfpack"),
         "Determine which type of direct solver should be used. All of them "
         "are implemented in external libraries. ",
         {"umfpack", "pardiso", "mumps"});
  
  db.add("iterative_solver_type", std::string("fgmres"),
         "Determine which type of iterative solver should be used.",
         {"jacobi", "sor", "ssor", "richardson", "cg", "cgs", "bi_cgstab", 
          "left_gmres", "right_gmres", "fgmres"});
  
  // the range is not so nice, maybe we have to store an integer interval with
  // min, max in the Parameter class 
  db.add<size_t>("max_n_iterations", 100,
                 "Maximum number of iterations of the iterative solver. This "
                 "is used as a stopping criterion of the iteration.",
                 {0, 1, 2, 3, 4, 5, 10, 100, 1000, 10000, 100000} );
  
  db.add<size_t>("min_n_iterations", 0,
                 "Minimum number of iterations of the iterative solver. This "
                 "enforces iterations even of other stopping criteria are "
                 "already satisfied.",
                 {0, 1, 2, 3, 4, 5, 10} );
  
  db.add("residual_tolerance", 1.0e-8,
         "The desired accuracy for the residual using an iterative solver. "
         "This is used as a stopping criterion of the iteration.",
         0.0, 1.0e10);
  
  db.add("residual_reduction", 0.0, 
         "The factor by which the residual should be reduced. This is used "
         "as a stopping criterion. Whenever the residual is smaller than "
         "the product of the initial residual and this parameter, the "
         "iteration is terminated. A value of 0.0 therefore effectively "
         "disables this stopping criterion.", 0.0, 1.0);
  
  db.add<size_t>("gmres_restart", 20, "The number of gmres iterations until "
                 "a restart is done. Larger numbers lead to more memory "
                 "consumption, smaller numbers typically mean more "
                 "iterations.", 1, 1e3);
  
  db.add("preconditioner", std::string("no_preconditioner"),
         "Determine the used preconditioner. Note that some of these are "
         "specific for some problem types.",
         {"no_preconditioner", "jacobi", "sor", "ssor", "multigrid", 
          "semi_implicit_method_for_pressure_linked_equations",
          "least_squares_commutator", "least_squares_commutator_boundary"});
  
  db.add("sor_omega", 1.5, "The overrelaxation parameter (typically called "
         "omega). This is only used for the (symmetric) successive "
         "overrelaxation method.", 0., 2.);
  
  db.add("saddle_point_preconditioner_direct_velocity_solve", true, 
         "During the application of a Saddle_point_preconditioner one has to "
         "solve a system involving only the velocity part of the matrix. Set "
         "this parameter to true if you want to solve this with a direct "
         "solver, otherwise some iterative scheme is used. Check out the class "
         "Saddle_point_preconditioner.");
  
  db.add("damping_factor", 1.0, "The damping in an iteration. A value of 1.0 "
         "means no damping while 0.0 would mean no progress. In general "
         "smaller values make iterations slower. This can still be necessary "
         "in cases where the iterations does not converge at all with larger "
         "values.", 0.0, 1.0);
  
  db.add("damping_factor_finest_grid", 1.0,
         "The damping of an iteration in case of a multigrid preconditioner. "
         "This only affects the update on the finest grid.", 0.0, 1.0);
  
  return db;
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
std::shared_ptr<Preconditioner<V>> get_preconditioner(
  std::string preconditioner_name, const L& matrix, const ParameterDatabase& db)
{
  if(preconditioner_name == "no_preconditioner")
  {
    return std::make_shared<NoPreconditioner<V>>();
  }
  else if(preconditioner_name == "jacobi")
  {
    return std::make_shared<Iteration_jacobi<L, V>>(matrix);
  }
  else if(preconditioner_name == "sor")
  {
    return std::make_shared<Iteration_sor<L, V>>(matrix, 0, db["sor_omega"]);
  }
  else if(preconditioner_name == "ssor")
  {
    return std::make_shared<Iteration_sor<L, V>>(matrix, 2, db["sor_omega"]);
  }
  else if(preconditioner_name == "least_squares_commutator")
  {
    return std::make_shared<Saddle_point_preconditioner>(
      matrix, Saddle_point_preconditioner::type::lsc,
      db["saddle_point_preconditioner_direct_velocity_solve"]);
  }
  else if(preconditioner_name == "least_squares_commutator_boundary")
  {
    return std::make_shared<Saddle_point_preconditioner>(
      matrix, Saddle_point_preconditioner::type::bd_lsc,
      db["saddle_point_preconditioner_direct_velocity_solve"]);
  }
  else if(preconditioner_name == 
          "semi_implicit_method_for_pressure_linked_equations")
  {
    return std::make_shared<Saddle_point_preconditioner>(
      matrix, Saddle_point_preconditioner::type::simple,
      db["saddle_point_preconditioner_direct_velocity_solve"]);
  }
  else
  {
    ErrThrow("preconditioner ", preconditioner_name, " not yet implemented");
  }
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
std::shared_ptr<IterativeMethod<L, V>> get_iterative_method(
  std::string iterative_solver_type, const ParameterDatabase& db,
  const L& matrix, std::shared_ptr<Preconditioner<V>> p)
{
  std::shared_ptr<IterativeMethod<L, V>> ret;
  if(iterative_solver_type == "jacobi")
  {
    ret = std::make_shared<Iteration_jacobi<L, V>>(matrix);
  }
  else if(iterative_solver_type == "sor")
  {
    ret = std::make_shared<Iteration_sor<L, V>>(matrix, 0, db["sor_omega"]);
  }
  else if(iterative_solver_type == "ssor")
  {
    ret = std::make_shared<Iteration_sor<L, V>>(matrix, 2, db["sor_omega"]);
  }
  else if(iterative_solver_type == "richardson")
  {
    ret = std::make_shared<Iteration_richardson<L, V>>(p);
  }
  else if(iterative_solver_type == "cg")
  {
    ret = std::make_shared<Iteration_cg<L, V>>(p);
  }
  else if(iterative_solver_type == "cgs")
  {
    ret = std::make_shared<Iteration_cgs<L, V>>(p);
  }
  else if(iterative_solver_type == "bi_cgstab")
  {
    ret = std::make_shared<Iteration_bicgstab<L, V>>(p);
  }
  else if(iterative_solver_type == "left_gmres")
  {
    ret = std::make_shared<Iteration_gmres<L, V>>(p, gmres_type::left);
  }
  else if(iterative_solver_type == "right_gmres")
  {
    ret = std::make_shared<Iteration_gmres<L, V>>(p, gmres_type::right);
  }
  else if(iterative_solver_type == "fgmres")
  {
    ret = std::make_shared<Iteration_gmres<L, V>>(p, gmres_type::flexible);
  }
  else
    ErrThrow("unknown iterative_solver_type: ", iterative_solver_type);
  
  return ret;
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
Solver<L, V>::Solver(const ParameterDatabase& param_db)
 : db(default_solver_database()), linear_operator(nullptr), direct_solver(), 
   iterative_method(), preconditioner(), multigrid(nullptr)
{
  this->db.merge(param_db, false);
  if(this->is_using_multigrid())
    multigrid = std::make_shared<Multigrid>(param_db);
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Solver<L, V>::update_matrix(const L& matrix)
{
  this->linear_operator = &matrix;
  if(db["solver_type"].is("direct")) // direct solver
  {
    DirectSolver::DirectSolverTypes t;
    if(db["direct_solver_type"].is("umfpack"))
      t = DirectSolver::DirectSolverTypes::umfpack;
    else if(db["direct_solver_type"].is("pardiso"))
      t = DirectSolver::DirectSolverTypes::pardiso;
    else //if(db["direct_solver_type"].is("mumps"))
      ErrThrow("currently the DirectSolver class only supports umfack and "
               "pardiso");
    this->direct_solver.reset(new DirectSolver(matrix, t));
  }
  else
  {
    // iterative solver
    std::string ist = db["iterative_solver_type"];
    std::string prec_name = db["preconditioner"];
    size_t max_it = db["max_n_iterations"];
    size_t min_it = db["min_n_iterations"];
    double tol = db["residual_tolerance"];
    double reduc = db["residual_reduction"];
    size_t restart = db["gmres_restart"]; // only for gmres
    double damping = db["damping_factor"];
    
    // find out if a preconditioner object already exists and if it is of type 
    // Saddle_point_preconditioner
    bool is_saddle_point_preconditioner = this->preconditioner 
      && (prec_name == "least_squares_commutator" 
          || prec_name == "least_squares_commutator_boundary"
          || prec_name == "semi_implicit_method_for_pressure_linked_equations");
    
    if(!is_saddle_point_preconditioner)
    {
      // create a new preconditioner
      if(this->is_using_multigrid())
      {
        this->preconditioner = std::make_shared<Iteration_multigrid<L, V>>(
          this->multigrid);
        this->multigrid->update();
      }
      else
      {
        this->preconditioner = get_preconditioner<L, V>(prec_name, matrix, 
                                                        this->db);
      }
    }
    else
    {
      // only update the matrix
      Saddle_point_preconditioner* spp = 
        dynamic_cast<Saddle_point_preconditioner*>(this->preconditioner.get());
      if(spp != nullptr)
        spp->update(matrix);
    }
    this->iterative_method = get_iterative_method<L, V>(ist, this->db, matrix, 
                                                        this->preconditioner);
    this->iterative_method->set_stopping_parameters(max_it, min_it, tol, reduc, 
                                                    2., damping, restart);
  }
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Solver<L, V>::solve(const V& rhs, V& solution)
{
  if(this->linear_operator == nullptr)
    ErrThrow("in order to use Solver::solve(rhs, solution), you have to call "
             "Solver::update_matrix first. Otherwise it's not clear which "
             "system is supposed to be solved");
  if(db["solver_type"].is("direct"))
  {
    this->direct_solver->solve(rhs, solution);
  }
  else
  {
    this->iterative_method->iterate(*this->linear_operator, rhs, solution);
  }
  //compute the residual by hand again.
  //V r(rhs);
  //linear_operator->apply_scaled_add(solution, r, -1.);
  //Output::print<2>("computed residual in Solver class: ", r.norm());
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Solver<L, V>::solve(const L& matrix, const V& rhs, V& solution)
{
  this->update_matrix(matrix);
  this->solve(rhs, solution);
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
const ParameterDatabase& Solver<LinearOperator, Vector>::get_db() const
{
  return this->db;
}
/* ************************************************************************** */

template <class L, class V>
bool Solver<L, V>::is_using_multigrid() const
{
  return db["solver_type"].is("iterative") 
      && db["preconditioner"].is("multigrid");
}

/* ************************************************************************** */
template <class L, class V>
std::shared_ptr<Multigrid> Solver<L, V>::get_multigrid()
{
  if(this->is_using_multigrid())
    return this->multigrid;
  ErrThrow("unable to return a multigrid object, because this solver object is "
           "not set to use multigrid");
}

/* ************************************************************************** */
template <class L, class V>
std::shared_ptr<const Multigrid> Solver<L, V>::get_multigrid() const
{
  if(this->is_using_multigrid())
    return this->multigrid;
  ErrThrow("unable to return a multigrid object, because this solver object is "
           "not set to use multigrid");
}

/* ************************************************************************** */

// explicit instantiations
template class Solver<BlockMatrix, BlockVector>;
template class Solver<BlockFEMatrix, BlockVector>;

#include <Solver.h>
#include <Iteration_bicgstab.h>
#include <Iteration_cg.h>
#include <Iteration_cgs.h>
#include <Iteration_gmres.h>
#include <Iteration_jacobi.h>
#include <Iteration_richardson.h>

ParameterDatabase get_default_solver_parameters(size_t solver_type)
{
  Output::print<3>("creating a default solver parameter database");
  ParameterDatabase db("default solver database");
  
  db.add("solver_type", solver_type,
         "Determine which kind of solver should be used. This can be an "
         "iterative (1) or a direct (2) solver", (size_t)1, (size_t)2);
  
  db.add("direct_solver_type", 1,
         "Determine which type of direct solver should be used. All of them "
         "are implemented in external libraries. The values have the "
         "following meanings: "
         "1 - umfpack, "
         "2 - pardiso, "
         "3 - mumps.",
         1, 3);
  
  db.add<size_t>("iterative_solver_type", 1,
                 "Determine which type of iterative solver should be used. "
                 "The values have the following meanings: "
                 "1 - (weighted) Jacobi iteration, "
                 "2 - successive over-relaxation iteration (sor), "
                 "3 - symmetric successive over-relaxation iteration (ssor), "
                 "4 - Richardson iteration, "
                 "5 - conjugate gradients, "
                 "6 - conjugate gradients squared, "
                 "7 - Biconjugate Gradient Stabilized, "
                 "8 - left generalized minimal residuals (gmres), "
                 "9 - right generalized minimal residuals (gmres), "
                 "10- flexible (right) generalized minimal residuals(gmres).",
                 1, 10);
  
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
  
  // more than 20 multigrid meshes is probably an error
  db.add("n_multigrid_levels", (size_t)3,
         "The number of different multigrid meshes.",
         (size_t)1, (size_t)20);
  
  db.add("preconditioner", std::string("no_preconditioner"),
         "Determine the used preconditioner. Note that some of these are "
         "specific for some problem types.",
         {"no_preconditioner", "jacobi", "gauss_seidel", "multigrid", 
           "least_square_commutator", "least_square_commutator_boundary"});
  
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
template <class LinearOperator, class Vector>
Solver<LinearOperator, Vector>::Solver(const ParameterDatabase& param_db)
 : db("empty database to be replaced"), direct_solver()
{
  size_t solver_type;
  try
  {
    solver_type = param_db["solver_type"];
  }
  catch(...)
  {
    solver_type = 2;
  }
  this->db = get_default_solver_parameters(solver_type);
  this->db.merge(param_db, false);
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
void Solver<LinearOperator, Vector>::update_matrix(const LinearOperator& matrix)
{
  if(db["solver_type"].is(2)) // direct solver
  {
    this->direct_solver.reset(
      new DirectSolver(matrix, DirectSolver::DirectSolverTypes::umfpack));
  }
  //else // we just ignore this call, nothing is done
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
void Solver<LinearOperator, Vector>::solve(const Vector& rhs, Vector& solution)
{
  if(!db["solver_type"].is(2))
    ErrThrow("Calling Solver::solve without the matrix as an argument only "
             "works for direct solvers. In such a case the method "
             "Solver::update_matrix has to be called in advance. For iterative "
             "solvers you should call the method Solver::solve with a "
             "BlockMatrix (and right hand side and solution).");
  if(this->direct_solver)
    this->direct_solver->solve(rhs, solution);
  else
    ErrThrow("Calling Solver::solve(rhs, solution) with a direct solver "
             "requires a preceeding call to Solver::update_matrix(matrix).");
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::shared_ptr<Preconditioner<Vector>> get_preconditioner(
  std::string preconditioner_name, const LinearOperator& matrix)
{
  if(preconditioner_name == "no_preconditioner")
  {
    return std::make_shared<NoPreconditioner<Vector>>();
  }
  else if(preconditioner_name == "jacobi")
  {
    return std::make_shared<Iteration_jacobi<LinearOperator, Vector>>(matrix);
  }
  else
  {
    ErrThrow("preconditioner ", preconditioner_name, " not yet implemented");
  }
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::shared_ptr<IterativeMethod<LinearOperator, Vector>> get_iterative_method(
  size_t iterative_solver_type, const LinearOperator& matrix,
  std::shared_ptr<Preconditioner<Vector>> p)
{
  switch(iterative_solver_type)
  {
    case 1:
      return std::make_shared<Iteration_jacobi<LinearOperator, Vector>>(matrix);
      break;
    case 4:
      return std::make_shared<Iteration_richardson<LinearOperator, Vector>>(p);
      break;
    case 5:
      return std::make_shared<Iteration_cg<LinearOperator, Vector>>(p);
      break;
    case 6:
      return std::make_shared<Iteration_cgs<LinearOperator, Vector>>(p);
      break;
    case 7:
      return std::make_shared<Iteration_bicgstab<LinearOperator, Vector>>(p);
      break;
    case 8:
      return std::make_shared<Iteration_gmres<LinearOperator, Vector>>(
        p, gmres_type::left);
      break;
    case 9:
      return std::make_shared<Iteration_gmres<LinearOperator, Vector>>(
        p, gmres_type::right);
      break;
    case 10:
      return std::make_shared<Iteration_gmres<LinearOperator, Vector>>(
        p, gmres_type::flexible);
      break;
    default:
      ErrThrow("unknown iterative_solver_type: ", iterative_solver_type);
      break;
  }
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
void Solver<LinearOperator, Vector>::solve(const LinearOperator& matrix,
                                           const Vector& rhs, Vector& solution)
{
  if(db["solver_type"].is(2)) // direct solver
  {
    this->update_matrix(matrix);
    this->solve(rhs, solution);
    return;
  }
  //else // iterative solver
  size_t ist = db["iterative_solver_type"];
  std::string prec_name = db["preconditioner"];
  size_t max_it = db["max_n_iterations"];
  size_t min_it = db["min_n_iterations"];
  double tol = db["residual_tolerance"];
  double reduc = db["residual_reduction"];
  size_t restart = db["gmres_restart"]; // only for gmres
  double damping = db["damping_factor"];
  
  auto prec = get_preconditioner<LinearOperator, Vector>(prec_name, matrix);
  auto it_method = get_iterative_method<LinearOperator, Vector>(ist, matrix,
                                                                prec);
  it_method->set_stopping_parameters(max_it, min_it, tol, reduc, 2., damping,
                                     restart);
  auto n_it_residual = it_method->iterate(matrix, rhs, solution);
  Output::print<2>(it_method->get_name(), " iterations: ", n_it_residual.first,
                   ", residual: ", n_it_residual.second);
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
const ParameterDatabase& Solver<LinearOperator, Vector>::get_db()
{
  return this->db;
}
/* ************************************************************************** */

// explicit instantiations
template class Solver<BlockMatrix, BlockVector>;
template class Solver<BlockFEMatrix, BlockVector>;

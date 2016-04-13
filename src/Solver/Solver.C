#include <Solver.h>

ParameterDatabase get_default_solver_parameters(size_t solver_type)
{
  Output::print<3>("creating a default solver parameter database");
  ParameterDatabase db("default solver database");
  
  db.add("solver_type", solver_type, 
         "Determine which kind of solver should be used. This can be an "
         "iterative (1) or a direct (2) solver", (size_t)1, (size_t)2);
  
  if(solver_type == 2)
  {
    
  }
  else
  {
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
           "means no daming while 0.0 would mean no progress. In general "
           "smaller values make iterations slower. This can still be necessary "
           "in cases where the iterations does not converge at all with larger "
           "values.", 0.0, 1.0);
    
    db.add("damping_factor_finest_grid", 1.0,
           "The damping of an iteration in case of a multigrid preconditioner. "
           "This only affects the update on the finest grid.", 0.0, 1.0);
  }
  
  return db;
}
/* ************************************************************************** */
Solver::Solver(const ParameterDatabase& param_db)
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
void Solver::update_matrix(const BlockMatrix& matrix)
{
  if(!db["solver_type"].is(2))
    ErrThrow("Solver class only implemented for direct solver currently");
  this->direct_solver.reset(
    new DirectSolver(matrix, DirectSolver::DirectSolverTypes::umfpack));
}

/* ************************************************************************** */
void Solver::solve(const BlockVector& rhs, BlockVector& solution)
{
  if(!db["solver_type"].is(2))
    ErrThrow("Solver class only implemented for direct solver currently");
  
  this->direct_solver->solve(rhs, solution);
}

/* ************************************************************************** */
const ParameterDatabase& Solver::get_db()
{
  return this->db;
}

/* ************************************************************************** */
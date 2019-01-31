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
#include <Preconditioner_vanka.h>
#include <Saddle_point_preconditioner.h>
#include <DirectSolver.h>
#include <PETScSolver.h>
#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

template <class L, class V>
ParameterDatabase Solver<L, V>::default_solver_database()
{
  Output::print<5>("creating a default solver parameter database");
  ParameterDatabase db("default solver database");
  
  db.add("solver_type", "direct",
         "Determine which kind of solver should be used. Set it to direct to "
         "call a direct solver. Set to iterative to call a ParMooN internal "
         "iterative solver. Furthermore it is possible to use the external "
         "library PETSc, see the class PETScSolver. This also includes direct "
         "and iterative solvers.",
         {"direct", "iterative", "petsc"});
  
  db.add("petsc_arguments", "",
		 "Setting PETSc arguments as if used from command line. "
		 "use -pc_type <pc_arg> for linear solver method "
		 "<pc_arg> = { lu, ilu, hypre, mg, ... }"
		 "(see https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html) "
		 "use -pc_factor_mat_solver_package <package> for the package that should be used "
		 "for the linear solver method "
		 "<package> = { umfpack, mumps, superlu, ... }"
		 "(see https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatSolverPackage.html) "
		 "for hypre use -pc_hypre_type <package>"
		 "(see https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCHYPRE.html)",
		 {""});

  db.add("direct_solver_type", "umfpack",
         "Determine which type of direct solver should be used. All of them "
         "are implemented in external libraries. ",
         {"umfpack", "pardiso", "mumps"});
  
  db.add("iterative_solver_type", "fgmres",
         "Determine which type of iterative solver should be used.",
         {"jacobi", "sor", "ssor", "richardson", "cg", "cgs", "bi_cgstab", 
          "left_gmres", "right_gmres", "fgmres"});
  
  // the range is not so nice, maybe we have to store an integer interval with
  // min, max in the Parameter class 
  db.add<size_t>("max_n_iterations", 100,
                 "Maximum number of iterations of the iterative solver. This "
                 "is used as a stopping criterion of the iteration.",
                  0, 100000 );
  
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
 
  db.add("divergence_factor", 1.e10,
	 "If the norm of the residual increases by a factor that is larger, " 
	 "the iteration is stopped and the execution of the code will be" 
	 "stopped, too.", 1.0, 1.e36);
 
  db.add<size_t>("gmres_restart", 20, "The number of gmres iterations until "
                 "a restart is done. Larger numbers lead to more memory "
                 "consumption, smaller numbers typically mean more "
                 "iterations.", 1, 1e3);
  
  db.add("preconditioner", "no_preconditioner",
         "Determine the used preconditioner. Note that some of these are "
         "specific for some problem types.",
         {"no_preconditioner", "jacobi", "sor", "ssor", "multigrid", 
          "semi_implicit_method_for_pressure_linked_equations",
          "least_squares_commutator", "least_squares_commutator_boundary",
          "vanka_cell", "vanka_nodal", "vanka_cell_jacobi",
					"augmented_Lagrangian_based",	"modified_augmented_Lagrangian_based"}); //TODO maybe these vanka preconditioner types should be controlled by another db parameter
  
  db.add("sor_omega", 1.5, "The overrelaxation parameter (typically called "
         "omega). This is only used for the (symmetric) successive "
         "overrelaxation method.", 0., 2.);
  
  db.add("damping_factor", 1.0, "The damping in an iteration. A value of 1.0 "
         "means no damping while 0.0 would mean no progress. In general "
         "smaller values make iterations slower. This can still be necessary "
         "in cases where the iterations does not converge at all with larger "
         "values.", 0.0, 1.0);

  db.add("gamma", 1.0, "This is the augmentation factor (scalar) for the augmented Lagrangian based preconditioner."
         "0.0 means no augmentation."
         "values.", -1.0e10, 1.0e10);
  

  return db;
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
std::shared_ptr<Preconditioner<V>> get_preconditioner(
  const std::string& preconditioner_name, const L& matrix,
  const ParameterDatabase& db)
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
      matrix, Saddle_point_preconditioner::type::lsc, db);
  }
  else if(preconditioner_name == "least_squares_commutator_boundary")
  {
    return std::make_shared<Saddle_point_preconditioner>(
      matrix, Saddle_point_preconditioner::type::bd_lsc, db);
  }
  else if(preconditioner_name == 
          "semi_implicit_method_for_pressure_linked_equations")
  {
    return std::make_shared<Saddle_point_preconditioner>(
      matrix, Saddle_point_preconditioner::type::simple, db);
  }
  else if(preconditioner_name == "vanka_cell")
  {
    return std::make_shared<Preconditioner_vanka<V>>(
        matrix, VankaType::CELL, 0.8);
  }
  else if(preconditioner_name == "vanka_nodal")
  {
    return std::make_shared<Preconditioner_vanka<V>>(
        matrix, VankaType::NODAL, 0.8);
  }
  else if(preconditioner_name == "vanka_cell_jacobi")
  {
    return std::make_shared<Preconditioner_vanka<V>>(
        matrix, VankaType::CELL_JACOBI, 0.8);
  }
  else if(preconditioner_name == "augmented_Lagrangian_based")
  {
      return std::make_shared<Saddle_point_preconditioner>(
      matrix, Saddle_point_preconditioner::type::AL, db);
  }
  else if(preconditioner_name == "modified_augmented_Lagrangian_based")
   {
       return std::make_shared<Saddle_point_preconditioner>(
       matrix, Saddle_point_preconditioner::type::mod_AL, db);
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
  const std::string& iterative_solver_type, const ParameterDatabase& db,
  const L& matrix, const std::shared_ptr<Preconditioner<V>>& p)
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
 : db(default_solver_database()), linear_operator(nullptr),
   direct_solver(nullptr), iterative_method(nullptr), preconditioner(nullptr),
   multigrid(nullptr), petsc_solver(nullptr)
{

this->db.merge(param_db, false, true);

 if(this->is_using_multigrid())
    multigrid = std::make_shared<Multigrid>(param_db);
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Solver<L, V>::update_matrix(const L& matrix)
{
  if(db["solver_type"].is("direct")) // direct solver
  {
    DirectSolver::DirectSolverTypes t;
    if(db["direct_solver_type"].is("umfpack"))
      t = DirectSolver::DirectSolverTypes::umfpack;
    else if(db["direct_solver_type"].is("pardiso"))
      t = DirectSolver::DirectSolverTypes::pardiso;
    else //if(db["direct_solver_type"].is("mumps"))
      ErrThrow("currently the DirectSolver class only supports umfpack and "
               "pardiso");
    this->direct_solver.reset(new DirectSolver(matrix, t));
  }
  else if(db["solver_type"].is("iterative")) // ParMooN internal iterative 
  {
    // iterative solver
    bool create_new = this->linear_operator == nullptr;
    if(create_new)
    {
      std::string ist = db["iterative_solver_type"];
      std::string prec_name = db["preconditioner"];
      size_t max_it = db["max_n_iterations"];
      size_t min_it = db["min_n_iterations"];
      double tol = db["residual_tolerance"];
      double reduc = db["residual_reduction"];
      size_t restart = db["gmres_restart"]; // only for gmres
      double damping = db["damping_factor"];
      double divergence_factor = db["divergence_factor"];

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
      this->iterative_method = get_iterative_method<L, V>(ist, this->db, 
                                                          matrix, 
                                                          this->preconditioner);
      this->iterative_method->set_stopping_parameters(max_it, min_it, tol, 
                                                      reduc, divergence_factor, damping, 
                                                      restart);
   }
    else
    {
      this->preconditioner->update();
      this->iterative_method->update(matrix);
    }
  }
  else if(db["solver_type"].is("petsc"))
  {
    this->petsc_solver = std::make_shared<PETScSolver>(matrix, db);
  }
  else
  {
    ErrThrow("unknown solver type ", db["solver_type"]);
  }
  
  // set the linear operator (needed for iterative solvers)
  this->linear_operator = &matrix;
}


/* NEW LB FOR AUGM LAGR PREC ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Solver<L, V>::solve_augmented(const V& rhs, V& solution)
{
  // lambda function to compute the residual, you can use it e.g. during debug
  // note that this changes the consistency level of the solution to 2
  auto compute_residual = [&]()
  {
    //compute the absolute residual by hand again.
#ifdef _MPI
    // update to consistency level 2, needed to properly call apply_scaled_add
    auto comms = linear_operator->get_communicators();
    for (size_t bl = 0; bl < comms.size() ; ++bl)
    {
      comms[bl]->consistency_update(solution.block(bl), 2);
    }
#endif
    V r(rhs);
    linear_operator->apply_scaled_add(solution, r, -1.);
#ifndef _MPI
    Output::info<4>("Iterative solver", "Absolute residual in Solver class, ",
                    setprecision(16), r.norm());
#elif _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    double norm = r.norm(comms);
    if(my_rank == 0)
      Output::info<4>("Iterative solver", "Absolute residual in Solver class, ",
                      setprecision(16), norm);
#endif
  };
  (void)compute_residual; // silence the unused variable warning when not used

  // compute_residual();

  if(this->linear_operator == nullptr)
	  ErrThrow("in order to use Solver::solve(rhs, solution), you have to call "
			  "Solver::update_matrix first. Otherwise it's not clear which "
			  "system is supposed to be solved");
  if(db["solver_type"].is("direct"))
  {
	  this->direct_solver->solve(rhs, solution);

	 /* //New LB 11.07.18 start
	  // augmented Lagrangian stuff:
	  		 		  auto p = std::dynamic_pointer_cast <  Saddle_point_preconditioner > ( this->preconditioner );
	  		      std::shared_ptr<IterativeMethod<BlockMatrix, BlockVector>> iterative_method = get_iterative_method<BlockMatrix, BlockVector>( db["iterative_solver_type"] ,
	  				  this->db,  *this->linear_operator,  this->preconditioner);
	  		  //use augmented matrix and rhs in case of AL
	  		  iterative_method->iterate(p->get_augmented_matrix(), p->get_augmented_blockvector(rhs), solution);
	   //New LB 11.07.18 end
	 */
  }
  else if(db["solver_type"].is("iterative")) // ParMooN internal iterative
  {
#ifndef _MPI
    // augmented Lagrangian stuff:
    auto p = std::dynamic_pointer_cast<Saddle_point_preconditioner>(
      this->preconditioner);
    std::shared_ptr<IterativeMethod<BlockMatrix, BlockVector>> 
    iterative_method = get_iterative_method<BlockMatrix, BlockVector>(
      db["iterative_solver_type"], this->db, *this->linear_operator,
      this->preconditioner);
    //use augmented matrix and rhs in case of AL
    iterative_method->iterate(p->get_augmented_matrix(),
                              p->get_augmented_blockvector(rhs), solution);
    this->iterative_method->iterate(*this->linear_operator, rhs, solution);
#else
    ErrThrow("no augmented Lagrangian implemented in parallel mode yet.");
#endif
  }
  else if(db["solver_type"].is("petsc"))
  {
	  this->petsc_solver->solve(rhs, solution);
  }
  else
  {
	  ErrThrow("unknown solver type ", db["solver_type"]);
  }
  compute_residual();
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Solver<L, V>::solve(const V& rhs, V& solution)
{
  // lambda function to compute the residual, you can use it e.g. during debug
  // note that this changes the consistency level of the solution to 2
  auto compute_residual = [&]()
  {
    //compute the absolute residual by hand again.
#ifdef _MPI
    // update to consistency level 2, needed to properly call apply_scaled_add
    auto comms = linear_operator->get_communicators();
    for (size_t bl = 0; bl < comms.size() ; ++bl)
    {
      comms[bl]->consistency_update(solution.block(bl), 2);
    }
#endif
    V r(rhs);
    linear_operator->apply_scaled_add(solution, r, -1.);
#ifndef _MPI
    Output::info<4>("Iterative solver", "Absolute residual in Solver class, ",
                    setprecision(16), r.norm());
#elif _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    double norm = r.norm(comms);
    if(my_rank == 0)
      Output::info<4>("Iterative solver", "Absolute residual in Solver class, ",
                      setprecision(16), norm);
#endif
  };
  (void)compute_residual; // silence the unused variable warning when not used

  // compute_residual();

  if(this->linear_operator == nullptr)
	  ErrThrow("in order to use Solver::solve(rhs, solution), you have to call "
			  "Solver::update_matrix first. Otherwise it's not clear which "
			  "system is supposed to be solved");
  if(db["solver_type"].is("direct"))
  {
	  this->direct_solver->solve(rhs, solution);
  }
  else if(db["solver_type"].is("iterative")) // ParMooN internal iterative 
  {
  	if (db["preconditioner"].is("augmented_Lagrangian_based"))
	  {
#ifndef _MPI
		 		  auto p = std::dynamic_pointer_cast <  Saddle_point_preconditioner > ( this->preconditioner );
		  //this->iterative_method->iterate(p->get_augmented_matrix(), p->get_augmented_blockvector(rhs), solution);
		 		  std::shared_ptr<IterativeMethod<BlockMatrix, BlockVector>> iterative_method = get_iterative_method<BlockMatrix, BlockVector>( db["iterative_solver_type"] ,
		 		  		this->db,  *this->linear_operator,  this->preconditioner);
		 		  //use augmented matrix and rhs in case of AL
		 		  iterative_method->iterate(p->get_augmented_matrix(), p->get_augmented_blockvector(rhs), solution);

		 		  this->iterative_method->iterate(*this->linear_operator, rhs, solution);
#else
          ErrThrow("augmented Lagrangian not available in parallel mode yet.")
#endif // not MPI
	  }
  	else
  	{
  		this->iterative_method->iterate(*this->linear_operator, rhs, solution);
  	}

  }
  else if(db["solver_type"].is("petsc"))
  {
	  this->petsc_solver->solve(rhs, solution);
  }
  else
  {
	  ErrThrow("unknown solver type ", db["solver_type"]);
  }
  compute_residual();
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

template <class L, class V>
const LoopInfo& Solver<L, V>::get_solver_loop_info() const
{
    if(!iterative_method)
      throw std::runtime_error("Non-iterative solver, this holds no loop info object!");

    return iterative_method->get_loop_info();
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
template class Solver<BlockFEMatrix, BlockVector>;
// In MPI case we are so dependent on the connection of Matrix and FESpace, that
// it does not make sense to instantiate the function for BlockMatrix.
#ifndef _MPI
template class Solver<BlockMatrix, BlockVector>;
#endif


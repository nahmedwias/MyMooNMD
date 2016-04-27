#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <ParameterDatabase.h>
#include <DirectSolver.h>
#include <IterativeMethod.h>
#include <Preconditioner.h>
#include <memory>

/** @brief Solve a linear system
 * 
 * This class can solve linear systems of the \f$Ax=b\f$. How exactly it is 
 * solved is controlled entirely through the ParameterDatabase object given to
 * the constructor. After construction there is no way of changing the behavior
 * of the class. Typically you only need to call `solve` with a matrix, right
 * hand side and a solution vector. 
 * 
 * Currently this class is instantiated for `BlockMatrix, BlockVector` and
 * `BlockFEMatrix, BlockVector`.
 */
template <class LinearOperator, class Vector>
class Solver
{
  public:
    /// @brief create a solver object
    ///
    /// All solver related parameters are set to default values except those
    /// which also exist in the given ParameterDatabase. 
    /// 
    /// If you want to solve a system only once, you can use the method `solve`
    /// which has the linear operator as well as the right hand side and 
    /// solution as arguments. In case of a direct solver you can use the 
    /// already computed factorization a second time through calling `solve` 
    /// with only the right hand side and the solution as arguments.
    Solver(const ParameterDatabase& param_db);
    
    /// @brief update the solver and preconditioner objects
    ///
    /// This methods prepares the members. It either creates new ones or calls
    /// an appropriate update method on the objects. For example a direct_solver
    /// object will be created (deleting the old one, if existing) in case a 
    /// direct solver is used at all. A Saddle_point_preconditioner on the other
    /// hand only needs to be created once, and then update is enough. This 
    /// method does this.
    void update_matrix(const LinearOperator& matrix);
    
    /// @brief solve after calling `Solver::update_matrix`
    ///
    /// This only makes sense for direct solvers, where the factorization is 
    /// stored. This way you can solve many times for different right hand 
    /// sides.
    ///
    /// Using iterative solvers one should use the other Solver::solve method
    /// which also has the matrix as an argument.
    void solve(const Vector& rhs, Vector& solution);
    
    /// @brief solve the sytem matrix*solution = rhs
    ///
    /// How exactly this is solved is determined by the ParameterDatabase.
    void solve(const LinearOperator& matrix, const Vector& rhs,
               Vector& solution);
    
    /// @brief return a constant reference to the local ParameterDatabase
    ///
    /// Note that you can not change the behavior of this class after 
    /// construction. This method only lets you inspect the solver parameters
    const ParameterDatabase& get_db();
    
  protected:
    
    /// @brief the ParameterDatabase which controls the entire solving process
    ParameterDatabase db;
    
    /// @brief this object is only created if needed.
    std::unique_ptr<DirectSolver> direct_solver;
    /// @brief this object is only created if needed.
    std::shared_ptr<IterativeMethod<LinearOperator, Vector>> iterative_method;
    /// @brief this object is only created if needed.
    std::shared_ptr<Preconditioner<Vector>> preconditioner;
};

#endif // __SOLVER_H__
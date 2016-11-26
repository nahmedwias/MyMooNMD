#ifndef PR_TIME_NSE2D_H
#define PR_TIME_NSE2D_H

#include <Time_NSE2D.h>


class PR_Time_NSE2D : public Time_NSE2D
{
  enum class matrix_{Type14, Type1, Type2, Type3, Type4};
protected:
 struct system_per_grid
 {
   /// @brief finite element projection space
   TFESpace2D projSpace_;
   /// @brief modified mass matrix which is the
   /// product of the projection matrix and the
   /// mass matrix constructed from the vector valued
   /// functions
   BlockFEMatrix modifedMassMatrix_;
   /// @brief constructor
   system_per_grid(const TFESpace2D& velsp, const TFESpace2D &prsp,
                   const Example_TimeNSE2D& ex,
                   int order);
 };
 /** @brief a local parameter database which controls this class
   *
   * The database given to the constructor will be merged into this one. Only
   * parameters which are of interest to this class are stored (and the
   * default ParMooN parameters). Note that this usually does not include
   * other parameters such as solver parameters. Those are only in the
   * Solver object.
 */
 ParameterDatabase db_;
 /** @brief a complete system on each grid
  *
  * Note that the size of this deque is at least one and larger only in case
  * of multigrid.
  */
 std::deque<system_per_grid> systems_;
 /** @brief Definition of the used example */
 Example_TimeNSE2D example_;

 /** @brief a solver object which will solve the linear system
   *
   * Storing it means that for a direct solver we also store the factorization
   * which is usually not necessary.
 */
 Solver<BlockFEMatrix, BlockVector> solver_;

 void setParameters();

 /// @brief set projection spaces
 void set_projection_space() const;

 /// @brief print useful information
 void print_info() const;
 /// @brief old_solution for next time steps
 BlockVector oldsol_;
 /// @brief old_rhs for the next time steps
 BlockVector oldrhs_;
public:
 /// @brief constructor
 PR_Time_NSE2D(const TDomain& domain,
               const ParameterDatabase& paramdb, int ref_id);
 /// @brief constructor
 PR_Time_NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
               const Example_TimeNSE2D& ex, int ref_id=-4711);
 /// @brief assemble the initial matrices and the right hand side
 /// this will assemble the mass matrices, all A's and B's matrices
 /// and the right hand side at the initial time
 void assemble_initial();
 /// @brief assemble the right hand side
 void rhs_assemble();
 /// @brief assemble the system matrix for solver
 void system_assemble();
 /// @brief solve the system with right hand side
 void system_solve();
 /// @brief descale
 /// subtract the mass matrix from the system matrix
 /// and recover the A-blocks by descaling
 void descale();
 
 /// @brief assembling of the nonlinear terms
 void nonlinear_assemble();
 
 /// @brief check stopping criterion
 int stop_iteration(unsigned int it_counter);
 /// getters and setters
 const TFESpace2D & get_projection_space() const
 { return this->systems_.front().projSpace_; }

 const Example_TimeNSE2D & get_example_() const
 { return example_; }


};

#endif // PR_TIME_NSE2D_H

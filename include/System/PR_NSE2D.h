#ifndef PR_NSE2D_H
#define PR_NSE2D_H

#include <NSE2D.h>

class PR_NSE2D : public NSE2D
{
  enum class matrix_{Type14, Type1, Type2, Type3, Type4};
protected:
  struct system_per_grid
  {
    /// @brief projection space
    TFESpace2D projSpace_;
    /// @brief constructor
    system_per_grid(const Example_NSE2D& ex_, TCollection& coll,
                    size_t order);
  };

  std::deque<system_per_grid> system_;

  /** @brief Definition of the used example */
  Example_NSE2D example_;


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

  /** @brief a solver object which will solve the linear system
     *
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
   */
  Solver<BlockFEMatrix, BlockVector> solver_;

  /// @brief set projection spaces
  void set_projection_space() const;

  /// @brief print useful information
  void print_info() const;


public:
  /// @brief constructor
    PR_NSE2D(const TDomain & domain, const ParameterDatabase& param_db,
             const Example_NSE2D _example, unsigned int reference_id = -4711);
  /// @brief assemble matrices and rhs
  /// in the case of Stokes: only rhs is needed to be assemble
  /// in the derived class, matrices will be assembled in the
  /// base class. In addition if the reconstruction is used also
  /// in the nonlinear term, then the matrices should also be assembled
  /// using the reconstruction concept
  void all_assemble();
  
  /// @brief assemble nonlinear matrices 
  void nonlinear_assemble();

  /// getters and setters
  const TFESpace2D & get_projection_space() const
  { return this->system_.front().projSpace_; }

  const Example_NSE2D & get_example_() const
  { return example_; }
};

#endif // PR_NSE2D_H

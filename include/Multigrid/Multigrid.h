/**
 * @file New declaration of a multigrid object, which holds the necessary
 * grid information for executing a multigrid iteration.
 *
 * @note When using multigrid you must take care of the matrices being
 * correctly assembled on all levels before calling "cycle()". This might,
 * should some nonlinearity or time-dependency be included
 * include a call like former RestrictToAllGrids which informs every grid about
 * a current approximate solution.
 *
 * TODO Write own little class for cycle control.
 * TODO Implement and enable construction and usage of smoothers
 * (esp.: keep a "no smoother" for debugging purpose.)
 * TODO Make MultigridLevel store a weak_ptr instead of a raw pointer.
 * TODO Complete comments.
 * TODO Reenable step length control.
 * TODO How about parallelization?
 *
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_MULTIGRID_H_
#define INCLUDE_MULTIGRID_MULTIGRID_H_

#include <MultigridLevel.h>

#include <list>
#include <vector>

//Forward declaration.
class ParameterDatabase;
class BlockVector;

enum class MGCycle{V, W, F};
enum class MGDirection{Down, Up, End, Dummy};

class Multigrid
{
  public:

    /**
     * Set up a multigrid object. Note that there is ABSOLUTELY NO CHECKS
     * performed whether the matrices found a reasonable Multigrid hierarchie.
     * The crucial point are the restriction and prolongation methods -
     * everything else will work for any set of Rubbish matrices.
     *
     * @param db A database object containing control parameters.
     * @param matrices A vector of the matrices per level, ordered from coarsest
     * (0) to finest level.
     */
    Multigrid(const ParameterDatabase& db,
              std::list<BlockFEMatrix*> matrices);

    /// @brief Apply one complete multigrid cycle.
    void cycle();

    const BlockVector& get_finest_sol();

    void set_finest_rhs(const BlockVector& bv);

    void set_finest_sol(const BlockVector& bv);

    /// Calling this method is the sign for updating the smoothers on all levels.
    /// It must be called before every solve which was preceded by a change in
    /// the matrices (e.g. assembling).
    void update();

    static ParameterDatabase default_multigrid_database();

  private:

    /// A list of the participating levels, ordered from coarsest (0)
    /// to finest level.
    std::vector<MultigridLevel> levels_;

    std::vector<double> damp_factors_;

    size_t n_pre_smooths_;

    size_t coarse_n_maxit;

    double coarse_epsilon;

    size_t n_post_smooths_;

    int cycle_step(size_t step, size_t level);

    /// Restrict defect on level lvl and store it as rhs in the next coarsest level.
    void update_rhs_in_coarser_grid(size_t lvl);

    /// Prolongate solution on level lvl and add it (damped!) to solution on the
    /// next finest level.
    void update_solution_in_finer_grid(size_t lvl);

    /// Nuke the solution on the next coarsest grid. This is done in order
    /// to ensure a zero start iterate when smoothing.
    void set_solution_in_coarser_grid_to_zero(size_t lvl);

    //TODO Put cycle control into its own little class
    MGCycle cycle_;

    void set_cycle_control();

    void fill_recursively(std::vector<int>& mg_recursions, int level);

    void print_cycle_control() const;

    /// This vector visibly controls moving up and down in the grid hierarchy.
    std::vector<MGDirection> cycle_control_;



};




#endif /* INCLUDE_MULTIGRID_MULTIGRID_H_ */

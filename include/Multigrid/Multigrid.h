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
 * Here are the (bigger) tasks and functionalities to regain. Most of the work
 * necessary will not amass in this class but elsewhere.
 *
 * @todo TODO Reenable step length control (work in VankaSmootherNew).
 * @todo TODO Reenable MDML (work in system classes )
 * @todo TODO Parallelize (work here and in different smoothers)
 *
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_MULTIGRID_H_
#define INCLUDE_MULTIGRID_MULTIGRID_H_

#include <CycleControl.h>
#include <MultigridLevel.h>

#include <list>
#include <vector>

//Forward declarations.
class BlockVector;
class ParameterDatabase;


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

    /// @brief Apply one complete multigrid cycle. Which kind of cycle that is
    /// is determined at the time of construction via a parameter. So far V,W
    /// and F cycle are implemented.
    void cycle();

    /// Get the solution on the finest grid. Call this after a cycle is complete
    /// to get the result.
    const BlockVector& get_finest_sol();

    /// Set the right hand side on the finest grid. It must of course fit the
    /// matrix stored on the finest grid.
    /// When right hand side and (initial) solution are set on the finest grid,
    /// and a call of "update" was made since the last change of the matrices,
    /// the multigrid object is ready for a cycle.
    void set_finest_rhs(const BlockVector& bv);

    /// Set the (initial) solution on the finest grid. It must of course fit the
    /// matrix stored on the finest grid.
    /// When right hand side and (initial) solution are set on the finest grid,
    /// and a call of "update" was made since the last change of the matrices,
    /// the multigrid object is ready for a cycle.
    void set_finest_sol(const BlockVector& bv);

    /// Calling this method is the sign for updating the smoothers on all levels.
    /// It must be called before every cycle which was preceded by a change in
    /// the matrices (e.g. assembling).
    void update();

    /// Set up a database which contains default values of all control parameters
    /// necessary to control a multigrid cycle.
    static ParameterDatabase default_multigrid_database();

  private:

    /// A list of the participating levels, ordered from coarsest (0)
    /// to finest level.
    std::vector<MultigridLevel> levels_;

    /// A list of damping factors which are used when updating the solution
    /// on a finer level by adding a coarser level's solution.
    std::vector<double> damp_factors_;

    /// The number of pre-smoothing steps to perform per level before
    /// descending to the next coarser level.
    size_t n_pre_smooths_;

    /// Maximal number of smoother iteration on the coarsest level.
    size_t coarse_n_maxit;

    /// The residual to reach in the coarsest level's equation by smoothing
    /// on the coarsest level.
    double coarse_epsilon;

    /// The number of pre-smoothing steps to perform per level after
    /// ascending from the next coarser level.
    size_t n_post_smooths_;

    /// An object taking care of the order of ascends and descends between the levels.
    CycleControl control_;

    /// Restrict defect on level lvl and store it as rhs in the next coarsest level.
    void update_rhs_in_coarser_grid(size_t lvl);

    /// Prolongate solution on level lvl and add it (damped!) to solution on the
    /// next finest level.
    void update_solution_in_finer_grid(size_t lvl);

    /// Nuke the solution on the next coarsest grid. This is done in order
    /// to ensure a zero start iterate when smoothing.
    void set_solution_in_coarser_grid_to_zero(size_t lvl);

    /// Perform the operations necessary on one grid.
    int cycle_step(size_t step, size_t level);





};




#endif /* INCLUDE_MULTIGRID_MULTIGRID_H_ */

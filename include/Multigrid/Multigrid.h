/**
 * @file New declaration of a multigrid object, which holds the necessary
 * grid information for executing a multigrid iteration.
 *
 * TODO Can't do anything so far, this is only a dummy.
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

enum class MGCycle{F, V, W};
enum class MGDirection{Down, Up, End};

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

    //void update_all_rhs(const BlockVector& bv);

    void set_finest_sol(const BlockVector& bv);

    void set_finest_rhs(const BlockVector& bv);

    //Apply one multigrid cycle.
    void cycle();

    // Copy solution into solution output.
    const BlockVector& get_finest_sol();

    static ParameterDatabase default_multigrid_database();

  private:

    /// A list of the participating levels, ordered from coarsest (0)
    /// to finest level.
    std::list<MultigridLevel> levels_;

    std::vector<double> damp_factors_;

    //TODO Cycle control into own little class

    //Let's hardcode the control parameters, but pick 'em from a database (for a start)
    MGCycle cycle_;

    void set_cycle_control();

    void cycle_step(size_t step);

    void fill_recursively(std::vector<int>& mg_recursions, int level);

    void print_cycle_control() const;

    /// This vector visibly controls moving up and down in the grid hierarchy.
    std::vector<MGDirection> cycle_control_;



};




#endif /* INCLUDE_MULTIGRID_MULTIGRID_H_ */

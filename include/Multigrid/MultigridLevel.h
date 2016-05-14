/**
 * @file New declaration of a level in a multigrid method.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_MULTIGRIDLEVEL_H_
#define INCLUDE_MULTIGRID_MULTIGRIDLEVEL_H_

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Smoother.h>

#include <memory>

class Smoother;

enum class SmootherCode{DIRECT_SOLVE, JACOBI, NODAL_VANKA, CELL_VANKA, BATCH_VANKA,
NO_SMOOTHER};

class MultigridLevel
{
    //make Multigrid, which manages data transfer between levels, a good friend
    friend class Multigrid;

  public:
    MultigridLevel(BlockFEMatrix* matrix, SmootherCode sm);

    void apply_smoother();

    /// Ask the level to compute and store its current defect and residual.
    void calculate_defect();

    /// Ask the level to update its smoother.
    /// This should be called whenever the matrix pointed to by "matrix_"
    /// has changed. The update method of the multigrid object does this for
    /// all its levels.
    void update_smoother();

  private:

    BlockFEMatrix* matrix_; // TODO This is a dangerous pointer so far -
                            // TODO change to weak_ptr as soon as the system classes store shared pointers to the matrices

    BlockVector defect_;

    double residual_;

    BlockVector rhs_;

    BlockVector solution_;

    /// The smoother object. Since Smoother is an abstract base class,
    /// we can only manage it as a pointer here.
    std::shared_ptr<Smoother> smoother_;

};



#endif /* INCLUDE_MULTIGRID_MULTIGRIDLEVEL_H_ */

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

enum class SmootherCode{DIRECT_SOLVE, NODAL_VANKA, CELL_VANKA, BATCH_VANKA};

class MultigridLevel
{
  public:
    MultigridLevel(BlockFEMatrix* matrix, double beta, SmootherCode sm);

    void set_rhs(const BlockVector& rhs);

    void apply_smoother(BlockVector& sol) const;

  private:

    BlockFEMatrix* matrix_; // TODO This is a dangerous pointer so far -
                            // TODO change to weak_ptr as soon as the system classes store shared pointers to the matrices

    BlockVector defect_;
    BlockVector rhs_;
    BlockVector solution_;


    std::shared_ptr<Smoother> smoother_;

    /// The update coefficient - some damping when adding to solution.
    double beta_;

};



#endif /* INCLUDE_MULTIGRID_MULTIGRIDLEVEL_H_ */

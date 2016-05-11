/**
 * @file Pure virtual base class for multigrid smoothers.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_SMOOTHER_H_
#define INCLUDE_MULTIGRID_SMOOTHER_H_

class Smoother
{
  public:
    /**
     * Constructor, takes the BlockFEMatrix which is supposed to be smoothed
     * by this smoother.
     */
    Smoother(const BlockFEMatrix&);

    /**
     * Pure virtual smooth method - applies the smoother once.
     * @param rhs
     * @param solution
     */
    virtual void smooth(const BlockVector& rhs, BlockVector& solution ) = 0;

    /// Default destructor - class does not manage resources.
    virtual ~Smoother() = default;
};


#endif /* INCLUDE_MULTIGRID_SMOOTHER_H_ */

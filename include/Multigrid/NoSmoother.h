/**
 * @file NoSmoother.h A dummy instantiation of abstract base class Smoother.
 * Can be aplied, but does nothing when applied.
 *
 * @date 2016/10/11
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_NOSMOOTHER_H_
#define INCLUDE_MULTIGRID_NOSMOOTHER_H_

#include <Smoother.h>

class NoSmoother : public Smoother
{
  public:
    NoSmoother(){};

    /**
     * This does nothing.
     */
    virtual void smooth(const BlockVector& rhs, BlockVector& solution)
    {
    };

    virtual void update(const BlockFEMatrix&){};

    virtual ~NoSmoother() = default;
};




#endif /* INCLUDE_MULTIGRID_NOSMOOTHER_H_ */

/**
 * @file DirectSmoother.h
 *
 * @date 2016/10/11
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_DIRECTSMOOTHER_H_
#define INCLUDE_MULTIGRID_DIRECTSMOOTHER_H_

#include <DirectSolver.h>
#include <Smoother.h>

#include <memory>

/**
 * A direct solver object wrapped up as a Smoother. Should be used on
 * the coarsest grid of a multigrid method. If used on a finer grid the method
 * might become pretty slow.
 */
class DirectSmoother : public Smoother
{
  public:
    virtual void smooth(const BlockVector& rhs, BlockVector& solution ) override;

    virtual void update(const BlockFEMatrix&) override;

    /* ************* *
     * Special member functions
     * ************* */
    //! Default constructor.
    DirectSmoother();

    //! Copy constructor.
    DirectSmoother( const DirectSmoother& );

    //! Move constructor.
    DirectSmoother( DirectSmoother&& );

    //! Copy assignment operator.
    DirectSmoother& operator=( const DirectSmoother& );

    //! Move assignment operator.
    DirectSmoother& operator=( DirectSmoother&& );


    ~DirectSmoother() = default;

  private:
    std::shared_ptr<DirectSolver> solver_;

};



#endif /* INCLUDE_MULTIGRID_DIRECTSMOOTHER_H_ */

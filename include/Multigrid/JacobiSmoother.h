/**
 * @file JacobiSmoother.h
 * Wrap up a Jacobi step as a smoother for multigrid method.
 *
 * @date 2016/05/13
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_JACOBISMOOTHER_H_
#define INCLUDE_MULTIGRID_JACOBISMOOTHER_H_

#include <Smoother.h>

#include <memory>

class BlockVector;
class BlockFEMatrix;
template <class LinearOperator, class Vector> class Iteration_jacobi;

class JacobiSmoother : public Smoother
{
  public:

    void smooth(const BlockVector& rhs, BlockVector& solution) override;

    void update(const BlockFEMatrix& matrix) override;

    /* ************* *
     * Special member functions
     * ************* */
    JacobiSmoother();

    //! Copy constructor.
    JacobiSmoother( const JacobiSmoother& );

    //! Move constructor.
    JacobiSmoother( JacobiSmoother&& );

    //! Copy assignment operator.
    JacobiSmoother& operator=( const JacobiSmoother& );

    //! Move assignment operator.
    JacobiSmoother& operator=( JacobiSmoother&& );

    ~JacobiSmoother() = default;


  private:
    std::shared_ptr<Iteration_jacobi<BlockFEMatrix, BlockVector>> jacobi;
};

#endif /* INCLUDE_MULTIGRID_JACOBISMOOTHER_H_ */

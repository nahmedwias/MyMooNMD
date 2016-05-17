/**
 * @file VankaSmootherNew.h
 *
 * @date 2016/05/16
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_VANKASMOOTHERNEW_H_
#define INCLUDE_MULTIGRID_VANKASMOOTHERNEW_H_

#include <DofBatch.h>
#include <Smoother.h>

enum class VankaType {NODAL, CELL, BATCH};

/**
 * A class for Vanka smoother as used in the multgrid method. Currently the
 * implementation and nomenclature are adapted to (Navier--)Stokes saddle
 * point problems - but extending it should be relatively simple.
 */
class VankaSmootherNew : public Smoother
{
  public:
    //! Default constructor.
    VankaSmootherNew(VankaType type, double damp_factor);

    void smooth(const BlockVector& rhs, BlockVector& solution ) override;

    void update(const BlockFEMatrix&) override;

    /* ************* *
     * Special member functions. Declared but not defined, since it
     * is not yet clear whether to shallow or deep copy here.
     * ************* */
    //! Default constructor.
    VankaSmootherNew() = delete;

    //! Copy constructor.
    VankaSmootherNew( const VankaSmootherNew& );

    //! Move constructor.
    VankaSmootherNew( VankaSmootherNew&& );

    //! Copy assignment operator.
    VankaSmootherNew& operator=( const VankaSmootherNew& );

    //! Move assignment operator.
    VankaSmootherNew& operator=( VankaSmootherNew&& );

    ~VankaSmootherNew() = default;

  private:
    VankaType type_;

    size_t dimension_;

    double damp_factor_;

    std::shared_ptr<TMatrix> matrix_global_;

    std::vector<DofBatch> press_dofs_local_;

    std::vector<DofBatch> velo_dofs_local_;

// TODO Change these to weak pointers as soon as FESpaces are stored in a
// shared_ptr fashion.
#ifdef __2D__
    const TFESpace2D* pressure_space_;

    const TFESpace2D* velocity_space_;
#elif __3D__
    const TFESpace3D* pressure_space_;

    const TFESpace3D* velocity_space_;
#endif


    void set_up_pressure_batches(const TFESpace& pressureSpace);

    void set_up_velocity_batches(const TMatrix& pressureVelocityMatrix,
                                 const TFESpace& velocitySpace);


};


#endif /* INCLUDE_MULTIGRID_VANKASMOOTHERNEW_H_ */

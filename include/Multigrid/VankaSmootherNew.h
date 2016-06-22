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

//forward declaration
class DenseMatrix;

/**
 * A class for Vanka smoother as used in the multgrid method. Currently the
 * implementation and nomenclature are adapted to (Navier--)Stokes saddle
 * point problems - but extending it should be relatively simple.
 */
class VankaSmootherNew : public Smoother
{
  public:
    /** Default constructor.
     * @param type The Vanka type to use. Currently nodal, cell and batch
     * Vanke are available.
     * @param damp_factor The damping factor to use when updating the solution
     * of the global system by the solution of the local system.
     * @param store Whether to store the local systems or not. This can considerably
     * speed up the computation (especially if the grid level is traversed very
     * often) but is very memory intensive. You should use it for test cases only.
     */
    VankaSmootherNew(VankaType type, double damp_factor,  bool store = false);

    /// Perform one step of Vanka smoothing (solve all local systems).
    void smooth(const BlockVector& rhs, BlockVector& solution ) override;

    /// Update the local matrices. Must be called whenever the global matrix has changed.
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
    /// The type of Vanka smoother (nodal, cell, batch)
    VankaType type_;

    /// The dimension of the velocity equation (usually 2 or 3).
    size_t dimension_;

    /// Damping factor to be used int the local-to-global updates.
    double damp_factor_;

    /// The corresponding global matrix.
    std::shared_ptr<TMatrix> matrix_global_;

    /// The collection of pressure dofs for the local systems.
    std::vector<DofBatch> press_dofs_local_;

    /// The collection of velocity dofs for the local systems.
    std::vector<DofBatch> velo_dofs_local_;

    //Store the matrices which were assembled in order to reuse their LU factorization.
    std::vector<DenseMatrix*> local_systems_;
    bool store_systems_;

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

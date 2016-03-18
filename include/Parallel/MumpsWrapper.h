/**
 * A Class that wraps up the MUMPS solver package and enables its use
 * for ParMooN's BlockFEMatrix and ParFECommunicator.
 * It can be used only when ParMooN is compiled for parallel type MPI.
 *
 * @note Replaces the older classes ParDirectSolver and MumpsSolver.h
 *
 * @date 2016/03/14
 * @author Clemens Bartsch
 */

#ifdef _MPI

#ifndef __MUMPSWRAPPER__
#define __MUMPSWRAPPER__

extern "C"
{
#  include "dmumps_c.h"
}
#include <vector>
#include <string>

//forward declarations
class BlockFEMatrix;
class BlockVector;
class TParFECommunicator3D;

class MumpsWrapper
{
  public:
    /** */
    MumpsWrapper(
        const BlockFEMatrix& bmatrix,
        std::vector<const TParFECommunicator3D*> comms);

    // Something like this would be nice
    // TODO Figure out how to reuse factorizations etc.!
    // int factorize_and_store();

    /** Solve method - maybe return an integer error code.
     * @param[in] rhs The right hand side vector - local to each process.
     * Must fit in dimension to the stored matrix and the given communicators,
     * obviously.
     * @param[out] solution The solution - local to each proccess.
     * Must fit in dimension to the stored matrix and the given communicators,
     * obviously.
     * @param[in] comms The communicators with which to interpret
     * the blocks of the BlockVectors - these should be the same as have been
     * given to the constructor, but this is not checked explicitly.
     */
    int solve(const BlockVector& rhs, BlockVector& solution,
              std::vector<TParFECommunicator3D*> comms);

    /**
     * Write the locally stored matrix to a file in MatrixMarket
     * coordinate format. It will automatically append the number of the
     * process to the filename.
     */
    void write_matrix_distributed(std::string filename) const;

    // Declare special member functions.

    /** @brief This class is not copy constructible */
    MumpsWrapper(const MumpsWrapper&) = delete;

    /** @brief move constructor */
    MumpsWrapper(MumpsWrapper&&);

    /** @brief This class is not copy assignable */
    MumpsWrapper& operator=(const MumpsWrapper&) = delete;

    /** @brief move assignment operator */
    MumpsWrapper& operator=(MumpsWrapper&&);

    /** @brief Destructor. Must be called before MPI_FINALIZE()! */
    ~MumpsWrapper();

  private:
    /**
     * Performs a couple of input checks for the solve method.
     * As a side effect, counts the global number of dofs
     * and the local number of masters across all blocks. (Yeah, that's pretty badass...)
     * @param[in] rhs
     * @param[in] solution
     * @param[in comms
     * @param[out] n_masters_local_comms
     * @param[out] n_dofs_global_comms
     */
    void check_input_solve(const BlockVector& rhs, const BlockVector& solution,
                           const std::vector<TParFECommunicator3D*>& comms,
                           int& n_masters_local_comms,
                           int& n_dofs_global_comms);

    /**
     * Wraps call and output error handling of a MUMPS job with the stored mumps
     * handler.
     * @param job The job to execute. Choose "analyze", "factorize" or "solve".
     */
    void kick_off_job(const std::string& job);



    /**
     * Set some hard coded parameters in the mumps solver object.
     * Is called during the constructor.
     * You should not change these parameters, unless you know what you are
     * doing. Take a look at the Mumps Solver Documentation, if you are interested.
     */
    void set_mumps_parameters();

    /**
     * A method in mumps solver wrapper which transforms the
     * BlockFEMatrix to the coordinate format which mumps requires.
     * On each process only those dofs which it is master of will be regarded.
     *
     * @param bmatrix The BlockFEMatrix which is to transform.
     * @param comms An array containing the ParFECommunicators which belong
     * to the fespaces of the matrix' block rows (in the same order).
     * So far they all must have dimension 1, this method cannot yet deal
     * with using the same communicator for more than one block row.
     *
     * Performs no input checks, the constructor who calls it does that.
     */
    void store_in_distributed_coordinate_form(
        const BlockFEMatrix& bmatrix,
        std::vector<const TParFECommunicator3D*> comms);



    /// An instance of the mumps solver. Naming it "id_" is common mumps style.
    DMUMPS_STRUC_C id_;

    /**
     * Structure holding a Matrix in coordinate formate. This is the matrix
     * format the Mumps solver deals with.
     *
     * The format works like this:
     * A[irn[k]][jcn[k]]=a[k]
     * (Matrix A holds the entry a[k] at row irn[k], column jcn[k].)
     * All other entries are treated as zeroes.
     *
     * The order is not of importance, entries occuring more than once will be
     * added. We should avoid this case.
     */
    struct CoordinateMatrix{

        /// The order of the matrix. Global!
        size_t n;

        /// " non-zero". The number of explicitely stored entries, local. This does not
        /// necessarily mean only non-zero entries, but it usually does.
        size_t nz_loc;

        /// "i - row number"
        std::vector<int> irn_loc;

        /// "j - column number"
        std::vector<int> jcn_loc;

        /// The vector of stored entry values.
        std::vector<double> a_loc;
    } matrix_;

    ///A vector which holds the globally gathered right hand side.
    /// Will be filled only on root process.
    std::vector<double> rhs_global_;


};

#endif
#endif

/**
 * A Class that wraps up the MUMPS solver package and enables its use
 * for ParMooN's BlockFEMatrix and ParFECommunicator.
 * It can be used only when ParMooN is compiled for parallel type MPI.
 *
 * Although the MUMPS solver can be configured in some ways, there is several
 * "invariants" if one wants to use it with MPI ParMooN. These are especially:
 * Matrix input is distributed, right hand side input dense, solution output
 * is centralized.
 * For some insight, read the MUMPS documentation
 *   http://mumps.enseeiht.fr/doc/userguide_5.0.1.pdf
 * If you want to hard-code changes to some parameters, the method
 * "set_mumps_parameters" is the place to go to.
 *
 * TODO 1) Try out for actual ParMooN problems.
 * TODO 2) Play with sets of parameters.
 * TODO 3) Make some reuse of a factorization possible.
 * TODO 4) Maybe store just a raw pointer to the DMUMPS_STRUC_C object,
 *         to make the class at least movable (although not copyable yet.)
 *
 * @note Replaces the older classes ParDirectSolver and MumpsSolver.
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
    /**
     * Set up a Mumps wrapper for a certain BlockFEMatrix and suitable
     * ParFECommunicators. The communicators have to belong to the
     * testspaces of the matrix and must be in the same order.
     * Also their dimension variable "N_Dim" must equal 1.
     *
     * @param[in] bmatrix The matrix to wrap a mumps solver around.
     * @param[in] comms The FE communicators for the testspaces. Must be of dimension 1.
     */
    MumpsWrapper(
        const BlockFEMatrix& bmatrix,
        std::vector<const TParFECommunicator3D*> comms);

    /**
     * Solve an equation system for the wrapped up matrix with the mumps solver.
     *
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
    void solve(const BlockVector& rhs, BlockVector& solution,
              std::vector<TParFECommunicator3D*> comms);

    /**
     * Write the locally stored matrix to a file in MatrixMarket
     * coordinate format. It will automatically append the number of the
     * process to the filename.
     * @param[in] filename The filename (and path) to write the matrix to.
     */
    void write_matrix_distributed(const std::string& filename) const;


    // Declare special member functions. Delete all copy and move operations,
    // because I have not yet figured out how the DMUMPS_STRUC_C object treats
    // its memory.
    /** @brief This class is not yet copy constructible */
    MumpsWrapper(const MumpsWrapper&) = delete;

    /** @brief This class is not yet movable. */
    MumpsWrapper(MumpsWrapper&&) = delete;

    /** @brief This class is not yet copy assignable */
    MumpsWrapper& operator=(const MumpsWrapper&) = delete;

    /** @brief This class is not yet movable.  */
    MumpsWrapper& operator=(MumpsWrapper&&) = delete;

    /** @brief Destructor. Must be called before MPI_FINALIZE()! */
    ~MumpsWrapper();

  private:
    /**
     * Performs a couple of input checks for the solve method.
     * As a side effect, counts the global number of dofs
     * and the local number of masters across all blocks.
     * (Yeah, that's pretty badass, but I wanted to double-use these loops...)
     *
     * @param[in] rhs The right hand side vector.
     * @param[in] solution The solution vector.
     * @param[in comms The communicators for the FE spaces.
     * @param[out] n_masters_local_comms The local number of masters, counted
     * across all blocks of the input vectors.
     * @param[out] n_dofs_global_comms The overall, global number of dofs.
     */
    void check_input_solve(const BlockVector& rhs, const BlockVector& solution,
                           const std::vector<TParFECommunicator3D*>& comms,
                           int& n_masters_local_comms,
                           int& n_dofs_global_comms);

    /**
     * Wraps call and output error handling of a MUMPS job with the stored mumps
     * handler.
     * Throws if you give an unknown input string.
     * @param job The job to execute. Choose "analyze", "factorize" or "solve".
     */
    void kick_off_job(const std::string& job);



    /**
     * Set some hard coded parameters in the mumps solver object.
     * Is called during the constructor.
     * You should not change these parameters, unless you know what you are
     * doing.
     * Take a look at the Mumps Solver Documentation, if you are interested.
     */
    void set_mumps_parameters();

    /**
     * A method in mumps solver wrapper which transforms the
     * BlockFEMatrix to the coordinate format which mumps requires.
     * On each process only those dofs which it is master of will be regarded.
     *
     * Performs no input checks, the constructor who calls it does that.
     *
     * @param bmatrix The BlockFEMatrix which is to transform.
     * @param comms An array containing the ParFECommunicators which belong
     * to the fespaces of the matrix' block rows (in the same order).
     * So far they all must have dimension 1, this method cannot yet deal
     * with using the same communicator for more than one block row.

     */
    void store_in_distributed_coordinate_form(
        const BlockFEMatrix& bmatrix,
        std::vector<const TParFECommunicator3D*> comms);

    /**
     * Wraps two MPI calls which are used to comunicate all local master values
     * at a certain block to a global vector present in root.
     * Will order the dofs as proposed by Sashi as
     *    global_dof_id(i,p) = \sum_{k=0}^{p-1} n_own_dofs(k) + local_dof_id(i,p)
     * where i,p is the i-th local master dof on process p.
     * Performs no input checks or anything, this is just pure raw pointers and
     * mpi routines.
     *
     * @note Use it blockwise - not for the whole thing at once!
     *
     * @param GlobalArray Start here to write the global vector - here should
     * be as much space for doubles allocated as there is dofs globally.
     * @param LocalArray Start here to read the local vector.
     * @param LocalSize The size of the local vector (should be number of
     * masters of the fitting communicator).
     * @param root The root process which will store the global vector.
     */
    void gather_vector(
        double* GlobalArray, double *LocalArray, int LocalSize, int root) const;

    /// This is the reverse operation of gather_vector. Read doc there, this
    /// works the same but the other way round.
    void scatter_vector(
        double* GlobalArray, double *LocalArray, int LocalSize, int root) const;


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


};

#endif
#endif

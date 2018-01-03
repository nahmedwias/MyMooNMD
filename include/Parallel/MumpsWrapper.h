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
 * TODO 1) Make some reuse of a factorization possible.
 * TODO 2) Maybe store just a raw pointer to the DMUMPS_STRUC_C object,
 *         to make the class at least movable (although not copyable yet.)
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
class BlockMatrix;
class BlockVector;
class TParFECommunicator3D;

class MumpsWrapper
{
  public:
    /**
     * Set up a Mumps wrapper for a certain BlockFEMatrix.
     * @param[in] bmatrix The matrix to wrap a mumps solver around.
     */
    MumpsWrapper( const BlockFEMatrix& bmatrix );

    /**
     * Set up a Mumps wrapper for a certain BlockMatrix, supplying additional
     * communicators for the interpretation of the matrix.
     *
     * Calling this will be necessary only in special cases - and
     * will need the caller to know what he or she is doing.
     *
     * @param[in] matrix The matrix to wrap a mumps solver around.
     * @param[in] comms The communicators with which to interpret
     *                  the parallel structure of the matrix.
     * @param loc_to_seq This optional parameter can specify a re-ordering
     *            of the d.o.f. ("local to sequential" was my use-case).
     *            This can be used for debugging purposes.
     */
    MumpsWrapper( const BlockMatrix& matrix,
                  std::vector<const TParFECommunicator3D*> comms,
                  std::vector<int> loc_to_seq = {});

    /**
     * Solve an equation system for the wrapped up matrix with the mumps solver.
     *
     * @param[in] rhs The right hand side vector - local to each process.
     * Must fit in dimension to the stored matrix and the given communicators,
     * obviously.
     * @param[out] solution The solution - local to each proccess.
     * Must fit in dimension to the stored matrix and the given communicators,
     * obviously.
     */
    void solve(const BlockVector& rhs, BlockVector& solution);

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
     *
     */
    static void check_input_matrix(const BlockMatrix& bmatrix);

    /**
     * Performs a couple of input checks for the solve method.
     * As a side effect, counts the global number of dofs
     * and the local number of masters across all blocks.
     * (Yeah, that's pretty badass, but I wanted to double-use these loops...)
     *
     * @param[in] rhs The right hand side vector.
     * @param[in] solution The solution vector.
     * @param[out] n_masters_local_comms The local number of masters, counted
     * across all blocks of the input vectors.
     * @param[out] n_dofs_global_comms The overall, global number of dofs.
     */
    void check_input_solve(const BlockVector& rhs, const BlockVector& solution,
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
     *
     * BlockMatrix to the coordinate format which Mumps requires.
     * On each process only those dofs which it is master of will be regarded.
     *
     * @param bmatrix The BlockMatrix which is to transform. In
     * case it is a BlockFEMatrix, some extra stuff on the non-actives etc.
     * is performed.
     */
    void store_in_distributed_coordinate_form(const BlockMatrix& bmatrix,
                                              std::vector<int> loc_to_seq = {});

    /**
     * This method is need for enclosed flows. Determines which row is the
     * first pressure row and sets it to 0 with a 1 on the diagonal.
     */
    void pressure_row_correction(std::vector<const TParFECommunicator3D*> comms);

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

    /// The vector of parallel communicators belonging to the matrix.
    std::vector< const TParFECommunicator3D* > comms_;

    /// Flag that checks if there was already a factorization of the matrix computed.
    bool analyzed_and_factorized;

};

#endif
#endif

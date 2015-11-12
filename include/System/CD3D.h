/** ***************************************************************************
 *
 * @name       CD3D
 * @brief      store everything needed to solve a convection-diffusion-reaction
 *             (cdr) problem
 *
 *             Store matrix, right hand side, FE spaces, FE functions and 
 *             the solution vector of a convection-diffusion problem. This 
 *             wraps up everything which is necessary to solve a convection 
 *             diffusion problem in 3D.
 *
 * @author     Ulrich Wilbrandt, Clemens Bartsch
 * @date       09.06.15, 2015/10/20
 *
 ******************************************************************************/

#ifndef __CD3D_H__
#define __CD3D_H__

#include <BlockMatrixCD3D.h>
#include <AssembleMat3D.h>
#include <Example_CD3D.h>
#include <FEFunction3D.h>
#include <vector>
#include <deque>

#ifdef _MPI
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#include <ParDirectSolver.h>
#endif

#ifdef _OMPONLY
#include <ParDirectSolver.h>
#endif

class CD3D
{
  protected:
    /** @brief store a complete system on a particular grid
     * 
     * This combines a matrix, rhs, solution, space and function needed to
     * describe one convection-diffusion-reaction problem in 3D.
     * In MPI case the parallel infrastructure is stored, too.
     */
    struct SystemPerGrid
    {
      /** @brief Finite Element space */
      TFESpace3D feSpace_;
      /** @brief the system matrix */
      BlockMatrixCD3D matrix_;
      /** @brief the right hand side vector */
      BlockVector rhs_;
      /** @brief solution vector with one component. */
      BlockVector solution_;
      /** @brief Finite Element function */
      TFEFunction3D feFunction_;

      #ifdef _MPI
      	/** @brief A parallel FE mapper storing parallel information for the one matrix of this grid.*/
        TParFEMapper3D parMapper_;
        /** @brief A parallel FE communicator taking care for the MPI communication on this grid. */
        TParFECommunicator3D parComm_;
      #endif



#ifdef _MPI
      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.*/
      SystemPerGrid(const Example_CD3D& example,
      	            TCollection& coll, int maxSubDomainPerDof);
#else
      /** @brief constructor
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.*/
      SystemPerGrid( const Example_CD3D& example, TCollection& coll );
#endif

      //Explicitely delete copy and move - member classes do not obey rule of 0/5 yet.
      //! Delete copy constructor. No copies allowed until TFESpace3D, BlockMatrixCD3D and TFEFunction3D obey rule of 0/5.
      SystemPerGrid(const SystemPerGrid&) = delete;

      //! Delete move constructor. No moves allowed.
      SystemPerGrid(SystemPerGrid&&) = delete;

      //! Delete copy assignment operator. No copies allowed.
      SystemPerGrid& operator=(const SystemPerGrid&) = delete;

      //! Default move assignment operator. No moves allowed.
      SystemPerGrid& operator=(SystemPerGrid&&) = delete;

      //! Default destructor.
      ~SystemPerGrid() = default;

    };
    
    /** @brief a complete system on each grid
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<SystemPerGrid> systems_;
    
    /** @brief Definition of the used example */
    const Example_CD3D& example_;
    
    /** @brief a multigrid object which is set to nullptr in case it is not 
     *         needed
     */
    std::shared_ptr<TMultiGrid3D> multigrid_;

    /** @brief set parameters in database
     * 
     * This functions checks if the parameters in the database are meaningful 
     * and resets them otherwise. The hope is that after calling this function
     * this class is fully functional. 
     * 
     * If some parameters are set to unsupported values, an error occurs and 
     * throws an exception.
     */
    void setParameters();
    
  public:

#ifdef _MPI
    /** @brief constructor
     *
     * The domain must have been refined a couple of times already. On the finest
     * level the finite element spaces and functions as well as matrices,
     * solution and right hand side vectors are initialized.
     */
    CD3D(const TDomain& domain, const Example_CD3D& _example, int maxSubDomainPerDof);
#else
    /** @brief constructor 
     * 
     * The domain must have been refined a couple of times already. On the finest
     * level the finite element spaces and functions as well as matrices, 
     * solution and right hand side vectors are initialized. 
     */
    CD3D(const TDomain& domain, const Example_CD3D& _example);
#endif
    
    /** @brief assemble matrix, 
     * 
     * depending on 'TDatabase::ParamDB->DISCTYPE' different (local) assembling 
     * routines are used.
     * TODO Also in case of multigrid the matrices on all grids are
     * assembled.
     */
    void assemble();
    
    /** @brief solve the system */
    void solve();
    
    /** 
     * @brief measure errors and write pictures 
     * 
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk files are created.
     * 
     * @param i suffix for output file name, -1 means no suffix
     */
    void output(int i = -1);
    
    // getters and setters
    const BlockMatrixCD3D & getMatrix() const
    { return systems_.front().matrix_; }
    BlockMatrixCD3D & getMatrix()
    { return systems_.front().matrix_; }
    const BlockVector & getRhs() const
    { return systems_.front().rhs_; }
    BlockVector & getRhs()
    { return systems_.front().rhs_; }
    const TFEFunction3D & getFunction() const
    { return systems_.front().feFunction_; }
    const TFESpace3D & getSpace() const
    { return systems_.front().feSpace_; }
    const BlockVector & getSolution() const
    { return systems_.front().solution_; }
    BlockVector & getSolution()
    { return systems_.front().solution_; }
    unsigned int getSize() const
    { return this->systems_.front().solution_.length(); }
    const Example_CD3D& getExample() const
    { return example_; }
    
    #ifdef _MPI
    const TParFECommunicator3D& getParComm(int level){ return systems_[level].parComm_; }
    unsigned int getNTotalDof() const;
    unsigned int getNTotalCells() const;
    #endif

    // Declaration of special member functions - delete all but destructor.
    // TODO CB Apply rule of zero/five as soon as underlying classes do.

    //! Delete copy constructor.
    CD3D(const CD3D&) = delete;

    //! delete move constructor.
    CD3D(CD3D&&) = delete;

    //! delete copy assignment operator.
    CD3D& operator=(const CD3D&) = delete;

    //! delete move assignment operator
    CD3D& operator=(CD3D&&) = delete;

    //! Default destructor. Does most likely cause memory leaks.
    ~CD3D() = default;
};

#endif // __CD3D_H__

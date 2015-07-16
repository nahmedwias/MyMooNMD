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
 * @author     Ulrich Wilbrandt
 * @date       09.06.15
 *
 ******************************************************************************/

#ifndef __CD3D_H__
#define __CD3D_H__

#include <FEFunction3D.h>
#include <Example_CD3D.h>
#include <SystemMatScalar3D.h>
#include <Example_CD3D.h>
#include <vector>

#ifdef _MPI
//#include "mpi.h"
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
    /** @brief the system matrix (here one block)
     * 
     * More entries in this vector only for multigrid.
     */
    std::vector<TSystemMatScalar3D*> matrix;
    
    /** @brief the right hand side vector 
     * 
     * More entries in this vector only for multigrid.
     */
    std::vector<double*> rhs;
    
    /** @brief Finite Element function 
     * 
     * The finite element function knows its space and finite element vector so
     * that those two are not explicitly stored in this class.
     * 
     * More entries in this vector only for multigrid. 
     */
    std::vector<TFEFunction3D*> function;
    
    /** @brief Definition of the used example */
    const Example_CD3D* example;
    
    /** @brief a multigrid object which is set to nullptr in case it is not 
     *         needed
     */
    TMultiGrid3D * multigrid;
    
    /** own fespace and parallel FE Communicator */
    #ifdef _MPI
    TParFEMapper3D **ParMapper;
    TParFECommunicator3D **ParComm;
    MPI_Comm Comm;
    TParDirectSolver *DS;
    #endif
        
    #ifdef _OMPONLY
    TParDirectSolver *DS;
    #endif
    
    /** @brief set parameters in database
     * 
     * This functions checks if the parameters in the database are meaningful 
     * and resets them otherwise. The hope is that after calling this function
     * this class is fully functional. 
     * 
     * If some parameters are set to unsupported values, an error occurs and 
     * throws an exception.
     */
    void set_parameters();
    
  public:
    /** @brief constructor 
     * 
     * The domain must have been refined a couple of times already. On the finest
     * level the finite element spaces and functions as well as matrices, 
     * solution and right hand side vectors are initialized. 
     */
    CD3D(TDomain *domain, const Example_CD3D* _example = NULL);
    
    /** @brief standard destructor */
    ~CD3D();
    
    /** @brief assemble matrix, 
     * 
     * depending on 'TDatabase::ParamDB->DISCTYPE' different (local) assembling 
     * routines are used. Also in case of multigrid the matrices on all grids are
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
    TSystemMatScalar3D* getMatrix() const
    { return matrix[0]; }
    double* getRhs() const
    { return rhs[0]; }
    TFEFunction3D *get_function() const
    { return function[0]; }
    TFESpace3D * getSpace() const
    { return function[0]->GetFESpace3D(); }
    double * getSolution() const
    { return function[0]->GetValues(); }
    unsigned int getSize() const
    { return function[0]->GetLength(); }
    const Example_CD3D* getExample() const
    { return example; }
    
    #ifdef _MPI
    TParFECommunicator3D *Get_ParComm(int level)
    { return ParComm[level]; }
    #endif
};

#endif // __CD3D_H__

/** ************************************************************************ 
*
* @class     Darcy2D
* @brief     stores the information of a 2D Darcy system matrix 
* @author    Ulrich Wilbrandt
* @date      15.03.15
 ************************************************************************  */


#ifndef __SYSTEMDARCY2D__
#define __SYSTEMDARCY2D__

#include <BlockMatrixDarcy2D.h>
#include <Example_Darcy2D.h>
#include <LocalAssembling2D.h>
#include <FEFunction2D.h>

#include <vector>

/**class for 2D scalar system matrix */
class Darcy2D
{
  protected:
    /** @brief the system matrix (here one block)
     * 
     * More entries in this vector only for multigrid.
     */
    std::vector<BlockMatrixDarcy2D*> matrix;
    
    /** @brief the right hand side vector 
     * 
     * More entries in this vector only for multigrid.
     */
    std::vector<double*> rhs;
    
    /** @brief Finite Element functions for velocity and pressure
     * 
     * The finite element functions know their spaces and finite element vectors
     * so that those two are not explicitly stored in this class.
     * 
     * More entries in these vectors only for multigrid. 
     */
    std::vector<TFEFunction2D *> u; // velocity
    std::vector<TFEFunction2D *> p; // pressure
    
    /** @brief Definition of the used example */
    Example_Darcy2D* example;
    
    /** @brief a multigrid object which is set to nullptr in case it is not 
     *         needed
     * 
     * @note multigrid for Darcy type problems is not yet implemented
     */
    TMultiGrid2D * multigrid;
    
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
     * The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and 
     * functions as well as matrices, solution and right hand side vectors are 
     * initialized. 
     */
    Darcy2D(TDomain *domain, Example_Darcy2D* _example = NULL);
    
    /** @brief standard destructor */
    ~Darcy2D();
    
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
    BlockMatrixDarcy2D* getMatrix() const
    { return matrix[0]; }
    double* getRhs() const
    { return rhs[0]; }
    TFEFunction2D *get_velocity() const
    { return u[0]; }
    TFEFunction2D *get_pressure() const
    { return p[0]; }
    TFESpace2D * get_velocity_space() const
    { return u[0]->GetFESpace2D(); }
    TFESpace2D * get_pressure_space() const
    { return p[0]->GetFESpace2D(); }
    double * get_solution() const
    { return u[0]->GetValues(); }
    unsigned int get_size() const
    { return u[0]->GetLength() + p[0]->GetLength(); }
    const Example_Darcy2D* getExample() const
    { return example; }
};

#endif // __SYSTEMMATDARCY2D__

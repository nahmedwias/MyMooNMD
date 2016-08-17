/** ***************************************************************************
 *
 * @name   2D_Assembler
 * @brief  Assemble
 *
 *
 *
 * @author     Alfonso Caiazzo and Laura Blank
 * @date       16.08.2016
 *
 ******************************************************************************/

#ifndef __ASSEMBLER2D__
#define __ASSEMBLER2D__

#include <Constants.h>
#include <FEDatabase2D.h>
#include <LocalAssembling2D.h>

#ifdef __3D__
  #include <Aux2D3D.h>
#endif

class Assembler{
    
protected:
    
public:
    
    /** @brief constructor
     *
     * This constructor calls the other constructor creating an ??
     * object for you. See there for more documentation.
     */
    Assembler(){};
    
    /** a function from a finite element space */
    void Assemble2D(int n_fespaces, const TFESpace2D **fespaces,
                    int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                    int n_matrices, TMatrix2D **matrices,
                    int n_rhs, double **rhs, const TFESpace2D **ferhs,
                    BoundCondFunct2D **BoundaryConditions,
                    BoundValueFunct2D **BoundaryValues,
                    LocalAssembling2D& la
#ifdef __3D__
                    , TAux2D3D *Aux2D3D
#endif
                    , int AssemblePhaseID = -1
                    );
    
};
#endif // __ASSEMBLE2D__

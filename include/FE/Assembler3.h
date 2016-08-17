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
#include <Example2D.h>

#ifdef __3D__
  #include <Aux2D3D.h>
#endif



class Assembler3{
    
protected:
    
public:
    
    /** @brief constructor
     *
     * This constructor calls the other constructor creating an ??
     * object for you. See there for more documentation.
     */
    Assembler3(LocalAssembling2D_type type, TFEFunction2D **fefunctions2d,
                                    CoeffFct2D *coeffs);
    
//    LocalAssembling2D la(LocalAssembling2D_type type, TFEFunction2D **fefunctions2d,
//                         CoeffFct2D *coeffs);
    
    
    LocalAssembling2D la;
    

    
    /** a function from a finite element space */
    void Assemble2D(int n_fespaces, const TFESpace2D **fespaces,
                    int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                    int n_matrices, TMatrix2D **matrices,
                    int n_rhs, double **rhs, const TFESpace2D **ferhs,
                    //BoundCondFunct2D **BoundaryConditions,
                    //BoundValueFunct2D **BoundaryValues,
                    //LocalAssembling2D& la,
                    const Example2D& example,
                    int AssemblePhaseID=-1
                    );
    
    
    
    void loop_over_cells(const TFESpace2D** fespaces, int AssemblePhaseID, int n_fespaces,
                         //LocalAssembling2D& la,
                         int n_sqmatrices, double ***LocMatrices,
                         TSquareMatrix2D **sqmatrices,int n_matrices,TMatrix2D **matrices,
                         int **TestGlobalNumbers, int **TestBeginIndex,
                         int **AnsatzGlobalNumbers, int **AnsatzBeginIndex, double **HangingEntries, double **HangingRhs,
                         double *righthand, int **RhsBeginIndex, int **RhsGlobalNumbers,
                         //BoundCondFunct2D** BoundaryConditions,
                         //BoundValueFunct2D** BoundaryValues,
                         int N_AllMatrices, double **Matrices, double **Param, int n_rhs, double** rhs,
                         double **AuxArray,double **LocRhs,const TFESpace2D** ferhs,
                         int *N_BaseFunct,  BaseFunct2D *BaseFuncts,int **GlobalNumbers, int **BeginIndex,
                         const Example2D& example);
    
    
    ~Assembler3();
    
   


    
    };
#endif // __ASSEMBLE2D__
//double *Param[MaxN_QuadPoints_2D]
//double *AuxArray[MaxN_QuadPoints_2D]
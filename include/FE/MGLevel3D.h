// =======================================================================
// %W% %G%
//
// Class:       TMGLevel3D
// Purpose:     store all data for one level in a multi grid method in 3d
//
// Author:      Gunar Matthies 26.06.2000
//
// History:     26.06.2000 start of implementation
//
// =======================================================================

#ifndef __MGLEVEL3D__
#define __MGLEVEL3D__

#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>

class TMGLevel3D
{
  protected:
    /** level number */
    int Level;

    /** FE space */
    TFESpace3D *FESpace;

    /** permutation vector */
    int *Permutation;

    /** number of active nodes */
    int N_Active;

    /** upper bound for hanging node number */
    int HangingNodeBound;

    /** number of Dirichlet nodes */
    int N_Dirichlet;

    /** number of all degrees of freedom */
    int N_DOF;

    /** used matrix */
    TSquareMatrix3D *A;

    /** structure of used matrix */
    TSquareStructure3D *MatrixStructure;

    /** row pointer for matrix */
    int *RowPtr;

    /** column number vector */
    int *KCol;

    /** matrix entries */
    double *Entries;

    /** array with right-hand sides */
    double *Rhs;

    /** array with approximate solution */
    double *X;

    /** number of auxiliary vectors */
    int N_Aux;

    /** array of auxiliary vectors */
    double **Aux;

    /** array for additional data, e.g. ILU decomposition */
    double *Additional;

    /** generate ILU decomposition */
    void ILUDecomposition();
    
#ifdef _MPI
     /** number of all degrees of freedom in own cells*/
     int OwnN_DOF;   
    
     TParFECommunicator3D *ParComm;

     /** FEFunction in own+Hallo fe space */
     TFEFunction3D *C;     
     
     /** Own FE space */
     TFESpace3D *OwnScalarSpace;

     /** FEFunction in own fe space */
     TFEFunction3D *OwnC;
 
     /** Own solution */
     double *OwnSolArray;
     
     /** Reorder of sol array */
     int *Reorder,N_Master,N_Int,N_Dept;
     
     /** Coloring variables */
     int N_CMaster, N_CDept, N_CInt, *ptrCMaster, *ptrCDept, *ptrCInt;
     
#endif    

  public:

    
    /** constructor */
    TMGLevel3D(int level, TSquareMatrix3D *A,
             double *rhs, double *sol, int n_aux,
             int *permutation);
    
#ifdef _MPI
/** constructor for parallel */
    TMGLevel3D(int level, TSquareMatrix3D *a, double *rhs, double *sol, 
                       TFEFunction3D *c, TParFECommunicator3D *parComm,TFESpace3D *ownScalarSpace, int n_aux,
                       int *permutation);
#endif  

    /** destructor */
    ~TMGLevel3D();

    /** return i-th auxiliary vector */
    double *GetAuxVector(int i);

    /** return FunctionVectors */
    double *GetSolution()
    { return X; }

    /** return Rhs */
    double *GetRhs()
    { return Rhs; }

    /** return AuxVectors */
    double **GetAuxVectors()
    { return Aux; }

    /** return system matrix */
    TSquareMatrix3D *GetMatrix()
    { return A; }

    /** return number of degrees of freedom */
    int GetN_DOF()
    { return N_DOF; }

    /** get HangingNodeBound */
    int GetHangingNodeBound()
    { return HangingNodeBound; }

    /** get number of Dirichlet nodes */
    int GetN_Dirichlet()
    { return N_Dirichlet; }

    /** calculate defect */
    void Defect(double *sol, double *f, double *d, double &res);

    /** update solution */
    void Update(double *sol, double *upd);

    /** correct Dirichlet and hanging nodes */
    void CorrectNodes(double *vect);

    /** correct defect */
    void CorrectDefect(double *vect);

    /** reset vector to zero */
    void Reset(double *vect);

    /** return FE space */
    TFESpace3D *GetFESpace()
    { return FESpace; }

    /** smoother */
    void ILU(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void SOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);
    
    void SOR_Re(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);
    
    void SOR_Color(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void SSOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void Jacobi(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void Block2x2(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** solve exact on this level */
    void SolveExact(double *u1, double *rhs1);

    /** step length control */
    double StepLengthControl(double *u, 
                         double *uold, 
                         double *def,
                         int N_Parameters, 
                         double *Parameters);
    
    
        
    #ifdef _MPI       
    TFESpace3D *GetOwnFESpace()
       { return OwnScalarSpace; }

    TParFECommunicator3D *GetParComm()
      { return ParComm; }

    double *GetOwnSolution()
      { return OwnSolArray; }  

     TFEFunction3D *GetFEFunction()
      { return C; }

     TFEFunction3D *GetOwnFEFunction()
      { return OwnC; }

    
#endif

};

#endif

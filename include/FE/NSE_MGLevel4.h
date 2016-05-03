// =======================================================================
// @(#)NSE_MGLevel4.h        1.6 07/03/00
//
// Class:       TNSE_MGLevel4
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system of
//              type 4 (all Aij, B1, B2, B3)
//
// Author:      Volker John 25.07.2000
//
// History:      25.07.2000start of implementation
//
// =======================================================================

#ifndef __NSE_MGLEVEL4__
#define __NSE_MGLEVEL4__

#include <NSE_MGLevel.h>

#include <VankaSmoother_NSE4.h>

class TNSE_MGLevel4 : public TNSE_MGLevel
{
  protected:
#ifdef __2D__
    /** matrix A11 */
    TSquareMatrix2D *A11;

    /** matrix A12 */
    TSquareMatrix2D *A12;

    /** matrix A21 */
    TSquareMatrix2D *A21;

    /** matrix A22 */
    TSquareMatrix2D *A22;

    /** structure of matrix A */
    const TStructure *StructureA;

    /** matrix B1 */
    TMatrix2D *B1;

    /** matrix B2 */
    TMatrix2D *B2;

    /** matrix B1T */
    TMatrix2D *B1T;

    /** matrix B2 */
    TMatrix2D *B2T;

    /** structure of matrix B */
    const TStructure *StructureB;

    /** structure of matrix BT */
    const TStructure *StructureBT;

    /** structure of matrix C */
    const TStructure *StructureC;

    /** matrix C */
    TMatrix2D *C;
#endif  

#ifdef __3D__
    /** matrix A11 */
    TSquareMatrix3D *A11;

    /** matrix A12 */
    TSquareMatrix3D *A12;

    /** matrix A13 */
    TSquareMatrix3D *A13;

    /** matrix A21 */
    TSquareMatrix3D *A21;

    /** matrix A22 */
    TSquareMatrix3D *A22;

    /** matrix A23 */
    TSquareMatrix3D *A23;

    /** matrix A31 */
    TSquareMatrix3D *A31;

    /** matrix A32 */
    TSquareMatrix3D *A32;

    /** matrix A33 */
    TSquareMatrix3D *A33;

    /** structure of matrix A */
    const TStructure *StructureA;

    /** matrix B1 */
    TMatrix3D *B1;

    /** matrix B2 */
    TMatrix3D *B2;

    /** matrix B3 */
    TMatrix3D *B3;

    /** matrix B1T */
    TMatrix3D *B1T;

    /** matrix B2 */
    TMatrix3D *B2T;

    /** matrix B3 */
    TMatrix3D *B3T;

    /** structure of matrix B */
    const TStructure *StructureB;
 
    /** structure of matrix BT */
    const TStructure *StructureBT;

    /** structure of matrix C */
    const TStructure *StructureC;

    /** matrix C */
    TMatrix3D *C;
#endif  

    /** row pointer for matrix A */
    const int *ARowPtr;

    /** column number vector for matrix A */
    const int *AKCol;

    /** matrix entries of matrix A */
    double *A11Entries;

    /** matrix entries of matrix A */
    double *A12Entries;

    /** matrix entries of matrix A */
    double *A21Entries;

    /** matrix entries of matrix A */
    double *A22Entries;

    /** matrix entries of matrix B1 */
    double *B1Entries;

    /** matrix entries of matrix B2 */
    double *B2Entries;

    /** matrix entries of matrix B1 */
    double *B1TEntries;

    /** matrix entries of matrix B2 */
    double *B2TEntries;
#ifdef __3D__
    /** matrix entries of matrix A */
    double *A13Entries;

    /** matrix entries of matrix A */
    double *A23Entries;

    /** matrix entries of matrix A */
    double *A31Entries;

    /** matrix entries of matrix A */
    double *A32Entries;

    /** matrix entries of matrix A */
    double *A33Entries;

    /** matrix entries of matrix B3 */
    double *B3Entries;

    /** matrix entries of matrix BT3 */
    double *B3TEntries;

#endif  

    /** row pointer for matrix B */
    const int *BRowPtr;

    /** column number vector for matrix B */
    const int *BKCol;

    /** row pointer for matrix BT */
    const int *BTRowPtr;

    /** column number vector for matrix BT */
    const int *BTKCol;

    /** row pointer for matrix C */
    const int *CRowPtr;

    /** column number vector for matrix C */
    const int *CKCol;

    /** matrix entries of matrix C */
    double *CEntries;
    
  public:
    /** constructor */
#ifdef __2D__
    TNSE_MGLevel4(int level, 
                  TSquareMatrix2D *A11, TSquareMatrix2D *A12, 
                  TSquareMatrix2D *A21, TSquareMatrix2D *A22, 
                  TMatrix2D *B1, TMatrix2D *B2,
                  TMatrix2D *B1T, TMatrix2D *B2T,
                  double *f1, double *u1,
                  int n_aux, double *al, int VelocitySpace, 
                  int PressureSpace, TCollection *coll, int *dw);

    TNSE_MGLevel4(int level,
                  TSquareMatrix2D *a11, TSquareMatrix2D *a12,
                  TSquareMatrix2D *a21, TSquareMatrix2D *a22,
                  TMatrix2D *b1, TMatrix2D *b2,
                  TMatrix2D *b1t, TMatrix2D *b2t,
                  TMatrix2D *c,
                  double *f1, double *u1,
                  int n_aux, double *al, int velocity_space,
                  int pressure_space, TCollection *Coll,
                  int *dw);
#endif  
#ifdef __3D__
    TNSE_MGLevel4(int level, 
                  TSquareMatrix3D *A11, TSquareMatrix3D *A12, 
                  TSquareMatrix3D *A13, 
                  TSquareMatrix3D *A21, TSquareMatrix3D *A22, 
                  TSquareMatrix3D *A23, 
                  TSquareMatrix3D *A31, TSquareMatrix3D *A32, 
                  TSquareMatrix3D *A33, 
                  TMatrix3D *B1, TMatrix3D *B2, TMatrix3D *B3,  
                  TMatrix3D *B1T, TMatrix3D *B2T, TMatrix3D *B3T,
                  double *f1, double *u1,
                  int n_aux, double *al, int VelocitySpace, 
                  int PressureSpace, TCollection *coll, int *dw);

TNSE_MGLevel4(int level,
                               TSquareMatrix3D *a11, TSquareMatrix3D *a12,
                               TSquareMatrix3D *a13, TSquareMatrix3D *a21,
                               TSquareMatrix3D *a22, TSquareMatrix3D *a23,
                               TSquareMatrix3D *a31, TSquareMatrix3D *a32,
                               TSquareMatrix3D *a33,
                               TMatrix3D *b1, TMatrix3D *b2, TMatrix3D *b3,
                               TMatrix3D *b1t, TMatrix3D *b2t, TMatrix3D *b3t,
                               TMatrix3D *c,
                               double *f1, double *u1,
                               int n_aux, double *al, int velocity_space,
                               int pressure_space, TCollection *Coll,
	      int *dw);

#endif  

    /** destructor */
    ~TNSE_MGLevel4();

    virtual void Defect(double *u1, double *f1, double *d1, double &res);

    /** correct Dirichlet and hanging nodes */
    virtual void CorrectNodes(double *u1);

    /** Vanka smoother */
    virtual void CellVanka(double *u1, double *rhs1, double *aux, 
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

    /** Vanka smoother */
    virtual void NodalVanka(double *u1, double *rhs1, double *aux, 
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

    /** solve exact on this level */
    virtual void SolveExact(double *u1, double *rhs1);

    /** solve exact on this level */
    virtual void SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag);

    /** Braess Sarazin smoother */
    virtual void BraessSarazin(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters,int N_Levels);

    /** step length control for Vanka */
    virtual double StepLengthControl(double *u1, double *u1old, double *def1,
                                     int N_Parameters, double *Parameter);

    /** print all matrices and both right hand sides */
    virtual void PrintAll();
    
    /****************************************************************
     * Getters for the matrices.
     ****************************************************************/
#ifdef __2D__
 /*!
     * @brief Return a pointer to matrix A11.
     *
     * @return A pointer to matrix A11 (2D case only!)
     */
    TSquareMatrix2D* getA11() const
    { return A11; }

    /*!
     * @brief Return a pointer to matrix A12.
     *
     * @return A pointer to matrix A12 (2D case only!)
     */
    TSquareMatrix2D* getA12() const
    { return A12; }

    /*!
     * @brief Return a pointer to matrix A21.
     *
     * @return A pointer to matrix A21 (2D case only!)
     */
    TSquareMatrix2D* getA21() const
    { return A21; }

    /*!
     * @brief Return a pointer to matrix A22.
     *
     * @return A pointer to matrix A22 (2D case only!)
     */
    TSquareMatrix2D* getA22() const
    { return A22; }

    /*!
     * @brief Return a pointer to matrix B1.
     *
     * @return A pointer to matrix B1 (2D case only!)
     */
    TMatrix2D* getB1() const
    { return B1; }

    /*!
     * @brief Return a pointer to matrix B2.
     *
     * @return A pointer to matrix B2 (2D case only!)
     */
    TMatrix2D* getB2() const
    { return B2; }

    /*!
     * @brief Return a pointer to matrix B1T.
     *
     * @return A pointer to matrix B1T (2D case only!)
     */
    TMatrix2D* getB1T() const
    { return B1T; }

    /*!
     * @brief Return a pointer to matrix B2T.
     *
     * @return A pointer to matrix B2T (2D case only!)
     */
    TMatrix2D* getB2T() const
    { return B2T; }

    /*!
     * @brief Return a pointer to matrix C.
     *
     * @return A pointer to matrix C (2D case only!)
     */
    TMatrix2D* getC() const
    { return C; }
#else
/*!
     * @brief Return a pointer to matrix A11.
     *
     * @return A pointer to matrix A11 (3D case only!)
     */
    TSquareMatrix3D* getA11() const
    { return A11; }

    /*!
     * @brief Return a pointer to matrix A12.
     *
     * @return A pointer to matrix A12 (3D case only!)
     */
    TSquareMatrix3D* getA12() const
    { return A12; }

    /*!
     * @brief Return a pointer to matrix A12.
     *
     * @return A pointer to matrix A12 (3D case only!)
     */
    TSquareMatrix3D* getA13() const
    { return A13; }

    /*!
     * @brief Return a pointer to matrix A21.
     *
     * @return A pointer to matrix A21 (3D case only!)
     */
    TSquareMatrix3D* getA21() const
    { return A21; }

    /*!
     * @brief Return a pointer to matrix A22.
     *
     * @return A pointer to matrix A22 (3D case only!)
     */
    TSquareMatrix3D* getA22() const
    { return A22; }

    /*!
     * @brief Return a pointer to matrix A22.
     *
     * @return A pointer to matrix A22 (3D case only!)
     */
    TSquareMatrix3D* getA23() const
    { return A23; }

    /*!
     * @brief Return a pointer to matrix A31.
     *
     * @return A pointer to matrix A31 (3D case only!)
     */
    TSquareMatrix3D* getA31() const
    { return A31; }

    /*!
     * @brief Return a pointer to matrix A32.
     *
     * @return A pointer to matrix A32 (3D case only!)
     */
    TSquareMatrix3D* getA32() const
    { return A32; }

    /*!
     * @brief Return a pointer to matrix A33.
     *
     * @return A pointer to matrix A33 (3D case only!)
     */
    TSquareMatrix3D* getA33() const
    { return A33; }

    /*!
     * @brief Return a pointer to matrix B1.
     *
     * @return A pointer to matrix B1 (3D case only!)
     */
    TMatrix3D* getB1() const
    { return B1; }

    /*!
     * @brief Return a pointer to matrix B2.
     *
     * @return A pointer to matrix B2 (3D case only!)
     */
    TMatrix3D* getB2() const
    { return B2; }

    /*!
     * @brief Return a pointer to matrix B2.
     *
     * @return A pointer to matrix B2 (3D case only!)
     */
    TMatrix3D* getB3() const
    { return B3; }

    /*!
     * @brief Return a pointer to matrix B1T.
     *
     * @return A pointer to matrix B1T (3D case only!)
     */
    TMatrix3D* getB1T() const
    { return B1T; }

    /*!
     * @brief Return a pointer to matrix B2T.
     *
     * @return A pointer to matrix B2T (3D case only!)
     */
    TMatrix3D* getB2T() const
    { return B2T; }

    /*!
     * @brief Return a pointer to matrix B2T.
     *
     * @return A pointer to matrix B2T (3D case only!)
     */
    TMatrix3D* getB3T() const
    { return B3T; }

    /*!
     * @brief Return a pointer to matrix C.
     *
     * @return A pointer to matrix C (3D case only!)
     */
    TMatrix3D* getC() const
    { return C; }
#endif
    
private:
  //! The Vanka smoother object used for Vanka
  VankaSmoother_NSE4 vankaSmoother_;
  
public:
  //! @brief Initialize the vankaSmoother_ object.
  //! Call only after everything global is initialized.
  void initializeSmoother() override;

  //! @brief Perform one step of Vanka smoothing on the level.
  //! Call only after initializer finished.
  void applySmoother(double *currentSolution, const double* const currentRHS,
                     double *storeOldSolution) override;


};

#endif


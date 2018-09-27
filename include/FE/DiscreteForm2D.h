// =======================================================================
// @(#)DiscreteForm2D.h        1.6 10/18/99
// 
// Class:       TDiscreteForm2D
// Purpose:     assemble on one cell a couple of bilinear and linear forms
//
// Author:      Gunar Matthies (02.10.98)
//
// History:     start of implementation 02.10.98 (Gunar Matthies)
//              add new data 14.04.99 (Gunar Matthies)
//
// =======================================================================

#ifndef __DISCRETEFORM2D__
#define __DISCRETEFORM2D__

#include <Enumerations.h>
#include <Constants.h>

/** a function from a finite element space */
class TDiscreteForm2D
{
  protected:
    /** name */
    char *Name;

    /** some more describing words */
    char *Description;

    /** number of terms */
    int N_Terms;

    /** number of involved spaces */
    int N_Spaces;

    /** Are second derivatives from space i needed */
    bool *Needs2ndDerivatives;

    /** multiindices for derivatives of ansatz functions */
    MultiIndex2D *Derivatives;

    /** number of FESpace2D which is used for a derivative */
    int *FESpaceNumber;

    /** number of matrices */
    int N_Matrices;

    /** number of right-hand sides */
    int N_Rhs;

    /** which FE space corresponds to each row */
    int *RowSpace;

    /** which FE space corresponds to each column */
    int *ColumnSpace;

    /** which FE space corresponds to each right-hand side */
    int *RhsSpace;

    /** function for calculating the coefficients */
    CoeffFct2D Coeffs;

    /** function doing the real assembling */
    AssembleFct2D *Assemble;

    /** function doing the real assembling using parameters from 
        argument list */
    AssembleFctParam AssembleParam;

    /** function for manipulating the coefficients */
    ManipulateFct2D *Manipulate;

    /** memory for storing the original value arrays */
    double ***AllOrigValues;

    /** memory for storing the original value arrays at one point */
    double **OrigValues;

  public:
    /** constructor */
    TDiscreteForm2D(char *name, char *description, int n_terms,
        MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFct2D *assemble, CoeffFct2D coeffs,
        ManipulateFct2D *manipulate);

    /** constructor with assembling using parameters */
    TDiscreteForm2D(char *name, char *description, int n_terms,
        MultiIndex2D *derivatives, int *fespacenumber,
        int n_matrices, int n_rhs,
        int *rowspace, int *columnspace, int *rhsspace,
        AssembleFctParam assembleparam, CoeffFct2D coeffs,
        ManipulateFct2D *manipulate);

    /** destructor */
    ~TDiscreteForm2D();

    /** return name */
    char *GetName() const
    { return Name; }

    /** return description */
    char *GetDescription() const
    { return Description; }

    /** return local stiffness matrix */
    void GetLocalForms(int N_Points, double *weights, double *AbsDetjk,
                       double hK, double *X, double *Y,
                       int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                       double **Parameters, double **AuxArray,
                       TBaseCell *Cell, int N_Matrices, int N_Rhs,
                       double ***LocMatrix, double **LocRhs,
                       double factor = 1.);
    
    /** assemble local matrices and right hand sides 
     * 
     * This is a simplified version of the above GetLocalForms(...).
     */
    void GetLocalForms(int N_Points, double *weights, double *AbsDetjk,
                        double hK, double *X, double *Y,
                        int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                        TBaseCell *Cell,double ***LocMatrix, double **LocRhs);

    /** return array Needs2ndDerivatives */
    bool *GetNeeds2ndDerivatives() const
    { return Needs2ndDerivatives; };

    /** function for calculating the coefficients */
    const CoeffFct2D& GetCoeffFct() const
    { return Coeffs; }
    
    /** return the index of the row space of the i-th matrix */
    int rowSpaceOfMat(int i) const
    { return RowSpace[i]; }
    
    /** return the index of the column space of the i-th matrix */
    int colSpaceOfMat(int i) const
    { return ColumnSpace[i]; }
    
    int Get_NTerms() const
    { return N_Terms; }
    
    int Get_N_Spaces() const
    { return N_Spaces; }
    
    MultiIndex2D get_derivative(int i) const // i should be smaller than N_Terms
    { return this->Derivatives[i]; }
    
    int get_FESpaceNumber(int i) const
    { return this->FESpaceNumber[i]; }
    
    int get_N_Matrices() const
    { return N_Matrices; }
    
    int get_N_Rhs() const
    { return N_Rhs; }
    
    int get_RhsSpace(int i) const
    { return RhsSpace[i]; }
    
    AssembleFctParam get_AssembleParam() const
    { return AssembleParam; }
    
    ManipulateFct2D * get_Manipulate() const
    { return Manipulate; }
};
#endif

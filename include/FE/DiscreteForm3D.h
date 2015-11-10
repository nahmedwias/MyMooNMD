// =======================================================================
// %W% %G%
// 
// Class:       TDiscreteForm3D
// Purpose:     assemble on one cell a couple of bilinear and linear forms
//
// Author:      Gunar Matthies (02.10.98)
//
// History:     start of implementation 02.10.98 (Gunar Matthies)
//              add new data 14.04.99 (Gunar Matthies)
//
// =======================================================================

#ifndef __DISCRETEFORM3D__
#define __DISCRETEFORM3D__

#include <Enumerations.h>
#include <Constants.h>

/** a function from a finite element space */
class TDiscreteForm3D
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
    MultiIndex3D *Derivatives;

    /** number of FESpace3D which is used for a derivative */
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
    CoeffFct3D *Coeffs;

    /** function doing the real assembling */
    AssembleFct3D *Assemble;

    /** function doing the real assembling using parameters from 
        argument list */
    AssembleFctParam3D *AssembleParam;

    /** function for manipulating the coefficients */
    ManipulateFct3D *Manipulate;

    /** memory for storing the original value arrays */
    double ***AllOrigValues;

    /** memory for storing the original value arrays at one point */
    double **OrigValues;

  public:
    /** constructor */
    TDiscreteForm3D(char *name, char *description, int n_terms,
                    MultiIndex3D *derivatives, int *fespacenumber,
                    int n_matrices, int n_rhs,
                    int *rowspace, int *columnspace, int *rhsspace,
                    AssembleFct3D *assemble, CoeffFct3D *coeffs,
                    ManipulateFct3D *manipulate);

    /** constructor with assembling using parameters */
    TDiscreteForm3D(char *name, char *description, int n_terms,
                    MultiIndex3D *derivatives, int *fespacenumber,
                    int n_matrices, int n_rhs,
                    int *rowspace, int *columnspace, int *rhsspace,
                    AssembleFctParam3D *assembleparam, CoeffFct3D *coeffs,
                    ManipulateFct3D *manipulate);

    /** destructor */
    ~TDiscreteForm3D();

    /** return local stiffness matrix */
    void GetLocalForms(int N_Points, double *weights, double *AbsDetjk,
                        double hK, double *X, double *Y, double *Z,
                        int *N_BaseFuncts, BaseFunct3D *BaseFuncts,
                        double **Parameters, double **AuxArray,
                        TBaseCell *Cell,
                        int N_Matrices, int N_Rhs,
                        double ***LocMatrix, double **LocRhs);

    //For backwards compatibility: getters which do not follow naming conventions.
    /** return name */
    [[deprecated]] char *GetName()
    { return Name; }

    /** return description */
    [[deprecated]] char *GetDescription()
    { return Description; }

    /** return array Needs2ndDerivatives */
    [[deprecated]] bool *GetNeeds2ndDerivatives()
    { return Needs2ndDerivatives; };

    /** function for calculating the coefficients */
    [[deprecated]] CoeffFct3D *GetCoeffFct()
    { return Coeffs; }

    //Getter methods.

	double*** getAllOrigValues() const {
		return AllOrigValues;
	}

	AssembleFct3D* getAssemble() const {
		return Assemble;
	}

	AssembleFctParam3D* getAssembleParam() const {
		return AssembleParam;
	}

	CoeffFct3D* getCoeffs() const {
		return Coeffs;
	}

	int* getColumnSpace() const {
		return ColumnSpace;
	}

    /** return the index of the column space of the i-th matrix */
    int colSpaceOfMat(int i) const
    {
    	if (i >= N_Matrices)
    	{
    		ErrMsg("Array out of bound: i >= N_Matrices.");
    	}
    	return ColumnSpace[i];
    }

	MultiIndex3D* getDerivatives() const {
		return Derivatives;
	}

    MultiIndex3D getDerivative(int i) const
    {
    	if (i >= N_Terms)
    	{
    		ErrMsg("Array out of bound: i >= N_Terms.");
    	}
    	return Derivatives[i];
    }



	char* getDescription() const {
		return Description;
	}

	int* getFeSpaceNumber() const {
		return FESpaceNumber;
	}

    int getFeSpaceNumber(int i) const
    {
    	if (i >= N_Terms)
    	{
    		ErrMsg("Array out of bound: i >= N_Terms.");
    	}
    	return this->FESpaceNumber[i];
    }

	ManipulateFct3D* getManipulate() const {
		return Manipulate;
	}

	int getNMatrices() const {
		return N_Matrices;
	}

	int getNRhs() const {
		return N_Rhs;
	}

	int getNSpaces() const {
		return N_Spaces;
	}

	int getNTerms() const {
		return N_Terms;
	}

	char* getName() const {
		return Name;
	}

	bool* getNeeds2ndDerivatives() const {
		return Needs2ndDerivatives;
	}

	double** getOrigValues() const {
		return OrigValues;
	}

	int* getRhsSpace() const {
		return RhsSpace;
	}

    int getRhsSpace(int i) const
    {
    	if (i >= N_Rhs)
    	{
    		ErrMsg("Array out of bound: i >= N_Rhs.");
    	}
    	return RhsSpace[i];
    }

	int* getRowSpace() const {
		return RowSpace;
	}

    /** return the index of the row space of the i-th matrix */
    int rowSpaceOfMat(int i) const
    {
    	if (i >= N_Matrices)
    	{
    		ErrMsg("Array out of bound: i >= N_Matrices.");
    	}
    	return RowSpace[i];
    }


};

/******************************************************************************/
//
// Routines for initializing discrete forms
//
/******************************************************************************/

void InitializeDiscreteForms(  
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormSDFEM,
  TDiscreteForm3D *&DiscreteFormUpwind,
  TDiscreteForm3D *&DiscreteFormSGS,
  TDiscreteForm3D *&DiscreteFormVMSProjection,
  TDiscreteForm3D *&DiscreteFormNLGalerkin,
  TDiscreteForm3D *&DiscreteFormNLSDFEM,
  TDiscreteForm3D *&DiscreteFormNLUpwind,
  TDiscreteForm3D *&DiscreteFormNLSGS,
  TDiscreteForm3D *&DiscreteFormNLVMSProjection,
  TDiscreteForm3D *&DiscreteFormNLSDFEM_DivDiv,
  TDiscreteForm3D *&DiscreteFormPressSep,
  TDiscreteForm3D *&DiscreteFormPressSepAuxProb,
  TDiscreteForm3D *&DiscreteFormNSRFBRhs,
  CoeffFct3D *LinCoeffs, int NSTYPE);

void InitializeDiscreteFormsFriction(  
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormSDFEM,
  TDiscreteForm3D *&DiscreteFormUpwind,
  TDiscreteForm3D *&DiscreteFormSGS,
  TDiscreteForm3D *&DiscreteFormNLGalerkin,
  TDiscreteForm3D *&DiscreteFormNLSDFEM,
  TDiscreteForm3D *&DiscreteFormNLUpwind,
  TDiscreteForm3D *&DiscreteFormNLSGS,
  TDiscreteForm3D *&DiscreteFormGalerkinFriction,
  TDiscreteForm3D *&DiscreteFormNLGalerkinFriction,
  CoeffFct3D *LinCoeffs, int NSTYPE);

void InitializeDiscreteForms(  
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormUpwind,
  TDiscreteForm3D *&DiscreteFormUpwindNC,
  TDiscreteForm3D *&DiscreteFormSmagorinsky,
  TDiscreteForm3D *&DiscreteFormClassicalLES,
  TDiscreteForm3D *&DiscreteFormGL00Convolution, 
  TDiscreteForm3D *&DiscreteFormGL00AuxProblem, 
  TDiscreteForm3D *&DiscreteFormVMS_Projection, 
  TDiscreteForm3D *&DiscreteFormVMS_SUPG, 
  TDiscreteForm3D *&DiscreteFormLerayAlpha,
  TDiscreteForm3D *&DiscreteFormNLGalerkin,
  TDiscreteForm3D *&DiscreteFormNLUpwind, 
  TDiscreteForm3D *&DiscreteFormNLUpwindNC, 
  TDiscreteForm3D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm3D *&DiscreteFormNLClassicalLES, 
  TDiscreteForm3D *&DiscreteFormNLGL00Convolution,
  TDiscreteForm3D *&DiscreteFormNLGL00AuxProblem,
  TDiscreteForm3D *&DiscreteFormNLVMS_Projection, 
  TDiscreteForm3D *&DiscreteFormNLVMS_ProjectionExpl, 
  TDiscreteForm3D *&DiscreteFormNSRFBExplRhs,
  TDiscreteForm3D *&DiscreteFormNLVMS_SUPG,
  TDiscreteForm3D *&DiscreteFormNLLerayAlpha,
  TDiscreteForm3D *&DiscreteFormRHS,
  TDiscreteForm3D *&DiscreteFormRHSClassicalLES,
  TDiscreteForm3D *&DiscreteFormRHSLES,
  TDiscreteForm3D *&DiscreteFormRHSSUPG,
  TDiscreteForm3D *&DiscreteFormMatrixGL00AuxProblem,
  TDiscreteForm3D *&DiscreteFormGL00AuxProblemRHS,
  TDiscreteForm3D *&DiscreteFormMatrixAuxProblemU,
  TDiscreteForm3D *&DiscreteFormRHSAuxProblemU,
  TDiscreteForm3D *&DiscreteFormRHSNewton,
  TDiscreteForm3D *&DiscreteFormRHSNewtonNL,			
  TDiscreteForm3D *&DiscreteFormC, 
  TDiscreteForm3D *&DiscreteFormJ,
  CoeffFct3D *LinCoeffs, int NSTYPE);

void InitializeDiscreteFormsVMS(  
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormUpwind,
  TDiscreteForm3D *&DiscreteFormSmagorinsky,
  TDiscreteForm3D *&DiscreteFormNLGalerkin,
  TDiscreteForm3D *&DiscreteFormNLUpwind, 
  TDiscreteForm3D *&DiscreteFormNLSmagorinsky,
  TDiscreteForm3D *&DiscreteFormRHS,
  TDiscreteForm3D *&DiscreteForm_ho_RHS,
  CoeffFct3D *LinCoeffs, int NSTYPE);
  
void InitializeDiscreteFormsOS_ST(
  TDiscreteForm3D *&DiscreteFormGalerkin,
  TDiscreteForm3D *&DiscreteFormStiffRhsOS_ST,
  CoeffFct3D *LinCoeffs, int NSTYPE);

void InitializeDiscreteForms
     ( TDiscreteForm3D *&DiscreteFormGalerkin,
       TDiscreteForm3D *&DiscreteFormNLGalerkin,
       CoeffFct3D *LinCoeffs, int NSTYPE);
       
void InitializeDiscreteForms
     ( TDiscreteForm3D *&DiscreteFormGalerkin,
       TDiscreteForm3D *&DiscreteFormNLGalerkin,
       TDiscreteForm3D *&DiscreteFormRHS,
       CoeffFct3D *LinCoeffs, int NSTYPE);
       
void InitializeDiscreteForms (TDiscreteForm3D *&DiscreteForm, CoeffFct3D *LinCoeff);

void InitializeDiscreteFormsScalar(TDiscreteForm3D *&DiscreteFormMRhs_Galerkin, TDiscreteForm3D *&DiscreteFormARhs_Galerkin, 
                                   TDiscreteForm3D *&DiscreteFormRhs, CoeffFct3D *LinCoeff);
#endif

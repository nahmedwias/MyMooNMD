/** =======================================================================
 * @(#)LocalAssembling2D.h        10.03.2015
 * 
 * @Class:    LocalAssembling2D
 * Purpose:   Assemble on one cell a couple of bilinear and linear forms. That
 *            means the loop over all quadrature points is done within this 
 *            class.
 * 
 *            Furthermore this class includes the computation of function values
 *            at quadrature points, where the function is given by some finite
 *            element function. 
 * 
 * @note This class replaces former TDiscreteForm2D and TAuxParam2D
 * 
 * @Author:      Ulrich Wilbrandt (10.03.2015)
 * 
 * History:     start of implementation 10.03.2015 (Ulrich Wilbrandt)
 * 
 * =======================================================================
 */
#ifndef __LOCAL_ASSEMBLING_2D__
#define __LOCAL_ASSEMBLING_2D__

#include <Enumerations.h>
#include <Constants.h>
#include <string>
#include <vector>
#include <AuxParam2D.h>
#include <DiscreteForm2D.h>

/**
 * A local assembling object has a type associated with it. The type
 * determines, which problem type the object can be used for.
 *
 * Assembling routines throughout ParMooN might check if they
 * get an assembling object of the correct type.
 *
 * So far there is nine built-in types in 3D and one custom type.
 */
enum LocalAssembling2D_type { ConvDiff,                              
                              TCD2D, // stiffness matrix and rhs
                              TCD2D_Mass,// mass matrix, (+ K matrix in case of SUPG)
                              NSE2D_Galerkin,
                              NSE2D_Galerkin_Nonlinear,
                              TNSE2D,
                              TNSE2D_NL,
                              TNSE2D_Rhs,
                              Darcy2D_Galerkin,
                              RECONSTR_GALERKIN,
                              RECONSTR_GALERKIN_Rhs,
                              Custom /// Customized local assembling object. To be used with non-standard problems.
};

/** a function from a finite element space */
class LocalAssembling2D
{
  protected:

    /** @brief The type of problem this assembling objects is made for. */
    const LocalAssembling2D_type type;
    
    /** name */
    std::string name;

    /** @brief number of terms */
    int N_Terms;

    /** @brief number of involved spaces (typically one or two) */
    int N_Spaces;

    /** @brief for each space we store a bool indicatin if second derivatives 
     *         are needed */
    bool *Needs2ndDerivatives;

    /** @brief multiindices for derivatives of ansatz and test functions 
     * 
     * This is an array of size N_Terms.
     */
    std::vector<MultiIndex2D> Derivatives;

    /** @brief for each term, there is one FESpace2D asociated with that term */
    std::vector<int> FESpaceNumber;

    /** @brief which FE space corresponds to each row */
    std::vector<int> RowSpace;

    /** @brief which FE space corresponds to each column */
    std::vector<int> ColumnSpace;

    /** @brief which FE space corresponds to each right-hand side */
    std::vector<int> RhsSpace;

    /** function for calculating the coefficients */
    CoeffFct2D *Coeffs;

    /** @brief function doing the real assembling using parameters from 
     *         argument list */
    AssembleFctParam2D *AssembleParam;

    /** function for manipulating the coefficients */
    ManipulateFct2D *Manipulate;

    /** memory for storing the original value arrays */
    double ***AllOrigValues;

    /** memory for storing the original value arrays at one point */
    double **OrigValues;
    
    int N_Matrices;
    int N_Rhs;
    
    
    /** number of stored parameter functions (ParamFct) */
    int N_ParamFct;
    
    /** array of stored parameter function */
    std::vector<ParamFct*> ParameterFct;
    
    /** index of first parameter produced by parameter function i */
    std::vector<int> BeginParameter;
    
    // number of parameters
    int N_Parameters;
    
    /** number of FE values */
    int N_FEValues;
    
    /** array of stored FEFunction2D */
    TFEFunction2D **FEFunctions2D;
        
    /** index of FEFunction2D used for FE value i */
    std::vector<int> FEValue_FctIndex;
    
    /** which multiindex is used for FE value i */
    std::vector<MultiIndex2D> FEValue_MultiIndex;

    /** Depending on the NSTYPE and the NSE_NONLINEAR_FORM all parameters are 
     * set within this function. This function is called from the constructor 
     * in case of Navier-Stokes problems. It only exists in order to not make 
     * the constructor huge. 
     * 
     * Basically this function implements four nested switches (discretization 
     * type, NSTYOE, Laplace type, nonlinear form type)
     */
    void set_parameters_for_nse(LocalAssembling2D_type type);
    /** 
     * setting every thing for the time dependent Navier-Stokes
     * problems
     */
    void set_parameters_for_tnse(LocalAssembling2D_type type);
    /**
     * parameters and local assembling functions for the 
     * Stokes and Navier-Stokes equations 
     */
    void set_parameters_for_Rec_nse(LocalAssembling2D_type type);
  public:
    /** constructor */
    LocalAssembling2D(LocalAssembling2D_type type, TFEFunction2D **fefunctions2d,
                      CoeffFct2D *coeffs);
    
    /** @brief constructor for backward compatibility
     * 
     * This uses the deprecated classes TAuxParam2D and TDiscreteForm2D to 
     * construct an object of this class.
     * 
     * The member variable 'type' is likely to be wrong using this constructor!
     * Please do not use it, unless you know what you are doing.
     * Type, aux and df will need proper tuning for this to work.
     *
     * @param[in] type The problem type this assemlbing object will be used for.
     * @param[in] aux A (deprecated) "aux object".
     * @param[in] df  A (deprecated) discrete form object.
     */
    LocalAssembling2D(LocalAssembling2D_type type,
                      const TAuxParam2D& aux, const TDiscreteForm2D& df);

    /*!
     * @brief Customized constructor.
     *
     * Set the data members individually. This is to be used when assembling a very problem specific
     * set of matrices and vectors, which is not (or not yet) covered by any of the types from LocalAssembling2D_type.
     *
     * //FROM HERE IT IS MEMBERS OF DEPRECATED TDiscreteForm2D.
     *
     * @param myN_Terms Number of terms to be assembled totally. Determines the size of AllOrigValues and OrigValues and should
     * 					equal the size of "myDerivatives" and "myFESpaceNumber".
     *
     * @param myDerivatives Stores which derivative of the ansatz functions is to be used in which term.
     * 						e.g. myDerivatives[0] = D10 - In the first term the first order x derivative
     * 						 of the ansatz function is to be used
     * 						Corresponds to the "Derivatives" of deprecated TDiscreteForm2D.
     *
     * @param myFESpaceNumber Determines which FE space from ? is to be used for the ansatz (?) function in which term.
     * 						  e.g. myFESpaceNumber[0] = 1 - Use fe space "1" (from where?) for the ansatz function in the first term.
     * 						  Corresponds to the "FESpaceNumber" of deprecated TDiscreteForm2D.
     * 						  NOTE: myFESpaceNumber[0] = 1 means: for first term, use " BaseFuncts[1]" from  "BaseFunct2D *BaseFuncts",
     * 						  which is handed as a parameter to LocalAssembling2D::GetLocalForms(...) by the assembling routine.
     * 						  TO DETERMINE THIS CORRECTLY REQUIRES KNOWLEDGE ABOUT THE BASEFUNCTIONS KNOWN IN THE ASSEMBLING ROUTINE -
     * 						  THOSE ARE TAKEN CELLWISE FROM THE SPACES GIVEN TO ASSEMBLE2D AS 2ND ARGUMENT!
     *
     * @param myRowSpace 	Stores which is the space used for the rows of matrix i.
     * 						e.g. myRowSpace[0]=1 - matrix 0 uses the space "1" (see above) as row (ansatz?) space
     *
     * @param myColumnSpace Stores which is the space used for the columns of matrix i.
     * 						e.g. myColumnSpace[0]=1 - matrix 0 uses the space "1" (see above) as column (test?) space
     *
     * @param myRhsSpace Stores which is the space used for the rhs.
     * 					 e.g. myRhsSpace[1] = 0 - rhs 1 uses the space "0" (see above)
     *
     * @param myCoeffs Pointer to the coefficient function used in "GetLocalForms" to (determine) the coefficients.
     *
     * @param myAssembleParam  Pointer to the assembling function used in "GetLocalForms" to do the entire assembling at one quad point.
     *
     * @param myManipulate	Manipulate function. Do not use this, pass nullptr.
     *
     * @param myN_Matrices How many matrices are to be assembled at once.
     * 					   NOTE: Actually this is never used, because the number of matrices (split into "square" and "rect")
     * 					   		 is given as a parameter to Assemble2D(...) and that values is used.
     * 					   		 This value could only be used to control if "myRowSpace" and "myColumnSpace" have the right length.
     *
     * @param myN_Rhs 	How many right hand sides are to be assembled at once.
     *				 	NOTE: Actually this is never used, because the number of matrices (split into "square" and "rect")
     * 					      is given as a parameter to Assemble2D(...) and that values is used.
     * 					   	  This value could only be used to control if "myRowSpace" and "myColumnSpace" have the right length.
     *
     * //FROM HERE IT IS MEMBERS OF DEPRECATED TAuxParam2D.
     * @param myN_ParamFct The number of Parameter functions working on the assembling.
     *
     * @param myParameterFct A list of pointers to parameter functions. Size should equal "myN_ParamFct".
     *
     * @param myBeginParameter Determine, where which parameter function starts working. Size should equal "myN_ParamFct".
     *
     * @param myN_Parameters	The number of parameters given out in the end. Could differ from myN_FEValues if
     * 							additionally e.g. space coordinates are given out.
     *
     * @param myFEFunctions2D An array of the FE functions which have to be evaluated to get the FE_Values.
     *
     * @param myN_FEValues	The number of parameters to be evaluated directly from FE functions.
     *
     * @param myFEValue_FctIndex Stores which FE function in "myFEFunctions2D" has to be evaluated in order to get a FE_Value.
     * 							 e.g. myFEValue_FctIndex[i]=j - evaluate "myFEFunctions2D[j]" for FE_Value i
     *
     * @param myFEValue_MultiIndex Stores which derivative of an FE function has to be evaluated in order to get a FE_Value.
     * 								e.g. myFEValue_MultiIndex[i]=D01 - for FE_Value i evaluate y-derivative of corresp. FE function
     *
     * @note Some data members get an extra treatment - "name" is set to CUSTOMIZED,
	 * The auxiliary arrays (All)OrigValues are dynamically allocated with size "N_Terms".
	 * "N_Spaces" is determined by finding the max in "FESpaceNumber" (+1).
	 * "Needs2ndDerivative" is dynamically allocated to the size "N_Spaces" and then filled
	 * according to the appearance of "D20", "D11" or "D02" in "Derivatives".
     */
    LocalAssembling2D(int myN_Terms,
    		 std::vector<MultiIndex2D> myDerivatives, std::vector<int> myFESpaceNumber,
			 std::vector<int> myRowSpace, std::vector<int> myColumnSpace, std::vector<int> myRhsSpace,
			 CoeffFct2D* myCoeffs, AssembleFctParam2D* myAssembleParam, ManipulateFct2D* myManipulate,
			 int myN_Matrices, int myN_Rhs,
			 int myN_ParamFct, std::vector<ParamFct*> myParameterFct, std::vector<int> myBeginParameter, int myN_Parameters,
			 TFEFunction2D **myFEFunctions2D,  int myN_FEValues,
			 std::vector<int> myFEValue_FctIndex, std::vector<MultiIndex2D> myFEValue_MultiIndex);

    /** destructor */
    ~LocalAssembling2D();
    
    

    /** return local stiffness matrix */
    void GetLocalForms(int N_Points, double *weights, double *AbsDetjk,
                       double *X, double *Y,
                       int *N_BaseFuncts, BaseFunct2D *BaseFuncts, 
                       double **Parameters, double **AuxArray,
                       TBaseCell *Cell, int N_Matrices, int N_Rhs,
                       double ***LocMatrix, double **LocRhs,
                       double factor = 1.);
    void GetLocalForms(int N_Points, double *weights, double *AbsDetjk,
                       double *X, double *Y, int *N_BaseFuncts,
                       BaseFunct2D *BaseFuncts, TBaseCell *Cell,
                       double ***LocMatrix, double **LocRhs,
                       double factor = 1.);
    
    
    /** return all parameters at all quadrature points */
    void GetParameters(int n_points, TCollection *Coll,
                       TBaseCell *cell, int cellnum,
                       double *x, double *y,
                       double **Parameters);

    /** return all parameters at boundary points */
    void GetParameters(int N_Points, TCollection *Coll,
                       TBaseCell *cell, int cellnum,
                       double *s, int joint,
                       double **Parameters);
    
    

    /** return name */
    const std::string& get_name() const
    { return name; }

    /** return array Needs2ndDerivatives */
    bool *GetNeeds2ndDerivatives() const
    { return Needs2ndDerivatives; }

    /** function for calculating the coefficients */
    CoeffFct2D *GetCoeffFct() const
    { return Coeffs; }
    
    /** return the index of the row space of the i-th matrix */
    int rowSpaceOfMat(int i) const
    { return RowSpace[i]; }
    
    /** return the index of the column space of the i-th matrix */
    int colSpaceOfMat(int i) const
    { return ColumnSpace[i]; }
    
    int GetN_ParamFct() const
    { return N_ParamFct; }
    
    int GetN_Parameters() const
    { return N_Parameters; }
    
    TFEFunction2D* get_fe_function(int i) const
    { return FEFunctions2D[i]; }
    
    LocalAssembling2D_type get_type() const
    { return type; }
};


#endif // __LOCAL_ASSEMBLING_2D__

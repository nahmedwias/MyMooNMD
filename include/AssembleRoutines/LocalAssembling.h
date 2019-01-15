/** =======================================================================
 * @(#)LocalAssembling.h        09.06.2015
 * 
 * @Class:    LocalAssembling
 * Purpose:   Assemble on one cell a couple of bilinear and linear forms. That
 *            means the loop over all quadrature points is done within this 
 *            class.
 * 
 *            Furthermore this class includes the computation of function values
 *            at quadrature points, where the function is given by some finite
 *            element function. 
 * 
 * @Author:      Naveed, Ulrich (09.06.2015)
 * 
 * History:     start of implementation 09.06.2015 (Naveed, Ulrich)
 * 
 * =======================================================================
 */
#ifndef __LOCAL_ASSEMBLING_3D__
#define __LOCAL_ASSEMBLING_3D__

#include <Enumerations.h>
#include <Constants.h>
#include <ParameterDatabase.h>
#include "templateNames.h"
#include <string>
#include <vector>

/**
 * A local assembling object has a type associated with it. The type
 * determines, which problem type the object can be used for.
 *
 * Assembling routines throughout ParMooN might check if they
 * get an assembling object of the correct type.
 *
 * So far there is three built-in types in 3D and one custom type.
 */
enum class LocalAssembling_type
 {
    Brinkman3D_Galerkin,
    ResidualStabPkPk_for_Brinkman3D_Galerkin1,
    GradDivStab_for_Brinkman3D_Galerkin1,
    Brinkman2D_Galerkin1,
    Brinkman2D_Galerkin2,
    Brinkman2D_Galerkin1ResidualStabP1,
    Brinkman2D_Galerkin1ResidualStabP2,
    Brinkman2D_GradDivStabilization,
    ConvDiff, /// Stationary convection diffusion reaction in 3D
    Darcy, // stationary Darcy problem (mixed form)
    TCDStiffMassRhs,
    TCDStiffRhs,
    NSE3D_Linear,    /// Linear part of stationary Navier--Stokes in 3D
    NSE3D_NonLinear, /// Non-linear part of stationary Navier--Stokes in 3D
    TNSE3D_LinGAL,   /// Linear part of time-dependent NS in 3D
    TNSE3D_NLGAL,    /// Non-linear part of time-dependant NS in 3D
    TNSE3D_Rhs,      /// Rhs part of time-dependent NS in 3D
    Custom /// Assembling object created with a custom constructor, probably for a non-standard proble
};
// this function allows you to write these (strongly typed) enums to stdout
// see http://stackoverflow.com/questions/3342726/c-print-out-enum-value-as-text
std::ostream& operator<<(std::ostream& out, const LocalAssembling_type value);

template <int d>
class LocalAssembling
{
  using FEFunction = typename Template_names<d>::FEFunction;
  using MultiIndex_vector = typename Template_names<d>::MultiIndex_vector;
  using CoeffFct = typename Template_names<d>::CoeffFct;
  using ManipulateFct = typename Template_names<d>::ManipulateFct;
  using BaseFunct = typename Template_names<d>::BaseFunct;
  
  protected:
    /// @brief a local parameter database
    ParameterDatabase db;
    
    /** @brief The type of problem this assembling objects is made for. */
    const LocalAssembling_type type;

    /** an integer to identify the space discretization_type **/
    int discretization_type;

    /** @brief number of terms */
    int N_Terms;

    /** @brief number of involved spaces (typically one or two) */
    int N_Spaces;

    /** @brief for each space we store a bool indicating if second derivatives
     *         are needed */
    bool *Needs2ndDerivatives;

    /** @brief multiindices for derivatives of ansatz and test functions 
     * 
     * This is an array of size N_Terms.
     */
    MultiIndex_vector Derivatives;

    
     /** @brief for each term, there is one FESpace3D asociated with that term */
    std::vector<int> FESpaceNumber;

    /** @brief which FE space corresponds to each row */
    std::vector<int> RowSpace;

    /** @brief which FE space corresponds to each column */
    std::vector<int> ColumnSpace;

    /** @brief which FE space corresponds to each right-hand side */
    std::vector<int> RhsSpace;

    /** function for calculating the coefficients */
    CoeffFct Coeffs;

    /** @brief function doing the real assembling using parameters from 
     *         argument list */
    std::vector<AssembleFctParam> local_assemblings_routines;

    /** function for manipulating the coefficients */
    ManipulateFct *Manipulate;

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
    
    /** array of stored finite element functions */
    FEFunction **FEFunctions3D;
        
    /** index of FEFunction3D used for FE value i */
    std::vector<int> FEValue_FctIndex;
    
     /** which multiindex is used for FE value i */
    MultiIndex_vector FEValue_MultiIndex;
    
    /** values of parameter functions stored for computing local forms */
    std::vector<double> parameter_functions_values;
    
    static constexpr int n_local_coefficients = 40;
    std::vector<std::array<double, n_local_coefficients>> local_coefficients;
    
    void set_parameters_for_tcd(LocalAssembling_type type);

    /** Depending on the NSTYPE and the 'nse_nonlinear_form' all parameters are 
     * set within this function. This function is called from the constructor 
     * in case of Navier-Stokes problems. It only exists in order to not make 
     * the constructor huge. 
     * 
     * Basically this function implements four nested switches (discretization 
     * type, NSTYOE, Laplace type, nonlinear form type)
     */
    void set_parameters_for_nse( LocalAssembling_type type);
    
    /** Depending on the NSTYPE and SC_NONLIN_ITE_TYPE_SADDLE all parameters are 
     * set within this function. 
     * 
     * For different discretization schemes: we tried to use different functions
     * in order to keep the function definition smaller.
     */
    /// standard case
    void set_parameters_for_tnse( LocalAssembling_type type);
    /// SMAGORINSKY model
    void set_parameters_for_tnse_smagorinsky( LocalAssembling_type type);
  public:
    /** Constructs a Local Assembling object of a certain type from an array
     *  of fe functions and coefficient functions.
     *
     * @param[in] param_db controlling various parts of this class
     * @param[in] type The type of problem this assembling object can be used
     *            for. Must not be "Custom" - program terminates.
     * @param fefunctions3d The fe  functions to be evaluated.
     * @param coeffs A function pointer to the problem coefficients. These must
     * be hard-coded somewhere, usually in the used example file.
     *
     */
    LocalAssembling(ParameterDatabase param_db, LocalAssembling_type type,
                    FEFunction **fefunctions3d, const CoeffFct& coeffs,
                    int disctype = 1);
    
    /** @brief custom constuctor setting all variables 
     * 
     * Only use this if you know what you are doing. 
     * 
     * @brief Customized constructor.
     *
     * Set the data members individually. This is to be used when assembling a very problem specific
     * set of matrices and vectors, which is not (or not yet) covered by any of the types from LocalAssembling2D_type.
     *
     * //FROM HERE IT IS MEMBERS OF DEPRECATED TDiscreteForm2D.
     *
     * @param myN_Terms Number of terms to be assembled totally. Determines the
     * size of AllOrigValues and OrigValues and should equal the size of
     * "myDerivatives" and "myFESpaceNumber".
     *
     * @param myDerivatives Stores which derivative of the ansatz functions is
     * to be used in which term. e.g. myDerivatives[0] = D10 - In the first term
     * the first order x derivative of the ansatz function is to be used.
     *
     * @param myFESpaceNumber Determines which FE space is to be used for the
     * function in each term. e.g. myFESpaceNumber[0] = 1 - Use fe space "1" 
     * for the ansatz function in the first term.
     * NOTE: myFESpaceNumber[0] = 1 means: for first term, use " BaseFuncts[1]"
     * from  "BaseFunct2D *BaseFuncts", which is handed as a parameter to 
     * LocalAssembling::GetLocalForms(...) by the assembling routine.
     * TO DETERMINE THIS CORRECTLY REQUIRES KNOWLEDGE ABOUT THE BASEFUNCTIONS
     * KNOWN IN THE ASSEMBLING ROUTINE - THOSE ARE TAKEN CELLWISE FROM THE
     * SPACES GIVEN TO ASSEMBLE2D AS 2ND ARGUMENT!
     *
     * @param myRowSpace Stores which is the space used for the rows of
     * matrix i. e.g. myRowSpace[0]=1 - matrix 0 uses the space "1" (see above)
     * as row space
     *
     * @param myColumnSpace Stores which is the space used for the columns of
     * matrix i. e.g. myColumnSpace[0]=1 - matrix 0 uses the space "1"
     * (see above) as column space
     *
     * @param myRhsSpace Stores which is the space used for the rhs.
     * e.g. myRhsSpace[1] = 0 - rhs 1 uses the space "0" (see above)
     *
     * @param myCoeffs Pointer to the coefficient function used in
     * "GetLocalForms" to (determine) the coefficients.
     *
     * @param myAssembleParam Pointer to the assembling function used in
     * "GetLocalForms" to do the entire assembling at one quad point.
     *
     * @param myManipulate Manipulate function. pass nullptr if you do not know
     * what this does.
     *
     * @param myN_Matrices How many matrices are to be assembled at once.
     * NOTE: Actually this is never used, because the number of matrices (split
     * into "square" and "rect") is given as a parameter to Assemble2D(...) and
     * that values is used. This value could only be used to control if
     * "myRowSpace" and "myColumnSpace" have the right length.
     *
     * @param myN_Rhs How many right hand sides are to be assembled at once.
     * NOTE: Actually this is never used, because the number of matrices (split
     * into "square" and "rect") is given as a parameter to Assemble2D(...) and
     * that values is used. This value could only be used to control if
     * "myRowSpace" and "myColumnSpace" have the right length.
     *
     * @param myN_ParamFct The number of Parameter functions working on the
     * assembling.
     *
     * @param myParameterFct A list of pointers to parameter functions. Size
     * should equal "myN_ParamFct".
     *
     * @param myBeginParameter Determine, where which parameter function starts
     * working. Size should equal "myN_ParamFct".
     *
     * @param myN_Parameters	The number of parameters given out in the end.
     * Could differ from myN_FEValues if additionally e.g. space coordinates are
     * given out.
     *
     * @param myFEFunctions2D An array of the FE functions which have to be
     * evaluated to get the FE_Values.
     *
     * @param myN_FEValues The number of parameters to be evaluated directly
     * from FE functions.
     *
     * @param myFEValue_FctIndex Stores which FE function in "myFEFunctions2D"
     * has to be evaluated in order to get a FE_Value. e.g.
     * myFEValue_FctIndex[i]=j - evaluate "myFEFunctions2D[j]" for FE_Value i
     *
     * @param myFEValue_MultiIndex Stores which derivative of an FE function has
     * to be evaluated in order to get a FE_Value. e.g.
     * myFEValue_MultiIndex[i]=D01 - for FE_Value i evaluate y-derivative of
     * corresp. FE function
     *
     * @note Some data members get an extra treatment - The auxiliary arrays
     * (All)OrigValues are dynamically allocated with size "N_Terms". "N_Spaces"
     * is determined by finding the max in "FESpaceNumber" (+1).
     * "Needs2ndDerivative" is dynamically allocated to the size "N_Spaces" and
     * then filled according to the appearance of "D20", "D11" or "D02" (or
     * D200, D110, ..., D002) in "Derivatives".
     */
    LocalAssembling(int myN_Terms,
                    MultiIndex_vector myDerivatives,
                    std::vector<int> myFESpaceNumber,
                    std::vector<int> myRowSpace,
                    std::vector<int> myColumnSpace,
                    std::vector<int> myRhsSpace,
                    CoeffFct myCoeffs, 
                    AssembleFctParam myAssembleParam,
                    ManipulateFct* myManipulate,
                    int myN_Matrices, int myN_Rhs,
                    int myN_ParamFct, std::vector<ParamFct*> myParameterFct,
                    std::vector<int> myBeginParameter, int myN_Parameters,
                    FEFunction** myFEFunctions3D,  int myN_FEValues,
                    std::vector<int> myFEValue_FctIndex, 
                    MultiIndex_vector myFEValue_MultiIndex,
                    int discretization_type = 1);
    
    LocalAssembling(LocalAssembling<d> &la) = delete;
    
    /** destructor */
    ~LocalAssembling();
    
    static ParameterDatabase default_local_assembling_database();
    
    /** return local stiffness matrix */
    void GetLocalForms(int N_Points, const double *weights, double *AbsDetjk,
                       std::array<double*, d> coordinates,
                       int *N_BaseFuncts, BaseFunct *BaseFuncts, 
                       TBaseCell *Cell, int cell_num, int N_Matrices, int N_Rhs,
                       double ***LocMatrix, double **LocRhs,
                       double factor = 1.);
    
     /** return all parameters at all quadrature points */
    void GetParameters(int n_points, TBaseCell *cell, int cellnum,
                       std::array<double*, d> coordinates);

    /** return array Needs2ndDerivatives */
    bool *GetNeeds2ndDerivatives() const
    { return Needs2ndDerivatives; }

    /** function for calculating the coefficients */
    const CoeffFct& GetCoeffFct() const
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
    
    FEFunction* get_fe_function(int i) const
    { return FEFunctions3D[i]; }

    /** @return The type of the local assembling object. */
    LocalAssembling_type get_type() const
    { return type; }

    int get_disctype() const
    { return discretization_type; }
    
    int get_n_rhs() const
    { return N_Rhs; }
    
    void ResetCoeffFct(CoeffFct f)
    { Coeffs = f; }
    
    void SetN_Parameters(int n)
    { N_Parameters = n; }
    
    void setBeginParameter(const std::vector<int>& beginParameter)
    { BeginParameter = beginParameter; }
    
    void setN_ParamFct(int N_paramFct)
    { N_ParamFct = N_paramFct; }
    
    void setParameterFct(const std::vector<ParamFct*>& parameterFct)
    { ParameterFct = parameterFct; }
    
    void setN_FeValues(int N_feValues)
    { N_FEValues = N_feValues; }
    
    void setFeValueFctIndex(const std::vector<int>& feValueFctIndex)
    { FEValue_FctIndex = feValueFctIndex; }
    
    void setFeValueMultiIndex(const MultiIndex_vector& feValueMultiIndex)
    { FEValue_MultiIndex = feValueMultiIndex; }
};

typedef LocalAssembling<3> LocalAssembling3D;
typedef LocalAssembling<2> LocalAssembling2D;

#endif


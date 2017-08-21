/** =======================================================================
 * @(#)LocalAssembling3D.h        09.06.2015
 * 
 * @Class:    LocalAssembling3D
 * Purpose:   Assemble on one cell a couple of bilinear and linear forms. That
 *            means the loop over all quadrature points is done within this 
 *            class.
 * 
 *            Furthermore this class includes the computation of function values
 *            at quadrature points, where the function is given by some finite
 *            element function. 
 * 
 * @note This class replaces former TDiscreteForm3D and TAuxParam3D
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
enum class LocalAssembling3D_type {
    Brinkman3D_Galerkin,
    ResidualStabPkPk_for_Brinkman3D_Galerkin1,
    GradDivStab_for_Brinkman3D_Galerkin1,
    CD3D, /// Stationary convection diffusion reaction in 3D
    TCD3D, // mass matrix (+ matrix comming from time discretization SUPG case)
    TCD3DStiffRhs, // stiffness matrix and right hand side
    NSE3D_Linear,    /// Linear part of stationary Navier--Stokes in 3D
    NSE3D_NonLinear, /// Non-linear part of stationary Navier--Stokes in 3D
    TNSE3D_LinGAL,   /// Linear part of time-dependent NS in 3D
    TNSE3D_NLGAL,    /// Non-linear part of time-dependant NS in 3D
    TNSE3D_Rhs,      /// Rhs part of time-dependent NS in 3D
    Custom /// Assembling object created with a custom constructor, probably for a non-standard proble
};

//Forward declarations
class TFEFunction3D;
class TAuxParam3D;
class TDiscreteForm3D;

/** a function from a finite element space */
class LocalAssembling3D
{
  protected:
    /** @brief The type of problem this assembling objects is made for. */
    const LocalAssembling3D_type type;

    /** an integer to identify the space discretization_type **/
    int discretization_type;

    /** name */
    std::string name;

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
    std::vector<MultiIndex3D> Derivatives;

    
     /** @brief for each term, there is one FESpace3D asociated with that term */
    std::vector<int> FESpaceNumber;

    /** @brief which FE space corresponds to each row */
    std::vector<int> RowSpace;

    /** @brief which FE space corresponds to each column */
    std::vector<int> ColumnSpace;

    /** @brief which FE space corresponds to each right-hand side */
    std::vector<int> RhsSpace;

    /** function for calculating the coefficients */
    CoeffFct3D *Coeffs;

    /** @brief function doing the real assembling using parameters from 
     *         argument list */
    AssembleFctParam3D *AssembleParam;

    /** function for manipulating the coefficients */
    ManipulateFct3D *Manipulate;

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
    
    /** array of stored FEFunction3D */
    TFEFunction3D **FEFunctions3D;
        
    /** index of FEFunction3D used for FE value i */
    std::vector<int> FEValue_FctIndex;
    
     /** which multiindex is used for FE value i */
    std::vector<MultiIndex3D> FEValue_MultiIndex;

    /** Depending on the NSTYPE and the NSE_NONLINEAR_FORM all parameters are 
     * set within this function. This function is called from the constructor 
     * in case of Navier-Stokes problems. It only exists in order to not make 
     * the constructor huge. 
     * 
     * Basically this function implements four nested switches (discretization 
     * type, NSTYOE, Laplace type, nonlinear form type)
     */
    void set_parameters_for_nse(LocalAssembling3D_type type);
    
    /** Depending on the NSTYPE and SC_NONLIN_ITE_TYPE_SADDLE all parameters are 
     * set within this function. 
     * 
     * For different discretization schemes: we tried to use different functions
     * in order to keep the function definition smaller.
     */
    /// standard case
    void set_parameters_for_tnse(LocalAssembling3D_type type);
    /// SMAGORINSKY model
    void set_parameters_for_tnse_smagorinsky(LocalAssembling3D_type type);
  public:
    /** Constructs a Local Assembling object of a certain type from an array
     *  of fe functions and coefficient functions.
     *
     * @param[in] type The type of problem this assembling object can be used
     *            for. Must not be "Custom" - program terminates.
     * @param fefunctions3d The fe  functions to be evaluated.
     * @param coeffs A function pointer to the problem coefficients. These must
     * be hard-coded somewhere, usually in the used example file.
     *
     */
    LocalAssembling3D(LocalAssembling3D_type type, TFEFunction3D **fefunctions3d,
                      CoeffFct3D *coeffs, int disctype = 1);
    
    /** @brief Constructor for backward compatibility
     * 
     * This uses the deprecated classes TAuxParam3D and TDiscreteForm3D to 
     * construct an object of this class. Use of this is discouraged, because the
     * user has to manually tune type, aux and df to make the created object funtioning.
     *
     * @param[in] type The type of problem this assembling object will be used
     *            for. It is the users responsibility to match the type and
     *            the deprecated TAuxParam3D and TDiscreteForm3D object.
     *
     */
    [[ deprecated ]] LocalAssembling3D(LocalAssembling3D_type la_type,
                      TAuxParam3D& aux, TDiscreteForm3D& df);

    /** @brief custom constuctor setting all variables 
     * 
     * Only use this if you know what you are doing. See the respective 
     * constructor in LocalAssembling2D for more documentation.
     */
    LocalAssembling3D(int myN_Terms,
                      std::vector<MultiIndex3D> myDerivatives,
                      std::vector<int> myFESpaceNumber,
                      std::vector<int> myRowSpace,
                      std::vector<int> myColumnSpace,
                      std::vector<int> myRhsSpace,
                      CoeffFct3D* myCoeffs, 
                      AssembleFctParam3D* myAssembleParam,
                      ManipulateFct3D* myManipulate,
                      int myN_Matrices, int myN_Rhs,
                      int myN_ParamFct, std::vector<ParamFct*> myParameterFct,
                      std::vector<int> myBeginParameter, int myN_Parameters,
                      TFEFunction3D** myFEFunctions3D,  int myN_FEValues,
                      std::vector<int> myFEValue_FctIndex, 
                      std::vector<MultiIndex3D> myFEValue_MultiIndex);
    
    /** destructor */
    ~LocalAssembling3D();
    
    /** return local stiffness matrix */
    void GetLocalForms(int N_Points, double *weights, double *AbsDetjk,
                       double *X, double *Y, double *Z,
                       int *N_BaseFuncts, BaseFunct3D *BaseFuncts, 
                       double **Parameters, double **AuxArray,
                       TBaseCell *Cell, int N_Matrices, int N_Rhs,
                       double ***LocMatrix, double **LocRhs,
                       double factor = 1.) const;
    
     /** return all parameters at all quadrature points */
    void GetParameters(int n_points, TCollection *Coll,
                       TBaseCell *cell, int cellnum,
                       double *x, double *y, double *z,
                       double **Parameters) const;

    /** return name */
    const std::string& get_name() const
    { return name; }

    /** return array Needs2ndDerivatives */
    bool *GetNeeds2ndDerivatives() const
    { return Needs2ndDerivatives; }

    /** function for calculating the coefficients */
    CoeffFct3D *GetCoeffFct() const
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
    
    TFEFunction3D* get_fe_function(int i) const
    { return FEFunctions3D[i]; }

    /** @return The type of the local assembling object. */
    LocalAssembling3D_type get_type() const
    { return type; }

    const int get_disctype() const
    { return discretization_type; }
};
#endif

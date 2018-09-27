// =======================================================================
// @(#)FEFunction2D.h        1.3 04/13/00
// 
// Class:       TFEFunction2D
// Purpose:     a function from a finite element space in 2D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
// =======================================================================

#ifndef __FEFUNCTION2D__
#define __FEFUNCTION2D__

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <functional>

class TAuxParam2D;
#include <FESpace2D.h>
//#include <AuxParam2D.h>

/** a function from a finite element space */
class TFEFunction2D
{
  protected:
    /** name of the function */
    std::string Name;

    /** some more words describing the function */
    std::string Description;

    /** space to which this function belongs to */
    const TFESpace2D *FESpace2D;

    /** double vector according to FE isomorphism */
    double *Values;

    /** length of vector */
    int Length;

  public:

    /// Default constructor. All empty object.
    TFEFunction2D();

    /** constructor with vector initialization */
    TFEFunction2D(const TFESpace2D *fespace2D, std::string name, std::string description,
        double *values, int length);

    /// Copy assignment operator. Shallow copy, as the
    /// FEFunction does not take any memory responsibility.
    TFEFunction2D & operator=(const TFEFunction2D & rhs);

    /** destructor */
    ~TFEFunction2D();

    /** return name */
    std::string GetName() const
    { return Name; }

    /** return description */
    std::string GetDescription() const
    { return Description; }

    /** return fe space */
    const TFESpace2D *GetFESpace2D() const
    { return FESpace2D; }
    
    /** return fe space */
    const TFESpace *GetFESpace() const
    { return FESpace2D; }

    /** return length */
    int GetLength() const
    { return Length; }

    /** return vector of data */
    double *GetValues() 
    { return Values; }

    const  double * GetValues() const 
    { return Values; }

    /** @brief calculate errors to given function 
     * 
     * @warning The array \p errors has to be of length at least `N_Errors+1` !!
     */
    void GetErrors(DoubleFunct2D *Exact, int N_Derivatives,
        MultiIndex2D *NeededDerivatives,
        int N_Errors, ErrorMethod2D *ErrorMeth, 
        CoeffFct2D Coeff, TAuxParam2D *Aux,
        int n_fespaces, const TFESpace2D **fespaces,
        double *errors, bool is_SDFEM = 0,
        std::function<bool(const TBaseCell*, int)>funct = [](const TBaseCell*, int){return false;}) const;


    /** @brief use this for vector valued basis functions (Raviart-Thomas (RT)
     *         or Brezzi-Douglas-Marini (BDM) elements) */
    void  GetErrorsForVectorValuedFunction(
        DoubleFunct2D * const * const Exact, 
        ErrorMethod2D * const ErrorMeth, 
        double * const errors) const;

    /** determine the value of function and its first derivatives at
      the given point */
    void FindGradient(double x, double y, double *values) const;

    /** determine the value of function and its first derivatives at
      the given point */
    void FindGradientLocal(const TBaseCell *cell, int cell_no, double x, double y,
        double *values) const ;

    /** determine the value of function at
      the given point */
    void FindValueLocal(const TBaseCell *cell, int cell_no, double x, double y, 
        double *values) const;

    /** calculate the interpolation of an exact function */
    void Interpolate(DoubleFunct2D *Exact);
    /** interpolate the old mesh fe function values to the new fe function
     * 
     * Note that this is rather slow, because no further information is 
     * required. The function 'OldFeFunction' could even live on a larger domain.
     */
    void Interpolate(TFEFunction2D *F);

    /**
     * @brief project this functions into the space L20 (having zero mean value)
     * 
     * After a call to this function the mean value (integral of this function
     * devided by the measure of its domain) has the value a. This is for 
     * example needed for the pressure in a Stokes problem with Dirichlet 
     * conditions on all boundaries.
     * 
     * @param a set mean value of this FEFunction2D to a
     */
    void project_into_L20(double a = 0.0);

    /**
     * @brief find the integral of this function and the measure of its domain
     * 
     * @param integral double value for the integral of this TFEFunction2D
     * @param measure double value for the measure of its domain 
     */
    void compute_integral_and_measure(double& integral, double& measure) const;

    /**
     * @brief compute the mean value of this TFEFunction2D
     * 
     * this functions uses 'compute_integral_and_measure'. Then the mean is the
     * integral divided by the measure.
     */
    double compute_mean() const;    

    /** calculate parameters which are connected to a mesh cell */
    void GetMeshCellParams(DoubleFunct2D *Exact, int N_Derivatives,
        MultiIndex2D *NeededDerivatives,
        int N_Errors, ErrorMethod2D *ErrorMeth, 
        CoeffFct2D Coeff, 
        TAuxParam2D *Aux,
        int n_fespaces, const TFESpace2D **fespaces,
        double *errors, double *parameters);

    /** calculate the super-convergence interpolation of an exact function */
    void InterpolateSuper(DoubleFunct2D *Exact);

    /** set Dirichlet values according to boundary conditions */
    void SetDirichletBC(BoundCondFunct2D *BoundaryCondition,
        BoundValueFunct2D *BoundaryValue);

    /** write the solution into a data file - written by Sashi **/
    void WriteSol(std::string directory=std::string("."),
        std::string basename=std::string("parmoon_solution"));

    /** Read the solution from a given data file - written by Sashi **/
    void ReadSol(std::string BaseName);


    /** sol will be correct to conserve the Old_Mass (remessing, temp, surfact, psd, etc) - added by sashi */   
    void CorrectMass(double OldMass);

    /** Retun the mass, domain volume and mean values of the function - added by sashi */
    void GetMassAndMean(double *OutVal);

    /** multiply function with a scalar alpha. Only non-Dirichlet dofs are 
      multiplied! */
    TFEFunction2D& operator*=(double alpha);

    /** add one TFEFunction2D to another one. Both have to be defined on the 
      same space. Only non-Dirichlet dofs are added!  */
    TFEFunction2D & operator+=(const TFEFunction2D & rhs);

    /** 
     * @brief find the smallest and largest value in the vector representing this 
     * FE-function
     * 
     * If nodal finite elements are used, this is indeed the minimum and maximum
     * of the FE-function. However in other cases this might be wrong (e.g.
     * nonconforming or discontiuous elements).
     */
    void MinMax(double & min, double & max) const;

    /** 
     * @brief print the largest and smallest element in the vector of this FE 
     * function 
     * 
     * This function calls TFEFunction2D::MinMax and prints the given \p name
     * together with the minimum and maximum value. If the string is empty 
     * (default) the used name will be TFEFunction2D::Name.
     */
    void PrintMinMax(std::string name = "") const;

    /** @brief calculate errors to a given function at the boundary (this is of interest, e.g., for the Nitsche method)
     * 
     */
    void GetL2BoundaryError(BoundValueFunct2D *Exact, 
        TAuxParam2D *Aux,
        int n_fespaces, const TFESpace2D **fespaces,
        double *final_boundary_error_l2,
        bool rescale_by_h_E = false);

    /// @brief compute the L^2 norm on a particular boundary component.
    /// Negative values mean that all components are considered.
    double get_L2_norm_on_boundary(int boundary_component = -1) const;
    
    /// @brief compute the L^2 norm.
    double get_L2_norm() const;
};
#endif

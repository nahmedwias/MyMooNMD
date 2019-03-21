// =======================================================================
// %W% %G%
// 
// Class:       TFEFunction3D
// Purpose:     a function from a finite element space in 3D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
// =======================================================================

#ifndef __FEFUNCTION3D__
#define __FEFUNCTION3D__

class TAuxParam3D;
#include <FESpace3D.h>
#include <AuxParam3D.h>
#include <Constants.h>
#include <vector>

/** a function from a finite element space */
class TFEFunction3D
{
  public:
    typedef void ErrorMethod(int, std::array<double*, 3>, double *,
                             const double *, double, double **, double **,
                             double **, double *);
  protected:
    /** name of the function */
    std::string Name;

    /** some more words describing the function */
    std::string Description;

    /** space to which this function belongs to */
    std::shared_ptr<const TFESpace3D> FESpace3D;

    /** double vector according to FE isomorphism */
    double *Values;

    /** length of vector */
    int Length;

  public:
    /// Default constructor. All empty object.
    TFEFunction3D();

    /** constructor with vector initialization */
    TFEFunction3D(std::shared_ptr<const TFESpace3D> fespace3D,
                  const std::string& name, const std::string& description,
                  double *values, int length);

    /// Copy assignment operator. Shallow copy, as the
    /// FEFunction does not take any memory responsibility.
    TFEFunction3D& operator=( const TFEFunction3D & );

    /** destructor */
    ~TFEFunction3D();

    /** return name */
    std::string GetName() const
    { return Name; }

    /** return description */
    std::string GetDescription()
    { return Description; }

    /** return fe space */
    std::shared_ptr<const TFESpace3D> GetFESpace3D() const
    { return FESpace3D; }
    
    /** return fe space */
    std::shared_ptr<const TFESpace> GetFESpace() const
    { return FESpace3D; }

    /** return length */
    int GetLength() const
    { return Length; }

    /** return vector of data */
    const double *GetValues() const
    { return Values; }
    
    /** return vector of data */
    double *GetValues()
    { return Values; }


    /** @brief calculate errors to given function 
     * 
     * @warning The array \p errors has to be of length at least `N_Errors+1` !!
     */
    void GetErrors(DoubleFunct3D *Exact, int N_Derivatives,
                   MultiIndex3D *NeededDerivatives,
                   int N_Errors, ErrorMethod *ErrorMeth, 
                   CoeffFct3D Coeff, TAuxParam3D *Aux,
                   int n_fespaces, const TFESpace3D **fespaces,
                   double *errors) const;
    
    void GetErrorsForVectorValuedFunction(DoubleFunct3D * const * const Exact,
                                          ErrorMethod * const ErrMeth,
                                          double * const errors);

    /** calculate errors to given function */
    void GetMeshCellParams(DoubleFunct3D *Exact, int N_Derivatives,
                   MultiIndex3D *NeededDerivatives,
                   int N_Errors, ErrorMethod *ErrorMeth, 
                   CoeffFct3D Coeff, TAuxParam3D *Aux,
                   int n_fespaces, const TFESpace3D **fespaces,
                   double *errors, double *cell_parameters);

    /**
     * Calculates value of this FEFunction at a given point and the
     * gradient recovery (arithmetic mean of the gradient value in
     *  all containing mesh cells).
     *
     * @param[in] x x value of the point at which to evaluate
     * @param[in] y y value of the point at which to evaluate
     * @param[in] z z value of the point at which to evaluate
     *
     * @param[out] values A vector of length 4 (checked)! Will be filled with
     *            function value
     *            gradient recovery dx
     *            gradient recovery dy
     *            gradient recovery dz
     *
     * @return True if the point was found in the mesh cell collection
     * of this function, false if not so.
     */
    bool FindGradient(double x, double y, double z, std::vector<double>& values) const;

    /** determine the value of function and its first derivatives at
        the given point */
    void FindGradientLocal(const TBaseCell *cell, int cell_no, 
                           double x, double y, double z, 
                           double *values);

    /** determine the value of function
        the given point */
    void FindValueLocal(const TBaseCell *cell, int cell_no, 
                        double x, double y, double z, 
                        double *values) const;

    /** calculate the interpolation of an exact function */
    void Interpolate(DoubleFunct3D *Exact);

   /** calculate the super-convergence interpolation of an exact function */
    void InterpolateSuper(DoubleFunct3D *Exact);
    
    /**  interpolation of an exact function with give FeFunction values */
    void Interpolate(int N_Coord, double *Coords, int N_AuxFeFcts, 
                     TFEFunction3D **AuxFeFcts, DoubleFunctND *Exact);
    
    /** @brief interpolate a vector valued function
     *
     * @param[in] Exact must be of length 3 (= space dimension)
     * @warning EvalAll must be correctly implemented for the used finite element
     */
    void Interpolate_vector_valued_function(
      const std::vector<DoubleFunct3D*>& Exact);
    
    typedef std::function<double(const TBaseCell* cell, int cell_index,
                                 std::array<double, 3> xyz)> AnalyticFunction;
    /**
     * @brief add a given function f to this fe function
     * 
     * The cell and cell_index are from the collection of this TFEFunction3D
     * and the point (x,y,z) is in that cell.
     * 
     * Note that this is similar to creating a second TFEFunction3D, 
     * interpolating f on it and then adding it via operator+=. Here, no second
     * TFEFunction3D is required.
     */
    void add(AnalyticFunction f);
    
    /**
     * @brief project this functions into the space L20 (having zero mean value)
     *
     * After a call to this function the mean value (integral of this function
     * devided by the measure of its domain) has the value a. This is for
     * example needed for the pressure in a Stokes problem with Dirichlet
     * conditions on all boundaries.
     *
     * @param[in] a set mean value of this FEFunction3D to a
     */
    void project_into_L20(double a = 0.0);

    /**
     * @brief find the integral of this function and the measure of its domain
     *
     * In MPI case it returns the global integral and measure, summed up over
     * the own domains of all processes.
     *
     * @param[out] integral double value for the integral of this TFEFunction3D
     * @param[out] measure double value for the measure of its domain
     */
    void compute_integral_and_measure(double& integral, double& measure) const;

    
    
  /** @brief Set Dirichlet values according to boundary conditions */
    void SetDirichletBC(BoundCondFunct3D *BoundaryCondition,
                                   BoundValueFunct3D *BoudaryValue);
    
   /** 
    * @brief find the smallest and largest value in the vector representing this 
    * FE-function
    * 
    * If nodal finite elements are used, this is indeed the minimum and maximum
    * of the FE-function. However in other cases this might be wrong (e.g.
    * nonconforming or discontiuous elements).
    *
    * @todo This kind of min/max computation is nothing but a source of
    * confusion - it would be much better, if this would compute the actual min/max
    * function values!!!
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
   void PrintMinMax(const std::string& name = "") const;
};

#endif

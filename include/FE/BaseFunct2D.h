#ifndef __BASEFUNCT2D__
#define __BASEFUNCT2D__

#include <QuadFormula1D.h>
#include <QuadFormula2D.h>

class TBaseCell;
class TGridCell;


/** set of all base function on the reference element for a finite 
    element in two dimensions */
class TBaseFunct2D
{
  protected:
    /** number of base functions = dimension of local space */
    int Dimension;

    /** Id for this set of base functions */
    BaseFunct2D BaseFunct;

    /** array for all functions and derivatives */
    DoubleFunct2D *Functions[N_MultiIndices2D];

    /** reference element used for this set of base functions */
    BF2DRefElements RefElement;

    /** polynomial degree */
    int PolynomialDegree;

    /** accuracy */
    int Accuracy;

    /** number of basis functions per joint where the sign
        has to be changed if needed */
    int N_BF2Change;

    /** indices of basis functions with changeable sign,
        sorted by joints */
    int const * const * BF2Change;    
 
    /** Dimension of the vector basis function */
    int BaseVectDim;
    
    /** Space dependent basis function, (LPS with Exp bubble) */
    bool SpaceDeptBasis;    
    
  public:

    /** constructor, fill in all information. */
    TBaseFunct2D(int dimension, BaseFunct2D basefunct,
                 BF2DRefElements refelement,
                 DoubleFunct2D* functions, 
                 DoubleFunct2D* derivativesxi,
                 DoubleFunct2D* derivativeseta,
                 DoubleFunct2D* derivativesxixi,
                 DoubleFunct2D* derivativesxieta,
                 DoubleFunct2D* derivativesetaeta,
                 int polynomialdegree,
                 int accuracy,
                 int n_bf2change,
                 int const * const * bf2change,
                 int baseVectDim = 1,
                 bool spaceDeptBasis = false);

    /** return the dimension of local space */
    int GetDimension() const
    { return Dimension; }

    /** return BaseFunct_ID */
    BaseFunct2D GetID() const
    { return BaseFunct; }

    /** return the values for derivative MultiIndex at (xi,eta) */
    void GetDerivatives(MultiIndex2D MultiIndex, double xi,
                        double eta, double *values) const
    { Functions[MultiIndex](xi, eta, values); };

    /** return the values for derivative MultiIndex at all
        quadrature points */
    void GetDerivatives(MultiIndex2D MultiIndex, 
                        TQuadFormula2D *formula, double **values) const;

    /** return values on joint i */
    void GetValues(int N_Points, const double *zeta, int i, double **Values)
      const;

    /** return derivatives on joint i */
    void GetDerivatives(MultiIndex2D MultiIndex, int N_Points,
                        const double *zeta, int i, double **Values) const;

   /** return values of derivative index on joint */
   void GetValues(int N_Points, const double *zeta, int i, 
                  MultiIndex2D index, double **Values) const;

    /** make date on reference element */
    void MakeRefElementData(QuadFormula1D QuadFormula) const;

    /** make date on reference element */
    void MakeRefElementData(QuadFormula2D QuadFormula) const;

    /** generate reference element */
    TGridCell *GenerateRefElement() const;

    /** return reference element */
    BF2DRefElements GetRefElement() const
      { return RefElement; };

    /** return polynomial degree */
    int GetPolynomialDegree() const
      { return PolynomialDegree; };

    /** return accuracy */
    int GetAccuracy() const
      { return Accuracy; };

    /** return number of changeable basis functions per joint */
    int GetN_BF2Change() const
      { return N_BF2Change; }

    /** return array with basis function indices */
    int const * const * GetBF2Change() const
      { return BF2Change; }

    /** change basis functions on cell if needed */
    void ChangeBF(const TCollection *Coll, const TBaseCell *Cell,
                  double *Values) 
      const;

    /** change basis functions on cell in all points if needed */
    void ChangeBF(const TCollection *Coll, const TBaseCell *Cell, int N_Points,
                  double **Values) const;
    
    /** return the dimension of the vector basis function */
    int GetBaseVectDim() const
    { return BaseVectDim; }
};

#endif

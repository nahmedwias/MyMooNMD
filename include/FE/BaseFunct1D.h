#ifndef __BASEFUNCT1D__
#define __BASEFUNCT1D__

#include <QuadFormula1D.h>
#include <Constants.h>
#include <Enumerations.h>

/** set of all base function on the reference element for a finite 
    element in two dimensions */
class TBaseFunct1D
{
  protected:
    /** number of base functions = dimension of local space */
    int Dimension;

    /** Id for this set of base functions */
    BaseFunct1D BaseFunct;

    /** array for all functions and derivatives */
    DoubleFunct1D *Functions[N_MultiIndices1D];

    /** polynomial degree */
    int PolynomialDegree;

    /** accuracy */
    int Accuracy;

  public:
    /** constructor, fill in all information */
    TBaseFunct1D(int dimension, BaseFunct1D basefunct,
                 DoubleFunct1D* functions, 
                 DoubleFunct1D* derivativesxi,
                 DoubleFunct1D* derivativesxixi,
                 int polynomialdegree,
                 int accuracy);

    /** return the dimension of local space */
    int GetDimension() const
    { return Dimension; }

    /** return the values for derivative MultiIndex at xi */
    void GetDerivatives(MultiIndex1D MultiIndex, double xi,
                        double *values) const
      { Functions[MultiIndex](xi, values); };

    /** return the values for derivative MultiIndex at all
        quadrature points */
    void GetDerivatives(MultiIndex1D MultiIndex, 
                        TQuadFormula1D *formula, double **values) const;

    /** make date on reference element */
    void MakeRefElementData(QuadFormula1D QuadFormula) const;

    /** return polynomial degree */
    int GetPolynomialDegree() const
      { return PolynomialDegree; };

    /** return accuracy */
    int GetAccuracy() const
      { return Accuracy; };

};

#endif

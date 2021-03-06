#ifndef __BASEFUNCT3D__
#define __BASEFUNCT3D__

#include <QuadFormula2D.h>
#include <QuadFormula3D.h>

class TBaseCell;
class TGridCell;

/** set of all base function on the reference element for a finite 
    element in two dimensions */
class TBaseFunct3D
{
  protected:
    /** number of base functions = dimension of local space */
    int Dimension;

    /** Id for this set of base functions */
    BaseFunct3D BaseFunct;

    /** array for all functions and derivatives */
    DoubleFunct3D *Functions[N_MultiIndices3D];

    /** reference element used for this set of base functions */
    BF3DRefElements RefElement;

    /** polynomial degree */
    int PolynomialDegree;

    /** accuracy */
    int Accuracy;

    /** number of basis functions per joint where the sign
        has to be changed if needed */
    int N_BF2Change;

    /** indices of basis functions with changeable sign,
        sorted by joints */
    int ***BF2Change;
    
    /** Dimension of the vector basis function */
    int BaseVectDim;

  public:
    /** constructor, fill in all information */
    TBaseFunct3D(int dimension, BaseFunct3D basefunct,
                 BF3DRefElements refelement,
                 DoubleFunct3D* functions, 
                 DoubleFunct3D* derivativesXi,
                 DoubleFunct3D* derivativesEta,
                 DoubleFunct3D* derivativesZeta,
                 DoubleFunct3D* derivativesXiXi,
                 DoubleFunct3D* derivativesXiEta,
                 DoubleFunct3D* derivativesXiZeta,
                 DoubleFunct3D* derivativesEtaEta,
                 DoubleFunct3D* derivativesEtaZeta,
                 DoubleFunct3D* derivativesZetaZeta,
                 int polynomialdegree,
                 int accuracy,
                 int n_bf2change,
                 int ***bf2change,
                 int baseVectDim = 1
                );

    /** constructor without filling data structure */
    explicit TBaseFunct3D(int dimension);

    /** return the dimension of local space */
    int GetDimension() const
    { return Dimension; }
    
    /** return the values for derivative MultiIndex at (xi,eta) */
    void GetDerivatives(MultiIndex3D MultiIndex, double xi,
                        double eta, double zeta, double *values) const
    { Functions[MultiIndex](xi, eta, zeta, values); };

    /** return the values for derivative MultiIndex at all
        quadrature points */
    void GetDerivatives(MultiIndex3D MultiIndex, 
                        TQuadFormula3D *formula, double **values) const;

    /** return values on joint i */
    void GetValues(int N_Points, const double *t, const double *s, 
                   int i, double **Values) const;

   /** return values of derivative index on joint */
   void GetValues(int N_Points, const double *t, const double *s, int i, 
                  MultiIndex3D index, double **Values) const;

    /** make date on reference element */
    void MakeRefElementData(QuadFormula2D QuadFormula) const;

    /** make date on reference element */
    void MakeRefElementData(QuadFormula3D QuadFormula) const;

    /** generate reference element */
    TGridCell *GenerateRefElement() const;

    /** return reference element */
    BF3DRefElements GetRefElement() const
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
    int ***GetBF2Change() const
      { return BF2Change; }

    /** change basis functions on cell if needed */
    void ChangeBF(const TCollection *Coll, const TBaseCell *Cell, double *Values)
      const;

    /** change basis functions on cell in all points if needed */
    void ChangeBF(const TCollection *Coll, const TBaseCell *Cell, int N_Points,
                  double **Values) const;

    /** return BaseFunct_ID */
    BaseFunct3D GetID() const
    { return BaseFunct; }
    
    /** return the dimension of the vector basis function */
    int GetBaseVectDim() const
    { return BaseVectDim; }
};

#endif

// =======================================================================
// @(#)TriaIsoparametric.h        1.5 04/13/00
//
// Class:      TTriaIsoparametric
//
// Purpose:    isoparametric reference transformations for triangle
//
// Author:     Gunar Matthies
//
// History:    29.04.98 start implementation
// 
// =======================================================================

#ifndef __TRIAISOPARAMETRIC__
#define __TRIAISOPARAMETRIC__

#include <RefTrans2D.h>
#include <Enumerations.h>

/** reference transformations for triangle */
class TTriaIsoparametric : public TRefTrans2D
{
  protected:
    /** x coordinate */
    double x[3];

    /** y coordinate */
    double y[3];

    /** x parameters for reference transformation */
    double xc0, xc1, xc2;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2;

    /** number of additional points */
    int N_AuxPoints;

    /** distance in x direction between real auxiliary point and 
        its position after a affine mapping */
    double XDistance[MaxN_BaseFunctions2D];

    /** distance in y direction between real auxiliary point and 
        its position after a affine mapping */
    double YDistance[MaxN_BaseFunctions2D];

    /** order of approximation */
    int ApproximationOrder;

    /** values of corresponding base function at quadpoints */
    double FctValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];

    /** xi-derivatives of corresponding base function at quadpoints */
    double XiDerValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];

    /** eta-derivatives of corresponding base function at quadpoints */
    double EtaDerValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];

    /** base function type for each order of approximation */
    static BaseFunct2D BaseFunctFromOrder[]; 

    /** base function type for each order of approximation */
    static FEDesc2D FEDescFromOrder[]; 

    /** auxiliary array */
    double DoubleAux[MaxN_BaseFunctions2D];

    /** auxiliary array */
    int IntAux[MaxN_BaseFunctions2D];

    /** used quadrature rule */
    QuadFormula2D QuadFormula;
  
    /** for data from quadrature formula */
    const double *XI, *ETA, *W;

    /** number of quadrature points */
    int N_QuadPoints;

  public:
    /** constuctor */
    TTriaIsoparametric();

    /** transfer form reference element to original element */
    void GetOrigFromRef(double eta, double xi, double &x, double &y) override;

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, const double *eta, const double *xi, 
                        double *x, double *y, double *absdetjk) override;

    /** transfer form reference element to original element */
    void GetOrigFromRef(const double *ref, double *orig) override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double &eta, double &xi) override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(const double *orig, double *ref) override;

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(BaseFunct2D BaseFunct,
                       int N_Points, const double *xi, const double *eta,
                       int N_Functs, QuadFormula2D formula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct2D *BaseFunct,
                       int N_Points, const double *xi, const double *eta,
                       QuadFormula2D formula,
                       bool *Needs2ndDer);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, int N_BaseFunct,
                       const double *uref, const double *uxiref,
                       const double *uetaref, 
                       double *uorig, double *uxorig, double *uyorig);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(int joint, double zeta, int N_BaseFunct,
                       const double *uref, const double *uxiref,
                       const double *uetaref,
                       double *uorig, double *uxorig, double *uyorig);

    /** set element to cell */
    void SetCell(const TBaseCell * cell) override;

    /** set order of approximation */
    void SetApproximationOrder(int order)
    {
      if(order <= 0)
        ApproximationOrder = 1;
       else
        ApproximationOrder = order;
    }

    /** set used quadrature formula */
    void SetQuadFormula(QuadFormula2D formula)
    { QuadFormula = formula; }

    /** return outer normal vector */
    void GetOuterNormal(int j, double zeta, double &n1, double &n2) const override;

    /** return tangent */
    void GetTangent(int j, double zeta, double &t1, double &t2) const override;

    /** return volume of cell according to reference transformation */
    double GetVolume() const override;

    /** return boundary vertices */
    void GetOrigBoundFromRef( int joint, int N_Points, double *zeta, double *X, double *Y);

    /** @brief Piola map, needed for vector values basis functions such as 
     *         Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM).
     */
    void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                             const double *refD00, double *origD00) override;
    
    /// @brief Piola map for derivatives
    void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                             const double *refD00, const double *refD10,
                             const double *refD01, double *origD10,
                             double *origD01) override;
};

#endif

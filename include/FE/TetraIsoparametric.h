// =======================================================================
// @(#)TetraIsoparametric.h        1.3 02/22/00
//
// Class:      TTetraIsoparametric
//
// Purpose:    Isoparametric reference transformations for Tetrahedron
//
// Author:     Gunar Matthies
//
// History:    2000/11/20 start implementation
// 
// =======================================================================

#ifndef __TetraIsoparametric__
#define __TetraIsoparametric__

#include <Enumerations.h>
#include <RefTrans3D.h>

/** reference transformations for Tetrahedron */
class TTetraIsoparametric : public TRefTrans3D
{
  protected:
    /** x coordinate */
    double x0, x1, x2, x3;

    /** y coordinate */
    double y0, y1, y2, y3;

    /** z coordinate */
    double z0, z1, z2, z3;

    /** x parameters for reference transformation */
    double xc0, xc1, xc2, xc3;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2, yc3;

    /** z parameters for reference transformation */
    double zc0, zc1, zc2, zc3;

    /** number of additional points */
    int N_AuxPoints;

    /** distance in x direction between real auxiliary point and
        its position after a trilinear mapping */
    double XDistance[MaxN_BaseFunctions3D];

    /** distance in y direction between real auxiliary point and 
        its position after a trilinear mapping */
    double YDistance[MaxN_BaseFunctions3D];

    /** distance in z direction between real auxiliary point and 
        its position after a trilinear mapping */
    double ZDistance[MaxN_BaseFunctions3D];

    /** order of approximation */
    int ApproximationOrder;

    /** values of corresponding base function at quadpoints */
    double FctValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];

    /** xi-derivatives of corresponding base function at quadpoints */
    double XiDerValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];

    /** eta-derivatives of corresponding base function at quadpoints */
    double EtaDerValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];

    /** zeta-derivatives of corresponding base function at quadpoints */
    double ZetaDerValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];

    /** base function type for each order of approximation */
    static BaseFunct3D BaseFunctFromOrder[]; 

    /** base function type for each order of approximation */
    static FEDesc3D FEDescFromOrder[]; 

    /** auxiliary array */
    double DoubleAux[MaxN_BaseFunctions3D];

    /** auxiliary array */
    int IntAux[MaxN_BaseFunctions3D];

    /** used quadrature rule */
    QuadFormula3D QuadFormula;

    /** for data from quadrature formula */
    const double *XI, *ETA, *ZETA, *W;

    /** number of quadrature points */
    int N_QuadPoints;

    /** detjk */
    double detjk;

    /** 1/detjk */
    double rec_detjk;

  public:
    /** constuctor */
    TTetraIsoparametric();

    /** transfer from reference face to original element face*/
    void GetOrigBoundFromRef(int Joint, double xi, double eta,
                             double &X, double &Y, double &Z);
    
    /** transfer from reference face to original element face*/
    void GetOrigBoundFromRef(int Joint, int N_Points, const double *xi,
                             const double *eta, double *X, double *Y,
                             double *Z);


    /** transfer from reference element to original element */
    void GetOrigFromRef(double eta, double xi, double zeta,
                        double &x, double &y, double &z) override;

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, const double *xi, const double *eta,
                        const double *zeta,  double *x, double *y, double *z,
                        double *absdetjk) override;

    /** transfer form reference element to original element */
    void GetOrigFromRef(const double *ref, double *orig) override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double z,
                        double &eta, double &xi, double &zeta) override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(const double *orig, double *ref) override;

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(BaseFunct3D BaseFunct, int N_Points, const double *xi,
                       const double *eta, const double *zeta,
                       int N_Functs, QuadFormula3D QuadFormula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct3D *BaseFunct, int N_Points,
                       const double *xi, const double *eta, const double *zeta,
                       QuadFormula3D QuadFormula,
                       bool *Needs2ndDer);

    void GetOrigValues(int JointNr, double p1, double p2, int N_BaseFunct,
                       const double *uref, const double *uxiref,
                       const double *uetaref, const double *uzetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       double *uzorig);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, double zeta, int N_BaseFunct,
                       const double *uref, const double *uxiref,
                       const double *uetaref, const double *uzetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       double *uzorig);

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
    void SetQuadFormula(QuadFormula3D formula)
    { QuadFormula = formula; }

    /** return outer normal unit vector */
    void GetOuterNormal(int j, double s, double t,
                        double &n1, double &n2, double &n3) const override;

    /** return two tangent vectors */
    void GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23) const override;

    /** @brief Piola map, needed for vector values basis functions such as 
     *         Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM).
     */
    virtual void PiolaMapOrigFromRef(double xi, double eta, double zeta,
                                     int N_Functs, const double *refD00, 
                                     double *origD00) override;
};

#endif

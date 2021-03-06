// =======================================================================
// @(#)HexaTrilinear.h        1.3 02/22/00
//
// Class:      THexaTrilinear
//
// Purpose:    trilinear reference transformations for Hexahedron
//
// Author:     Daniel Quoos  
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#ifndef __HexaTrilinear__
#define __HexaTrilinear__

#include <Enumerations.h>
#include <RefTrans3D.h>

/** reference transformations for Hexahedron */
class THexaTrilinear : public TRefTrans3D
{
  protected:
    /** x coordinate */
    double x0, x1, x2, x3, x4, x5, x6, x7;

    /** y coordinate */
    double y0, y1, y2, y3, y4, y5, y6, y7;

    /** z coordinate */
    double z0, z1, z2, z3, z4, z5, z6, z7;

    /** x parameters for reference transformation */
    double xc0, xc1, xc2, xc3, xc4, xc5, xc6, xc7;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2, yc3, yc4, yc5, yc6, yc7;

    /** z parameters for reference transformation */
    double zc0, zc1, zc2, zc3, zc4, zc5, zc6, zc7;

    /** detjk */
    double detjk;

    /** 1/detjk */
    double rec_detjk;

  public:
    /** constuctor */
    THexaTrilinear();

    /** transfer from reference element to original element */
    void GetOrigFromRef(double eta, double xi, double zeta, double &x,
                        double &y, double &z) override;

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, const double *eta, const double *xi,
                        const double *zeta,  double *x, double *y, double *z,
                        double *absdetjk) override;

    /** transfer form reference element to original element */
    void GetOrigFromRef(const double *ref, double *orig) override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double z, double &eta, double &xi,
                        double &zeta) override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(const double *orig, double *ref) override;

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(BaseFunct3D BaseFunct, int N_Points, const double *xi,
                       const double *eta, const double *zeta,
                       int N_Functs, QuadFormula3D HexaFormula);

    /** calculate functions and derivatives from reference element
        to original element, for all given elements */
    void GetOrigValues(int N_Sets, BaseFunct3D *BaseFunct,
                       int N_Points, const double *xi, const double *eta,
                       const double *zeta, QuadFormula3D HexaFormula,
                       bool *Needs2ndDer);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, double zeta, int N_BaseFunct,
                const double *uref, const double *uxiref, const double *uetaref,
                const double *uzetaref, double *uorig, double *uxorig,
                double *uyorig, double *uzorig, int _BaseVectDim = 1);

    /** set element to cell */
    void SetCell(const TBaseCell * cell) override;

    /** return outer normal unit vector */
    void GetOuterNormal(int j, double s, double t,
                        double &n1, double &n2, double &n3) const override;

    /** return two tangent vectors */
    void GetTangentVectors(int j, double p1, double p2,
        double &t11, double &t12, double &t13,
        double &t21, double &t22, double &t23) const override;
    
    /** @brief Piola transformation for vector valued basis functions **/
    void PiolaMapOrigFromRef(double xi, double eta, double zeta,
                             int N_Functs, const double *refD000,
                             double *origD000) override;
};

#endif
